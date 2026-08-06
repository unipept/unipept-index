use std::error::Error;
use std::fs::{File, OpenOptions, create_dir_all};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rayon::prelude::*;
use sa_index::kmer_table::AMINO_ACID_COUNT;
use sa_index::{sa_searcher::{SearchAllSuffixesResult, Searcher}, KmerTable, SuffixArray, SuffixArrayBackend};
use sa_server::{load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Serialize;
use sysinfo::{Pid, System};
use sa_index::ProteinsBackend as _;
use sa_index::suffix_to_protein_index::SuffixToProteinMappingBackend as _;
use text_compression::ProteinTextBackend as _;

/// Schema version — increment when the output JSON format changes.
/// v2: matrix records aggregate `runs` reps into one line and carry a `stats` spread.
const SCHEMA_VERSION: u32 = 2;

/// Canonical 20 amino acids used for random peptide generation
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// ---------------------------------------------------------------------------
// Warmup mode
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum WarmupMode {
    /// Touch every byte of every mmap-backed region to fully populate the page cache.
    All,
    /// Search the first N peptides from warmup.txt to partially warm the page cache.
    Count(usize),
    /// Touch all pages, then run N peptides through the search + retrieval pipeline.
    /// Recommended for mmap: page-cache warmup alone leaves CPU caches and TLB cold.
    AllThenCount(usize),
}

impl std::str::FromStr for WarmupMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "all" { return Ok(WarmupMode::All); }
        if let Some(rest) = s.strip_prefix("all:") {
            return rest.parse::<usize>()
                .map(WarmupMode::AllThenCount)
                .map_err(|_| format!("expected 'all:<count>', got '{}'", s));
        }
        s.parse::<usize>()
            .map(WarmupMode::Count)
            .map_err(|_| format!("expected 'all', 'all:<N>', or a non-negative integer, got '{}'", s))
    }
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

/// Benchmark the suffix array index by searching randomly generated peptides
/// and writing timing/statistics as JSONL records.
#[derive(Parser, Debug)]
#[command(about = "Benchmark the suffix array index")]
struct Args {
    /// Folder containing sa.bin, proteins.bin, mapping.bin (and warmup.txt if --warmup is set)
    #[arg(short, long)]
    index_dir: PathBuf,

    /// Output folder where .jsonl files will be stored
    #[arg(short, long)]
    output: PathBuf,

    /// Label for this run series — used as the output file name and JSON label field
    #[arg(long, default_value = "benchmark")]
    label: String,

    /// Equate I and L during search
    #[arg(long, action = clap::ArgAction::Set, default_value_t = true)]
    equate_il: bool,

    /// Only return tryptic matches
    #[arg(long, action = clap::ArgAction::Set, default_value_t = false)]
    tryptic: bool,

    /// Maximum number of suffix matches per peptide before the cutoff is applied
    #[arg(long, default_value_t = 10_000)]
    max_matches: usize,

    /// Number of peptides per run.
    /// In random mode: how many random peptides to generate.
    /// In file mode (--peptide-file): how many lines to consume per run.
    #[arg(long, default_value_t = 10_000)]
    amount_of_peptides: usize,

    /// Minimum length of a randomly generated peptide (ignored when --peptide-file is set)
    #[arg(long, default_value_t = 5)]
    peptide_length_min: usize,

    /// Maximum length of a randomly generated peptide (ignored when --peptide-file is set)
    #[arg(long, default_value_t = 50)]
    peptide_length_max: usize,

    /// Path to a pre-generated peptide file (one peptide per line).
    /// When set, peptides are consumed sequentially: run 1 reads lines 0..N, run 2 reads N..2N, etc.
    /// The file must contain at least amount_of_peptides × runs lines.
    /// --peptide-length-min and --peptide-length-max are ignored.
    #[arg(long)]
    peptide_file: Option<PathBuf>,

    /// Number of timed benchmark runs to perform
    #[arg(long, default_value_t = 100)]
    runs: u32,

    /// Seed for random peptide generation (random mode only — no effect with --peptide-file
    /// or --matrix, which read a fixed query stream from disk). Fixes the generated stream so
    /// runs and separate invocations are reproducible. Omit for a fresh nondeterministic stream.
    #[arg(long)]
    seed: Option<u64>,

    /// Warm up the index before timing.
    /// "all": touch every page of every mmap-backed region (populates page cache only).
    /// "all:N": touch every page then push N peptides from warmup.txt through the pipeline
    ///          (recommended for mmap — also warms CPU caches and TLB).
    /// N: push N peptides from warmup.txt through the pipeline (used for preloaded).
    #[arg(long, num_args = 0..=1, default_missing_value = "all")]
    warmup: Option<WarmupMode>,

    /// Optional path to a pre-built k-mer bounds table file.
    /// When provided, binary search is accelerated by narrowing the initial search window.
    #[arg(long)]
    kmer_table_file: Option<PathBuf>,

    /// Build an in-memory k-mer bounds table of size k from the loaded index instead of
    /// loading one from a file. Handy for A/B testing the k-mer acceleration.
    /// Ignored when --kmer-table-file is also given.
    #[arg(long)]
    build_kmer_table: Option<usize>,

    /// Run the full parameter matrix in one process: loads the index once, then iterates
    /// equate_il × tryptic × {scalar,batched} × {none,5-mer,6-mer} for each --matrix-files
    /// entry. Writes one record per (config × run) to <output>/<label>.jsonl.
    #[arg(long)]
    matrix: bool,

    /// Matrix mode: comma-separated peptide files; each becomes one "file" dimension.
    #[arg(long, value_delimiter = ',')]
    matrix_files: Vec<PathBuf>,

    /// Matrix mode: pre-built 5-mer table file (falls back to building one if omitted).
    #[arg(long)]
    kmer5_file: Option<PathBuf>,

    /// Matrix mode: pre-built 6-mer table file (falls back to building one if omitted).
    #[arg(long)]
    kmer6_file: Option<PathBuf>,

    /// Matrix mode: MLP batch sizes to sweep, comma-separated (1 = scalar). e.g. 1,8,16,32.
    #[arg(long, value_delimiter = ',', default_values_t = vec![1usize, 16])]
    matrix_batches: Vec<usize>,
}

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

#[derive(Clone, Serialize)]
struct BenchmarkConfig {
    sa_type: String,
    mapping_type: String,
    use_mmap: bool,
    sample_rate: u8,
    bits_per_value: usize,
    equate_il: bool,
    tryptic: bool,
    max_matches: usize,
    /// MLP batch size for this run (1 = scalar, one peptide per task).
    batch_size: usize,
    /// k of the attached k-mer table (0 = no table).
    kmer_k: usize,
    amount_of_peptides: usize,
    peptide_length_min: usize,
    peptide_length_max: usize,
    peptide_source: String,
}

#[derive(Serialize)]
struct BenchmarkResult {
    search_duration_ns: u64,
    retrieval_duration_ns: u64,
    total_duration_ns: u64,
    throughput_qps: f64,
    amount_of_queries: usize,
    query_hit_count: usize,
    suffix_hit_count: usize,
    protein_hit_count: usize,
    cutoff_reached: bool,
    total_memory: u64,
    theoretical_max_memory: u64,
    /// Nanoseconds spent exclusively inside `search_bounds()` (binary search + k-mer lookup).
    search_bounds_ns: u64,
    /// Nanoseconds spent iterating over the matched suffix range after `search_bounds()` returns.
    match_iter_ns: u64,
}

/// Spread of throughput across the `runs` timed reps of a single config. Emitted in matrix
/// mode so a ±band is visible and cell-to-cell noise isn't read as signal. `result` above still
/// carries the detailed timing of the representative (median-qps) rep.
#[derive(Serialize)]
struct RunStats {
    runs: u32,
    qps_min: f64,
    qps_p10: f64,
    qps_p50: f64,
    qps_p90: f64,
    qps_max: f64,
}

#[derive(Serialize)]
struct BenchmarkRecord {
    version: u32,
    label: String,
    commit: String,
    config: BenchmarkConfig,
    result: BenchmarkResult,
    /// Per-config throughput spread over all reps (matrix mode only; omitted otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    stats: Option<RunStats>,
}

/// Linear-interpolated percentile of an already-sorted (ascending) slice. `p` in [0, 1].
fn percentile(sorted: &[f64], p: f64) -> f64 {
    match sorted.len() {
        0 => 0.0,
        1 => sorted[0],
        n => {
            let rank = p * (n - 1) as f64;
            let lo = rank.floor() as usize;
            let frac = rank - lo as f64;
            sorted[lo] + (sorted[(lo + 1).min(n - 1)] - sorted[lo]) * frac
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Reads the first byte of `path` to detect the mapping type without fully loading the file.
fn first_byte_of(path: &Path) -> Result<u8, Box<dyn Error>> {
    let mut f = File::open(path)?;
    let mut buf = [0u8; 1];
    f.read_exact(&mut buf)?;
    Ok(buf[0])
}

/// Generates `count` random peptides whose length is in `[min_len, max_len]`, drawing from
/// `rng` so the caller controls determinism (seeded → reproducible stream).
fn generate_peptides(rng: &mut impl Rng, count: usize, min_len: usize, max_len: usize) -> Vec<String> {
    (0..count)
        .map(|_| {
            let len = rng.gen_range(min_len..=max_len);
            (0..len)
                .map(|_| AMINO_ACIDS[rng.gen_range(0..AMINO_ACIDS.len())] as char)
                .collect()
        })
        .collect()
}

/// Computes the theoretical in-memory footprint of the loaded index structures.
///
/// This is derived from the actual data sizes, **not** from disk file sizes, so it remains
/// accurate when new structures are added to the `Searcher`. When you add a new structure,
/// extend this function with its memory calculation.
fn theoretical_memory(searcher: &Searcher<SuffixArray>, mapping_type: &str, use_mmap: bool) -> u64 {
    let text_len = searcher.proteins.text().len() as u64;
    let protein_count = searcher.proteins.len() as u64;

    // Suffix array: one entry per SA item at bits_per_value bits each
    let sa_bytes = (searcher.sa.len() as u64 * searcher.sa.bits_per_value() as u64).div_ceil(8);

    // ProteinText: 5 bits per character (BitArray), rounded up to whole bytes
    let text_bytes = (text_len * 5).div_ceil(8);

    // Protein metadata
    let string_bytes: u64 = (0..searcher.proteins.len())
        .map(|i| {
            let p = searcher.proteins.get(i);
            p.uniprot_id.len() as u64 + p.functional_annotations.len() as u64
        })
        .sum();
    let metadata_bytes = if use_mmap {
        // MmapBacked: 16-byte fixed table entry per protein + concatenated string blobs
        protein_count * 16 + string_bytes
    } else {
        // InMemory: Vec<Protein> on heap — each Protein struct is 56 bytes
        // (String=24 + u32=4 + padding=4 + Vec<u8>=24) plus the heap-allocated string data
        protein_count * 56 + string_bytes
    };

    // Suffix-to-protein mapping (see suffix_to_protein_index/{dense,sparse,bitvec}.rs)
    let mapping_bytes = match mapping_type {
        "dense" => text_len * 4, // Vec<u32> with one u32 per text character
        "sparse" => (protein_count + 2) * 8, // Vec<i64> with one i64 per protein boundary
        "bitvec" => {
            // One bit per text character + Rank9 superblock overhead (~16 bytes per 512 bits)
            let bits_bytes = text_len.div_ceil(8);
            let superblock_count = text_len / 512 + 1;
            let rank9_bytes = superblock_count * 16;
            bits_bytes + rank9_bytes
        }
        _ => 0,
    };

    // k-mer table size: 16 bytes per entry, AMINO_ACID_COUNT^k entries total
    let kmer_table_bytes = searcher.kmer_table.as_ref().map_or(0, |t| {
        (AMINO_ACID_COUNT as u64).pow(t.k as u32) * 16
    });

    sa_bytes + text_bytes + metadata_bytes + mapping_bytes + kmer_table_bytes
}

/// Returns the current process's resident set size in bytes via sysinfo.
fn measure_process_memory() -> u64 {
    let pid = Pid::from(std::process::id() as usize);
    let mut sys = System::new();
    sys.refresh_process(pid);
    sys.process(pid).map(|p| p.memory()).unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Benchmark run
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn run_benchmark(searcher: &Searcher<SuffixArray>, peptides: &[String], max_matches: usize, equate_il: bool, tryptic: bool, mlp_batch: usize, theoretical_max_memory: u64, baseline_memory: u64) -> BenchmarkResult {
    // Memory snapshot before any timing starts — captures index-resident pages only
    let index_memory = measure_process_memory().saturating_sub(baseline_memory);

    // Reset per-run timing accumulators before the search phase.
    searcher.drain_timing_ns();

    // Phase 1: suffix array search (parallel), via the same orchestrator production uses.
    // mlp_batch > 1 interleaves that many peptides per rayon task for memory-level parallelism;
    // 1 = scalar one-peptide-at-a-time.
    let refs: Vec<&[u8]> = peptides.iter().map(|p| p.as_bytes()).collect();
    let search_start = Instant::now();
    let suffix_results = searcher.search_all_matching_suffixes(&refs, max_matches, equate_il, tryptic, mlp_batch);
    let search_duration_ns = search_start.elapsed().as_nanos() as u64;

    // Read internal timing breakdown (bounds lookup vs match iteration).
    let (search_bounds_ns, match_iter_ns) = searcher.drain_timing_ns();

    // Phase 2: protein retrieval (parallel) — returns (protein_count, cutoff_reached)
    let retrieval_start = Instant::now();
    let retrieval_stats: Vec<(usize, bool)> = suffix_results
        .par_iter()
        .map(|r| match r {
            SearchAllSuffixesResult::MaxMatches(suf) => {
                (searcher.retrieve_proteins(suf).len(), true)
            }
            SearchAllSuffixesResult::SearchResult(suf) => {
                (searcher.retrieve_proteins(suf).len(), false)
            }
            SearchAllSuffixesResult::NoMatches => (0, false),
        })
        .collect();
    let retrieval_duration_ns = retrieval_start.elapsed().as_nanos() as u64;

    // Aggregate stats
    let query_hit_count = suffix_results
        .iter()
        .filter(|r| !matches!(r, SearchAllSuffixesResult::NoMatches))
        .count();

    let suffix_hit_count: usize = suffix_results
        .iter()
        .map(|r| match r {
            SearchAllSuffixesResult::MaxMatches(suf)
            | SearchAllSuffixesResult::SearchResult(suf) => suf.len(),
            SearchAllSuffixesResult::NoMatches => 0,
        })
        .sum();

    let protein_hit_count: usize = retrieval_stats.iter().map(|(c, _)| *c).sum();
    let cutoff_reached = retrieval_stats.iter().any(|(_, c)| *c);

    let total_duration_ns = search_duration_ns + retrieval_duration_ns;
    let throughput_qps = if total_duration_ns > 0 {
        peptides.len() as f64 / (total_duration_ns as f64 / 1e9)
    } else {
        0.0
    };

    BenchmarkResult {
        search_duration_ns,
        retrieval_duration_ns,
        total_duration_ns,
        throughput_qps,
        amount_of_queries: peptides.len(),
        query_hit_count,
        suffix_hit_count,
        protein_hit_count,
        cutoff_reached,
        total_memory: index_memory,
        theoretical_max_memory,
        search_bounds_ns,
        match_iter_ns,
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Matrix mode: load the index once, then sweep the full config grid for each peptide
/// file. The k-mer table is swapped in/out (moved, not cloned) so the big 6-mer table is
/// built/loaded only once. Each config runs `runs` reps but writes a single aggregated record
/// (median-qps rep as `result`, plus a `stats` spread) to <output>/<label>.jsonl.
#[allow(clippy::too_many_arguments)]
fn run_matrix(
    mut searcher: Searcher<SuffixArray>,
    args: &Args,
    mapping_type: &str,
    sa_type: &str,
    sample_rate: u8,
    bits_per_value: usize,
    baseline_memory: u64,
) -> Result<(), Box<dyn Error>> {
    if args.matrix_files.is_empty() {
        return Err("--matrix requires at least one --matrix-files entry".into());
    }
    let use_mmap = cfg!(feature = "mmap");

    // Build/load the 5-mer and 6-mer tables once (swapped in/out per config below).
    let mut table5: Option<KmerTable> = Some(match &args.kmer5_file {
        Some(p) => { eprintln!("Loading 5-mer table from {}...", p.display()); load_kmer_table_file(p.to_str().unwrap())? }
        None => { eprintln!("Building 5-mer table..."); KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 5) }
    });
    let mut table6: Option<KmerTable> = Some(match &args.kmer6_file {
        Some(p) => { eprintln!("Loading 6-mer table from {}...", p.display()); load_kmer_table_file(p.to_str().unwrap())? }
        None => { eprintln!("Building 6-mer table..."); KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 6) }
    });

    // Warm the page cache once (matters for mmap); CPU caches warm over the repeated runs.
    eprintln!("Warming up (touching all pages)...");
    rayon::scope(|s| {
        s.spawn(|_| searcher.sa.touch_all_pages());
        s.spawn(|_| searcher.proteins.touch_all_pages());
        s.spawn(|_| searcher.suffix_index_to_protein.touch_all_pages());
    });

    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new().create(true).write(true).truncate(true).open(&output_path)?;
    let commit = env!("GIT_COMMIT_HASH").to_string();

    for pep_path in &args.matrix_files {
        let peptides: Vec<String> = BufReader::new(File::open(pep_path)?)
            .lines().take(args.amount_of_peptides).collect::<Result<_, _>>()?;
        if peptides.len() < args.amount_of_peptides {
            return Err(format!("{} has only {} peptides, need {}", pep_path.display(), peptides.len(), args.amount_of_peptides).into());
        }
        let (p_min, p_max) = peptides.iter().fold((usize::MAX, 0usize), |(lo, hi), p| (lo.min(p.len()), hi.max(p.len())));
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?").to_string();

        for kmer_k in [0usize, 5, 6] {
            searcher.kmer_table = match kmer_k { 5 => table5.take(), 6 => table6.take(), _ => None };
            let theoretical_max = theoretical_memory(&searcher, mapping_type, use_mmap);

            for equate_il in [true, false] {
                for tryptic in [true, false] {
                    for &batch in &args.matrix_batches {
                        // Run every rep, then summarise: one record per config with a spread,
                        // and the median-qps rep kept as the representative detailed `result`.
                        let mut results: Vec<BenchmarkResult> = (0..args.runs)
                            .map(|_| run_benchmark(&searcher, &peptides, args.max_matches, equate_il, tryptic, batch, theoretical_max, baseline_memory))
                            .collect();

                        let mut qps: Vec<f64> = results.iter().map(|r| r.throughput_qps).collect();
                        qps.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        let stats = RunStats {
                            runs: args.runs,
                            qps_min: qps[0],
                            qps_p10: percentile(&qps, 0.10),
                            qps_p50: percentile(&qps, 0.50),
                            qps_p90: percentile(&qps, 0.90),
                            qps_max: *qps.last().unwrap(),
                        };
                        let band = if stats.qps_p50 > 0.0 {
                            (stats.qps_p90 - stats.qps_p10) / 2.0 / stats.qps_p50 * 100.0
                        } else { 0.0 };
                        eprintln!("  {} il={} tr={} batch={} kmer={}  ->  {:.0} qps  (±{:.1}%, p10 {:.0} .. p90 {:.0})",
                            source, equate_il, tryptic, batch, kmer_k, stats.qps_p50, band, stats.qps_p10, stats.qps_p90);

                        // Representative rep = the one nearest the median throughput.
                        results.sort_by(|a, b| a.throughput_qps.partial_cmp(&b.throughput_qps).unwrap());
                        let representative = results.remove(results.len() / 2);

                        let record = BenchmarkRecord {
                            version: SCHEMA_VERSION,
                            label: args.label.clone(),
                            commit: commit.clone(),
                            config: BenchmarkConfig {
                                sa_type: sa_type.to_string(),
                                mapping_type: mapping_type.to_string(),
                                use_mmap,
                                sample_rate,
                                bits_per_value,
                                equate_il,
                                tryptic,
                                max_matches: args.max_matches,
                                batch_size: batch,
                                kmer_k,
                                amount_of_peptides: peptides.len(),
                                peptide_length_min: p_min,
                                peptide_length_max: p_max,
                                peptide_source: source.clone(),
                            },
                            result: representative,
                            stats: Some(stats),
                        };
                        writeln!(output_file, "{}", serde_json::to_string(&record)?)?;
                    }
                }
            }
            match kmer_k { 5 => table5 = searcher.kmer_table.take(), 6 => table6 = searcher.kmer_table.take(), _ => searcher.kmer_table = None };
        }
    }
    eprintln!("Matrix complete → {}", output_path.display());
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    // Validate length range only for random mode
    if args.peptide_file.is_none() && args.peptide_length_min > args.peptide_length_max {
        return Err("--peptide-length-min must be <= --peptide-length-max".into());
    }

    let sa_path = args.index_dir.join("sa.bin");
    let proteins_path = args.index_dir.join("proteins.bin");
    let mapping_path = args.index_dir.join("mapping.bin");

    // Detect mapping type before loading (peek at the type byte)
    let mapping_type_str = match first_byte_of(&mapping_path)? {
        0 => "dense",
        1 => "sparse",
        2 => "bitvec",
        _ => "unknown",
    };

    // Warm up the Rayon thread pool so its thread stacks are included in the baseline.
    // Without this, thread stacks (~2 MB × N threads) would be attributed to index memory
    // because Rayon initialises lazily on the first par_iter() call inside run_benchmark.
    let _ = (0..rayon::current_num_threads()).into_par_iter().sum::<usize>();

    // Baseline RSS after Rayon is warm but before any index data is loaded
    let baseline_memory = measure_process_memory();

    // Load index
    eprintln!("Loading suffix array from {}...", sa_path.display());
    let suffix_array = load_suffix_array_file(sa_path.to_str().unwrap())?;
    eprintln!("  {} items, {} bits/value, sample rate {}",
        suffix_array.len(), suffix_array.bits_per_value(), suffix_array.sample_rate());

    eprintln!("Loading proteins from {}...", proteins_path.display());
    let proteins = load_proteins_file(proteins_path.to_str().unwrap())?;

    eprintln!("Loading mapping from {} (type: {})...", mapping_path.display(), mapping_type_str);
    let mapping = load_mapping_file(mapping_path.to_str().unwrap())?;

    let sa_type = if suffix_array.bits_per_value() == 64 { "original" } else { "compressed" };
    let sample_rate = suffix_array.sample_rate();
    let bits_per_value = suffix_array.bits_per_value();

    let mut searcher = Searcher::new(suffix_array, proteins, mapping);

    if args.matrix {
        return run_matrix(searcher, &args, mapping_type_str, sa_type, sample_rate, bits_per_value, baseline_memory);
    }

    if let Some(ref path) = args.kmer_table_file {
        eprintln!("Loading k-mer table from {}...", path.display());
        let table = load_kmer_table_file(path.to_str().unwrap())?;
        eprintln!("  k={}", table.k);
        searcher = searcher.with_kmer_table(table);
    } else if let Some(k) = args.build_kmer_table {
        eprintln!("Building in-memory k-mer table (k={})...", k);
        searcher.build_kmer_table(k);
        eprintln!("  done.");
    }

    let theoretical_max = theoretical_memory(&searcher, mapping_type_str, cfg!(feature = "mmap"));
    eprintln!("Theoretical max memory: {} bytes ({:.1} MB)", theoretical_max, theoretical_max as f64 / 1_048_576.0);

    // Load peptides: either all at once from a file (for sequential chunking across runs)
    // or generate fresh random peptides for each run.
    let all_peptides: Option<Vec<String>> = if let Some(ref path) = args.peptide_file {
        let peptides: Vec<String> = BufReader::new(File::open(path)?)
            .lines()
            .collect::<Result<_, _>>()?;
        if peptides.is_empty() {
            return Err(format!("peptide file '{}' is empty", path.display()).into());
        }
        let required = args.amount_of_peptides * args.runs as usize;
        if peptides.len() < required {
            return Err(format!(
                "peptide file '{}' has {} lines, but {} runs × {} peptides/run = {} lines are required",
                path.display(),
                peptides.len(),
                args.runs,
                args.amount_of_peptides,
                required,
            ).into());
        }
        eprintln!("Loaded {} peptides from {}", peptides.len(), path.display());
        Some(peptides)
    } else {
        None
    };

    let peptide_source = match &args.peptide_file {
        Some(p) => p.display().to_string(),
        None => "random".to_string(),
    };

    // Optional warmup pass (not timed)
    const WARMUP_BATCH_SIZE: usize = 100_000;
    match &args.warmup {
        None => {}
        Some(WarmupMode::All) => {
            eprintln!("Warming up: touching all mmap pages...");
            rayon::scope(|s| {
                s.spawn(|_| searcher.sa.touch_all_pages());
                s.spawn(|_| searcher.proteins.touch_all_pages());
                s.spawn(|_| searcher.suffix_index_to_protein.touch_all_pages());
            });
            eprintln!("Warmup complete.");
        }
        Some(WarmupMode::Count(warmup_count)) => {
            let warmup_path = args.index_dir.join("warmup.txt");
            eprintln!("Warming up with {} peptides from {} (batch size {})...", warmup_count, warmup_path.display(), WARMUP_BATCH_SIZE);
            let mut lines = BufReader::new(File::open(&warmup_path)?).lines();
            let mut remaining = *warmup_count;
            while remaining > 0 {
                let batch_size = remaining.min(WARMUP_BATCH_SIZE);
                let batch: Vec<String> = lines.by_ref().take(batch_size).collect::<Result<_, _>>()?;
                if batch.is_empty() {
                    break;
                }
                remaining -= batch.len();
                batch.par_iter().for_each(|peptide| {
                    let result = searcher.search_matching_suffixes(
                        peptide.trim_end().to_uppercase().as_bytes(),
                        args.max_matches,
                        args.equate_il,
                        args.tryptic,
                    );
                    match result {
                        SearchAllSuffixesResult::SearchResult(ref suf)
                        | SearchAllSuffixesResult::MaxMatches(ref suf) => {
                            let _ = searcher.retrieve_proteins(suf);
                        }
                        SearchAllSuffixesResult::NoMatches => {}
                    }
                });
            }
            eprintln!("Warmup complete.");
        }
        Some(WarmupMode::AllThenCount(warmup_count)) => {
            eprintln!("Warming up: touching all mmap pages...");
            rayon::scope(|s| {
                s.spawn(|_| searcher.sa.touch_all_pages());
                s.spawn(|_| searcher.proteins.touch_all_pages());
                s.spawn(|_| searcher.suffix_index_to_protein.touch_all_pages());
            });
            eprintln!("Page warmup complete. Running pipeline warmup with {} peptides (batch size {})...", warmup_count, WARMUP_BATCH_SIZE);
            let warmup_path = args.index_dir.join("warmup.txt");
            let mut lines = BufReader::new(File::open(&warmup_path)?).lines();
            let mut remaining = *warmup_count;
            while remaining > 0 {
                let batch_size = remaining.min(WARMUP_BATCH_SIZE);
                let batch: Vec<String> = lines.by_ref().take(batch_size).collect::<Result<_, _>>()?;
                if batch.is_empty() {
                    break;
                }
                remaining -= batch.len();
                batch.par_iter().for_each(|peptide| {
                    let result = searcher.search_matching_suffixes(
                        peptide.trim_end().to_uppercase().as_bytes(),
                        args.max_matches,
                        args.equate_il,
                        args.tryptic,
                    );
                    match result {
                        SearchAllSuffixesResult::SearchResult(ref suf)
                        | SearchAllSuffixesResult::MaxMatches(ref suf) => {
                            let _ = searcher.retrieve_proteins(suf);
                        }
                        SearchAllSuffixesResult::NoMatches => {}
                    }
                });
            }
            eprintln!("Warmup complete.");
        }
    }

    // Prepare output file
    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&output_path)?;

    eprintln!();
    eprintln!("Starting {} benchmark run(s) — results → {}", args.runs, output_path.display());
    eprintln!();

    // SA_MLP_BATCH=B (>1) selects the batched searcher for the whole run.
    let mlp_batch: usize = std::env::var("SA_MLP_BATCH")
        .ok().and_then(|v| v.parse().ok()).filter(|&b| b > 1).unwrap_or(1);

    // RNG for random mode. Seeded → reproducible stream across runs and invocations; otherwise
    // seeded from OS entropy. Unused in file mode (peptides come from disk).
    let mut rng = match args.seed {
        Some(s) => { eprintln!("Random peptide seed: {}", s); StdRng::seed_from_u64(s) }
        None => StdRng::from_entropy(),
    };

    for run in 1..=args.runs {
        // In file mode: consume the next chunk sequentially.
        // In random mode: generate a fresh set of random peptides.
        let run_peptides: Vec<String> = match &all_peptides {
            Some(all) => {
                let start = (run as usize - 1) * args.amount_of_peptides;
                all[start..start + args.amount_of_peptides].to_vec()
            }
            None => generate_peptides(
                &mut rng,
                args.amount_of_peptides,
                args.peptide_length_min,
                args.peptide_length_max,
            ),
        };

        let (p_min, p_max) = run_peptides.iter().fold(
            (usize::MAX, 0usize),
            |(lo, hi), p| (lo.min(p.len()), hi.max(p.len())),
        );

        let config = BenchmarkConfig {
            sa_type: sa_type.to_string(),
            mapping_type: mapping_type_str.to_string(),
            use_mmap: cfg!(feature = "mmap"),
            sample_rate,
            bits_per_value,
            equate_il: args.equate_il,
            tryptic: args.tryptic,
            max_matches: args.max_matches,
            batch_size: mlp_batch,
            kmer_k: searcher.kmer_table.as_ref().map_or(0, |t| t.k),
            amount_of_peptides: run_peptides.len(),
            peptide_length_min: p_min,
            peptide_length_max: p_max,
            peptide_source: peptide_source.clone(),
        };

        let result = run_benchmark(&searcher, &run_peptides, args.max_matches, args.equate_il, args.tryptic, mlp_batch, theoretical_max, baseline_memory);

        let record = BenchmarkRecord {
            version: SCHEMA_VERSION,
            label: args.label.clone(),
            commit: env!("GIT_COMMIT_HASH").to_string(),
            config,
            result,
            stats: None,
        };

        let line = serde_json::to_string(&record)?;
        writeln!(output_file, "{}", line)?;

        eprintln!(
            "Run {:>3}/{}: {:.1} qps  |  total {:.1} ms  (search {:.1} ms [bounds {:.1} ms + iter {:.1} ms] + retrieval {:.1} ms)  \
             |  hits: {} queries, {} suffixes, {} proteins{}",
            run,
            args.runs,
            record.result.throughput_qps,
            record.result.total_duration_ns as f64 / 1e6,
            record.result.search_duration_ns as f64 / 1e6,
            record.result.search_bounds_ns as f64 / 1e6,
            record.result.match_iter_ns as f64 / 1e6,
            record.result.retrieval_duration_ns as f64 / 1e6,
            record.result.query_hit_count,
            record.result.suffix_hit_count,
            record.result.protein_hit_count,
            if record.result.cutoff_reached { "  [cutoff reached]" } else { "" },
        );
    }

    eprintln!();
    eprintln!("Done. {} lines written to {}", args.runs, output_path.display());

    Ok(())
}
