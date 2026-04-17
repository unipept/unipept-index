use std::error::Error;
use std::fs::{File, OpenOptions, create_dir_all};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use rand::Rng;
use rayon::prelude::*;
use sa_index::sa_searcher::{SearchAllSuffixesResult, Searcher};
use sa_index::suffix_to_protein_index::SuffixToProteinMapping;
use sa_server::{load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Serialize;
use sysinfo::{Pid, System};
use text_compression::WriteBinary;

/// Schema version — increment when the output JSON format changes
const SCHEMA_VERSION: u32 = 1;

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
}

impl std::str::FromStr for WarmupMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "all" { return Ok(WarmupMode::All); }
        s.parse::<usize>()
            .map(WarmupMode::Count)
            .map_err(|_| format!("expected a non-negative integer, got '{}'", s))
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
    #[arg(long, default_value_t = true)]
    equate_il: bool,

    /// Only return tryptic matches
    #[arg(long, default_value_t = false)]
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

    /// Warm up the index before timing.
    /// Without a value: touch every page of every mmap-backed region (fully populates page cache).
    /// With N: search the first N peptides from {index-dir}/warmup.txt.
    #[arg(long, num_args = 0..=1, default_missing_value = "all")]
    warmup: Option<WarmupMode>,

    /// Use memory-mapped I/O when loading the index files
    #[arg(long, default_value_t = false)]
    mmap: bool,
    /// Optional path to a pre-built k-mer bounds table file.
    /// When provided, binary search is accelerated by narrowing the initial search window.
    #[arg(long)]
    kmer_table_file: Option<PathBuf>,
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

#[derive(Serialize)]
struct BenchmarkRecord {
    version: u32,
    label: String,
    commit: String,
    config: BenchmarkConfig,
    result: BenchmarkResult,
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

/// Generates `count` random peptides whose length is in `[min_len, max_len]`.
fn generate_peptides(count: usize, min_len: usize, max_len: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
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
fn theoretical_memory(searcher: &Searcher, mapping_type: &str, use_mmap: bool) -> u64 {
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

    // k-mer table size
    let kmer_table_bytes = 25_u64.pow(6) * 16;

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

fn run_benchmark(searcher: &Searcher, args: &Args, peptides: &[String], theoretical_max_memory: u64, baseline_memory: u64) -> BenchmarkResult {
    // Memory snapshot before any timing starts — captures index-resident pages only
    let index_memory = measure_process_memory().saturating_sub(baseline_memory);

    // Reset per-run timing accumulators before the search phase.
    searcher.drain_timing_ns();

    // Phase 1: suffix array search (parallel)
    let search_start = Instant::now();
    let suffix_results: Vec<SearchAllSuffixesResult> = peptides
        .par_iter()
        .map(|p| {
            searcher.search_matching_suffixes(
                p.as_bytes(),
                args.max_matches,
                args.equate_il,
                args.tryptic,
            )
        })
        .collect();
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
    let suffix_array = load_suffix_array_file(sa_path.to_str().unwrap(), args.mmap)?;
    eprintln!("  {} items, {} bits/value, sample rate {}",
        suffix_array.len(), suffix_array.bits_per_value(), suffix_array.sample_rate());

    eprintln!("Loading proteins from {}...", proteins_path.display());
    let proteins = load_proteins_file(proteins_path.to_str().unwrap(), args.mmap)?;

    eprintln!("Loading mapping from {} (type: {})...", mapping_path.display(), mapping_type_str);
    let SuffixToProteinMapping(mapping) = load_mapping_file(mapping_path.to_str().unwrap(), args.mmap)?;

    let sa_type = if suffix_array.bits_per_value() == 64 { "original" } else { "compressed" };
    let sample_rate = suffix_array.sample_rate();
    let bits_per_value = suffix_array.bits_per_value();

    let mut searcher = Searcher::new(suffix_array, proteins, mapping);

    if let Some(ref path) = args.kmer_table_file {
        eprintln!("Loading k-mer table from {}...", path.display());
        let table = load_kmer_table_file(path.to_str().unwrap())?;
        eprintln!("  k={}", table.k);
        searcher = searcher.with_kmer_table(table);
    } else {
        eprintln!("Building k-mer table with k=6...");
        searcher.build_kmer_table(6);
    }

    let theoretical_max = theoretical_memory(&searcher, mapping_type_str, args.mmap);
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
            searcher.sa.touch_all_pages();
            searcher.proteins.touch_all_pages();
            searcher.suffix_index_to_protein.touch_all_pages();
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

    for run in 1..=args.runs {
        // In file mode: consume the next chunk sequentially.
        // In random mode: generate a fresh set of random peptides.
        let run_peptides: Vec<String> = match &all_peptides {
            Some(all) => {
                let start = (run as usize - 1) * args.amount_of_peptides;
                all[start..start + args.amount_of_peptides].to_vec()
            }
            None => generate_peptides(
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
            use_mmap: args.mmap,
            sample_rate,
            bits_per_value,
            equate_il: args.equate_il,
            tryptic: args.tryptic,
            max_matches: args.max_matches,
            amount_of_peptides: run_peptides.len(),
            peptide_length_min: p_min,
            peptide_length_max: p_max,
            peptide_source: peptide_source.clone(),
        };

        let result = run_benchmark(&searcher, &args, &run_peptides, theoretical_max, baseline_memory);

        let record = BenchmarkRecord {
            version: SCHEMA_VERSION,
            label: args.label.clone(),
            commit: env!("GIT_COMMIT_HASH").to_string(),
            config,
            result,
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
