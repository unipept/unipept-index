use std::error::Error;
use std::fs::{File, OpenOptions, create_dir_all};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use rand::Rng;
use rayon::prelude::*;
use sa_index::peptide_search::search_all_peptides;
use sa_index::sa_searcher::{SearchAllSuffixesResult, Searcher};
use sa_index::suffix_to_protein_index::SuffixToProteinMapping;
use sa_server::{load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Serialize;
use sysinfo::{Pid, System};

/// Schema version — increment when the output JSON format changes
const SCHEMA_VERSION: u32 = 1;

/// Canonical 20 amino acids used for random peptide generation
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

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

    /// Number of randomly generated peptides per run
    #[arg(long, default_value_t = 10_000)]
    amount_of_peptides: usize,

    /// Minimum length of a randomly generated peptide
    #[arg(long, default_value_t = 5)]
    peptide_length_min: usize,

    /// Maximum length of a randomly generated peptide
    #[arg(long, default_value_t = 50)]
    peptide_length_max: usize,

    /// Number of timed benchmark runs to perform
    #[arg(long, default_value_t = 100)]
    runs: u32,

    /// Warm up the index with the first N peptides from {index-dir}/warmup.txt before timing
    #[arg(long)]
    warmup: Option<usize>,

    /// Use memory-mapped I/O when loading the index files
    #[arg(long, default_value_t = false)]
    mmap: bool,
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
fn theoretical_memory(searcher: &Searcher, mapping_type: &str) -> u64 {
    let text_len = searcher.proteins.text().len() as u64;
    let protein_count = searcher.proteins.len() as u64;

    // Suffix array: one entry per SA item at bits_per_value bits each
    let sa_bytes = (searcher.sa.len() as u64 * searcher.sa.bits_per_value() as u64).div_ceil(8);

    // ProteinText: 5 bits per character (BitArray), rounded up to whole bytes
    let text_bytes = (text_len * 5).div_ceil(8);

    // Protein metadata: compact 16-byte table entries + raw string bytes
    let string_bytes: u64 = (0..searcher.proteins.len())
        .map(|i| {
            let p = searcher.proteins.get(i);
            p.uniprot_id.len() as u64 + p.functional_annotations.len() as u64
        })
        .sum();
    let metadata_bytes = protein_count * 16 + string_bytes;

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

    sa_bytes + text_bytes + metadata_bytes + mapping_bytes
    // ← add new Searcher structures here
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

fn run_benchmark(searcher: &Searcher, args: &Args, theoretical_max_memory: u64) -> BenchmarkResult {
    let peptides = generate_peptides(
        args.amount_of_peptides,
        args.peptide_length_min,
        args.peptide_length_max,
    );

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
        args.amount_of_peptides as f64 / (total_duration_ns as f64 / 1e9)
    } else {
        0.0
    };

    BenchmarkResult {
        search_duration_ns,
        retrieval_duration_ns,
        total_duration_ns,
        throughput_qps,
        amount_of_queries: args.amount_of_peptides,
        query_hit_count,
        suffix_hit_count,
        protein_hit_count,
        cutoff_reached,
        total_memory: measure_process_memory(),
        theoretical_max_memory,
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    if args.peptide_length_min > args.peptide_length_max {
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

    let searcher = Searcher::new(suffix_array, proteins, mapping);

    let theoretical_max = theoretical_memory(&searcher, mapping_type_str);
    eprintln!("Theoretical max memory: {} bytes ({:.1} MB)", theoretical_max, theoretical_max as f64 / 1_048_576.0);

    let config = BenchmarkConfig {
        sa_type: sa_type.to_string(),
        mapping_type: mapping_type_str.to_string(),
        use_mmap: args.mmap,
        sample_rate,
        bits_per_value,
        equate_il: args.equate_il,
        tryptic: args.tryptic,
        max_matches: args.max_matches,
        amount_of_peptides: args.amount_of_peptides,
        peptide_length_min: args.peptide_length_min,
        peptide_length_max: args.peptide_length_max,
    };

    // Optional warmup pass (not timed)
    if let Some(warmup_count) = args.warmup {
        let warmup_path = args.index_dir.join("warmup.txt");
        let warmup_peptides: Vec<String> = BufReader::new(File::open(&warmup_path)?)
            .lines()
            .take(warmup_count)
            .collect::<Result<_, _>>()?;
        eprintln!("Warming up with {} peptides from {}...", warmup_peptides.len(), warmup_path.display());
        let _ = search_all_peptides(&searcher, &warmup_peptides, args.max_matches, args.equate_il, args.tryptic);
        eprintln!("Warmup complete.");
    }

    // Prepare output file
    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&output_path)?;

    eprintln!();
    eprintln!("Starting {} benchmark run(s) — results → {}", args.runs, output_path.display());
    eprintln!();

    for run in 1..=args.runs {
        let result = run_benchmark(&searcher, &args, theoretical_max);

        let record = BenchmarkRecord {
            version: SCHEMA_VERSION,
            label: args.label.clone(),
            commit: env!("GIT_COMMIT_HASH").to_string(),
            config: config.clone(),
            result,
        };

        let line = serde_json::to_string(&record)?;
        writeln!(output_file, "{}", line)?;

        eprintln!(
            "Run {:>3}/{}: {:.1} qps  |  total {:.1} ms  (search {:.1} ms + retrieval {:.1} ms)  \
             |  hits: {} queries, {} suffixes, {} proteins{}",
            run,
            args.runs,
            record.result.throughput_qps,
            record.result.total_duration_ns as f64 / 1e6,
            record.result.search_duration_ns as f64 / 1e6,
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
