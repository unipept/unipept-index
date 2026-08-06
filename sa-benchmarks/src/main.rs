use std::collections::HashMap;
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
use sa_index::{sa_searcher::{SearchAllSuffixesResult, Searcher}, KmerTable, SearchTuning, SuffixArray, SuffixArrayBackend};
use sa_server::{load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Serialize;
use sysinfo::{Pid, System};
use sa_index::ProteinsBackend as _;
use sa_index::suffix_to_protein_index::SuffixToProteinMappingBackend as _;
use text_compression::ProteinTextBackend as _;

/// Schema version — increment when the output JSON format changes.
/// v2: matrix records aggregate `runs` reps into one line and carry a `stats` spread.
/// v3: every `SearchTuning` field is recorded in `config` (was previously implicit/default),
///     plus `phase` / `ofat_baseline` / `ofat_knob` so grid, OFAT, and confirm records are
///     groupable in one jsonl file; `result` gains `candidates_examined` / `candidates_accepted`.
const SCHEMA_VERSION: u32 = 3;

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

    // -- SearchTuning knobs (applied to `searcher.tuning`). In single-run mode these come
    // straight from the CLI; in matrix mode they're the defaults for the "grid"/"ofat" phases
    // and the explicit combination under test for the "confirm" phase (Task 3). Defaults match
    // `SearchTuning::default()` so omitting them is a no-op.
    /// Candidates per two-pass validation batch in `iterate_sa_range` (only reachable on the
    /// tryptic / non-fast-path route — see sa_index::SearchTuning::validate_batch). Clamped to
    /// 1..=256 internally.
    #[arg(long, default_value_t = SearchTuning::default().validate_batch)]
    validate_batch: usize,

    /// Minimum SA range size before `iterate_sa_range` switches from a straight loop to
    /// two-pass validation.
    #[arg(long, default_value_t = SearchTuning::default().validate_prefetch_threshold)]
    validate_prefetch_threshold: usize,

    /// Prefetch look-ahead distance (in suffixes) inside protein retrieval.
    #[arg(long, default_value_t = SearchTuning::default().retrieval_prefetch_distance)]
    retrieval_prefetch_distance: usize,

    /// Queries whose retrieval is interleaved in one cross-query group.
    #[arg(long, default_value_t = SearchTuning::default().retrieval_batch)]
    retrieval_batch: usize,

    /// Whether the scalar search path issues a redundant `prefetch_kmer_range` before each
    /// skip's `search_bounds`. Only affects the scalar (mlp_batch=1) path.
    #[arg(long, action = clap::ArgAction::Set, default_value_t = SearchTuning::default().scalar_kmer_prefetch)]
    scalar_kmer_prefetch: bool,

    /// Run the full parameter matrix in one process: loads the index once, then iterates
    /// the selected `--matrix-phases` for each `--matrix-files` entry. Writes one record per
    /// (config × run) to <output>/<label>.jsonl.
    #[arg(long)]
    matrix: bool,

    /// Matrix mode: comma-separated peptide files; each becomes one "file" dimension.
    /// File stems must be "small" / "medium" / "large" for the tryptic-collapse rule in
    /// `grid_for_file` to apply; unrecognised stems get the full (uncollapsed) grid.
    #[arg(long, value_delimiter = ',')]
    matrix_files: Vec<PathBuf>,

    /// Matrix mode: pre-built 5-mer table file (falls back to building one if omitted).
    #[arg(long)]
    kmer5_file: Option<PathBuf>,

    /// Matrix mode: pre-built 6-mer table file (falls back to building one if omitted).
    /// Only loaded/built at all when --matrix-kmer6 is set.
    #[arg(long)]
    kmer6_file: Option<PathBuf>,

    /// Matrix mode: MLP batch sizes to sweep, comma-separated (1 = scalar). e.g. 1,8,16,32.
    #[arg(long, value_delimiter = ',', default_values_t = vec![1usize, 16])]
    matrix_batches: Vec<usize>,

    /// Matrix mode: include the 6-mer table in the "grid" and "confirm" phases. Off by
    /// default — the full-DB sweep showed 6-mer vs 5-mer inside the noise floor (p90 3.9%)
    /// on medium/small and only +4.1% on large, for 3.06 GB vs 127 MB resident. Also gates
    /// whether the (expensive) 6-mer table is built/loaded at all.
    #[arg(long)]
    matrix_kmer6: bool,

    /// Matrix mode: comma-separated phases to run in this one process/index-load.
    ///   grid    - Task 1 trimmed default grid (kmer x equate_il x tryptic x mlp_batch, with
    ///             small/medium tryptic collapsed to one representative cell).
    ///   ofat    - Task 2 one-factor-at-a-time sweep of the 5 SearchTuning knobs around two
    ///             fixed baselines (B1 = fast path, B2 = iterate_sa_range path).
    ///   confirm - Task 3: the --validate-batch/--validate-prefetch-threshold/
    ///             --retrieval-prefetch-distance/--retrieval-batch/--scalar-kmer-prefetch CLI
    ///             values applied as one combined tuning across the same grid as "grid", so a
    ///             best-of-each-knob combination can be checked against the OFAT baselines
    ///             (knobs may not be separable — validate_batch and retrieval_prefetch_distance
    ///             both consume line-fill buffers).
    #[arg(long, value_delimiter = ',', default_values_t = vec!["grid".to_string()])]
    matrix_phases: Vec<String>,

    /// Print the planned config list for the selected --matrix-phases and exit, without
    /// loading the index. Use this to eyeball a sweep before committing a multi-hour run.
    #[arg(long)]
    dry_run: bool,
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
    // -- SearchTuning, recorded in full so runs are groupable/reproducible from the jsonl alone
    // (Task 4) without cross-referencing the CLI invocation that produced them.
    validate_batch: usize,
    validate_prefetch_threshold: usize,
    retrieval_prefetch_distance: usize,
    retrieval_batch: usize,
    scalar_kmer_prefetch: bool,
    /// Which sweep produced this record: "single" (non-matrix CLI run), "grid" (Task 1
    /// trimmed grid), "ofat" (Task 2 one-factor-at-a-time), or "confirm" (Task 3 combo check).
    phase: String,
    /// OFAT only: which baseline (B1/B2) this cell was swept around.
    #[serde(skip_serializing_if = "Option::is_none")]
    ofat_baseline: Option<String>,
    /// OFAT only: which SearchTuning field this cell varies.
    #[serde(skip_serializing_if = "Option::is_none")]
    ofat_knob: Option<String>,
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
    /// Candidate suffixes `iterate_sa_range` examined during the search phase (0 without the
    /// `metrics` feature). Settles whether tryptic's ~12.5x slowdown is a low acceptance rate
    /// (work already minimal) or unbounded exhaustive scanning (needs a scan cap) — see
    /// `candidates_accepted` and `Searcher::candidates_accepted`'s doc comment.
    candidates_examined: u64,
    /// Candidate suffixes `iterate_sa_range` accepted as real matches (0 without `metrics`).
    candidates_accepted: u64,
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

    // Reset per-run timing/candidate accumulators before the search phase.
    searcher.drain_timing_ns();
    searcher.drain_candidate_counts();

    // Phase 1: suffix array search (parallel), via the same orchestrator production uses.
    // mlp_batch > 1 interleaves that many peptides per rayon task for memory-level parallelism;
    // 1 = scalar one-peptide-at-a-time.
    let refs: Vec<&[u8]> = peptides.iter().map(|p| p.as_bytes()).collect();
    let search_start = Instant::now();
    let suffix_results = searcher.search_all_matching_suffixes(&refs, max_matches, equate_il, tryptic, mlp_batch);
    let search_duration_ns = search_start.elapsed().as_nanos() as u64;

    // Read internal timing/candidate breakdown accumulated during the search phase above.
    let (search_bounds_ns, match_iter_ns) = searcher.drain_timing_ns();
    let (candidates_examined, candidates_accepted) = searcher.drain_candidate_counts();

    // Phase 2: protein retrieval — the exact function production's `search_all_peptides` calls
    // (`retrieve_all_proteins`), which splits into `tuning.retrieval_batch`-sized cross-query
    // groups run in parallel across rayon, interleaving the random reads *within* each group
    // (`retrieve_proteins_batched`) so the prefetch look-ahead keeps firing even when
    // individual result lists are far shorter than `tuning.retrieval_prefetch_distance`.
    // Calling the same function production calls — instead of reimplementing its chunking
    // here — is what makes this benchmark measure what ships; the old per-query
    // `retrieve_proteins` loop it replaces was invisible to the batched-retrieval change.
    //
    // NoMatches queries are dropped before retrieval, exactly as `search_all_peptides` does —
    // there is nothing to look up for them.
    let matched_suffixes: Vec<&[i64]> = suffix_results
        .iter()
        .filter_map(|r| match r {
            SearchAllSuffixesResult::MaxMatches(suf) | SearchAllSuffixesResult::SearchResult(suf) => {
                Some(suf.as_slice())
            }
            SearchAllSuffixesResult::NoMatches => None,
        })
        .collect();

    let retrieval_start = Instant::now();
    let retrieved = searcher.retrieve_all_proteins(&matched_suffixes);
    let retrieval_duration_ns = retrieval_start.elapsed().as_nanos() as u64;

    // Aggregate stats
    let query_hit_count = matched_suffixes.len();
    let suffix_hit_count: usize = matched_suffixes.iter().map(|s| s.len()).sum();
    let protein_hit_count: usize = retrieved.iter().map(|p| p.len()).sum();
    // Cutoff status comes from the search result, not retrieval — an empty-but-cutoff result
    // (max_matches == 0) is a real cutoff even though retrieval finds nothing to look up.
    let cutoff_reached = suffix_results.iter().any(|r| matches!(r, SearchAllSuffixesResult::MaxMatches(_)));

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
        candidates_examined,
        candidates_accepted,
    }
}

// ---------------------------------------------------------------------------
// Matrix mode: grid / OFAT / confirm sweep generation
// ---------------------------------------------------------------------------

/// One cell of the Task 1 trimmed default grid, or the Task 3 confirm sweep (which reuses the
/// same cell shape with a different tuning).
#[derive(Clone, Copy, Debug)]
struct GridCell {
    kmer_k: usize,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize,
}

/// Builds the trimmed default grid for one peptide-file bucket (Task 1).
///
/// Full sweep is kmer × equate_il × mlp_batch, doubled for tryptic=true/false — except tryptic
/// on the small/medium buckets, which collapses to one representative cell: the last full run
/// showed all 30 small/tryptic cells landing at 653-684 qps (a flat line, not a grid), so
/// sweeping kmer/batch/equate_il there just re-measures a constant at ~1/6 of the whole
/// matrix's wall time. `large` keeps the full sweep since tryptic there is retrieval/search
/// volume bound, not constant.
fn grid_for_file(file_bucket: &str, kmers: &[usize], batches: &[usize]) -> Vec<GridCell> {
    let sweep = |tryptic: bool| -> Vec<GridCell> {
        let mut v = Vec::new();
        for &kmer_k in kmers {
            for equate_il in [true, false] {
                for &mlp_batch in batches {
                    v.push(GridCell { kmer_k, equate_il, tryptic, mlp_batch });
                }
            }
        }
        v
    };

    let mut cells = sweep(false);
    if matches!(file_bucket, "small" | "medium") {
        // Representative cell only, at production defaults (5-mer table, MLP batch 16,
        // equate_il on) — enough to catch a gross regression without re-measuring a constant.
        cells.push(GridCell { kmer_k: 5, equate_il: true, tryptic: true, mlp_batch: 16 });
    } else {
        cells.extend(sweep(true));
    }
    cells
}

/// One cell of the Task 2 OFAT knob sweep: one `SearchTuning` field varied around a fixed
/// baseline, everything else held at that baseline's defaults.
#[derive(Clone, Debug)]
struct OfatCell {
    /// "B1" (fast path, bypasses `iterate_sa_range` entirely) or "B2" (forces it).
    baseline: &'static str,
    knob: &'static str,
    tuning: SearchTuning,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize,
    kmer_k: usize,
}

/// k-mer table used throughout the OFAT sweep — fixed at the production default so knob
/// effects aren't confounded with kmer-table effects.
const OFAT_KMER: usize = 5;
/// MLP batch used throughout the OFAT sweep except the scalar-only `scalar_kmer_prefetch` cell.
const OFAT_MLP_BATCH: usize = 16;

/// Builds the Task 2 OFAT knob sweep: one factor at a time around two fixed baselines, instead
/// of the 46,080-config full cross-product.
///
/// B1 = equate_il=true, tryptic=false — the bulk fast path. `validate_batch` and
/// `validate_prefetch_threshold` only take effect inside `iterate_sa_range`, and B1 never calls
/// it, so sweeping those two knobs at B1 would measure pure noise — their baseline must be B2.
/// B2 = equate_il=true, tryptic=true — forces `iterate_sa_range` on every query.
///
/// `retrieval_prefetch_distance` and `retrieval_batch` affect retrieval regardless of which
/// search path ran, so both are swept at both baselines. `scalar_kmer_prefetch` only affects
/// the scalar search path, so it's swept at B1 with `mlp_batch` forced to 1 (the batched path
/// at mlp_batch=16 never reaches the code this knob guards).
fn build_ofat_cells() -> Vec<OfatCell> {
    let base = SearchTuning::default();
    let mut cells = Vec::new();

    let mut push = |baseline: &'static str, knob: &'static str, tuning: SearchTuning, tryptic: bool, mlp_batch: usize| {
        cells.push(OfatCell { baseline, knob, tuning, equate_il: true, tryptic, mlp_batch, kmer_k: OFAT_KMER });
    };

    // validate_batch / validate_prefetch_threshold: iterate_sa_range only, so baseline = B2.
    for v in [16, 32, 64, 128] {
        push("B2", "validate_batch", SearchTuning { validate_batch: v, ..base }, true, OFAT_MLP_BATCH);
    }
    for v in [8, 16, 32, 64] {
        push("B2", "validate_prefetch_threshold", SearchTuning { validate_prefetch_threshold: v, ..base }, true, OFAT_MLP_BATCH);
    }

    // retrieval_prefetch_distance / retrieval_batch: retrieval always runs, so both baselines.
    for &(label, tryptic) in &[("B1", false), ("B2", true)] {
        for v in [8, 16, 32, 64] {
            push(label, "retrieval_prefetch_distance", SearchTuning { retrieval_prefetch_distance: v, ..base }, tryptic, OFAT_MLP_BATCH);
        }
        for v in [1, 4, 16, 64] {
            push(label, "retrieval_batch", SearchTuning { retrieval_batch: v, ..base }, tryptic, OFAT_MLP_BATCH);
        }
    }

    // scalar_kmer_prefetch: scalar search path only, at B1, forced to mlp_batch=1.
    for v in [true, false] {
        push("B1", "scalar_kmer_prefetch", SearchTuning { scalar_kmer_prefetch: v, ..base }, false, 1);
    }

    cells
}

/// Validates and returns `--matrix-phases`.
fn parse_phases(raw: &[String]) -> Result<Vec<String>, Box<dyn Error>> {
    const VALID: [&str; 3] = ["grid", "ofat", "confirm"];
    if raw.is_empty() {
        return Err("--matrix-phases must name at least one phase".into());
    }
    for p in raw {
        if !VALID.contains(&p.as_str()) {
            return Err(format!("unknown --matrix-phases entry '{}': expected one of {:?}", p, VALID).into());
        }
    }
    Ok(raw.to_vec())
}

/// Swaps the k-mer table of size `k` into `searcher.kmer_table`, returning whatever was
/// previously attached to its owning slot (`table5`/`table6`) first. All swaps are `Option`
/// moves (pointer-sized), so this is cheap to call once per cell regardless of sweep order.
fn ensure_kmer_table(searcher: &mut Searcher<SuffixArray>, table5: &mut Option<KmerTable>, table6: &mut Option<KmerTable>, k: usize) {
    if let Some(t) = searcher.kmer_table.take() {
        match t.k {
            5 => *table5 = Some(t),
            6 => *table6 = Some(t),
            _ => {}
        }
    }
    searcher.kmer_table = match k {
        5 => table5.take(),
        6 => table6.take(),
        _ => None,
    };
}

/// Runs one config for `args.runs` reps, prints the summary line, and appends one aggregated
/// record to `output_file`. Shared by the grid, ofat, and confirm phases of `run_matrix`.
#[allow(clippy::too_many_arguments)]
fn run_cell(
    searcher: &Searcher<SuffixArray>,
    peptides: &[String],
    args: &Args,
    mapping_type: &str,
    sa_type: &str,
    sample_rate: u8,
    bits_per_value: usize,
    use_mmap: bool,
    baseline_memory: u64,
    theoretical_max: u64,
    source: &str,
    p_min: usize,
    p_max: usize,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize,
    kmer_k: usize,
    phase: &str,
    ofat_baseline: Option<&str>,
    ofat_knob: Option<&str>,
    commit: &str,
    output_file: &mut File,
) -> Result<(), Box<dyn Error>> {
    // Run every rep, then summarise: one record per config with a spread, and the median-qps
    // rep kept as the representative detailed `result`.
    let mut results: Vec<BenchmarkResult> = (0..args.runs)
        .map(|_| run_benchmark(searcher, peptides, args.max_matches, equate_il, tryptic, mlp_batch, theoretical_max, baseline_memory))
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

    // Representative rep = the one nearest the median throughput; also carries the detailed
    // per-phase metrics (search_bounds_ns / match_iter_ns / candidate counts) printed below.
    results.sort_by(|a, b| a.throughput_qps.partial_cmp(&b.throughput_qps).unwrap());
    let representative = results.remove(results.len() / 2);

    // Acceptance rate settles whether tryptic's slowdown is "already minimal work, just a low
    // hit rate" or "unbounded exhaustive scanning" (see BenchmarkResult::candidates_examined).
    // Only meaningful with `metrics` — without it both counters are always 0.
    let accept_note = if cfg!(feature = "metrics") {
        let (ex, ac) = (representative.candidates_examined, representative.candidates_accepted);
        let rate = if ex > 0 { ac as f64 / ex as f64 * 100.0 } else { 0.0 };
        format!("  |  candidates: {} examined, {} accepted ({:.1}% accept)", ex, ac, rate)
    } else {
        String::new()
    };

    let tag = match (ofat_baseline, ofat_knob) {
        (Some(b), Some(k)) => format!(" [{} {}]", b, k),
        _ => String::new(),
    };
    eprintln!(
        "  {} {}{} il={} tr={} batch={} kmer={} tuning{{vb={} vpt={} rpd={} rb={} skp={}}}  ->  {:.0} qps  (±{:.1}%, p10 {:.0} .. p90 {:.0}){}",
        source, phase, tag, equate_il, tryptic, mlp_batch, kmer_k,
        searcher.tuning.validate_batch, searcher.tuning.validate_prefetch_threshold,
        searcher.tuning.retrieval_prefetch_distance, searcher.tuning.retrieval_batch, searcher.tuning.scalar_kmer_prefetch,
        stats.qps_p50, band, stats.qps_p10, stats.qps_p90, accept_note,
    );

    let record = BenchmarkRecord {
        version: SCHEMA_VERSION,
        label: args.label.clone(),
        commit: commit.to_string(),
        config: BenchmarkConfig {
            sa_type: sa_type.to_string(),
            mapping_type: mapping_type.to_string(),
            use_mmap,
            sample_rate,
            bits_per_value,
            equate_il,
            tryptic,
            max_matches: args.max_matches,
            batch_size: mlp_batch,
            kmer_k,
            amount_of_peptides: peptides.len(),
            peptide_length_min: p_min,
            peptide_length_max: p_max,
            peptide_source: source.to_string(),
            validate_batch: searcher.tuning.validate_batch,
            validate_prefetch_threshold: searcher.tuning.validate_prefetch_threshold,
            retrieval_prefetch_distance: searcher.tuning.retrieval_prefetch_distance,
            retrieval_batch: searcher.tuning.retrieval_batch,
            scalar_kmer_prefetch: searcher.tuning.scalar_kmer_prefetch,
            phase: phase.to_string(),
            ofat_baseline: ofat_baseline.map(|s| s.to_string()),
            ofat_knob: ofat_knob.map(|s| s.to_string()),
        },
        result: representative,
        stats: Some(stats),
    };
    writeln!(output_file, "{}", serde_json::to_string(&record)?)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Dry run: print the planned config list without loading the index
// ---------------------------------------------------------------------------

fn print_dry_run(args: &Args) -> Result<(), Box<dyn Error>> {
    if args.matrix_files.is_empty() {
        return Err("--matrix requires at least one --matrix-files entry".into());
    }
    let phases = parse_phases(&args.matrix_phases)?;
    let kmers: Vec<usize> = if args.matrix_kmer6 { vec![0, 5, 6] } else { vec![0, 5] };
    let confirm_tuning = SearchTuning {
        validate_batch: args.validate_batch,
        validate_prefetch_threshold: args.validate_prefetch_threshold,
        retrieval_prefetch_distance: args.retrieval_prefetch_distance,
        retrieval_batch: args.retrieval_batch,
        scalar_kmer_prefetch: args.scalar_kmer_prefetch,
    };

    println!("DRY RUN — planned matrix config, no index will be loaded");
    println!("backend        : {}", if cfg!(feature = "mmap") { "mmap" } else { "preloaded" });
    println!("phases         : {:?}", phases);
    println!("kmer sizes     : {:?}", kmers);
    println!("mlp batches    : {:?}", args.matrix_batches);
    if phases.iter().any(|p| p == "confirm") {
        println!("confirm tuning : {:?}", confirm_tuning);
    }
    println!("runs/config    : {}", args.runs);
    println!();

    let mut grand_total = 0usize;

    for pep_path in &args.matrix_files {
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?");
        println!("== {} ==", source);
        let mut file_total = 0usize;

        if phases.iter().any(|p| p == "grid") {
            let cells = grid_for_file(source, &kmers, &args.matrix_batches);
            println!("  grid: {} configs", cells.len());
            for c in &cells {
                println!("    kmer={:<2} il={:<5} tr={:<5} batch={:<3}", c.kmer_k, c.equate_il, c.tryptic, c.mlp_batch);
            }
            file_total += cells.len();
        }

        if phases.iter().any(|p| p == "ofat") {
            let cells = build_ofat_cells();
            println!("  ofat: {} configs", cells.len());
            for c in &cells {
                println!(
                    "    baseline={} knob={:<28} il={} tr={:<5} batch={:<2} kmer={} tuning={:?}",
                    c.baseline, c.knob, c.equate_il, c.tryptic, c.mlp_batch, c.kmer_k, c.tuning
                );
            }
            file_total += cells.len();
        }

        if phases.iter().any(|p| p == "confirm") {
            let cells = grid_for_file(source, &kmers, &args.matrix_batches);
            println!("  confirm: {} configs (tuning: {:?})", cells.len(), confirm_tuning);
            for c in &cells {
                println!("    kmer={:<2} il={:<5} tr={:<5} batch={:<3}", c.kmer_k, c.equate_il, c.tryptic, c.mlp_batch);
            }
            file_total += cells.len();
        }

        println!("  -- {} total: {} configs", source, file_total);
        grand_total += file_total;
        println!();
    }

    println!("TOTAL this backend: {} configs x {} runs = {} timed executions", grand_total, args.runs, grand_total * args.runs as usize);
    println!("Run once per backend (preloaded, mmap) for the full sweep.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Matrix mode: load the index once, then run the selected `--matrix-phases` for each peptide
/// file. The k-mer tables are swapped in/out (moved, not cloned) so the big 6-mer table is
/// built/loaded at most once (and only at all if `--matrix-kmer6` is set). Each config runs
/// `runs` reps but writes a single aggregated record (median-qps rep as `result`, plus a
/// `stats` spread) to <output>/<label>.jsonl.
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
    let phases = parse_phases(&args.matrix_phases)?;
    let use_mmap = cfg!(feature = "mmap");
    let kmers: Vec<usize> = if args.matrix_kmer6 { vec![0, 5, 6] } else { vec![0, 5] };
    let confirm_tuning = SearchTuning {
        validate_batch: args.validate_batch,
        validate_prefetch_threshold: args.validate_prefetch_threshold,
        retrieval_prefetch_distance: args.retrieval_prefetch_distance,
        retrieval_batch: args.retrieval_batch,
        scalar_kmer_prefetch: args.scalar_kmer_prefetch,
    };

    // Build/load the 5-mer table (always in scope — the grid's default kmer set includes it)
    // and, only if requested, the 6-mer table: at 3.06 GB vs 127 MB for a sub-noise-floor
    // difference (see grid_for_file / --matrix-kmer6), it isn't worth the build/load cost
    // unless someone explicitly opts in.
    let mut table5: Option<KmerTable> = Some(match &args.kmer5_file {
        Some(p) => { eprintln!("Loading 5-mer table from {}...", p.display()); load_kmer_table_file(p.to_str().unwrap())? }
        None => { eprintln!("Building 5-mer table..."); KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 5) }
    });
    let mut table6: Option<KmerTable> = if args.matrix_kmer6 {
        Some(match &args.kmer6_file {
            Some(p) => { eprintln!("Loading 6-mer table from {}...", p.display()); load_kmer_table_file(p.to_str().unwrap())? }
            None => { eprintln!("Building 6-mer table..."); KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 6) }
        })
    } else {
        None
    };

    // Warm the page cache once (matters for mmap); CPU caches warm over the repeated runs.
    eprintln!("Warming up (touching all pages)...");
    rayon::scope(|s| {
        s.spawn(|_| searcher.sa.touch_all_pages());
        s.spawn(|_| searcher.proteins.touch_all_pages());
        s.spawn(|_| searcher.suffix_index_to_protein.touch_all_pages());
    });

    // Theoretical memory footprint per k-mer table size — index-wide, independent of the
    // peptide file, so this is computed once per k value up front rather than per cell (it
    // walks every protein's metadata length, which is not free at full-DB scale). `kmers` is
    // always {0,5} or {0,5,6}, so it already covers OFAT_KMER (5) — the ofat phase's lookups
    // below never miss.
    debug_assert!(kmers.contains(&OFAT_KMER), "kmers must include OFAT_KMER for the ofat phase's theoretical_by_k lookups");
    let mut theoretical_by_k: HashMap<usize, u64> = HashMap::new();
    for &k in &kmers {
        ensure_kmer_table(&mut searcher, &mut table5, &mut table6, k);
        theoretical_by_k.insert(k, theoretical_memory(&searcher, mapping_type, use_mmap));
    }

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

        if phases.iter().any(|p| p == "grid") {
            eprintln!("-- {} : grid --", source);
            for cell in grid_for_file(&source, &kmers, &args.matrix_batches) {
                ensure_kmer_table(&mut searcher, &mut table5, &mut table6, cell.kmer_k);
                searcher.tuning = SearchTuning::default();
                run_cell(
                    &searcher, &peptides, args, mapping_type, sa_type, sample_rate, bits_per_value,
                    use_mmap, baseline_memory, theoretical_by_k[&cell.kmer_k], &source, p_min, p_max,
                    cell.equate_il, cell.tryptic, cell.mlp_batch, cell.kmer_k, "grid", None, None,
                    &commit, &mut output_file,
                )?;
            }
        }

        if phases.iter().any(|p| p == "ofat") {
            eprintln!("-- {} : ofat --", source);
            for cell in build_ofat_cells() {
                ensure_kmer_table(&mut searcher, &mut table5, &mut table6, cell.kmer_k);
                searcher.tuning = cell.tuning;
                run_cell(
                    &searcher, &peptides, args, mapping_type, sa_type, sample_rate, bits_per_value,
                    use_mmap, baseline_memory, theoretical_by_k[&cell.kmer_k], &source, p_min, p_max,
                    cell.equate_il, cell.tryptic, cell.mlp_batch, cell.kmer_k, "ofat",
                    Some(cell.baseline), Some(cell.knob),
                    &commit, &mut output_file,
                )?;
            }
        }

        if phases.iter().any(|p| p == "confirm") {
            eprintln!("-- {} : confirm (tuning: {:?}) --", source, confirm_tuning);
            for cell in grid_for_file(&source, &kmers, &args.matrix_batches) {
                ensure_kmer_table(&mut searcher, &mut table5, &mut table6, cell.kmer_k);
                searcher.tuning = confirm_tuning;
                run_cell(
                    &searcher, &peptides, args, mapping_type, sa_type, sample_rate, bits_per_value,
                    use_mmap, baseline_memory, theoretical_by_k[&cell.kmer_k], &source, p_min, p_max,
                    cell.equate_il, cell.tryptic, cell.mlp_batch, cell.kmer_k, "confirm", None, None,
                    &commit, &mut output_file,
                )?;
            }
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

    // --dry-run never touches the index (which may not even exist locally) — it only expands
    // and prints the config list the requested matrix phases would run.
    if args.dry_run {
        if !args.matrix {
            return Err("--dry-run only applies to --matrix mode".into());
        }
        return print_dry_run(&args);
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

    // Apply the SearchTuning knobs from the CLI (defaults match SearchTuning::default(), so
    // this is a no-op unless the caller overrides one of --validate-batch/etc).
    searcher.tuning = SearchTuning {
        validate_batch: args.validate_batch,
        validate_prefetch_threshold: args.validate_prefetch_threshold,
        retrieval_prefetch_distance: args.retrieval_prefetch_distance,
        retrieval_batch: args.retrieval_batch,
        scalar_kmer_prefetch: args.scalar_kmer_prefetch,
    };

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
            validate_batch: searcher.tuning.validate_batch,
            validate_prefetch_threshold: searcher.tuning.validate_prefetch_threshold,
            retrieval_prefetch_distance: searcher.tuning.retrieval_prefetch_distance,
            retrieval_batch: searcher.tuning.retrieval_batch,
            scalar_kmer_prefetch: searcher.tuning.scalar_kmer_prefetch,
            phase: "single".to_string(),
            ofat_baseline: None,
            ofat_knob: None,
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

        // Acceptance rate is only meaningful with `metrics` (both counters are always 0
        // without it) — see BenchmarkResult::candidates_examined for what it settles.
        let accept_note = if cfg!(feature = "metrics") {
            let (ex, ac) = (record.result.candidates_examined, record.result.candidates_accepted);
            let rate = if ex > 0 { ac as f64 / ex as f64 * 100.0 } else { 0.0 };
            format!("  |  candidates: {} examined, {} accepted ({:.1}% accept)", ex, ac, rate)
        } else {
            String::new()
        };

        eprintln!(
            "Run {:>3}/{}: {:.1} qps  |  total {:.1} ms  (search {:.1} ms [bounds {:.1} ms + iter {:.1} ms] + retrieval {:.1} ms)  \
             |  hits: {} queries, {} suffixes, {} proteins{}{}",
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
            accept_note,
        );
    }

    eprintln!();
    eprintln!("Done. {} lines written to {}", args.runs, output_path.display());

    Ok(())
}
