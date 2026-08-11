//! Throughput/memory measurement harness for the suffix-array index.
//!
//! Loads a built index (`sa.bin` / `proteins.bin` / `mapping.bin`, plus an optional k-mer
//! table), pushes peptides through the same search + retrieval pipeline production uses, and
//! writes one JSONL record per measured config to `<output>/<label>.jsonl`.
//!
//! Two modes:
//!   * single run  — one config from the CLI flags, `--runs` reps, one record per rep;
//!   * `--matrix`  — load the index once and sweep the `grid` of configs
//!     (k-mer table × equate_il × tryptic × MLP batch) across several peptide files,
//!     writing one aggregated record (median rep + spread) per config. `--dry-run` prints
//!     the planned config list without touching the index.
//!
//! Dev-only: this crate is a workspace member but is excluded from `default-members`, so a
//! plain `cargo build` skips it. Build and run it explicitly:
//!
//! ```text
//! cargo build --release -p sa-benchmarks                  # preloaded backend
//! cargo build --release -p sa-benchmarks --features mmap  # mmap backend
//! ./target/release/sa-benchmarks --index-dir <idx> --output /tmp/bench --label smoke \
//!     --peptide-file <peptides.txt> --amount-of-peptides 10000 --runs 20 --warmup all
//! ```
//!
//! The `metrics` feature adds the per-candidate counters and the internal phase breakdown, at
//! the cost of perturbing what it measures — keep it off for timing runs.
//! See `matrix_bench.sh` / `mlp_sweep.sh` next to this crate for the driver scripts.

use std::{
    collections::HashMap,
    error::Error,
    fs::{File, OpenOptions, create_dir_all},
    io::{BufRead, BufReader, Read, Write},
    path::{Path, PathBuf},
    time::Instant
};

use clap::Parser;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use sa_index::{
    KmerTable, ProteinsBackend as _, SearchTuning, SuffixArray, SuffixArrayBackend,
    kmer_table::AMINO_ACID_COUNT,
    sa_searcher::{DEFAULT_MLP_BATCH, SearchAllSuffixesResult, Searcher},
    suffix_to_protein_index::SuffixToProteinMappingBackend as _
};
use sa_server::{load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Serialize;
use sysinfo::{Pid, System};
use text_compression::ProteinTextBackend as _;

/// Schema version — increment when the output JSON format changes.
/// v2: matrix records aggregate `runs` reps into one line and carry a `stats` spread.
/// v3: every `SearchTuning` field is recorded in `config` (was previously implicit/default),
///     plus a `phase` tag so records from different sweeps are groupable in one jsonl file;
///     `result` gains `candidates_examined` / `candidates_accepted`.
/// v4: the OFAT and confirm sweeps were retired once they had settled which knobs matter, so
///     `config.ofat_baseline` / `config.ofat_knob` are gone and `config.phase` is now only
///     "single" (non-matrix CLI run) or "grid" (matrix sweep).
/// v5: storage is chosen per structure, so the single `use_mmap` bool is replaced by
///     `sa_storage` / `text_storage` / `proteins_storage` / `mapping_storage`, each
///     "mmap" or "preloaded".
/// v6: records carry a `startup` section timing each structure's load and the warmup pass.
///     Preloading a structure moves cost from steady-state page faults to startup, and that
///     trade was previously invisible.
/// v7: `result` gains `major_faults` / `minor_faults`, counted across the timed region. They are
///     what separates "slow because it is waiting on I/O" from "slow for some other reason" when
///     the index does not fit in RAM.
const SCHEMA_VERSION: u32 = 7;

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
    AllThenCount(usize)
}

impl std::str::FromStr for WarmupMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "all" {
            return Ok(WarmupMode::All);
        }
        if let Some(rest) = s.strip_prefix("all:") {
            return rest
                .parse::<usize>()
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

    /// Cross-query MLP batch size: how many independent peptide searches are interleaved per
    /// rayon task to hide random-access DRAM latency. 1 = scalar (one peptide per task).
    /// Defaults to the production value so an unqualified run measures what ships.
    /// Single-run mode only — matrix mode sweeps `--matrix-batches` instead.
    #[arg(long, default_value_t = DEFAULT_MLP_BATCH)]
    mlp_batch: usize,

    // -- SearchTuning knobs (applied to `searcher.tuning`), in both single-run and matrix mode.
    // Defaults match `SearchTuning::default()` so omitting them is a no-op.
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

    /// Run the full parameter matrix in one process: loads the index once, then sweeps the
    /// grid (see `expand_cells`) for each `--matrix-files` entry. Writes one aggregated
    /// record per config to <output>/<label>.jsonl.
    #[arg(long)]
    matrix: bool,

    /// Matrix mode: comma-separated peptide files; each becomes one "file" dimension.
    /// File stems must be "small" / "medium" / "large" for the tryptic-collapse rule in
    /// `expand_cells` to apply; unrecognised stems get the full (uncollapsed) grid.
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

    /// Matrix mode: include the 6-mer table in the grid. Off by
    /// default — the full-DB sweep showed 6-mer vs 5-mer inside the noise floor (p90 3.9%)
    /// on medium/small and only +4.1% on large, for 3.06 GB vs 127 MB resident. Also gates
    /// whether the (expensive) 6-mer table is built/loaded at all.
    ///
    /// That verdict holds only while the index is resident. Under a memory ceiling the 6-mer is
    /// +18.4% and -27.9% major faults against no table, where a 5-mer is +3.2% / -6.2% — barely
    /// distinguishable from nothing. The table's value there is working-set size (~1 SA page per
    /// query at k=6 vs ~7 at k=5), not probe count, and that only matters once pages can be
    /// evicted. See the `sa-index` crate docs.
    #[arg(long)]
    matrix_kmer6: bool,

    /// Print the planned config list for the matrix sweep and exit, without
    /// loading the index. Use this to eyeball a sweep before committing a multi-hour run.
    #[arg(long)]
    dry_run: bool,

    /// Issue madvise(MADV_WILLNEED) over each SA range before scanning it.
    ///
    /// Off by default and it has regressed before — -16.8% qps with the index resident, from
    /// mmap_lock contention across rayon threads (see MmapBackedSA::advise_willneed_range). It is
    /// exposed to test the opposite regime: under a memory ceiling a CPU prefetch hint cannot help
    /// (it cannot fault), so the syscall may replace a real disk stall rather than nothing.
    #[arg(long)]
    willneed: bool,

    /// Skip the theoretical memory calculation, reporting `theoretical_max_memory: 0`.
    ///
    /// That calculation walks *every* protein's metadata (see `theoretical_memory`), which on an
    /// mmap backend faults the entire metadata section in before anything is timed. Harmless when
    /// the index fits in RAM, but under a cgroup memory cap it both pre-warms what the run is
    /// supposed to be faulting on demand and spends the budget being measured. Off by default, so
    /// ordinary runs still report the figure.
    #[arg(long)]
    no_theoretical_memory: bool
}

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

/// What the process paid before it could answer the first query.
///
/// Recorded per invocation (so every rep in a file repeats it) because it is a property of the
/// build and the index, not of a rep. It exists because preloading a structure does not make its
/// cost go away — it moves it from page faults spread across the steady state to one bulk copy at
/// startup, and comparing configurations on throughput alone hides that entirely.
#[derive(Clone, Copy, Serialize)]
struct StartupTiming {
    /// Reading `sa.bin`.
    load_sa_ms: u64,
    /// Reading `proteins.bin` — both the metadata table and the text, which may be stored
    /// differently from each other.
    load_proteins_ms: u64,
    /// Reading `mapping.bin`.
    load_mapping_ms: u64,
    /// Loading or building the k-mer table, 0 when there is none. Also 0 in matrix mode, where
    /// tables are swapped per cell rather than once at startup.
    kmer_table_ms: u64,
    /// The `--warmup` pass: touching pages and/or pushing peptides through. 0 without `--warmup`.
    warmup_ms: u64,
    /// The three index loads plus the k-mer table. Excludes warmup, which is opt-in.
    load_total_ms: u64
}

#[derive(Clone, Serialize)]
struct BenchmarkConfig {
    sa_type: String,
    mapping_type: String,
    /// Where each structure was stored, "mmap" or "preloaded". Independent per structure, and
    /// fixed at compile time — two binaries from the same commit can differ here, so a record
    /// without these four fields cannot be attributed to a configuration.
    sa_storage: &'static str,
    text_storage: &'static str,
    proteins_storage: &'static str,
    mapping_storage: &'static str,
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
    /// Whether madvise(MADV_WILLNEED) was issued over each SA range before scanning.
    willneed: bool,
    /// Which sweep produced this record: "single" (non-matrix CLI run) or "grid" (the trimmed
    /// default matrix grid).
    phase: String
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
    /// Page faults taken across the timed region (search + retrieval), not for the whole process.
    ///
    /// `major_faults` are the ones that reached disk. On a box where the index fits they stay near
    /// zero in both backends; when it does not fit they are the whole story, and a run that slows
    /// down *without* them is slow for some other reason.
    major_faults: u64,
    minor_faults: u64
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
    qps_max: f64
}

#[derive(Serialize)]
struct BenchmarkRecord {
    version: u32,
    label: String,
    commit: String,
    config: BenchmarkConfig,
    /// Absent in records written before schema v6.
    startup: StartupTiming,
    result: BenchmarkResult,
    /// Per-config throughput spread over all reps (matrix mode only; omitted otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    stats: Option<RunStats>
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
            (0..len).map(|_| AMINO_ACIDS[rng.gen_range(0..AMINO_ACIDS.len())] as char).collect()
        })
        .collect()
}

/// Whether the protein metadata table is stored in the mapping rather than on the heap.
///
/// Only the metadata layout differs between the two, so this is the single axis
/// [`theoretical_memory`] needs — the text is accounted for identically either way.
fn proteins_mapped() -> bool {
    sa_server::PROTEINS_BACKEND == "mmap"
}

/// Computes the theoretical in-memory footprint of the loaded index structures.
///
/// This is derived from the actual data sizes, **not** from disk file sizes, so it remains
/// accurate when new structures are added to the `Searcher`. When you add a new structure,
/// extend this function with its memory calculation.
fn theoretical_memory(searcher: &Searcher<SuffixArray>, mapping_type: &str, proteins_mapped: bool) -> u64 {
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
    let metadata_bytes = if proteins_mapped {
        // MmapBacked: 16-byte fixed table entry per protein + concatenated string blobs
        protein_count * 16 + string_bytes
    } else {
        // InMemory: Vec<Protein> on heap — each Protein struct is 56 bytes
        // (String=24 + u32=4 + padding=4 + Vec<u8>=24) plus the heap-allocated string data
        protein_count * 56 + string_bytes
    };

    // Suffix-to-protein mapping (see suffix_to_protein_index/{dense,sparse,bitvec}.rs)
    let mapping_bytes = match mapping_type {
        "dense" => text_len * 4,             // Vec<u32> with one u32 per text character
        "sparse" => (protein_count + 2) * 8, // Vec<i64> with one i64 per protein boundary
        "bitvec" => {
            // One bit per text character + Rank9 superblock overhead (~16 bytes per 512 bits)
            let bits_bytes = text_len.div_ceil(8);
            let superblock_count = text_len / 512 + 1;
            let rank9_bytes = superblock_count * 16;
            bits_bytes + rank9_bytes
        }
        _ => 0
    };

    // k-mer table size: 16 bytes per entry, AMINO_ACID_COUNT^k entries total
    let kmer_table_bytes = searcher.kmer_table.as_ref().map_or(0, |t| (AMINO_ACID_COUNT as u64).pow(t.k as u32) * 16);

    sa_bytes + text_bytes + metadata_bytes + mapping_bytes + kmer_table_bytes
}

/// Page faults this process has taken so far, as `(minor, major)`.
///
/// Major faults are the ones that went to disk; on an mmap backend whose index does not fit in
/// RAM they are the dominant cost, and their absence in a slow run means the slowness is *not*
/// residency. Counted around the timed region rather than for the whole process so warmup and
/// index loading do not swamp the measurement.
///
/// Read straight from `/proc/self/stat` — `sysinfo` does not expose fault counts portably. Linux
/// only; every other platform reports zeros, which keeps the crate building on macOS where these
/// numbers are not the point.
#[cfg(target_os = "linux")]
fn page_faults() -> (u64, u64) {
    std::fs::read_to_string("/proc/self/stat").map(|s| parse_proc_stat_faults(&s)).unwrap_or((0, 0))
}

/// Extracts `(minflt, majflt)` from the contents of `/proc/self/stat`.
///
/// Separate from the file read so it can be tested off Linux — the field arithmetic is the part
/// that breaks, and it breaks by returning zeros, which is indistinguishable from a run that
/// genuinely took no faults. That failure mode is why this is a pure function with tests rather
/// than four lines inline.
///
/// Field 2 is the executable name in parentheses and may itself contain spaces and parentheses,
/// so parsing starts after the **last** `)`, never at the first space.
#[allow(dead_code)] // Only called on Linux; the tests below exercise it everywhere.
fn parse_proc_stat_faults(stat: &str) -> (u64, u64) {
    let Some((_, rest)) = stat.rsplit_once(')') else {
        return (0, 0);
    };
    // `rest` begins at field 3 (state), so field 10 (minflt) is index 7 and field 12 (majflt) 9.
    let fields: Vec<&str> = rest.split_whitespace().collect();
    let get = |i: usize| fields.get(i).and_then(|v| v.parse().ok()).unwrap_or(0);
    (get(7), get(9))
}

/// Non-Linux stub — see the Linux implementation for why this is not worth emulating.
#[cfg(not(target_os = "linux"))]
fn page_faults() -> (u64, u64) {
    (0, 0)
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
fn run_benchmark(
    searcher: &Searcher<SuffixArray>,
    peptides: &[String],
    max_matches: usize,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize,
    theoretical_max_memory: u64,
    baseline_memory: u64
) -> BenchmarkResult {
    // Memory snapshot before any timing starts — captures index-resident pages only
    let index_memory = measure_process_memory().saturating_sub(baseline_memory);

    // Reset per-run timing/candidate accumulators before the search phase.
    searcher.drain_timing_ns();
    searcher.drain_candidate_counts();

    // Fault counters bracket the timed region only, so index loading and warmup are excluded.
    let (minflt_before, majflt_before) = page_faults();

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

    // Phase 2: protein retrieval — per query via `retrieve_proteins`, which is exactly what
    // production's `search_all_peptides` does. Keeping these two in step is what makes this
    // benchmark measure what ships: an earlier revision called a batched retrieval here while
    // production called the per-query one, which made a whole change invisible to measurement.
    // If production's retrieval shape changes, change it here too.
    //
    // NoMatches queries are dropped before retrieval, exactly as `search_all_peptides` does —
    // there is nothing to look up for them.
    let matched_suffixes: Vec<&[i64]> = suffix_results
        .iter()
        .filter_map(|r| match r {
            SearchAllSuffixesResult::MaxMatches(suf) | SearchAllSuffixesResult::SearchResult(suf) => {
                Some(suf.as_slice())
            }
            SearchAllSuffixesResult::NoMatches => None
        })
        .collect();

    let retrieval_start = Instant::now();
    let retrieved: Vec<Vec<_>> = matched_suffixes.par_iter().map(|suf| searcher.retrieve_proteins(suf)).collect();
    let retrieval_duration_ns = retrieval_start.elapsed().as_nanos() as u64;

    // Aggregate stats
    let query_hit_count = matched_suffixes.len();
    let suffix_hit_count: usize = matched_suffixes.iter().map(|s| s.len()).sum();
    let protein_hit_count: usize = retrieved.iter().map(|p| p.len()).sum();
    // Cutoff status comes from the search result, not retrieval — an empty-but-cutoff result
    // (max_matches == 0) is a real cutoff even though retrieval finds nothing to look up.
    let cutoff_reached = suffix_results.iter().any(|r| matches!(r, SearchAllSuffixesResult::MaxMatches(_)));

    let total_duration_ns = search_duration_ns + retrieval_duration_ns;
    let (minflt_after, majflt_after) = page_faults();

    let throughput_qps =
        if total_duration_ns > 0 { peptides.len() as f64 / (total_duration_ns as f64 / 1e9) } else { 0.0 };

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
        major_faults: majflt_after.saturating_sub(majflt_before),
        minor_faults: minflt_after.saturating_sub(minflt_before)
    }
}

// ---------------------------------------------------------------------------
// Matrix mode: grid generation
// ---------------------------------------------------------------------------

/// One cell of the trimmed default grid.
#[derive(Clone, Copy, Debug)]
struct GridCell {
    kmer_k: usize,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize
}

/// The `SearchTuning` the CLI asks for. Defaults match `SearchTuning::default()`, so this is
/// that default unless one of the `--validate-*` / `--retrieval-*` flags was passed.
fn tuning_from(args: &Args) -> SearchTuning {
    SearchTuning {
        validate_batch: args.validate_batch,
        validate_prefetch_threshold: args.validate_prefetch_threshold,
        retrieval_prefetch_distance: args.retrieval_prefetch_distance,
        willneed: args.willneed
    }
}

/// The k-mer table sizes the matrix sweeps. `0` means "no table attached".
fn matrix_kmers(args: &Args) -> Vec<usize> {
    if args.matrix_kmer6 { vec![0, 5, 6] } else { vec![0, 5] }
}

/// Expands the grid for one peptide-file bucket. This is the single source of truth for the
/// planned cell list: both `run_matrix` and `print_dry_run` go through it, so `--dry-run`
/// cannot drift from what a real run would execute (`matrix_bench.sh` shells out to
/// `--dry-run` for its expected-config count precisely to rely on that).
///
/// Full sweep is kmer × equate_il × mlp_batch, doubled for tryptic=true/false — except tryptic
/// on the small/medium buckets, which collapses to one representative cell: the last full run
/// showed all 30 small/tryptic cells landing at 653-684 qps (a flat line, not a grid), so
/// sweeping kmer/batch/equate_il there just re-measures a constant at ~1/6 of the whole
/// matrix's wall time. `large` keeps the full sweep since tryptic there is retrieval/search
/// volume bound, not constant.
fn expand_cells(args: &Args, file_bucket: &str) -> Vec<GridCell> {
    let kmers = matrix_kmers(args);
    let batches = &args.matrix_batches;
    let sweep = |tryptic: bool| -> Vec<GridCell> {
        let mut v = Vec::new();
        for &kmer_k in &kmers {
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

/// Swaps the k-mer table of size `k` into `searcher.kmer_table`, returning whatever was
/// previously attached to its owning slot (`table5`/`table6`) first. All swaps are `Option`
/// moves (pointer-sized), so this is cheap to call once per cell regardless of sweep order.
fn ensure_kmer_table(
    searcher: &mut Searcher<SuffixArray>,
    table5: &mut Option<KmerTable>,
    table6: &mut Option<KmerTable>,
    k: usize
) {
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
        _ => None
    };
}

/// Everything `run_cell` needs about one matrix cell that isn't the searcher, the peptides,
/// the CLI args, or the output file: the index-wide facts (fixed for a whole run), the
/// peptide-file facts (fixed per file), and the cell's own coordinates.
#[derive(Clone, Copy)]
struct CellSpec<'a> {
    // -- index-wide
    startup: StartupTiming,
    mapping_type: &'a str,
    sa_type: &'a str,
    sample_rate: u8,
    bits_per_value: usize,
    baseline_memory: u64,
    commit: &'a str,
    // -- per peptide file
    source: &'a str,
    p_min: usize,
    p_max: usize,
    // -- per cell
    theoretical_max: u64,
    equate_il: bool,
    tryptic: bool,
    mlp_batch: usize,
    kmer_k: usize,
    phase: &'a str
}

/// Runs one config for `args.runs` reps, prints the summary line, and appends one aggregated
/// record to `output_file`.
fn run_cell(
    searcher: &Searcher<SuffixArray>,
    peptides: &[String],
    args: &Args,
    spec: CellSpec,
    output_file: &mut File
) -> Result<(), Box<dyn Error>> {
    // Run every rep, then summarise: one record per config with a spread, and the median-qps
    // rep kept as the representative detailed `result`.
    let mut results: Vec<BenchmarkResult> = (0..args.runs)
        .map(|_| {
            run_benchmark(
                searcher,
                peptides,
                args.max_matches,
                spec.equate_il,
                spec.tryptic,
                spec.mlp_batch,
                spec.theoretical_max,
                spec.baseline_memory
            )
        })
        .collect();

    let mut qps: Vec<f64> = results.iter().map(|r| r.throughput_qps).collect();
    qps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let stats = RunStats {
        runs: args.runs,
        qps_min: qps[0],
        qps_p10: percentile(&qps, 0.10),
        qps_p50: percentile(&qps, 0.50),
        qps_p90: percentile(&qps, 0.90),
        qps_max: *qps.last().unwrap()
    };
    let band = if stats.qps_p50 > 0.0 { (stats.qps_p90 - stats.qps_p10) / 2.0 / stats.qps_p50 * 100.0 } else { 0.0 };

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

    eprintln!(
        "  {} {} il={} tr={} batch={} kmer={} tuning{{vb={} vpt={} rpd={}}}  ->  {:.0} qps  (±{:.1}%, p10 {:.0} .. p90 {:.0}){}",
        spec.source,
        spec.phase,
        spec.equate_il,
        spec.tryptic,
        spec.mlp_batch,
        spec.kmer_k,
        searcher.tuning.validate_batch,
        searcher.tuning.validate_prefetch_threshold,
        searcher.tuning.retrieval_prefetch_distance,
        stats.qps_p50,
        band,
        stats.qps_p10,
        stats.qps_p90,
        accept_note,
    );

    let record = BenchmarkRecord {
        version: SCHEMA_VERSION,
        label: args.label.clone(),
        commit: spec.commit.to_string(),
        startup: spec.startup,
        config: BenchmarkConfig {
            sa_type: spec.sa_type.to_string(),
            mapping_type: spec.mapping_type.to_string(),
            sa_storage: sa_server::SA_BACKEND,
            text_storage: sa_server::TEXT_BACKEND,
            proteins_storage: sa_server::PROTEINS_BACKEND,
            mapping_storage: sa_server::MAPPING_BACKEND,
            sample_rate: spec.sample_rate,
            bits_per_value: spec.bits_per_value,
            equate_il: spec.equate_il,
            tryptic: spec.tryptic,
            max_matches: args.max_matches,
            batch_size: spec.mlp_batch,
            kmer_k: spec.kmer_k,
            amount_of_peptides: peptides.len(),
            peptide_length_min: spec.p_min,
            peptide_length_max: spec.p_max,
            peptide_source: spec.source.to_string(),
            validate_batch: searcher.tuning.validate_batch,
            validate_prefetch_threshold: searcher.tuning.validate_prefetch_threshold,
            retrieval_prefetch_distance: searcher.tuning.retrieval_prefetch_distance,
            willneed: searcher.tuning.willneed,
            phase: spec.phase.to_string()
        },
        result: representative,
        stats: Some(stats)
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

    println!("DRY RUN — planned matrix config, no index will be loaded");
    println!("backend        : {}", if cfg!(feature = "mmap") { "mmap" } else { "preloaded" });
    println!("kmer sizes     : {:?}", matrix_kmers(args));
    println!("mlp batches    : {:?}", args.matrix_batches);
    println!("tuning         : {:?}", tuning_from(args));
    println!("runs/config    : {}", args.runs);
    println!();

    let mut grand_total = 0usize;

    for pep_path in &args.matrix_files {
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?");
        // Same expansion the real run uses — a dry run cannot diverge from it.
        let cells = expand_cells(args, source);
        println!("== {} ==", source);
        println!("  grid: {} configs", cells.len());
        for c in &cells {
            println!("    kmer={:<2} il={:<5} tr={:<5} batch={:<3}", c.kmer_k, c.equate_il, c.tryptic, c.mlp_batch);
        }
        grand_total += cells.len();
        println!();
    }

    println!(
        "TOTAL this backend: {} configs x {} runs = {} timed executions",
        grand_total,
        args.runs,
        grand_total * args.runs as usize
    );
    println!("Run once per backend (preloaded, mmap) for the full sweep.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Warmup (never timed)
// ---------------------------------------------------------------------------

/// Peptides per warmup batch — bounds the memory held for the peptide strings themselves
/// while still giving rayon a big enough chunk to parallelise over.
const WARMUP_BATCH_SIZE: usize = 100_000;

/// Touches every page of every mmap-backed region, populating the page cache. Leaves CPU
/// caches and the TLB cold — pair with `warmup_pipeline` for those.
fn warmup_touch_pages(searcher: &Searcher<SuffixArray>) {
    rayon::scope(|s| {
        s.spawn(|_| searcher.sa.touch_all_pages());
        s.spawn(|_| searcher.proteins.touch_all_pages());
        s.spawn(|_| searcher.suffix_index_to_protein.touch_all_pages());
    });
}

/// Pushes `count` peptides from `<index-dir>/warmup.txt` through the full search + retrieval
/// pipeline, in batches, discarding the results. Stops early if the file runs out.
fn warmup_pipeline(searcher: &Searcher<SuffixArray>, args: &Args, count: usize) -> Result<(), Box<dyn Error>> {
    let warmup_path = args.index_dir.join("warmup.txt");
    eprintln!(
        "Warming up with {} peptides from {} (batch size {})...",
        count,
        warmup_path.display(),
        WARMUP_BATCH_SIZE
    );
    let mut lines = BufReader::new(File::open(&warmup_path)?).lines();
    let mut remaining = count;
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
                args.tryptic
            );
            match result {
                SearchAllSuffixesResult::SearchResult(ref suf) | SearchAllSuffixesResult::MaxMatches(ref suf) => {
                    let _ = searcher.retrieve_proteins(suf);
                }
                SearchAllSuffixesResult::NoMatches => {}
            }
        });
    }
    eprintln!("Warmup complete.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Matrix mode: load the index once, then sweep the grid (see `expand_cells`) for each
/// peptide file. The k-mer tables are swapped in/out (moved, not cloned) so the big 6-mer
/// table is built/loaded at most once (and only at all if `--matrix-kmer6` is set). Each
/// config runs `runs` reps but writes a single aggregated record (median-qps rep as `result`,
/// plus a `stats` spread) to <output>/<label>.jsonl.
// Index-wide facts the caller already computed while loading; they are threaded through rather
// than recomputed here, which is what pushes this past the argument limit.
#[allow(clippy::too_many_arguments)]
fn run_matrix(
    mut searcher: Searcher<SuffixArray>,
    args: &Args,
    mapping_type: &str,
    sa_type: &str,
    sample_rate: u8,
    bits_per_value: usize,
    baseline_memory: u64,
    mut startup: StartupTiming
) -> Result<(), Box<dyn Error>> {
    if args.matrix_files.is_empty() {
        return Err("--matrix requires at least one --matrix-files entry".into());
    }
    let kmers = matrix_kmers(args);

    // Build/load the 5-mer table (always in scope — the grid's default kmer set includes it)
    // and, only if requested, the 6-mer table: at 3.06 GB vs 127 MB for a sub-noise-floor
    // difference (see expand_cells / --matrix-kmer6), it isn't worth the build/load cost
    // unless someone explicitly opts in.
    let mut table5: Option<KmerTable> = Some(match &args.kmer5_file {
        Some(p) => {
            eprintln!("Loading 5-mer table from {}...", p.display());
            load_kmer_table_file(p.to_str().unwrap())?
        }
        None => {
            eprintln!("Building 5-mer table...");
            KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 5)
        }
    });
    let mut table6: Option<KmerTable> = if args.matrix_kmer6 {
        Some(match &args.kmer6_file {
            Some(p) => {
                eprintln!("Loading 6-mer table from {}...", p.display());
                load_kmer_table_file(p.to_str().unwrap())?
            }
            None => {
                eprintln!("Building 6-mer table...");
                KmerTable::build_from_sa(&searcher.sa, searcher.proteins.text(), 6)
            }
        })
    } else {
        None
    };

    // Warm the page cache once (matters for mmap); CPU caches warm over the repeated runs.
    eprintln!("Warming up (touching all pages)...");
    let warmup_start = Instant::now();
    warmup_touch_pages(&searcher);
    startup.warmup_ms = warmup_start.elapsed().as_millis() as u64;
    startup.load_total_ms =
        startup.load_sa_ms + startup.load_proteins_ms + startup.load_mapping_ms + startup.kmer_table_ms;
    eprintln!(
        "Startup: load {} ms (sa {} / proteins {} / mapping {}), warmup {} ms",
        startup.load_total_ms, startup.load_sa_ms, startup.load_proteins_ms, startup.load_mapping_ms, startup.warmup_ms
    );

    // Theoretical memory footprint per k-mer table size — index-wide, independent of the
    // peptide file, so this is computed once per k value up front rather than per cell (it
    // walks every protein's metadata length, which is not free at full-DB scale).
    let mut theoretical_by_k: HashMap<usize, u64> = HashMap::new();
    for &k in &kmers {
        ensure_kmer_table(&mut searcher, &mut table5, &mut table6, k);
        let m =
            if args.no_theoretical_memory { 0 } else { theoretical_memory(&searcher, mapping_type, proteins_mapped()) };
        theoretical_by_k.insert(k, m);
    }

    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new().create(true).write(true).truncate(true).open(&output_path)?;
    let commit = env!("GIT_COMMIT_HASH").to_string();
    let tuning = tuning_from(args);

    for pep_path in &args.matrix_files {
        let peptides: Vec<String> = BufReader::new(File::open(pep_path)?)
            .lines()
            .take(args.amount_of_peptides)
            .collect::<Result<_, _>>()?;
        if peptides.len() < args.amount_of_peptides {
            return Err(format!(
                "{} has only {} peptides, need {}",
                pep_path.display(),
                peptides.len(),
                args.amount_of_peptides
            )
            .into());
        }
        let (p_min, p_max) =
            peptides.iter().fold((usize::MAX, 0usize), |(lo, hi), p| (lo.min(p.len()), hi.max(p.len())));
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?").to_string();

        // Everything the cells of this file share; the per-cell fields are filled in below.
        let base_spec = CellSpec {
            startup,
            mapping_type,
            sa_type,
            sample_rate,
            bits_per_value,
            baseline_memory,
            commit: &commit,
            source: &source,
            p_min,
            p_max,
            theoretical_max: 0,
            equate_il: false,
            tryptic: false,
            mlp_batch: 1,
            kmer_k: 0,
            phase: "grid"
        };

        eprintln!("-- {} : grid --", source);
        for cell in expand_cells(args, &source) {
            ensure_kmer_table(&mut searcher, &mut table5, &mut table6, cell.kmer_k);
            searcher.tuning = tuning;
            let spec = CellSpec {
                theoretical_max: theoretical_by_k[&cell.kmer_k],
                equate_il: cell.equate_il,
                tryptic: cell.tryptic,
                mlp_batch: cell.mlp_batch,
                kmer_k: cell.kmer_k,
                ..base_spec
            };
            run_cell(&searcher, &peptides, args, spec, &mut output_file)?;
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

    if args.mlp_batch == 0 {
        return Err("--mlp-batch must be >= 1 (1 = scalar)".into());
    }

    // --dry-run never touches the index (which may not even exist locally) — it only expands
    // and prints the config list the matrix sweep would run.
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
        _ => "unknown"
    };

    // Warm up the Rayon thread pool so its thread stacks are included in the baseline.
    // Without this, thread stacks (~2 MB × N threads) would be attributed to index memory
    // because Rayon initialises lazily on the first par_iter() call inside run_benchmark.
    let _ = (0..rayon::current_num_threads()).into_par_iter().sum::<usize>();

    // Baseline RSS after Rayon is warm but before any index data is loaded
    let baseline_memory = measure_process_memory();

    // Load index. Each load is timed separately: a `preloaded-*` feature only moves the cost of
    // *its own* structure, so a single total would not say which one paid.
    let mut startup = StartupTiming {
        load_sa_ms: 0,
        load_proteins_ms: 0,
        load_mapping_ms: 0,
        kmer_table_ms: 0,
        warmup_ms: 0,
        load_total_ms: 0
    };

    eprintln!("Loading suffix array from {}...", sa_path.display());
    let t0 = Instant::now();
    let suffix_array = load_suffix_array_file(sa_path.to_str().unwrap())?;
    startup.load_sa_ms = t0.elapsed().as_millis() as u64;
    eprintln!(
        "  {} items, {} bits/value, sample rate {} ({} ms)",
        suffix_array.len(),
        suffix_array.bits_per_value(),
        suffix_array.sample_rate(),
        startup.load_sa_ms
    );

    eprintln!("Loading proteins from {}...", proteins_path.display());
    let t0 = Instant::now();
    let proteins = load_proteins_file(proteins_path.to_str().unwrap())?;
    startup.load_proteins_ms = t0.elapsed().as_millis() as u64;
    eprintln!("  {} ms", startup.load_proteins_ms);

    eprintln!("Loading mapping from {} (type: {})...", mapping_path.display(), mapping_type_str);
    let t0 = Instant::now();
    let mapping = load_mapping_file(mapping_path.to_str().unwrap())?;
    startup.load_mapping_ms = t0.elapsed().as_millis() as u64;
    eprintln!("  {} ms", startup.load_mapping_ms);

    let sa_type = if suffix_array.bits_per_value() == 64 { "original" } else { "compressed" };
    let sample_rate = suffix_array.sample_rate();
    let bits_per_value = suffix_array.bits_per_value();

    let mut searcher = Searcher::new(suffix_array, proteins, mapping);

    if args.matrix {
        return run_matrix(
            searcher,
            &args,
            mapping_type_str,
            sa_type,
            sample_rate,
            bits_per_value,
            baseline_memory,
            startup
        );
    }

    let t0 = Instant::now();
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
    startup.kmer_table_ms = t0.elapsed().as_millis() as u64;
    startup.load_total_ms =
        startup.load_sa_ms + startup.load_proteins_ms + startup.load_mapping_ms + startup.kmer_table_ms;

    // Apply the SearchTuning knobs from the CLI (defaults match SearchTuning::default(), so
    // this is a no-op unless the caller overrides one of --validate-batch/etc).
    searcher.tuning = tuning_from(&args);

    let theoretical_max = if args.no_theoretical_memory {
        eprintln!("Theoretical max memory: skipped (--no-theoretical-memory)");
        0
    } else {
        let m = theoretical_memory(&searcher, mapping_type_str, proteins_mapped());
        eprintln!("Theoretical max memory: {} bytes ({:.1} MB)", m, m as f64 / 1_048_576.0);
        m
    };

    // Load peptides: either all at once from a file (for sequential chunking across runs)
    // or generate fresh random peptides for each run.
    let all_peptides: Option<Vec<String>> = if let Some(ref path) = args.peptide_file {
        let peptides: Vec<String> = BufReader::new(File::open(path)?).lines().collect::<Result<_, _>>()?;
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
            )
            .into());
        }
        eprintln!("Loaded {} peptides from {}", peptides.len(), path.display());
        Some(peptides)
    } else {
        None
    };

    let peptide_source = match &args.peptide_file {
        Some(p) => p.display().to_string(),
        None => "random".to_string()
    };

    // Optional warmup pass. Not part of the timed benchmark, but recorded: for mmap builds
    // `--warmup all` faults in the whole index, and preloading a structure removes it from that
    // sweep, so this is where part of the load cost reappears.
    let warmup_start = Instant::now();
    match &args.warmup {
        None => {}
        Some(WarmupMode::All) => {
            eprintln!("Warming up: touching all mmap pages...");
            warmup_touch_pages(&searcher);
            eprintln!("Warmup complete.");
        }
        Some(WarmupMode::Count(warmup_count)) => {
            warmup_pipeline(&searcher, &args, *warmup_count)?;
        }
        Some(WarmupMode::AllThenCount(warmup_count)) => {
            eprintln!("Warming up: touching all mmap pages...");
            warmup_touch_pages(&searcher);
            eprintln!("Page warmup complete.");
            warmup_pipeline(&searcher, &args, *warmup_count)?;
        }
    }
    if args.warmup.is_some() {
        startup.warmup_ms = warmup_start.elapsed().as_millis() as u64;
        eprintln!("Warmup took {} ms.", startup.warmup_ms);
    }
    eprintln!(
        "Startup: load {} ms (sa {} / proteins {} / mapping {} / kmer {}), warmup {} ms",
        startup.load_total_ms,
        startup.load_sa_ms,
        startup.load_proteins_ms,
        startup.load_mapping_ms,
        startup.kmer_table_ms,
        startup.warmup_ms
    );

    // Prepare output file
    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new().create(true).write(true).truncate(true).open(&output_path)?;

    eprintln!();
    eprintln!("Starting {} benchmark run(s) — results → {}", args.runs, output_path.display());
    eprintln!();

    // RNG for random mode. Seeded → reproducible stream across runs and invocations; otherwise
    // seeded from OS entropy. Unused in file mode (peptides come from disk).
    let mut rng = match args.seed {
        Some(s) => {
            eprintln!("Random peptide seed: {}", s);
            StdRng::seed_from_u64(s)
        }
        None => StdRng::from_entropy()
    };

    for run in 1..=args.runs {
        // In file mode: consume the next chunk sequentially.
        // In random mode: generate a fresh set of random peptides.
        let run_peptides: Vec<String> = match &all_peptides {
            Some(all) => {
                let start = (run as usize - 1) * args.amount_of_peptides;
                all[start..start + args.amount_of_peptides].to_vec()
            }
            None => {
                generate_peptides(&mut rng, args.amount_of_peptides, args.peptide_length_min, args.peptide_length_max)
            }
        };

        let (p_min, p_max) =
            run_peptides.iter().fold((usize::MAX, 0usize), |(lo, hi), p| (lo.min(p.len()), hi.max(p.len())));

        let config = BenchmarkConfig {
            sa_type: sa_type.to_string(),
            mapping_type: mapping_type_str.to_string(),
            sa_storage: sa_server::SA_BACKEND,
            text_storage: sa_server::TEXT_BACKEND,
            proteins_storage: sa_server::PROTEINS_BACKEND,
            mapping_storage: sa_server::MAPPING_BACKEND,
            sample_rate,
            bits_per_value,
            equate_il: args.equate_il,
            tryptic: args.tryptic,
            max_matches: args.max_matches,
            batch_size: args.mlp_batch,
            kmer_k: searcher.kmer_table.as_ref().map_or(0, |t| t.k),
            amount_of_peptides: run_peptides.len(),
            peptide_length_min: p_min,
            peptide_length_max: p_max,
            peptide_source: peptide_source.clone(),
            validate_batch: searcher.tuning.validate_batch,
            validate_prefetch_threshold: searcher.tuning.validate_prefetch_threshold,
            retrieval_prefetch_distance: searcher.tuning.retrieval_prefetch_distance,
            willneed: searcher.tuning.willneed,
            phase: "single".to_string()
        };

        let result = run_benchmark(
            &searcher,
            &run_peptides,
            args.max_matches,
            args.equate_il,
            args.tryptic,
            args.mlp_batch,
            theoretical_max,
            baseline_memory
        );

        let record = BenchmarkRecord {
            version: SCHEMA_VERSION,
            label: args.label.clone(),
            commit: env!("GIT_COMMIT_HASH").to_string(),
            config,
            startup,
            result,
            stats: None
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
             |  hits: {} queries, {} suffixes, {} proteins{}{}{}",
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
            if record.result.major_faults > 0 {
                format!("  |  {} major faults", record.result.major_faults)
            } else {
                String::new()
            },
        );
    }

    eprintln!();
    eprintln!("Done. {} lines written to {}", args.runs, output_path.display());

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// A real `/proc/self/stat` prefix. The fault counters are the whole instrument of the
    /// RAM-scaling sweep, and a mis-indexed parse returns zeros — which reads as "this run took
    /// no page faults", a perfectly plausible result. Nothing else would catch that.
    #[test]
    fn parses_fault_counts_from_proc_stat() {
        let stat = "1234 (sa-benchmarks) R 1200 1234 1234 0 -1 4194304 6789 0 12 0 5 3 0 0 20 0 8 0 123456";
        assert_eq!(parse_proc_stat_faults(stat), (6789, 12));
    }

    /// The process name can contain spaces and parentheses, which is why the parse starts after
    /// the last `)` rather than at the first space. Splitting on whitespace first would shift
    /// every field and silently report the wrong numbers.
    #[test]
    fn comm_field_with_spaces_and_parens_does_not_shift_the_fields() {
        let stat = "42 (my (weird) proc) S 1 42 42 0 -1 4194304 111 0 222 0 5 3 0 0 20 0 8 0 99";
        assert_eq!(parse_proc_stat_faults(stat), (111, 222));
    }

    /// Truncated or unexpected input reports zeros rather than panicking mid-benchmark. The
    /// counters are diagnostic, so losing them must not take a multi-hour run down with them.
    #[test]
    fn malformed_input_reports_zeros() {
        assert_eq!(parse_proc_stat_faults(""), (0, 0));
        assert_eq!(parse_proc_stat_faults("1234 (no-close-paren 1 2 3"), (0, 0));
        assert_eq!(parse_proc_stat_faults("1234 (short) R 1 2"), (0, 0));
    }
}
