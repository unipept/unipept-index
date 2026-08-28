//! Throughput/memory measurement harness for the suffix-array index, at commit fbc9328.
//!
//! # What this is
//!
//! A port of `sa-benchmarks/src/main.rs` from the `feature/preloaded-sa-improvements` branch, cut
//! down to what this commit can actually do, for one purpose: producing a BASELINE that the
//! branch's report machinery can read. It exists because fbc9328 predates the harness entirely, so
//! the only way to say what the branch changed is to measure the commit before it with the same
//! instrument.
//!
//! It is not a general-purpose harness and is not meant to grow into one. It measures one
//! configuration of one commit.
//!
//! # The schema is the contract
//!
//! `SCHEMA_VERSION`, every field of `BenchmarkRecord`, and the JSONL layout are IDENTICAL to the
//! branch's. That is the whole point — the driver's `bench` package joins cells on
//! `(peptide_source, equate_il, tryptic, kmer_k, amount_of_peptides)` and splits them by the `arm`
//! dim, so records from here drop into a branch session directory and land beside the branch's own
//! arms in one table. Do not "improve" a field. If the branch's schema moves, this moves with it or
//! the baseline stops being readable.
//!
//! # What this commit does not have, and what was done about it
//!
//! * **No `sa-server` library.** At fbc9328 `sa-server` is a binary with no lib target, so the
//!   loaders and the storage constants the branch's harness imports do not exist. They are inlined
//!   below, and the four `*_BACKEND` constants are pinned to `"preloaded"` — which is not a
//!   convenience, it is what this tree does: every structure is read into owned memory.
//!
//! * **No k-mer bound table.** The whole axis is gone, `kmer_k` is pinned to 0, and a grid cell
//!   asking for a table is REJECTED rather than run. Compare against the branch's `kmer_k = 0`
//!   cells (`suites/kmer.toml` sweeps `[0, 5, 6]`); the branch's `defaults` suite pins `kmer = [5]`
//!   and those cells have no counterpart here.
//!
//! * **No batched search.** `search_all_matching_suffixes_batched` and its MLP interleaving are the
//!   branch's. Phase 1 here is the flat `par_iter` over peptides that this commit's own
//!   `search_all_peptides` uses, so the measurement reflects this commit rather than a hybrid.
//!
//! * **No storage arms.** There is exactly one build, so a suite run against this binary has one
//!   arm and the driver's feature-forwarding checks have nothing to check.
//!
//! * **No `proteins.bin` / `mapping.bin`.** Both are branch artefacts. The proteins are parsed from
//!   the `proteins.tsv` the index was built from and the mapping is rebuilt from the text, both at
//!   startup — which is what this commit's server does. `mapping.bin` is still read for its type
//!   byte alone, so the mapping rebuilt here is the same KIND the branch's runs walked.
//!
//! Two `startup` fields are therefore NOT comparable across the two trees, and are recorded anyway
//! because "what did each pay before its first query" is a real question: `load_proteins_ms` (a TSV
//! parse and a bit-packing here, a file read there) and `load_mapping_ms` (a build here, a read
//! there). Throughput, the phase split, and the fault counters are comparable.
//!
//! # Running it
//!
//! ```text
//! cargo build --release -p sa-benchmarks
//! ./sa-benchmarks/run.sh legacy --out <session-dir>
//! ```
//!
//! See `sa-benchmarks/README.md` in this worktree for the full procedure and how the resulting
//! records are merged into a branch session.

use std::{
    collections::BTreeMap,
    error::Error,
    fs::{File, OpenOptions, create_dir_all},
    io::{BufRead, BufReader, Read, Write},
    ops::Deref,
    path::{Path, PathBuf},
    time::Instant
};

use clap::Parser;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use sa_compression::load_compressed_suffix_array;
use sa_index::{
    SuffixArray,
    binary::load_suffix_array,
    peptide_search::{ProteinInfo, SearchResult},
    sa_searcher::{BitVecSearcher, DenseSearcher, SearchAllSuffixesResult, Searcher, SparseSearcher}
};
use sa_mappings::proteins::Proteins;
use serde::{Deserialize, Serialize};
use sysinfo::{Pid, System};

/// Schema version — increment when the output JSON format changes.
/// v2: matrix records aggregate `runs` reps into one line and carry a `stats` spread.
/// v3: every `SearchTuning` field is recorded in `config` (was previously implicit/default),
///   plus a `phase` tag so records from different sweeps are groupable in one jsonl file;
///   `result` gains `candidates_examined` / `candidates_accepted` (removed again in v15).
/// v12: `SearchTuning` is gone — its fields are compile-time constants now — so `config.tuning`
///   and `config.tuning_defaults` are no longer written. A reader that wants to know what a
///   record ran at reads the commit, not the record.
/// v4: the OFAT and confirm sweeps were retired once they had settled which knobs matter, so
///   `config.ofat_baseline` / `config.ofat_knob` are gone and `config.phase` is now only
///   "single" (non-matrix CLI run) or "grid" (matrix sweep).
/// v5: storage is chosen per structure, so the single `use_mmap` bool is replaced by
///   `sa_storage` / `text_storage` / `proteins_storage` / `mapping_storage`, each
///   "mmap" or "preloaded".
/// v6: records carry a `startup` section timing each structure's load and the warmup pass.
///   Preloading a structure moves cost from steady-state page faults to startup, and that
///   trade was previously invisible.
/// v7: `result` gains `major_faults` / `minor_faults`, counted across the timed region. They are
///   what separates "slow because it is waiting on I/O" from "slow for some other reason" when
///   the index does not fit in RAM.
/// v8: records carry `suite` and a free-form `dims` map, both supplied by the caller. A sweep
///   dimension this binary knows nothing about — a cgroup ceiling, a `RAYON_NUM_THREADS` value,
///   which storage arm was built — used to survive only in the file name (`c167_pprot_t96.jsonl`),
///   so every driver script needed its own file-name parser and every one of them could disagree
///   about what a run meant. Now the coordinates travel inside the record.
/// v10: matrix grids can be supplied by the driver (`--grid-file`), so `config` gains `sweep` and
///   `grid_slot`. `sweep` names the block of the suite a cell came from; two blocks may measure
///   the same configuration at different rep counts, and comparing across them would compare
///   precisions rather than configurations. `grid_slot` distinguishes repeats of one
///   configuration inside a single process — the drift cadence writes the same tuning point
///   several times, and keyed on `config` alone those records would collapse into one.
/// v11: `peptide_source` is the name the SUITE gave the file, not the file's stem, when a grid cell
///   supplies one (`bucket`). A profile maps `mixed` to whatever that machine calls the file, so
///   keying a report on the stem meant the same bucket was called `mixed` on one box and
///   `peptides_5_50` on another — and a baseline comparison across two such boxes silently
///   matched nothing. Cells without a `bucket` still record the stem.
/// v12: the two phases production runs after retrieval are measured too, when a cell asks for them
///   (`response`): turning each `ProteinRef` into a `ProteinInfo` — an fa-compression decode plus
///   an accession `String` — and serialising the lot to JSON, which is what `sa-server` returns.
///   They are recorded BESIDE `total_duration_ns` and never inside it: widening the timed region
///   would change what every suite's throughput means and invalidate every baseline, including
///   the regression gate. `throughput_qps` is still search plus retrieval, exactly as in v11.
/// v13: v12's `decode_duration_ns` and `serialise_duration_ns` are replaced by a single
///   `response_duration_ns`, because the pair they formed was not measurable against the rest of
///   the record. v12 timed both of them SERIALLY while it timed search and retrieval in parallel,
///   so on a 12-core box the decode share was overstated by up to the core count, and the
///   "measured share" heatmap built on that comparison was correspondingly wrong. v12 also
///   serialised a bare `Vec<Vec<ProteinInfo>>` rather than the `Vec<SearchResult>` the server
///   returns, so `response_bytes` omitted every peptide's `sequence`, `cutoff_used` and object
///   framing and under-reported the response.
///
///   v13 measures the whole proteins-to-bytes phase the way production runs it — the decode
///   parallel across peptides, the JSON in the shape `sa-server` actually returns — and records it
///   as one number. **v12 and v13 shares are not comparable**; a baseline comparison across the
///   boundary will show a large apparent shift in the response phase that is the fix, not a
///   regression. Expect `response_bytes` to grow for the same reason.
///
///   The isolated decoder cost is no longer a field here: it is not a request phase, and
///   `cargo bench -p fa-compression` measures it properly. `total_duration_ns` and
///   `throughput_qps` are unchanged and still mean search plus retrieval.
/// v14: `stats` gains `search_p50_ns` / `retrieval_p50_ns` / `response_p50_ns`, the phase times
///   pooled over the same reps as `qps_p50`. Readers previously took the phase split off
///   `result`, which is a single rep, and printed it beside a throughput pooled over twenty —
///   so the two could disagree, and on the mixed file they disagreed in direction. Records
///   written before v14 have no such fields and readers must fall back to `result`, which is
///   what they always did.
/// v15: `result` loses `search_bounds_ns`, `match_iter_ns`, `candidates_examined` and
///   `candidates_accepted`. They came from `sa-index`'s `measure` feature, which has been
///   removed: the counters were atomics on cache lines every rayon worker shared, so the only
///   builds that carried them measured themselves ~2% slower than the ones that shipped. The
///   questions they settled — the binary-search/range-scan split, and tryptic's candidate
///   acceptance rate — are recorded in `sa-index`'s crate docs and in this crate's README.
const SCHEMA_VERSION: u32 = 15;

/// Canonical 20 amino acids used for random peptide generation
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// ---------------------------------------------------------------------------
// What sa-server's library provided on the branch, and does not provide here
// ---------------------------------------------------------------------------

// At this commit `sa-server` is a binary with no library target, and there is no per-structure
// storage selection to report: every structure is read into owned memory at startup. These four
// are the constants the branch's `sa-server` exports, pinned to what this tree actually does, so a
// record from here lands in the same `preloaded` column as the branch's fully-preloaded arm rather
// than in a column of its own.
const SA_BACKEND: &str = "preloaded";
const TEXT_BACKEND: &str = "preloaded";
const PROTEINS_BACKEND: &str = "preloaded";
const MAPPING_BACKEND: &str = "preloaded";

fn backend_summary() -> String {
    format!(
        "sa={} text={} proteins={} mapping={}",
        SA_BACKEND, TEXT_BACKEND, PROTEINS_BACKEND, MAPPING_BACKEND
    )
}

/// The searcher this run is using, selected to match the mapping type of the index on disk.
///
/// The branch has one `ActiveSearcher` whose mapping is chosen by `mapping.bin`; here the mapping
/// is rebuilt from the text at startup, so the choice is ours to make and we make it the same way —
/// off `mapping.bin`'s type byte. Comparing a run that walks a bitvec mapping against one that
/// walks a sparse mapping would be measuring the mapping, not the commit.
enum ActiveSearcher {
    Sparse(SparseSearcher),
    BitVec(BitVecSearcher),
    Dense(DenseSearcher)
}

impl Deref for ActiveSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Searcher {
        match self {
            ActiveSearcher::Sparse(searcher) => searcher,
            ActiveSearcher::BitVec(searcher) => searcher,
            ActiveSearcher::Dense(searcher) => searcher
        }
    }
}

impl ActiveSearcher {
    fn new(mapping_type: &str, sa: SuffixArray, proteins: Proteins) -> Result<Self, Box<dyn Error>> {
        Ok(match mapping_type {
            "sparse" => ActiveSearcher::Sparse(SparseSearcher::new(sa, proteins)),
            "bitvec" => ActiveSearcher::BitVec(BitVecSearcher::new(sa, proteins)),
            "dense" => ActiveSearcher::Dense(DenseSearcher::new(sa, proteins)),
            other => return Err(format!("unknown mapping type '{}'", other).into())
        })
    }
}

/// Reads `sa.bin`, dispatching on the packing declared in its first byte.
///
/// Lifted verbatim from this commit's `sa-server/src/main.rs` — the header is the same one the
/// branch documents in `sa-index/src/array/mod.rs`, which is what lets both trees read one index
/// file with no rebuild.
fn load_suffix_array_file(file: &str) -> Result<SuffixArray, Box<dyn Error>> {
    let mut sa_file = File::open(file)?;
    let mut reader = BufReader::new(&mut sa_file);

    let mut bits_per_value_buffer = [0_u8; 1];
    reader
        .read_exact(&mut bits_per_value_buffer)
        .map_err(|_| "Could not read the flags from the binary file")?;
    let bits_per_value = bits_per_value_buffer[0];

    if bits_per_value == 64 {
        load_suffix_array(&mut reader)
    } else {
        load_compressed_suffix_array(&mut reader, bits_per_value as usize)
    }
}

/// The branch's `json_chunk`, reimplemented against this commit's owned `SearchResult`.
///
/// Byte-for-byte the same framing — a leading `,` per chunk, closed into an array by
/// [`frame_chunks`] — so `response_bytes` from the two trees describe the same response. What it
/// cannot mirror is the allocation: this commit's `ProteinInfo` clones the accession and decodes
/// the annotations into an owned `String` per hit, where the branch borrows both. That difference
/// is a change the branch made, so it belongs in the measurement rather than being papered over.
fn json_chunk(result: &SearchResult) -> Vec<u8> {
    let estimate = 64
        + result.sequence.len()
        + result
            .proteins
            .iter()
            .map(|protein| 64 + protein.uniprot_accession.len() + protein.functional_annotations.len())
            .sum::<usize>()
        + 1;
    let mut chunk = Vec::with_capacity(estimate);
    chunk.push(b',');
    serde_json::to_writer(&mut chunk, result).expect("serialising a SearchResult into a Vec cannot fail");
    chunk
}

/// The branch's `frame_chunks`, unchanged.
fn frame_chunks(chunks: &mut Vec<Vec<u8>>) {
    match chunks.first_mut() {
        None => chunks.push(b"[]".to_vec()),
        Some(first) => {
            assert_eq!(
                first.first(),
                Some(&b','),
                "frame_chunks expects json_chunk output: every chunk must be non-empty and start with ','"
            );
            first[0] = b'[';
            chunks.last_mut().expect("non-empty").push(b']');
        }
    }
}

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

    // There is no `--kmer-table-file` / `--build-kmer-table` here: the k-mer bound table does not
    // exist at this commit. A grid cell asking for one is rejected in `load_grid_file` rather than
    // ignored, so a suite cannot quietly report table-accelerated cells that ran without a table.
    /// Run a parameter matrix in one process: loads the index once, then sweeps the cell list from
    /// `--grid-file` for each `--matrix-files` entry. Writes one aggregated record per cell to
    /// `<output>/<label>.jsonl`.
    #[arg(long)]
    matrix: bool,

    /// Matrix mode: comma-separated peptide files; each becomes one "file" dimension. A grid cell
    /// may name a file stem to restrict itself to one of them.
    #[arg(long, value_delimiter = ',')]
    matrix_files: Vec<PathBuf>,

    /// Matrix mode: the cell list to sweep, as JSONL — one `GridCell` per line. Required by
    /// `--matrix`.
    ///
    /// The grid is the driver's, not the harness's. Every cell names its k-mer size and its search
    /// options, and may carry its own `runs` and `amount`: a screening sweep looking for 10%
    /// effects does not need the rep count that resolves 4%, and paying for it on every cell is
    /// most of what makes a wide sweep slow.
    #[arg(long)]
    grid_file: Option<PathBuf>,

    /// Stop a matrix cell early once its own p10..p90 half-spread is under this percentage.
    ///
    /// A fixed rep count spends the worst cell's budget on every cell, and most cells are quiet.
    /// With a target band each cell runs `--min-runs` reps, then keeps going only while it is still
    /// too noisy to read, up to `--runs`. The reps actually taken are recorded in `stats.runs`, so
    /// a cell that hit the cap is visible rather than merely noisy. Off by default: `defaults` is
    /// the regression gate and wants a fixed rep count, so that two sessions are comparable rather
    /// than each stopping wherever its own noise let it. (The 3.9% floor was measured at 20 reps;
    /// `defaults.toml` now runs 100 fixed, and says why.)
    #[arg(long)]
    runs_target_band: Option<f64>,

    /// Matrix mode: reps to run before `--runs-target-band` may stop a cell. No effect without it.
    #[arg(long, default_value_t = 5)]
    min_runs: u32,

    /// Print the planned config list for the matrix sweep and exit, without
    /// loading the index. Use this to eyeball a sweep before committing a multi-hour run.
    #[arg(long)]
    dry_run: bool,

    /// Skip the theoretical memory calculation, reporting `theoretical_max_memory: 0`.
    ///
    /// That calculation walks *every* protein's metadata (see `theoretical_memory`), which on an
    /// mmap backend faults the entire metadata section in before anything is timed. Harmless when
    /// the index fits in RAM, but under a cgroup memory cap it both pre-warms what the run is
    /// supposed to be faulting on demand and spends the budget being measured. Off by default, so
    /// ordinary runs still report the figure.
    #[arg(long)]
    no_theoretical_memory: bool,

    /// Name of the suite this invocation belongs to, copied verbatim into every record.
    ///
    /// Set by the driver (`sa-benchmarks/run.sh <suite>`); a hand-run invocation leaves it at
    /// "adhoc" so its records are still groupable but obviously not part of a sweep.
    #[arg(long, default_value = "adhoc")]
    suite: String,

    /// Sweep coordinate for this invocation, as `key=value`. Repeatable.
    ///
    /// These are the facts about a cell that this binary cannot observe about itself — the cgroup
    /// ceiling it was launched under, its `RAYON_NUM_THREADS`, which storage arm was built, which
    /// slot of an ABBA ordering it occupies. The driver knows them, so the driver passes them, and
    /// they end up in `dims` on every record rather than encoded in the output file name.
    #[arg(long = "dim", value_parser = parse_dim, value_name = "KEY=VALUE")]
    dims: Vec<(String, String)>
}

/// Parses one `--dim key=value`. Splits on the FIRST `=` so values may contain `=` themselves
/// (a feature list, a path). An empty key is rejected: it would silently collide in the map.
fn parse_dim(s: &str) -> Result<(String, String), String> {
    let (key, value) = s.split_once('=').ok_or_else(|| format!("expected KEY=VALUE, got '{}'", s))?;
    if key.is_empty() {
        return Err(format!("empty key in '{}'", s));
    }
    Ok((key.to_string(), value.to_string()))
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
    load_total_ms: u64,

    // What the page sweep actually did
    //
    // Elapsed milliseconds alone cannot tell a sweep that read from disk apart from one that hit
    // the page cache, and the two differ by an order of magnitude on the same bytes. Since the
    // suites share a page cache and run their arms in a fixed order, a later arm can sweep the
    // same structure the earlier one left resident and look eight times faster for it. Bytes and
    // faults are what make that visible instead of being read as an arm difference.
    /// Bytes swept per structure by `--warmup all`. 0 for a structure with nothing mapped, which
    /// is itself the answer to "why was this arm's warmup cheap".
    warmup_sa_bytes: u64,
    warmup_proteins_bytes: u64,
    warmup_mapping_bytes: u64,
    /// Wall time of each structure's sweep. They run concurrently, so these do not sum to
    /// `warmup_ms` — the slowest one sets it.
    warmup_sa_ms: u64,
    warmup_proteins_ms: u64,
    warmup_mapping_ms: u64,
    /// Faults taken across the load phase and across the warmup phase. The timed region has its
    /// own counters and deliberately excludes both; these exist so a cold sweep is legible as a
    /// cold sweep. A major fault here is a read from the device.
    load_major_faults: u64,
    load_minor_faults: u64,
    warmup_major_faults: u64,
    warmup_minor_faults: u64
}

impl StartupTiming {
    /// Bytes the page sweep touched, across every structure.
    fn warmup_bytes(&self) -> u64 {
        self.warmup_sa_bytes + self.warmup_proteins_bytes + self.warmup_mapping_bytes
    }
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
    /// k of the attached k-mer table (0 = no table).
    kmer_k: usize,
    amount_of_peptides: usize,
    peptide_length_min: usize,
    peptide_length_max: usize,
    peptide_source: String,
    /// Which sweep produced this record: "single" (non-matrix CLI run) or "grid" (the trimmed
    /// default matrix grid).
    phase: String,
    /// The suite block this cell came from, verbatim from the grid file ("" for the built-in grid).
    ///
    /// Blocks carry their own rep and query counts, so two of them may measure the same
    /// configuration at different precision. Comparing across them would compare how carefully each
    /// was measured rather than what it measured, which is why the block travels with the record
    /// instead of being reconstructed from the cell's other coordinates.
    sweep: String,
    /// Distinguishes repeats of one configuration inside a single process ("a" unless set).
    ///
    /// The drift cadence re-measures the reference configuration every few cells; keyed on `config`
    /// alone those records would be indistinguishable and collapse into one. With the slot they
    /// stay separate and their trend over the process becomes readable.
    grid_slot: String
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
    /// Nanoseconds for everything between "we have the proteins" and "the client has bytes":
    /// decoding each hit's functional annotations, allocating its accession, and serialising the
    /// whole answer to JSON. 0 unless the cell set `response`.
    ///
    /// Measured the way production runs it, which is the whole point of the field — the decode
    /// parallel across peptides and the JSON in the shape `sa-server` returns. v12 split this in two
    /// and timed both serially, which inflated the decode share by up to the core count; see the
    /// v13 note on `SCHEMA_VERSION`.
    response_duration_ns: u64,
    /// Bytes of JSON the request would have returned, including per-peptide `sequence`,
    /// `cutoff_used` and object framing. 0 unless the cell set `response`.
    response_bytes: u64,
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
    qps_max: f64,

    /// Median phase times across the same reps `qps_p50` is taken over.
    ///
    /// These exist because the alternative was reading the phase split off `result`, which is one
    /// rep — the representative one. A phase decomposition drawn from a single rep and printed
    /// beside a throughput pooled over twenty of them can disagree with it, and on the mixed file
    /// it did: the phase times made `pprot` 11% slower than `mmap` where the pooled qps made it
    /// 0.3-2.8% faster. Whichever was right, a table that answers "where does this configuration
    /// spend itself" has to rest on the same reps as the table that says how fast it is.
    ///
    /// The candidate counters stay on `result` and are deliberately not aggregated: they count work
    /// the search did, which is a property of the configuration rather than of a rep, so every rep
    /// reports the same value and a median of them says nothing new.
    search_p50_ns: u64,
    retrieval_p50_ns: u64,
    /// 0 unless the cell set `response`, matching [`BenchmarkResult::response_duration_ns`].
    response_p50_ns: u64
}

#[derive(Serialize)]
struct BenchmarkRecord {
    version: u32,
    label: String,
    commit: String,
    /// Which suite produced this record ("adhoc" for a hand-run invocation).
    suite: String,
    /// The driver's sweep coordinates for this cell — see `Args::dims`. `BTreeMap` rather than
    /// `HashMap` so the serialized key order is stable and two records diff cleanly.
    dims: BTreeMap<String, String>,
    config: BenchmarkConfig,
    /// Absent in records written before schema v6.
    startup: StartupTiming,
    result: BenchmarkResult,
    /// Per-config throughput spread over all reps (matrix mode only; omitted otherwise).
    #[serde(skip_serializing_if = "Option::is_none")]
    stats: Option<RunStats>
}

/// Half the p10..p90 spread of these reps, as a percent of their median.
///
/// The same statistic the driver calls `band`, computed here so a cell can decide whether it has
/// run enough reps yet. Costs one sort per rep, against a rep that pushes tens of thousands of
/// peptides through the index — not a cost worth avoiding.
fn band_of(results: &[BenchmarkResult]) -> f64 {
    let mut qps: Vec<f64> = results.iter().map(|result| result.throughput_qps).collect();
    qps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = percentile(&qps, 0.50);
    if median <= 0.0 {
        return f64::INFINITY;
    }
    (percentile(&qps, 0.90) - percentile(&qps, 0.10)) / 2.0 / median * 100.0
}

/// Median of one field across a cell's reps, in nanoseconds.
///
/// Sorts a copy rather than the results themselves: the caller still needs them in their original
/// order to pick the representative rep, and a helper that reorders its input behind the caller's
/// back is the kind of thing that works until someone adds a second call to it.
fn median_ns(results: &[BenchmarkResult], field: impl Fn(&BenchmarkResult) -> u64) -> u64 {
    if results.is_empty() {
        return 0;
    }
    let mut values: Vec<f64> = results.iter().map(|r| field(r) as f64).collect();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    percentile(&values, 0.50) as u64
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
    false
}

/// Computes the theoretical in-memory footprint of the loaded index structures.
///
/// This is derived from the actual data sizes, **not** from disk file sizes, so it remains
/// accurate when new structures are added to the `Searcher`. When you add a new structure,
/// extend this function with its memory calculation.
fn theoretical_memory(searcher: &ActiveSearcher, mapping_type: &str, proteins_mapped: bool) -> u64 {
    let text_len = searcher.proteins.text.len() as u64;
    let protein_count = searcher.proteins.proteins.len() as u64;

    // Suffix array: one entry per SA item at bits_per_value bits each
    let sa_bytes = (searcher.sa.len() as u64 * searcher.sa.bits_per_value() as u64).div_ceil(8);

    // ProteinText: 5 bits per character (BitArray), rounded up to whole bytes
    let text_bytes = (text_len * 5).div_ceil(8);

    // Protein metadata
    let string_bytes: u64 = searcher
        .proteins
        .proteins
        .iter()
        .map(|p| p.uniprot_id.len() as u64 + p.functional_annotations.len() as u64)
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
            // One bit per text character + the two-level rank cells (16 bytes per 512 bits)
            let bits_bytes = text_len.div_ceil(8);
            let superblock_count = text_len / 512 + 1;
            let rank_bytes = superblock_count * 16;
            bits_bytes + rank_bytes
        }
        _ => 0
    };

    // No k-mer table exists at this commit, so there is nothing to add for one. The `kmer_k` field
    // on the record stays, pinned to 0, because it is part of the cell key the report joins on.
    sa_bytes + text_bytes + metadata_bytes + mapping_bytes
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
    searcher: &ActiveSearcher,
    peptides: &[String],
    max_matches: usize,
    equate_il: bool,
    tryptic: bool,
    response: bool,
    theoretical_max_memory: u64,
    baseline_memory: u64
) -> BenchmarkResult {
    // Memory snapshot before any timing starts — captures index-resident pages only
    let index_memory = measure_process_memory().saturating_sub(baseline_memory);

    // Fault counters bracket the timed region only, so index loading and warmup are excluded.
    let (minflt_before, majflt_before) = page_faults();

    // Phase 1: suffix array search (parallel). There is no batched orchestrator at this commit —
    // `search_all_matching_suffixes_batched` and the `par_chunks(MLP_BATCH)` interleaving are
    // things the branch added — so this is the decomposition this tree's own production path uses:
    // `search_all_peptides` is a flat `par_iter` over peptides, one independent search each.
    //
    // That is the point rather than a shortcut. Measuring this commit with the branch's batching
    // would attribute the branch's win to something other than the branch.
    let refs: Vec<&[u8]> = peptides.iter().map(|p| p.as_bytes()).collect();
    let search_start = Instant::now();
    let suffix_results: Vec<SearchAllSuffixesResult> = refs
        .par_iter()
        .map(|peptide| searcher.search_matching_suffixes(peptide, max_matches, equate_il, tryptic))
        .collect();
    let search_duration_ns = search_start.elapsed().as_nanos() as u64;

    // Phase 2: protein retrieval — per query via `retrieve_proteins`, which is exactly what
    // production's `search_all_peptides` does. Keeping these two in step is what makes this
    // benchmark measure what ships: an earlier revision called a batched retrieval here while
    // production called the per-query one, which made a whole change invisible to measurement.
    // If production's retrieval shape changes, change it here too.
    //
    // NoMatches queries are dropped before retrieval, exactly as `search_all_peptides` does —
    // there is nothing to look up for them.
    //
    // The peptide index and the cutoff flag are carried alongside because phase 3 needs them to
    // rebuild the `SearchResult` the server actually returns; dropping NoMatches loses the
    // correspondence with `peptides` otherwise.
    //
    // `&Vec<i64>` rather than `&[i64]`: this commit's `retrieve_proteins` takes a `&Vec<i64>`, and
    // going through a slice would mean rebuilding a Vec per query inside the timed region.
    let matched: Vec<(usize, &Vec<i64>, bool)> = suffix_results
        .iter()
        .enumerate()
        .filter_map(|(index, r)| match r {
            SearchAllSuffixesResult::MaxMatches(suf) => Some((index, suf, true)),
            SearchAllSuffixesResult::SearchResult(suf) => Some((index, suf, false)),
            SearchAllSuffixesResult::NoMatches => None
        })
        .collect();
    let matched_suffixes: Vec<&Vec<i64>> = matched.iter().map(|(_, suffixes, _)| *suffixes).collect();

    let retrieval_start = Instant::now();
    let retrieved: Vec<Vec<_>> = matched_suffixes.par_iter().map(|suf| searcher.retrieve_proteins(suf)).collect();
    let retrieval_duration_ns = retrieval_start.elapsed().as_nanos() as u64;

    // Phase 3, only when the cell asks for it: everything production does between "we have the
    // proteins" and "the client has bytes" — `search_all_peptides` builds a `ProteinInfo` per hit
    // (an fa-compression decode plus an accession `String`) and `sa-server` serialises the lot.
    //
    // Measured in production's shape, which is the entire point of the field and what v12 got
    // wrong. This calls the same `json_chunk` / `frame_chunks` the server does, on the same rayon
    // fan-out, so the two cannot drift — the contract phase 2 already carries. v12 instead timed a
    // serial decode and a serial `serde_json::to_vec` against parallel search and retrieval
    // numbers, which overstated the phase by up to the core count, and it serialised a bare
    // `Vec<Vec<ProteinInfo>>` carrying none of the per-peptide `sequence` / `cutoff_used` / object
    // framing the server returns.
    //
    // Opt-in because it is potentially the largest cost in a run: at a 10,000 cutoff it decodes up
    // to that many annotations per peptide. The knob suites are not measuring it and must not pay
    // for it.
    //
    // NOT added to `total_duration_ns`. See the v12 schema note: throughput keeps meaning search
    // plus retrieval, or every baseline in every suite silently changes meaning.
    let (response_duration_ns, response_bytes) = if response {
        let response_start = Instant::now();
        let mut chunks: Vec<Vec<u8>> = retrieved
            .par_iter()
            .zip(matched.par_iter())
            .map(|(proteins, &(index, _, cutoff_used))| {
                json_chunk(&SearchResult {
                    sequence: peptides[index].clone(),
                    proteins: proteins.iter().map(|protein| ProteinInfo::from(*protein)).collect(),
                    cutoff_used
                })
            })
            .collect();
        frame_chunks(&mut chunks);
        let bytes: usize = chunks.iter().map(Vec::len).sum();
        (response_start.elapsed().as_nanos() as u64, bytes as u64)
    } else {
        (0, 0)
    };

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
        response_duration_ns,
        response_bytes,
        major_faults: majflt_after.saturating_sub(majflt_before),
        minor_faults: minflt_after.saturating_sub(minflt_before)
    }
}

// ---------------------------------------------------------------------------
// Matrix mode: grid generation
// ---------------------------------------------------------------------------

/// One line of a `--grid-file`, before its tuning overrides are resolved.
///
/// `deny_unknown_fields` is what makes a stale suite file an error: a driver that writes a key this
/// binary does not know is told so, rather than having it dropped and reporting the shipped
/// configuration under the swept one's name. A suite still carrying `tune` fails here, which is
/// the intended way to find out that axis is gone.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct GridSpec {
    /// Peptide-file stem this cell belongs to. Absent = run it against every `--matrix-files` entry.
    #[serde(default)]
    file: Option<String>,
    /// What the suite calls that file (`mixed`, `small`). Recorded as `peptide_source` so a report
    /// reads in the suite's vocabulary rather than in one machine's file names. Absent = the stem.
    #[serde(default)]
    bucket: Option<String>,
    #[serde(default)]
    kmer_k: usize,
    #[serde(default = "yes")]
    equate_il: bool,
    #[serde(default)]
    tryptic: bool,
    /// Reps and queries for this cell. Absent = the invocation's `--runs` / `--amount-of-peptides`.
    #[serde(default)]
    runs: Option<u32>,
    #[serde(default)]
    amount: Option<usize>,
    /// Match cutoff for this cell. Absent = the invocation's `--max-matches`.
    ///
    /// This is NOT a pure performance parameter: a cutoff that binds truncates the result and sets
    /// `cutoff_used`. A suite that sweeps it is trading answers for time, and its report has to
    /// say so.
    #[serde(default)]
    max_matches: Option<usize>,
    /// Time the two phases production runs after retrieval — the `ProteinInfo` decode and the JSON
    /// serialisation. Off by default because it is expensive and most suites are not measuring it.
    #[serde(default)]
    response: bool,
    /// The suite block this cell came from, and its slot among repeats of the same configuration.
    #[serde(default)]
    sweep: String,
    #[serde(default = "slot_a")]
    grid_slot: String
}

fn yes() -> bool {
    true
}

fn slot_a() -> String {
    "a".to_string()
}

/// One cell of a matrix sweep, fully resolved.
///
/// `--dry-run` and the real sweep both go through this same resolved form, so the planned cell
/// list cannot diverge from the one that runs — which is what makes a dry run worth eyeballing
/// before committing to a multi-hour sweep.
#[derive(Clone, Debug)]
struct GridCell {
    bucket: Option<String>,
    max_matches: Option<usize>,
    response: bool,
    kmer_k: usize,
    equate_il: bool,
    tryptic: bool,
    runs: u32,
    amount: usize,
    sweep: String,
    grid_slot: String
}

impl GridCell {
    /// How this cell reads in the dry run and in the per-cell progress line.
    fn describe(&self) -> String {
        format!(
            "kmer={:<2} il={:<5} tr={:<5}  [{} runs x {} peptides]{}",
            self.kmer_k,
            self.equate_il,
            self.tryptic,
            self.runs,
            self.amount,
            if self.sweep.is_empty() { String::new() } else { format!("  <{}:{}>", self.sweep, self.grid_slot) }
        )
    }
}

/// A resolved grid: cells paired with the peptide file each is restricted to (None = every file).
type Grid = Vec<(Option<String>, GridCell)>;

/// `--grid-file`, or the error explaining that matrix mode has no grid of its own.
fn require_grid_file(args: &Args) -> Result<&Path, Box<dyn Error>> {
    args.grid_file
        .as_deref()
        .ok_or_else(|| "--matrix requires --grid-file: the grid is the driver's, not the harness's".into())
}

/// Reads a `--grid-file` into resolved cells.
///
/// Parsing happens here, up front, rather than per cell during the sweep: a typo in the last
/// line of a grid file should fail before the index is loaded, not four hours into the run.
fn load_grid_file(path: &Path, args: &Args) -> Result<Grid, Box<dyn Error>> {
    let mut cells = Vec::new();

    for (lineno, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let spec: GridSpec =
            serde_json::from_str(&line).map_err(|error| format!("{}:{}: {}", path.display(), lineno + 1, error))?;

        // Refused, not ignored. The k-mer bound table arrived with the branch; a cell asking for one
        // here can only be a suite file copied across without being cut down, and running it anyway
        // would emit a record claiming `kmer_k = 5` for a search that had no table.
        if spec.kmer_k != 0 {
            return Err(format!(
                "{}:{}: this commit has no k-mer bound table, so a cell cannot ask for kmer_k = {}. \
                 Compare against the branch's kmer_k = 0 cells.",
                path.display(),
                lineno + 1,
                spec.kmer_k
            )
            .into());
        }

        cells.push((spec.file.clone(), GridCell {
            bucket: spec.bucket,
            max_matches: spec.max_matches,
            response: spec.response,
            kmer_k: spec.kmer_k,
            equate_il: spec.equate_il,
            tryptic: spec.tryptic,
            runs: spec.runs.unwrap_or(args.runs),
            amount: spec.amount.unwrap_or(args.amount_of_peptides),
            sweep: spec.sweep,
            grid_slot: spec.grid_slot
        }));
    }

    if cells.is_empty() {
        return Err(format!("{} contains no grid cells", path.display()).into());
    }
    Ok(cells)
}

/// The cells to run for one peptide file: those the grid did not restrict to a different one.
fn cells_for(grid: &Grid, file_bucket: &str) -> Vec<GridCell> {
    grid.iter()
        .filter(|(file, _)| file.as_deref().is_none_or(|name| name == file_bucket))
        .map(|(_, cell)| cell.clone())
        .collect()
}

/// Everything `run_cell` needs about one matrix cell that isn't the searcher, the peptides,
/// the CLI args, the cell itself, or the output file: the index-wide facts (fixed for a whole
/// run) and the peptide-file facts (fixed per file).
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
    phase: &'a str
}

/// Runs one config for its rep count, prints the summary line, and appends one aggregated
/// record to `output_file`.
fn run_cell(
    searcher: &ActiveSearcher,
    peptides: &[String],
    args: &Args,
    spec: CellSpec,
    cell: &GridCell,
    output_file: &mut File
) -> Result<(), Box<dyn Error>> {
    // Run reps, then summarise: one record per config with a spread, and the median-qps rep kept
    // as the representative detailed `result`. With `--runs-target-band` the loop stops as soon as
    // the spread is tight enough to read, which is most cells; without it every cell runs the full
    // count, which is what `defaults` needs.
    let mut results: Vec<BenchmarkResult> = Vec::with_capacity(cell.runs as usize);
    while (results.len() as u32) < cell.runs {
        results.push(run_benchmark(
            searcher,
            peptides,
            cell.max_matches.unwrap_or(args.max_matches),
            cell.equate_il,
            cell.tryptic,
            cell.response,
            spec.theoretical_max,
            spec.baseline_memory
        ));
        if let Some(target) = args.runs_target_band {
            if results.len() as u32 >= args.min_runs && band_of(&results) <= target {
                break;
            }
        }
    }

    // A cell with no reps has nothing to summarise, and every statistic below indexes into `qps`.
    // Reachable from a grid cell (or a `--runs` override) of zero, which is a mistake in a suite
    // file rather than a measurement — so it is named as one here, before the sweep has spent
    // hours reaching whichever cell carries it.
    if results.is_empty() {
        return Err(format!(
            "{} {}: runs = 0, so this cell measures nothing — give it a positive rep count",
            spec.source,
            cell.describe()
        )
        .into());
    }

    let mut qps: Vec<f64> = results.iter().map(|r| r.throughput_qps).collect();
    qps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // Each phase is pooled over the same reps as `qps_p50`, and independently of the others: the
    // median search time and the median retrieval time need not come from the same rep, and forcing
    // them to would be picking one rep again under another name. They are read as "what this phase
    // typically cost", not as a decomposition of any single request.
    let stats = RunStats {
        runs: results.len() as u32,
        qps_min: qps[0],
        qps_p10: percentile(&qps, 0.10),
        qps_p50: percentile(&qps, 0.50),
        qps_p90: percentile(&qps, 0.90),
        qps_max: *qps.last().unwrap(),
        search_p50_ns: median_ns(&results, |r| r.search_duration_ns),
        retrieval_p50_ns: median_ns(&results, |r| r.retrieval_duration_ns),
        response_p50_ns: median_ns(&results, |r| r.response_duration_ns)
    };
    let band = if stats.qps_p50 > 0.0 { (stats.qps_p90 - stats.qps_p10) / 2.0 / stats.qps_p50 * 100.0 } else { 0.0 };

    // Representative rep = the one nearest the median throughput.
    results.sort_by(|a, b| a.throughput_qps.partial_cmp(&b.throughput_qps).unwrap());
    let representative = results.remove(results.len() / 2);

    eprintln!(
        "  {} {} {}  ->  {:.0} qps  (±{:.1}%, p10 {:.0} .. p90 {:.0}, {} reps)",
        spec.source,
        spec.phase,
        cell.describe(),
        stats.qps_p50,
        band,
        stats.qps_p10,
        stats.qps_p90,
        stats.runs,
    );

    let record = BenchmarkRecord {
        version: SCHEMA_VERSION,
        label: args.label.clone(),
        commit: spec.commit.to_string(),
        suite: args.suite.clone(),
        dims: args.dims.iter().cloned().collect(),
        startup: spec.startup,
        config: BenchmarkConfig {
            sa_type: spec.sa_type.to_string(),
            mapping_type: spec.mapping_type.to_string(),
            sa_storage: SA_BACKEND,
            text_storage: TEXT_BACKEND,
            proteins_storage: PROTEINS_BACKEND,
            mapping_storage: MAPPING_BACKEND,
            sample_rate: spec.sample_rate,
            bits_per_value: spec.bits_per_value,
            equate_il: cell.equate_il,
            tryptic: cell.tryptic,
            // The cell's cutoff, not the invocation's — a suite that sweeps `max_matches` runs
            // every cell at a different one, and recording the CLI value made all of them look
            // like the same measurement.
            max_matches: cell.max_matches.unwrap_or(args.max_matches),

            kmer_k: cell.kmer_k,
            amount_of_peptides: peptides.len(),
            peptide_length_min: spec.p_min,
            peptide_length_max: spec.p_max,
            peptide_source: spec.source.to_string(),
            phase: spec.phase.to_string(),
            sweep: cell.sweep.clone(),
            grid_slot: cell.grid_slot.clone()
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
    println!("suite          : {}", args.suite);
    if !args.dims.is_empty() {
        let dims: Vec<String> = args.dims.iter().map(|(k, v)| format!("{}={}", k, v)).collect();
        println!("dims           : {}", dims.join(" "));
    }
    // Same load-and-resolve the real run does, so a typo in the grid file surfaces here rather
    // than after the index is mapped.
    let path = require_grid_file(args)?;
    let grid = load_grid_file(path, args)?;

    println!("backend        : {}", backend_summary());
    println!("grid           : {}", path.display());
    println!("runs/config    : {}{}", args.runs, match args.runs_target_band {
        Some(target) => format!(" (max; stops at a ±{:.1}% band after {} reps)", target, args.min_runs),
        None => " (fixed)".to_string()
    });
    println!();

    let mut grand_total = 0usize;
    let mut grand_queries = 0usize;

    for pep_path in &args.matrix_files {
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?");
        // Same expansion the real run uses — a dry run cannot diverge from it.
        let cells = cells_for(&grid, source);
        println!("== {} ==", source);
        println!("  grid: {} configs", cells.len());
        for cell in &cells {
            println!("    {}", cell.describe());
            grand_queries += cell.runs as usize * cell.amount;
        }
        grand_total += cells.len();
        println!();
    }

    // Queries, not configs, is what the wall clock tracks once cells may differ in size: a tryptic
    // cell at a tenth the query count is a tenth of the cost, and a config count hides that.
    println!("TOTAL this backend: {} configs, {} timed queries", grand_total, grand_queries);
    println!("This commit has one backend, so this is the whole sweep.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Warmup (never timed)
// ---------------------------------------------------------------------------

/// Peptides per warmup batch — bounds the memory held for the peptide strings themselves
/// while still giving rayon a big enough chunk to parallelise over.
const WARMUP_BATCH_SIZE: usize = 100_000;

/// What one structure's page sweep did: bytes touched and how long it took.
///
/// A structure that is preloaded reports `(0, ~0)` because its sweep is the trait default, and
/// that zero is the useful part — it says the arm skipped the section rather than warming it
/// quickly.
#[derive(Clone, Copy, Default)]
struct SweepStat {
    bytes: u64,
    ms: u64
}

/// The page sweep, which has nothing to sweep at this commit.
///
/// Every structure is owned memory here — there are no mapped regions, so there are no pages to
/// fault in ahead of time. Reporting `(0, 0)` per structure is exactly what the branch's fully
/// preloaded arm reports (its `touch_all_pages` is the trait default), so the `warmup_*_bytes`
/// fields stay comparable: a zero means "this arm had nothing mapped", in both trees.
///
/// Kept as a function rather than inlined so the startup record is assembled by the same code path
/// in both trees, and so what warms this build is unambiguous: only `warmup_pipeline` does.
fn warmup_touch_pages(_searcher: &ActiveSearcher) -> (SweepStat, SweepStat, SweepStat) {
    (SweepStat::default(), SweepStat::default(), SweepStat::default())
}

/// Folds a page sweep's per-structure results into the record, and says what it found.
fn record_sweep(startup: &mut StartupTiming, sweeps: (SweepStat, SweepStat, SweepStat)) {
    let (sa, proteins, mapping) = sweeps;
    startup.warmup_sa_bytes = sa.bytes;
    startup.warmup_proteins_bytes = proteins.bytes;
    startup.warmup_mapping_bytes = mapping.bytes;
    startup.warmup_sa_ms = sa.ms;
    startup.warmup_proteins_ms = proteins.ms;
    startup.warmup_mapping_ms = mapping.ms;
    for (what, sweep) in [("sa", sa), ("proteins", proteins), ("mapping", mapping)] {
        eprintln!(
            "  swept {}: {:.2} GB in {} ms ({})",
            what,
            sweep.bytes as f64 / 2f64.powi(30),
            sweep.ms,
            rate(sweep.bytes, sweep.ms)
        );
    }
}

/// Bytes per second as a human-readable rate, or `-` when there was nothing to sweep.
fn rate(bytes: u64, ms: u64) -> String {
    if bytes == 0 {
        return "nothing mapped".to_string();
    }
    if ms == 0 {
        return "instant".to_string();
    }
    format!("{:.2} GB/s", bytes as f64 / 2f64.powi(30) / (ms as f64 / 1000.0))
}

/// Pushes `count` peptides from `<index-dir>/warmup.txt` through the full search + retrieval
/// pipeline, in batches, discarding the results. Stops early if the file runs out.
fn warmup_pipeline(searcher: &ActiveSearcher, args: &Args, count: usize) -> Result<(), Box<dyn Error>> {
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

/// Matrix mode: load the index once, then sweep the driver's `--grid-file` for each peptide file.
/// The k-mer tables are swapped in/out (moved, not cloned) so each is built or loaded at most once,
/// and only the sizes the grid actually names are touched at all. Each config runs its own rep
/// count but writes a single aggregated record (median-qps rep as `result`, plus a `stats` spread)
/// to `<output>/<label>.jsonl`.
// Index-wide facts the caller already computed while loading; they are threaded through rather
// than recomputed here, which is what pushes this past the argument limit.
#[allow(clippy::too_many_arguments)]
fn run_matrix(
    searcher: ActiveSearcher,
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
    let grid = load_grid_file(require_grid_file(args)?, args)?;

    // Warm the page cache once (matters for mmap); CPU caches warm over the repeated runs.
    eprintln!("Warming up (touching all pages)...");
    let (minflt_before_warmup, majflt_before_warmup) = page_faults();
    let warmup_start = Instant::now();
    let sweeps = warmup_touch_pages(&searcher);
    record_sweep(&mut startup, sweeps);

    // Then the pipeline half, if the suite asked for one. Matrix mode used to ignore `--warmup`
    // entirely and page-sweep only, which is not the same warmup the single-mode suites get — and
    // the difference falls on one arm. A page sweep touches mapped pages; a build that PRELOADS a
    // structure has it in anonymous memory, where no sweep reaches it, so the only thing that warms
    // it is running real queries. `ram` and `threads` do that (`warmup = "all:1000000"`) and put
    // `pprot` +14-16% over `mmap` on the mixed file; the matrix suites did not and put the same
    // comparison on the same file at +0.3-9.6%. Whether that is the whole explanation is what the
    // re-run this was added for will say, but the two modes could not be compared until they warmed
    // the same way.
    //
    // The pipeline warms at the top-level `--equate-il` / `--tryptic` / `--max-matches`, not at each
    // cell's: a grid has many cells and one warmup, so there is no per-cell value to use. That is
    // fine for warming — the structures touched are the same whatever the search options say.
    if let Some(WarmupMode::Count(count) | WarmupMode::AllThenCount(count)) = &args.warmup {
        warmup_pipeline(&searcher, args, *count)?;
    }

    startup.warmup_ms = warmup_start.elapsed().as_millis() as u64;
    let (minflt_after_warmup, majflt_after_warmup) = page_faults();
    startup.warmup_minor_faults = minflt_after_warmup.saturating_sub(minflt_before_warmup);
    startup.warmup_major_faults = majflt_after_warmup.saturating_sub(majflt_before_warmup);
    startup.load_total_ms =
        startup.load_sa_ms + startup.load_proteins_ms + startup.load_mapping_ms + startup.kmer_table_ms;
    eprintln!(
        "Startup: load {} ms (sa {} / proteins {} / mapping {}), warmup {} ms",
        startup.load_total_ms, startup.load_sa_ms, startup.load_proteins_ms, startup.load_mapping_ms, startup.warmup_ms
    );
    eprintln!(
        "         warmup swept {:.2} GB overall ({}), {} major faults",
        startup.warmup_bytes() as f64 / 2f64.powi(30),
        rate(startup.warmup_bytes(), startup.warmup_ms),
        startup.warmup_major_faults
    );

    // Index-wide and independent of the peptide file, so it is computed once up front rather than
    // per cell — it walks every protein's metadata length, which is not free at full-DB scale.
    // The branch computes one of these per k-mer size; with no table there is a single value.
    let theoretical_max =
        if args.no_theoretical_memory { 0 } else { theoretical_memory(&searcher, mapping_type, proteins_mapped()) };

    create_dir_all(&args.output)?;
    let output_path = args.output.join(format!("{}.jsonl", args.label));
    let mut output_file = OpenOptions::new().create(true).write(true).truncate(true).open(&output_path)?;
    let commit = env!("GIT_COMMIT_HASH").to_string();

    for pep_path in &args.matrix_files {
        let source = pep_path.file_stem().and_then(|s| s.to_str()).unwrap_or("?").to_string();
        let cells = cells_for(&grid, &source);
        if cells.is_empty() {
            eprintln!("-- {} : no cells in the grid, skipping --", source);
            continue;
        }

        // Cells may ask for different query counts, so read the largest once and let each cell
        // take the prefix it needs. Every cell then queries the same peptides in the same order,
        // which is what keeps a cheap screening cell comparable with the cells around it.
        let wanted = cells.iter().map(|cell| cell.amount).max().unwrap_or(args.amount_of_peptides);
        let peptides: Vec<String> =
            BufReader::new(File::open(pep_path)?).lines().take(wanted).collect::<Result<_, _>>()?;
        if peptides.len() < wanted {
            return Err(format!("{} has only {} peptides, need {}", pep_path.display(), peptides.len(), wanted).into());
        }
        let bucket = cells.iter().find_map(|cell| cell.bucket.clone()).unwrap_or_else(|| source.clone());

        // Everything the cells of this file share; the per-cell fields are filled in below.
        let base_spec = CellSpec {
            startup,
            mapping_type,
            sa_type,
            sample_rate,
            bits_per_value,
            baseline_memory,
            commit: &commit,
            // What the suite calls this file, not what the filesystem does. Every cell of one file
            // carries the same bucket, so the first one that names it settles it for the file.
            source: &bucket,
            p_min: 0,
            p_max: 0,
            theoretical_max: 0,
            phase: "grid"
        };

        eprintln!("-- {} : {} cells --", bucket, cells.len());
        for cell in &cells {
            // Per cell, not per file: cells may take different prefixes of the same stream, and a
            // length range describing lines the cell never queried would misreport the workload.
            let queried = &peptides[..cell.amount];
            let (p_min, p_max) =
                queried.iter().fold((usize::MAX, 0usize), |(lo, hi), p| (lo.min(p.len()), hi.max(p.len())));
            let spec = CellSpec {
                theoretical_max,
                p_min,
                p_max,
                ..base_spec
            };
            run_cell(&searcher, queried, args, spec, cell, &mut output_file)?;
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
    // and prints the config list the matrix sweep would run.
    if args.dry_run {
        if !args.matrix {
            return Err("--dry-run only applies to --matrix mode".into());
        }
        return print_dry_run(&args);
    }

    let sa_path = args.index_dir.join("sa.bin");
    // `proteins.tsv`, not `proteins.bin`: the prebuilt protein store and suffix-to-protein mapping
    // are artefacts the branch introduced. At this commit the server reads the same database TSV
    // the index was built from and reconstructs both at startup, so that is what is measured here.
    // The branch's index directory still ships the TSV beside the .bin files, which is what lets
    // both trees run off one index with nothing rebuilt.
    let proteins_path = args.index_dir.join("proteins.tsv");
    let mapping_path = args.index_dir.join("mapping.bin");

    // Detect mapping type before loading (peek at the type byte). `mapping.bin` is not loaded here
    // — nothing at this commit can read it — but its type byte says which mapping the branch's runs
    // walked, and this run rebuilds the same one so the two are comparable. That is the only use
    // this file has here, and if it is absent the run has no way to know what to match.
    let mapping_type_str = match first_byte_of(&mapping_path)? {
        0 => "dense",
        1 => "sparse",
        2 => "bitvec",
        other => return Err(format!("{} declares unknown mapping type {}", mapping_path.display(), other).into())
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
        load_total_ms: 0,
        warmup_sa_bytes: 0,
        warmup_proteins_bytes: 0,
        warmup_mapping_bytes: 0,
        warmup_sa_ms: 0,
        warmup_proteins_ms: 0,
        warmup_mapping_ms: 0,
        load_major_faults: 0,
        load_minor_faults: 0,
        warmup_major_faults: 0,
        warmup_minor_faults: 0
    };
    let (minflt_at_start, majflt_at_start) = page_faults();

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

    // NOT comparable with the branch's `load_proteins_ms`, and deliberately recorded anyway. The
    // branch reads a prebuilt `proteins.bin`; this parses a ~293 MB TSV, uppercases and concatenates
    // every sequence, fa-compresses every annotation and bit-packs the text. It is a different
    // operation, not a slower one — read the two as "what each tree pays before its first query",
    // never as a ratio.
    eprintln!("Building proteins from {}...", proteins_path.display());
    let t0 = Instant::now();
    let proteins = Proteins::try_from_database_file(proteins_path.to_str().unwrap())?;
    startup.load_proteins_ms = t0.elapsed().as_millis() as u64;
    eprintln!("  {} proteins, {} text characters ({} ms)", proteins.proteins.len(), proteins.text.len(), startup.load_proteins_ms);

    let sa_type = if suffix_array.bits_per_value() == 64 { "original" } else { "compressed" };
    let sample_rate = suffix_array.sample_rate();
    let bits_per_value = suffix_array.bits_per_value();

    // Same caveat as the proteins: the branch loads `mapping.bin`, this builds the same mapping from
    // the text. Charged to `load_mapping_ms` because it is the same structure arriving by a
    // different route, which is exactly what a startup comparison should show.
    eprintln!("Building the {} mapping...", mapping_type_str);
    let t0 = Instant::now();
    let searcher = ActiveSearcher::new(mapping_type_str, suffix_array, proteins)?;
    startup.load_mapping_ms = t0.elapsed().as_millis() as u64;
    eprintln!("  {} ms", startup.load_mapping_ms);

    // The three index loads, charged before the paths diverge: matrix mode loads its k-mer tables
    // inside `run_matrix` and would otherwise report no load faults at all. Single mode widens this
    // to include the k-mer table below.
    let (minflt_after_index, majflt_after_index) = page_faults();
    startup.load_minor_faults = minflt_after_index.saturating_sub(minflt_at_start);
    startup.load_major_faults = majflt_after_index.saturating_sub(majflt_at_start);

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

    // `kmer_table_ms` stays in the record and stays 0: there is no table at this commit.
    startup.load_total_ms =
        startup.load_sa_ms + startup.load_proteins_ms + startup.load_mapping_ms + startup.kmer_table_ms;
    let (minflt_after_load, majflt_after_load) = page_faults();
    startup.load_minor_faults = minflt_after_load.saturating_sub(minflt_at_start);
    startup.load_major_faults = majflt_after_load.saturating_sub(majflt_at_start);

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
            let sweeps = warmup_touch_pages(&searcher);
            record_sweep(&mut startup, sweeps);
            eprintln!("Warmup complete.");
        }
        Some(WarmupMode::Count(warmup_count)) => {
            warmup_pipeline(&searcher, &args, *warmup_count)?;
        }
        Some(WarmupMode::AllThenCount(warmup_count)) => {
            eprintln!("Warming up: touching all mmap pages...");
            let sweeps = warmup_touch_pages(&searcher);
            record_sweep(&mut startup, sweeps);
            eprintln!("Page warmup complete.");
            warmup_pipeline(&searcher, &args, *warmup_count)?;
        }
    }
    if args.warmup.is_some() {
        startup.warmup_ms = warmup_start.elapsed().as_millis() as u64;
        let (minflt_after_warmup, majflt_after_warmup) = page_faults();
        startup.warmup_minor_faults = minflt_after_warmup.saturating_sub(minflt_after_load);
        startup.warmup_major_faults = majflt_after_warmup.saturating_sub(majflt_after_load);
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
    eprintln!(
        "         warmup swept {:.2} GB overall ({}), {} major / {} minor faults; load took {} major faults",
        startup.warmup_bytes() as f64 / 2f64.powi(30),
        rate(startup.warmup_bytes(), startup.warmup_ms),
        startup.warmup_major_faults,
        startup.warmup_minor_faults,
        startup.load_major_faults
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
            sa_storage: SA_BACKEND,
            text_storage: TEXT_BACKEND,
            proteins_storage: PROTEINS_BACKEND,
            mapping_storage: MAPPING_BACKEND,
            sample_rate,
            bits_per_value,
            equate_il: args.equate_il,
            tryptic: args.tryptic,
            max_matches: args.max_matches,

            kmer_k: 0,
            amount_of_peptides: run_peptides.len(),
            peptide_length_min: p_min,
            peptide_length_max: p_max,
            peptide_source: peptide_source.clone(),
            phase: "single".to_string(),
            // A single run is one configuration measured once: it belongs to no suite block, and
            // there is nothing for it to be a repeat of.
            sweep: String::new(),
            grid_slot: "a".to_string()
        };

        let result = run_benchmark(
            &searcher,
            &run_peptides,
            args.max_matches,
            args.equate_il,
            args.tryptic,
            // A single CLI run has no grid to opt in with; `--matrix` is where the response phase
            // lives, because only a suite can decide it is worth what it costs.
            false,
            theoretical_max,
            baseline_memory
        );

        let record = BenchmarkRecord {
            version: SCHEMA_VERSION,
            label: args.label.clone(),
            commit: env!("GIT_COMMIT_HASH").to_string(),
            suite: args.suite.clone(),
            dims: args.dims.iter().cloned().collect(),
            config,
            startup,
            result,
            stats: None
        };

        let line = serde_json::to_string(&record)?;
        writeln!(output_file, "{}", line)?;

        eprintln!(
            "Run {:>3}/{}: {:.1} qps  |  total {:.1} ms  (search {:.1} ms + retrieval {:.1} ms)  \
             |  hits: {} queries, {} suffixes, {} proteins{}{}",
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
