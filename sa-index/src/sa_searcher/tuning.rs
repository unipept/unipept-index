//! The constants the search path is fixed at, and what measured them.
//!
//! None of these changes an answer, only how long it takes to produce one.
//!
//! # Why they are constants and not knobs
//!
//! They were four fields of a `SearchTuning` struct, settable at runtime and swept by the benchmark
//! harness. A full-database sweep — eleven suites across five storage arms, against a measured
//! noise floor — could not distinguish three of them from noise anywhere it looked. The run itself
//! is in `sa-index/BENCHMARKS.md`; what it found:
//!
//! * `validate_batch`, swept 16..256 over 40 contexts — **0 of 40** cleared their own floor. The
//!   largest gain any value showed was +5.6%, against that context's own floor of ±11.2%; the
//!   floors across the 40 ran ±3.9% to ±19.0%, and on 14 of them the shipped 64 was itself the
//!   peak.
//! * `prefetch_threshold` × `retrieval_prefetch_distance`, swept as a full 4×4 cross on five arms
//!   — **0 of 80 pairs** cleared their own floor, and the best pair differed on every arm
//!   (64/32, 32/8, 64/16, 32/32, 8/16), which is what noise looks like rather than a plateau.
//! * `combos`, crossing all three accelerators, found **0 of 10** — so nothing was hiding in an
//!   interaction between them either.
//!
//! The fourth, `mlp_batch`, was the one that moved: 64 beat the shipped 16 in 8 of 40 contexts,
//! every one of them on the long-peptide file. It stays at 16 anyway, and the reasoning is under
//! [`MLP_BATCH`] — no value wins everywhere, so there was nothing for a runtime knob to select
//! that a constant could not.
//!
//! Two knobs before all of them (`retrieval_batch`, `scalar_kmer_prefetch`) were removed the same
//! way. The values below are exactly the defaults those sweeps ran against, so fixing them changed
//! no measurement.
//!
//! # What this evidence does not license
//!
//! Deleting the machinery. Every swept value leaves two-pass validation on for ranges above 64 and
//! off for ranges below 8, so the sweeps priced the *crossover* and never the mechanism. Turning
//! the two-pass paths off entirely is an unmeasured change, not a corollary of this one.
//!
//! # Re-opening one of these
//!
//! There is no longer a runtime path to sweep, so re-tuning means restoring one: give the searcher
//! a field, thread it to the constant's use, and teach the benchmark harness to vary it. That is
//! deliberately more work than editing a number here — a knob that no measurement can move is a
//! knob that costs more to carry than to rebuild.

/// Peptides interleaved per rayon task in `search_all_matching_suffixes_batched`, for cross-query
/// memory-level parallelism.
///
/// A suffix-array probe is a random read whose address depends on the previous probe's result, so
/// one search is a dependent chain of cache misses with nothing to overlap. Batching B searches per
/// task hands the memory system B independent chains at once, and the win is however much of that
/// latency the hardware can then hide.
///
/// 16 is the value that hurts nothing, which is not the same as the fastest. The full-database
/// sweep found the curve inverts with peptide length — larger batches help the search-bound long
/// peptides and hurt the retrieval-bound short and mixed ones — so no single value wins everywhere:
///
/// | batch | resolved gains | resolved regressions | median vs 16 |
/// |---|---|---|---|
/// | scalar | 0 | 13 | -2.1% |
/// | 4 | 0 | 6 | -0.9% |
/// | 8 | 0 | 0 | -0.5% |
/// | **16** | — | — | — |
/// | 32 | 1 | 0 | +1.7% |
/// | 64 | 8 | 1 | +1.1% |
/// | 128 | 4 | 4 | -2.1% |
///
/// 64's eight wins are all on the long-peptide file, at +8.8% to +13.5%. It pays for them on
/// `mixed` — the 5..50 mix a server actually sees — where it reads -5.6% to -7.7% and the mapped
/// arm's cell resolves as a regression. Anything below 8 is a resolved loss on long peptides.
///
/// Revisit this only with a per-request or per-length decision; a single global value cannot take
/// the long-peptide gain without paying on the mix.
pub(crate) const MLP_BATCH: usize = 16;

/// Candidates per two-pass validation batch in `iterate_sa_range` and
/// `iterate_extended_sa_range`, and the size of the on-stack buffer that holds them.
///
/// 64 `i64`s is 512 bytes of stack per call. Tuned on x86_64 Zen4/Intel Sapphire Rapids: DRAM
/// latency ~80–100 ns, one SA entry read per ~2–3 ns at that cache level → 64 entries ≈ 192 ns
/// gap, comfortably above the latency floor.
///
/// Fixed rather than tunable because the sweep could not separate any value from any other: see
/// the module docs. One thing the sweep *did* show consistently is the direction below this value
/// — `validate_batch = 16` read −3.2% to −9.0% in all ten mixed cells. None of those resolved, but
/// the sign never flipped, so this is a floor to stay above rather than a free parameter.
pub(crate) const VALIDATE_BATCH: usize = 64;

/// Minimum SA range size, in entries, before the two-pass validation paths run instead of a
/// straight loop. Below it the two-pass overhead exceeds the latency it hides.
///
/// This is the gate on the *candidate-validation* path and the partner of [`VALIDATE_BATCH`]:
/// that one sets the batch size, this one decides whether a range is big enough to batch at all.
/// It does not reach retrieval, which prefetches unconditionally at
/// [`RETRIEVAL_PREFETCH_DISTANCE`].
pub(crate) const PREFETCH_THRESHOLD: usize = 32;

/// Prefetch look-ahead distance, in suffixes, inside protein retrieval.
///
/// D/2 iterations × ~5 ns ≈ 80–100 ns before the protein read in `proteins.get()`, giving the
/// hint time to complete on most DRAM configurations.
///
/// Worth knowing when reading a retrieval measurement: a query matching fewer suffixes than this
/// distance issues no prefetches at all. Tryptic queries match almost nothing and are mostly
/// already in that case.
pub(crate) const RETRIEVAL_PREFETCH_DISTANCE: usize = 32;

/// Cap on how much result space is pre-allocated per peptide, in suffixes.
///
/// Callers pass a `max_matches` cutoff that is an upper bound, not an estimate — the server's
/// default is 10 000 — while the overwhelming majority of peptides match a handful of times.
/// Reserving the full cutoff would allocate 80 KB per peptide and touch none of it; capping at
/// 4096 entries (32 KB) keeps the common case to one allocation while letting rare high-hit
/// peptides grow normally.
pub(crate) const MAX_RESULT_PREALLOC: usize = 4096;
