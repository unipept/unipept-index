//! Suffix-array search over the concatenated protein text.
//!
//! Given a peptide, find every protein containing it. The pipeline is:
//!
//! 1. [`sa_searcher`] binary-searches the suffix array for the range of suffixes sharing the
//!    peptide as a prefix, optionally starting from bounds looked up in a [`kmer_table`];
//! 2. it validates each candidate in that range against the text, since a sparse suffix array
//!    only indexes every n-th position and I/L equating and tryptic filtering add further
//!    conditions;
//! 3. [`suffix_to_protein_index`] maps surviving text positions to protein indices;
//! 4. [`peptide_search`] turns those into results with accessions and annotations.
//!
//! # Storage features: one build, one API, four independent choices
//!
//! Every storage structure has two implementations — one holding owned memory, one borrowing a
//! memory mapping — and features pick which. The selection is a *type alias*, resolved at compile
//! time, so there is no runtime branch and no dispatch anywhere in the search path.
//!
//! `mmap` maps everything. Each `preloaded-*` feature then pulls **one** structure back into owned
//! memory, leaving the rest mapped:
//!
//! | alias | owned | mapped | mapped when |
//! |---|---|---|---|
//! | [`SuffixArray`] | `InMemorySA` | `MmapBackedSA` | `mmap` |
//! | `text_compression::ProteinText` | `InMemoryProteinText` | `MmapBackedProteinText` | `mmap` and not `preloaded-text` |
//! | `sa_mappings::proteins::Proteins` | `InMemoryProteins<T>` | `MmapBackedProteins<T>` | `mmap` and not `preloaded-proteins` |
//! | [`suffix_to_protein_index::SuffixToProteinMapping`] | `InMemorySuffixToProteinMapping` | `MmapBackedSuffixToProteinMapping` | `mmap` and not `preloaded-mapping` |
//!
//! Nine configurations in all: everything preloaded, everything mapped, and the seven mixtures.
//! The point is that the best place for one structure is not the best place for another — the
//! text is the hottest and the metadata table the biggest — so, for instance
//! `--features mmap,preloaded-text` keeps the multi-gigabyte index mapped while the ~190 MB text
//! the search reads once per character compared sits in owned RAM.
//!
//! Consequences worth knowing before reading further:
//!
//! * **No crate declares `default = [...]`**, so a plain `cargo build` or `cargo test` gives you
//!   the fully *preloaded* configuration. The production server is built `--features mmap`.
//! * The suffix array follows `mmap` and has no override; there is no `preloaded-sa`.
//! * A `preloaded-*` feature without `mmap` is a no-op — everything is preloaded already. Cargo
//!   features are additive and cannot be negated by a dependent crate, so they only ever *remove*
//!   mapping, never add it.
//! * The protein text and the protein metadata share one file (`proteins.bin`) but are separate
//!   axes, which is why `Proteins` is generic over its text type. All four pairings load from the
//!   same file; see `sa_mappings::proteins::mmap`.
//! * The preloaded half is compiled in **every** configuration, because it owns the `WriteBinary`
//!   implementations that produce the files the mmap half reads.
//! * The two halves are kept deliberately separate rather than sharing code, so that a tuning
//!   change to one cannot perturb the other. Where that produces near-duplicate code it is marked
//!   as intentional at the site.
//!
//! # The `metrics` feature
//!
//! Off by default. It swaps the counters on [`sa_searcher::Searcher`] from zero-sized no-ops to
//! real atomics, so the benchmark harness can attribute time and count candidates. Enabling it
//! costs throughput; see `sa_searcher::metrics`.
//!
//! # Why this crate is written the way it is
//!
//! Search is dominated by DRAM latency, not by instruction count: the suffix array and the
//! suffix-to-protein mapping are far larger than any cache and are walked in an order the
//! hardware prefetcher cannot predict. Most of the non-obvious code here — two-pass batching,
//! software prefetch hints, cross-query batching, the k-mer table — exists to overlap those
//! misses rather than to do less work.
//!
//! One consequence shows up everywhere and is easy to undo by accident: because the workspace
//! sets no `[profile.release]`, there is **no cross-crate LTO**, so a call into `bitarray`,
//! `text-compression`, `sa-mappings` or `prefetch` is a real cross-crate call unless the callee's
//! body reaches the caller's codegen unit — which happens only when the callee is `#[inline]` or
//! generic (both export their MIR). Those attributes on the small getters are load-bearing, not
//! decoration.
//!
//! **Enabling LTO was measured and rejected** (2026-08-09, at c00cc53, full UniProt index on the
//! benchmark server; n=200 timed reps per cell, ABBA-interleaved, both backends). Adding
//! `[profile.release] lto = "thin", codegen-units = 1` moved median throughput by **-0.1% on
//! mmap and -3.0% on preloaded** — no gain on either backend, in exchange for a materially
//! slower clean release build.
//!
//! Neither delta is an effect. On mmap the two arms interleave and the base arm's own two
//! invocations spread wider (1.7%) than the difference between arms. On preloaded the base
//! arm's two invocations differ by **8.4%** — far more than the 3.0% gap — and all four
//! invocations decline monotonically with position, i.e. the machine drifted over the ~50
//! minutes the block ran. That drift decays rather than being linear, which the ABBA ordering
//! does not cancel: base holds the two endpoint slots and LTO the two middle ones, and for a
//! convex-decaying curve the endpoint average is the higher one. The preloaded arm therefore
//! cannot resolve anything below roughly 8%, and there is no evidence of a real regression
//! either.
//!
//! The likely reason for the null result is that this crate has already hand-annotated the hot
//! path, so the bodies LTO would have inlined were reaching the caller anyway. Do not re-add the
//! profile block without re-measuring on the full database.
//!
//! Decisions that were measured and *rejected* are recorded next to the code they would have
//! touched, so they are not rediscovered and retried.
//!
//! # When the index does not fit in RAM
//!
//! Everything above was tuned with the whole index resident. That is not the regime the mmap
//! backend exists for, and **several of those decisions do not survive the regime it does exist
//! for.** Measured 2026-08-10/11 at 3259427: full UniProt index (223 GB total) on a 295 GB /
//! 12-core server, ceilings imposed with cgroup v2 `MemoryMax` and swap off, page cache dropped
//! before every cell, 40-100 timed reps per cell, 6-mer table attached unless stated.
//!
//! For scale, the index divides roughly as SA ~149 GB, text ~40 GB, protein metadata ~24 GB,
//! mapping ~10 GB — derived from the entry widths and sample rate rather than measured, except
//! the metadata, which is the RSS delta observed when `preloaded-proteins` moves it to the heap.
//! The ratios are what matter: the SA is two thirds of the working set, so it dominates residency.
//!
//! **The degradation is a concurrency limit, not a bandwidth one.** A major fault blocks its
//! thread, and `prefetch_read` cannot help: a CPU hint instruction cannot fault, so every
//! prefetch in this crate is inert against an absent page. With rayon at the core count, each
//! faulting thread idles a core. Raising `RAYON_NUM_THREADS` leaves the fault *count* unchanged
//! to within 0.24% and still buys:
//!
//! | ceiling | major faults/rep | default threads | best | gain |
//! |---|---|---|---|---|
//! | none | 0 | 35,710 | 35,046 @ 48 | **-1.9%** |
//! | 167 GB (75%) | 24,190 | 15,739 | 26,071 @ 48 | **+65.6%** |
//! | 112 GB (50%) | 46,690 | 10,561 | 19,654 @ 96 | **+86.1%** |
//!
//! The gain tracks the fault rate and the cost is real when resident (-7.8% at 96 threads), so
//! this is a deployment knob, not a default. It is also the largest single effect anywhere in
//! this investigation — larger than any storage-backend choice.
//!
//! **A 6-mer k-mer table is worth its 3.06 GB here**, against the note in `sa-benchmarks` that
//! measured it inside the noise floor when resident. At a 167 GB ceiling it is +18.4% and -27.9%
//! faults versus no table; a 5-mer table is +3.2% and -6.2%, i.e. barely distinguishable from
//! nothing. The difference is working-set size, not probe count: a 5-mer narrows the search to
//! ~7 SA pages per query, a 6-mer to ~1.
//!
//! **All of the loss is in the search phase.** Retrieval is flat at ~147 ms per rep across every
//! ceiling while search goes 135 ms → 1127 ms, so the two-pass prefetch pipeline in
//! `sa_searcher::retrieval` keeps working under paging and needs nothing. Within search, the split
//! is roughly even between the dependent binary-search chain and the contiguous SA range scan
//! (52% / 48% of thread-time, `metrics` build at a 167 GB ceiling).
//!
//! Two further ideas were measured and rejected; see `array::mmap::MmapBackedSA::advise_willneed_range`
//! and `sa_searcher::SearchTuning::willneed` for `MADV_WILLNEED`, and note here that **sorting
//! queries by k-mer prefix to create page locality does not work**: it changed the fault count by
//! -0.1% and cost 4.4% throughput. With 10,000 queries per rep drawn against 20^6 possible
//! 6-mers, the expected number of queries sharing a prefix is under one, so there is no page reuse
//! for sorting to expose. Locality needs reuse, and this workload has none.

pub use text_compression::{ReadBinary, ReadBinaryMmap, WriteBinary};

pub mod array;
pub mod kmer_table;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

pub use array::{SuffixArray, SuffixArrayBackend};
pub use kmer_table::KmerTable;
pub use sa_mappings::proteins::ProteinsBackend;
pub use sa_searcher::{MAX_VALIDATE_BATCH, SearchTuning};

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    /// Returns whether the value is NULL.
    ///
    /// # Returns
    ///
    /// True if the value is NULL, false otherwise.
    fn is_null(&self) -> bool;
}

/// Implementation of the `Nullable` trait for the `u32` type.
impl Nullable<u32> for u32 {
    const NULL: u32 = u32::MAX;

    fn is_null(&self) -> bool {
        *self == Self::NULL
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nullable_is_null() {
        assert!(u32::NULL.is_null());
        assert!(!0u32.is_null());
    }
}
