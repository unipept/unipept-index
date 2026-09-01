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
//! # Storage: two backends per structure, chosen by the caller
//!
//! Every storage structure has two implementations — one holding owned memory, one borrowing a
//! memory mapping — and **this crate has no opinion about which**. Both are always compiled, the
//! searcher is generic over all three of them, and nothing here names a concrete one:
//!
//! | structure | owned | mapped |
//! |---|---|---|
//! | suffix array | [`array::InMemorySA`] | [`array::MmapBackedSA`] |
//! | protein text | `protein_text::InMemoryProteinText` | `protein_text::MmapBackedProteinText` |
//! | protein metadata | `protein_metadata::InMemoryProteins<T>` | `protein_metadata::MmapBackedProteins<T>` |
//! | suffix→protein | [`suffix_to_protein_index::InMemorySuffixToProteinMapping`] | [`suffix_to_protein_index::MmapBackedSuffixToProteinMapping`] |
//!
//! The choice is made once per build, by the binary: `sa-server`'s `backends` module resolves four
//! Cargo features into one concrete type per structure. That is the *only* place in the workspace
//! a storage feature is read. Selection is by type, so there is still no runtime branch and no
//! dispatch anywhere in the search path.
//!
//! Sixteen combinations are constructible, of which the binaries expose nine. The point is that
//! the best place for one structure is not the best place for another: the text is read once per
//! character compared and is the hottest thing in the index, while the metadata table is read once
//! per reported result and is the one that grows most when preloaded — roughly tripling, since it
//! becomes a `Vec` of owned strings rather than bytes in a file. So, for instance,
//! `--features mmap,preloaded-text` keeps the index mapped while the text sits in owned RAM. The
//! text is the cheapest structure to preload — 5 bits per residue, roughly a fifth of the index —
//! and the hottest, which is what makes that pairing worth having.
//!
//! Sizes in these docs are given relative to the index, not in bytes, because the two databases
//! this crate is run against differ by more than two orders of magnitude. Absolute figures for a
//! specific run are in [`BENCHMARKS.md`](../BENCHMARKS.md), which names its index and host.
//!
//! Consequences worth knowing before reading further:
//!
//! * The protein text and the protein metadata share one file (`proteins.bin`) but are separate
//!   axes, which is why the protein structs are generic over their text type. All four pairings
//!   load from the same file; see `protein_metadata::mmap`.
//! * Which reader a structure needs is not a decision any caller makes — it is the `LoadIndex`
//!   implementation on that concrete type. That is what lets a test build all sixteen combinations
//!   without a single `#[cfg]`; see `sa_searcher::tests::every_backend_combination_returns_identical_results`,
//!   which asserts they answer identically.
//! * The owned half owns the `WriteBinary` implementations that produce the files the mapped half
//!   reads, which is why `sa-builder` never mentions a backend.
//! * The two halves are kept deliberately separate rather than sharing code, so that a tuning
//!   change to one cannot perturb the other. Where that produces near-duplicate code it is marked
//!   as intentional at the site.
//!
//! # Measurement code
//!
//! **This crate declares no Cargo features at all and carries no instrumentation**: nothing in the
//! search path reads a clock or bumps a counter, and there is no build configuration that makes it
//! do so. Instrumentation here perturbs the very numbers it produces — the atomics and clock reads
//! of the feature that once gated it cost a couple of percent at the smallest batch size — so the
//! two questions it existed to answer were settled and it was removed. Anything in these docs
//! attributed to instrumentation is a historical measurement; it cannot be reproduced by
//! rebuilding this crate. See [`BENCHMARKS.md`](../BENCHMARKS.md).
//!
//! Measurement that is a property of a whole run rather than of the hot path — load timings,
//! page-fault counts — lives in the `sa-benchmarks` crate, which is excluded from the workspace's
//! `default-members` and never ships.
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
//! `protein-text`, `protein-metadata` or `memory-hints` is a real cross-crate call unless the callee's
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
//! for.** The full sweeps, the host they ran on and the raw figures are in
//! [`BENCHMARKS.md`](../BENCHMARKS.md); the conclusions are here.
//!
//! For scale, the SA is roughly 72% of the index by size, so it dominates residency — and nothing
//! can preload it, since there is no `preloaded-sa`. The protein metadata roughly triples when
//! preloaded, which makes it the most expensive structure to move to the heap and the least
//! rewarding.
//!
//! **The degradation is a concurrency limit, not a bandwidth one.** A major fault blocks its
//! thread, and `prefetch_read` cannot help: a CPU hint instruction cannot fault, so every prefetch
//! in this crate is inert against an absent page. With rayon at the core count, each faulting
//! thread idles a core. Raising `RAYON_NUM_THREADS` leaves the fault *count* essentially unchanged
//! and still buys a large fraction of throughput back under a ceiling — and costs about a tenth of
//! it when the index is resident, so it is a deployment knob rather than a default.
//!
//! **Under a ceiling, plain `mmap` wins.** From the first ceiling that binds, no preloading arm is
//! ahead of it by more than the noise floor, and `preloaded-proteins` is behind it at every one.
//! Below roughly a third of the index the preloaded arms do not degrade gracefully — they collapse,
//! at tens of times `mmap`'s fault rate — and the fully preloaded build is OOM-killed at every
//! ceiling. So preloading is worth having exactly when the whole index is guaranteed resident, and
//! is the wrong default anywhere the ceiling might move.
//!
//! **Even when residency is guaranteed, `preloaded-proteins` is not the arm to reach for.** The
//! sweep is its best case: at `tryptic=false` the search accepts 9-13% of the candidates it
//! examines, so retrieval is a third of the work and the structure it preloads is hot. Under
//! `tryptic` acceptance falls to ~0.5% and retrieval to 1-7% of the work, and
//! `mmap,preloaded-text,preloaded-proteins` is then ahead of it in every cell measured, at a lower
//! resident footprint. Preloading the metadata alone closes the retrieval gap and none of the
//! search one, which is the larger half.
//!
//! **A 6-mer k-mer table is worth its 3.06 GB under a ceiling**, where the resident-case
//! measurement cannot separate it from a 5-mer: it removes far more of the fault load, because the
//! difference is working-set size rather than probe count — a 5-mer narrows the search to several
//! SA pages per query, a 6-mer to about one — and that only matters once pages can be evicted.
//! This is why `sa-builder --kmer-size` documents 6 as the tuning step for constrained
//! deployments; the default stays 5, which is the size the resident measurements support.
//!
//! **All of the loss is in the search phase.** Retrieval is flat across every ceiling, so the
//! two-pass prefetch pipeline in `sa_searcher::retrieval` keeps working under paging and needs
//! nothing. Within search the split is roughly even between the dependent binary-search chain and
//! the contiguous SA range scan.
//!
//! Two further ideas were measured and rejected. **`MADV_WILLNEED` over the SA range about to be
//! scanned does not pay**: the advice lands, but the throughput it buys decays to nothing as
//! threads rise, since oversubscription already overlaps those faults, and it costs resident
//! performance. The comment in `array::mmap` carries the reasoning in full. And **sorting queries
//! by k-mer prefix to create page locality does not work either**: it left the fault count
//! unchanged and cost throughput. With ten thousand queries per rep drawn against 20^6 possible
//! 6-mers, the expected number sharing a prefix is under one, so there is no page reuse for
//! sorting to expose. Locality needs reuse, and this workload has none.

pub mod array;
pub mod kmer_table;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

pub use array::SuffixArrayBackend;
pub use kmer_table::KmerTable;
pub use protein_metadata::ProteinsBackend;

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    /// The sentinel value standing for "no value".
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
