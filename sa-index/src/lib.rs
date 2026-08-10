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
