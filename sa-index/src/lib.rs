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
//! # The `mmap` feature: two builds, one API
//!
//! Every storage structure has two implementations — one holding owned memory, one borrowing a
//! memory mapping — and the `mmap` feature picks which. The selection is a *type alias*, resolved
//! at compile time, not a runtime branch:
//!
//! | alias | preloaded (default) | with `mmap` |
//! |---|---|---|
//! | [`SuffixArray`] | `InMemorySA` | `MmapBackedSA` |
//! | [`suffix_to_protein_index::SuffixToProteinMapping`] | `PreloadedSuffixToProteinMapping` | `MmapBackedSuffixToProteinMapping` |
//! | `sa_mappings::proteins::Proteins` | `InMemoryProteins` | `MmapBackedProteins` |
//! | `text_compression::ProteinText` | `InMemoryProteinText` | `MmapBackedProteinText` |
//!
//! Consequences worth knowing before reading further:
//!
//! * **No crate declares `default = [...]`**, so a plain `cargo build` or `cargo test` gives you
//!   the *preloaded* configuration. The production server is built `--features mmap`.
//! * The preloaded half is compiled in **both** configurations, because it owns the `WriteBinary`
//!   implementations that produce the files the mmap half reads.
//! * The two are kept deliberately separate rather than sharing code, so that a tuning change to
//!   one cannot perturb the other. Where that produces near-duplicate code it is marked as
//!   intentional at the site.
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
//! sets no `[profile.release]`, there is **no LTO**, so a call into `bitarray`,
//! `text-compression`, `sa-mappings` or `prefetch` is a real cross-crate call unless the callee
//! is `#[inline]`. Those attributes on the small getters are load-bearing, not decoration.
//! Enabling `lto = "thin"` would be the principled fix and would likely let many of them go, but
//! it changes codegen and therefore needs its own measurement.
//!
//! Decisions that were measured and *rejected* are recorded next to the code they would have
//! touched, so they are not rediscovered and retried. See `docs/design/` for the long form.

pub use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};

pub mod array;
pub mod kmer_table;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

pub use array::{SuffixArray, SuffixArrayBackend};
pub use kmer_table::KmerTable;
pub use sa_mappings::proteins::ProteinsBackend;
pub use sa_searcher::{SearchTuning, MAX_VALIDATE_BATCH};

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
