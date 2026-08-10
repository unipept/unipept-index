//! Turning a text position into a protein index.
//!
//! Search yields positions in the concatenated text; results are about proteins. This module
//! answers "which protein contains position i?", and returns `u32::NULL` for the separator and
//! terminator bytes, which belong to no protein.
//!
//! Three representations, trading space against lookup cost, chosen at build time
//! (`sa-builder --mapping-style`) and recorded in the file:
//!
//! * **Dense** — one `u32` per text position. One load per lookup, but ~4 bytes per residue,
//!   which at UniProt scale is over a gigabyte.
//! * **Sparse** — the start position of each protein, binary-searched. Smallest, but O(log n)
//!   dependent loads per lookup, each likely a cache miss.
//! * **BitVec** — a bit per text position marking separators, with a rank structure over it.
//!   Near-dense speed at a fraction of the size; the default.
//!
//! Each has a preloaded and an mmap implementation, selected through the
//! [`SuffixToProteinMapping`] alias: `mmap` maps it, and `preloaded-mapping` pulls this one
//! structure back into owned memory while the rest of the index stays mapped. See the crate docs.

#[cfg(feature = "mmap")]
pub mod mmap;
pub mod preloaded;

#[cfg(feature = "mmap")]
pub use mmap::{
    MmapBackedSuffixToProteinMapping, MmapBitVecSuffixToProtein, MmapDenseSuffixToProtein, MmapSparseSuffixToProtein
};
pub use preloaded::{
    BitVecSuffixToProtein, DenseSuffixToProtein, InMemorySuffixToProteinMapping, SparseSuffixToProtein
};

/// Type alias — resolves to the active mapping backend for this build.
#[cfg(all(feature = "mmap", not(feature = "preloaded-mapping")))]
pub type SuffixToProteinMapping = MmapBackedSuffixToProteinMapping;
/// Type alias — resolves to the active mapping backend for this build.
///
/// Owned memory either because this is a preloaded build, or because `preloaded-mapping` asked
/// for this structure specifically.
#[cfg(any(not(feature = "mmap"), feature = "preloaded-mapping"))]
pub type SuffixToProteinMapping = InMemorySuffixToProteinMapping;

/// Trait implemented by the SuffixToProtein mappings
pub trait SuffixToProteinMappingBackend: Send + Sync {
    /// Returns the index of the protein in the protein list for the given suffix
    fn suffix_to_protein(&self, suffix: i64) -> u32;

    /// Non-blocking hardware prefetch hint for the data that `suffix_to_protein(suffix)` will access.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn prefetch_for_suffix(&self, _suffix: i64) {}

    /// Reads at least one byte from every OS page in the mmap backing this mapping,
    /// ensuring all pages are resident in the page cache.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn touch_all_pages(&self) {}
}
