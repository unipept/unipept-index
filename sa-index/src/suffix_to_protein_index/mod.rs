pub mod preloaded;
#[cfg(feature = "mmap")]
pub mod mmap;

pub use preloaded::{DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein, InMemorySuffixToProteinMapping};
#[cfg(feature = "mmap")]
pub use mmap::{MmapDenseSuffixToProtein, MmapSparseSuffixToProtein, MmapBitVecSuffixToProtein, MmapBackedSuffixToProteinMapping};

#[cfg(feature = "mmap")]
pub type SuffixToProteinMapping = MmapBackedSuffixToProteinMapping;
#[cfg(not(feature = "mmap"))]
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
