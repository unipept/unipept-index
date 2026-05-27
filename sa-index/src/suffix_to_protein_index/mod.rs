use clap::ValueEnum;

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

/// Enum used to define the commandline arguments and choose which index style is used
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SuffixToProteinMappingStyle {
    Dense,
    Sparse,
    BitVec
}

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

#[cfg(test)]
mod tests {
    use clap::ValueEnum;
    use crate::suffix_to_protein_index::SuffixToProteinMappingStyle;

    #[test]
    fn test_suffix_to_protein_mapping_style() {
        assert_eq!(SuffixToProteinMappingStyle::Dense, SuffixToProteinMappingStyle::from_str("dense", false).unwrap());
        assert_eq!(
            SuffixToProteinMappingStyle::Sparse,
            SuffixToProteinMappingStyle::from_str("sparse", false).unwrap()
        );
        assert_eq!(
            SuffixToProteinMappingStyle::BitVec,
            SuffixToProteinMappingStyle::from_str("bit-vec", false).unwrap()
        );
    }
}
