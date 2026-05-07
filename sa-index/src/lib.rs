pub use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};

pub mod array;
pub mod kmer_table;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

pub use array::{SuffixArray, SuffixArrayBackend};
pub use kmer_table::KmerTable;

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
