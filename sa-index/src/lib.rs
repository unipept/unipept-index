use bitarray::BitArray;

pub mod binary;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

/// Represents a suffix array.
pub enum SuffixArray {
    /// The original suffix array.
    Original(Vec<i64>),
    /// The compressed suffix array.
    Compressed(BitArray)
}

impl SuffixArray {
    /// Returns the length of the suffix array.
    ///
    /// # Returns
    ///
    /// The length of the suffix array.
    pub fn len(&self) -> usize {
        match self {
            SuffixArray::Original(sa) => sa.len(),
            SuffixArray::Compressed(sa) => sa.len()
        }
    }

    /// Returns the suffix array at the given index.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the suffix array.
    ///
    /// # Returns
    ///
    /// The suffix array at the given index.
    pub fn get(&self, index: usize) -> i64 {
        match self {
            SuffixArray::Original(sa) => sa[index],
            SuffixArray::Compressed(sa) => sa.get(index) as i64
        }
    }

    /// Returns whether the suffix array is empty.
    ///
    /// # Returns
    ///
    /// True if the suffix array is empty, false otherwise.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

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
    fn test_suffix_array_original() {
        let sa = SuffixArray::Original(vec![1, 2, 3, 4, 5]);
        assert_eq!(sa.len(), 5);
        assert_eq!(sa.get(0), 1);
        assert_eq!(sa.get(1), 2);
        assert_eq!(sa.get(2), 3);
        assert_eq!(sa.get(3), 4);
        assert_eq!(sa.get(4), 5);
    }

    #[test]
    fn test_suffix_array_compressed() {
        let mut bitarray = BitArray::with_capacity(5, 40);
        bitarray.set(0, 1);
        bitarray.set(1, 2);
        bitarray.set(2, 3);
        bitarray.set(3, 4);
        bitarray.set(4, 5);

        let sa = SuffixArray::Compressed(bitarray);
        assert_eq!(sa.len(), 5);
        assert_eq!(sa.get(0), 1);
        assert_eq!(sa.get(1), 2);
        assert_eq!(sa.get(2), 3);
        assert_eq!(sa.get(3), 4);
        assert_eq!(sa.get(4), 5);
    }

    #[test]
    fn test_nullable_is_null() {
        assert_eq!(u32::NULL.is_null(), true);
        assert_eq!(0u32.is_null(), false);
    }
}
