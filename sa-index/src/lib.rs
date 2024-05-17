use bitarray::BitArray;

pub mod binary;
pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;

pub enum SuffixArray {
    Original(Vec<i64>),
    Compressed(BitArray<37>)
}

impl SuffixArray {
    pub fn len(&self) -> usize {
        match self {
            SuffixArray::Original(sa) => sa.len(),
            SuffixArray::Compressed(sa) => sa.len()
        }
    }

    pub fn get(&self, index: usize) -> i64 {
        match self {
            SuffixArray::Original(sa) => sa[index],
            SuffixArray::Compressed(sa) => sa.get(index) as i64
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    fn is_null(&self) -> bool;
}

impl Nullable<u32> for u32 {
    const NULL: u32 = u32::MAX;

    fn is_null(&self) -> bool {
        *self == Self::NULL
    }
}
