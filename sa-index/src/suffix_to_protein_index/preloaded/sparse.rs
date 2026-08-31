use std::{
    error::Error,
    io::{Read, Write}
};

use binary_traits::WriteBinary;
use protein_metadata::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use protein_text::ProteinTextBackend;

use super::super::SuffixToProteinMappingBackend;
use crate::Nullable;

/// Mapping that uses O(m) memory with m the number of proteins, but retrieval of the protein is
/// O(log m)
#[derive(Debug, PartialEq)]
pub struct SparseSuffixToProtein {
    mapping: Vec<i64>
}

impl SuffixToProteinMappingBackend for SparseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let protein_index = self.mapping.binary_search(&suffix).unwrap_or_else(|index| index - 1);
        // if the next value in the mapping is 1 larger than the current suffix, that means that the
        // current suffix starts with a SEPARATION_CHARACTER or TERMINATION_CHARACTER
        // this means it does not belong to a protein
        if self.mapping[protein_index + 1] == suffix + 1 {
            return u32::NULL;
        }
        protein_index as u32
    }
}

impl SparseSuffixToProtein {
    /// Creates a new SparseSuffixToProtein mapping
    pub fn new<T: ProteinTextBackend>(text: &T) -> Self {
        Self::from_text_parts(text.len(), |i| text.get(i))
    }

    /// Closure-based constructor — works with any text type that exposes `len()` + `get()`.
    ///
    /// The result holds the start position of every protein — the leading `0` is the first one —
    /// and then, from the terminator, one past-the-end entry. That trailing entry is what lets
    /// `suffix_to_protein` read the entry after the one it found for any position in the text.
    pub fn from_text_parts(text_len: usize, get_char: impl Fn(usize) -> u8) -> Self {
        let mut suffix_index_to_protein: Vec<i64> = vec![0];
        for i in 0..text_len {
            if get_char(i) == SEPARATION_CHARACTER || get_char(i) == TERMINATION_CHARACTER {
                suffix_index_to_protein.push(i as i64 + 1);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        SparseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

impl WriteBinary for SparseSuffixToProtein {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        writer.write_all(&[1u8])?;
        writer.write_all(&(self.mapping.len() as u64).to_le_bytes())?;
        for &val in &self.mapping {
            writer.write_all(&val.to_le_bytes())?;
        }
        Ok(())
    }
}

/// Reads the body of a sparse mapping, after the type byte
/// [`InMemorySuffixToProteinMapping::read_binary`](super::InMemorySuffixToProteinMapping) consumed.
pub(super) fn read_sparse_mapping<R: Read>(reader: &mut R) -> Result<SparseSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let count = u64::from_le_bytes(buf8) as usize;
    let mut mapping = super::try_alloc_exact(count, "sparse")?;
    for _ in 0..count {
        reader.read_exact(&mut buf8)?;
        mapping.push(i64::from_le_bytes(buf8));
    }
    Ok(SparseSuffixToProtein { mapping })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{SparseSuffixToProtein, read_sparse_mapping};
    use crate::suffix_to_protein_index::test_utils::{assert_sample_lookups, sample_text, to_binary};

    #[test]
    fn test_sparse_build() {
        let index = SparseSuffixToProtein::new(&sample_text());
        let expected = SparseSuffixToProtein { mapping: vec![0, 4, 7, 11] };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_sparse() {
        assert_sample_lookups(&SparseSuffixToProtein::new(&sample_text()));
    }

    #[test]
    fn test_sparse_roundtrip() {
        let text = sample_text();
        let buf = to_binary(SparseSuffixToProtein::new(&text));
        assert_eq!(buf[0], 1u8);
        let restored = read_sparse_mapping(&mut Cursor::new(&buf[1..])).unwrap();
        assert_eq!(SparseSuffixToProtein::new(&text), restored);
    }
}
