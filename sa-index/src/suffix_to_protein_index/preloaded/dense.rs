use std::{
    error::Error,
    io::{Read, Write}
};

use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinTextBackend;

use super::super::SuffixToProteinMappingBackend;
use crate::{Nullable, WriteBinary};

/// Mapping that uses O(n) memory with n the size of the input text, but retrieval of the protein is
/// in O(1)
#[derive(Debug, PartialEq)]
pub struct DenseSuffixToProtein {
    // UniProtKB does not have more that u32::MAX proteins, so a larger type is not needed
    mapping: Vec<u32>
}

impl SuffixToProteinMappingBackend for DenseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        self.mapping[suffix as usize]
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let idx = suffix as usize;
        if idx < self.mapping.len() {
            prefetch::prefetch_read(&self.mapping[idx] as *const u32);
        }
    }
}

impl DenseSuffixToProtein {
    /// Creates a new DenseSuffixToProtein mapping
    pub fn new<T: ProteinTextBackend>(text: &T) -> Self {
        Self::from_text_parts(text.len(), |i| text.get(i))
    }

    /// Closure-based constructor — works with any text type that exposes `len()` + `get()`.
    pub fn from_text_parts(text_len: usize, get_char: impl Fn(usize) -> u8) -> Self {
        let mut current_protein_index: u32 = 0;
        let mut suffix_index_to_protein: Vec<u32> = vec![];
        for i in 0..text_len {
            let char = get_char(i);
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                current_protein_index += 1;
                suffix_index_to_protein.push(u32::NULL);
            } else {
                assert_ne!(current_protein_index, u32::NULL);
                suffix_index_to_protein.push(current_protein_index);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        DenseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

impl WriteBinary for DenseSuffixToProtein {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        writer.write_all(&[0u8])?;
        writer.write_all(&(self.mapping.len() as u64).to_le_bytes())?;
        for &val in &self.mapping {
            writer.write_all(&val.to_le_bytes())?;
        }
        Ok(())
    }
}

/// Reads the body of a dense mapping. Unlike the mmap readers, which are handed the whole file,
/// this starts after the type byte that [`InMemorySuffixToProteinMapping::read_binary`] has
/// already consumed to get here.
///
/// [`InMemorySuffixToProteinMapping::read_binary`]: super::InMemorySuffixToProteinMapping
pub(super) fn read_dense_mapping<R: Read>(reader: &mut R) -> Result<DenseSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let count = u64::from_le_bytes(buf8) as usize;
    let mut mapping = Vec::with_capacity(count);
    for _ in 0..count {
        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        mapping.push(u32::from_le_bytes(buf4));
    }
    Ok(DenseSuffixToProtein { mapping })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{DenseSuffixToProtein, read_dense_mapping};
    use crate::{
        Nullable,
        suffix_to_protein_index::test_utils::{assert_sample_lookups, sample_text, to_binary}
    };

    #[test]
    fn test_dense_build() {
        let index = DenseSuffixToProtein::new(&sample_text());
        let expected = DenseSuffixToProtein {
            mapping: vec![0, 0, 0, u32::NULL, 1, 1, u32::NULL, 2, 2, 2, u32::NULL]
        };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_dense() {
        assert_sample_lookups(&DenseSuffixToProtein::new(&sample_text()));
    }

    #[test]
    fn test_dense_roundtrip() {
        let text = sample_text();
        let buf = to_binary(DenseSuffixToProtein::new(&text));
        assert_eq!(buf[0], 0u8);
        let restored = read_dense_mapping(&mut Cursor::new(&buf[1..])).unwrap();
        assert_eq!(DenseSuffixToProtein::new(&text), restored);
    }
}
