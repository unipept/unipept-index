use std::{
    error::Error,
    io::{Read, Write}
};

use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinTextBackend;

use super::super::SuffixToProteinMappingBackend;
use crate::{Nullable, WriteBinary};

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

pub(super) fn read_sparse_mapping<R: Read>(reader: &mut R) -> Result<SparseSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let count = u64::from_le_bytes(buf8) as usize;
    let mut mapping = Vec::with_capacity(count);
    for _ in 0..count {
        reader.read_exact(&mut buf8)?;
        mapping.push(i64::from_le_bytes(buf8));
    }
    Ok(SparseSuffixToProtein { mapping })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::InMemoryProteinText;

    use super::{SparseSuffixToProtein, read_sparse_mapping};
    use crate::{Nullable, WriteBinary, suffix_to_protein_index::SuffixToProteinMappingBackend};

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    #[test]
    fn test_sparse_build() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        let expected = SparseSuffixToProtein { mapping: vec![0, 4, 7, 11] };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_sparse() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_sparse_roundtrip() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let mut cursor = Cursor::new(&buf[1..]);
        let restored = read_sparse_mapping(&mut cursor).unwrap();
        assert_eq!(SparseSuffixToProtein::new(&text), restored);
    }
}
