use std::{
    error::Error,
    io::{Read, Write}
};

use memmap2::Mmap;
use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinText;

use super::SuffixToProteinIndex;
use crate::Nullable;

/// Mapping that uses O(m) memory with m the number of proteins, but retrieval of the protein is
/// O(log m)
#[derive(Debug, PartialEq)]
pub struct SparseSuffixToProtein {
    mapping: Vec<i64>
}

/// Mapping backed by a memory-mapped Sparse binary file.
/// Format: [1 byte type=0x01] [8 bytes count (u64 LE)] [count × 8 bytes (i64 LE)]
pub struct MmapSparseSuffixToProtein {
    pub(super) mmap: Mmap,
    pub(super) data_offset: usize, // 9 = 1 (type) + 8 (count)
    pub(super) count: usize
}

impl SuffixToProteinIndex for SparseSuffixToProtein {
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

impl SuffixToProteinIndex for MmapSparseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let read_val = |i: usize| -> i64 {
            let off = self.data_offset + i * 8;
            i64::from_le_bytes(self.mmap[off..off + 8].try_into().unwrap())
        };

        // Binary search for the largest index where mapping[index] <= suffix
        let mut lo = 0usize;
        let mut hi = self.count;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if read_val(mid) <= suffix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let protein_index = lo - 1;

        // If the next boundary == suffix + 1, this suffix is a separator/terminator
        if read_val(protein_index + 1) == suffix + 1 {
            return u32::NULL;
        }
        protein_index as u32
    }
}

impl SparseSuffixToProtein {
    /// Creates a new SparseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new SparseSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
        let mut suffix_index_to_protein: Vec<i64> = vec![0];
        for (index, char) in text.iter().enumerate() {
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                suffix_index_to_protein.push(index as i64 + 1);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        SparseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

pub(super) fn write_sparse_mapping<W: Write>(
    mapping: &SparseSuffixToProtein,
    writer: &mut W
) -> Result<(), Box<dyn Error>> {
    let count = mapping.mapping.len() as u64;
    writer.write_all(&count.to_le_bytes())?;
    for &val in &mapping.mapping {
        writer.write_all(&val.to_le_bytes())?;
    }
    Ok(())
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
    use std::io::{Cursor, Write as IoWrite};

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::ProteinText;

    use super::{SparseSuffixToProtein, read_sparse_mapping, write_sparse_mapping};
    use crate::{
        Nullable, ReadBinaryMmap,
        suffix_to_protein_index::legacy::{SuffixToProteinIndex, SuffixToProteinMapping}
    };

    fn build_text() -> ProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        ProteinText::from_string(&text)
    }

    fn write_to_tempfile(buf: &[u8]) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(buf).unwrap();
        tmp.flush().unwrap();
        tmp
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
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_sparse_roundtrip() {
        let text = build_text();
        let original = SparseSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_sparse_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_sparse_mapping(&mut cursor).unwrap();
        assert_eq!(original, restored);
    }

    #[test]
    fn test_mmap_sparse_roundtrip() {
        let text = build_text();
        let original = SparseSuffixToProtein::new(&text);

        let mut buf = Vec::new();
        buf.push(1u8); // type byte
        write_sparse_mapping(&original, &mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;

        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), loaded.suffix_to_protein(i), "mismatch at suffix {}", i);
        }
    }

    #[test]
    fn test_search_mmap_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        buf.push(1u8);
        write_sparse_mapping(&SparseSuffixToProtein::new(&text), &mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL); // SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL); // TERMINATION_CHARACTER
    }
}
