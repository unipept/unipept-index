use std::io::{Read, Write};
use std::error::Error;

use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinText;

use crate::Nullable;
use super::SuffixToProteinIndex;
#[cfg(feature = "mmap")]
use memmap2::Mmap;

/// Mapping that uses O(n) memory with n the size of the input text, but retrieval of the protein is
/// in O(1)
#[derive(Debug, PartialEq)]
pub struct DenseSuffixToProtein {
    // UniProtKB does not have more that u32::MAX proteins, so a larger type is not needed
    mapping: Vec<u32>
}

/// Mapping backed by a memory-mapped Dense binary file.
/// Format: [1 byte type=0x00] [8 bytes count (u64 LE)] [count × 4 bytes (u32 LE)]
#[cfg(feature = "mmap")]
pub struct MmapDenseSuffixToProtein {
    pub(super) mmap: Mmap,
    pub(super) data_offset: usize, // 9 = 1 (type) + 8 (count)
}

impl SuffixToProteinIndex for DenseSuffixToProtein {
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

#[cfg(feature = "mmap")]
impl SuffixToProteinIndex for MmapDenseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let off = self.data_offset + suffix as usize * 4;
        u32::from_le_bytes(self.mmap[off..off + 4].try_into().unwrap())
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let off = self.data_offset + suffix as usize * 4;
        if off < self.mmap.len() {
            // safe: Mmap: Deref<Target=[u8]>, bounds checked above
            prefetch::prefetch_read(&self.mmap[off] as *const u8);
        }
    }

    fn touch_all_pages(&self) {
        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Sequential);

        for chunk in self.mmap[self.data_offset..].chunks(4096) {
            std::hint::black_box(chunk[0]);
        }

        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Random);
    }
}

impl DenseSuffixToProtein {
    /// Creates a new DenseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new DenseSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
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

pub(super) fn write_dense_mapping<W: Write>(mapping: &DenseSuffixToProtein, writer: &mut W) -> Result<(), Box<dyn Error>> {
    let count = mapping.mapping.len() as u64;
    writer.write_all(&count.to_le_bytes())?;
    for &val in &mapping.mapping {
        writer.write_all(&val.to_le_bytes())?;
    }
    Ok(())
}

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

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use std::io::Cursor;
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::ProteinText;

    use crate::Nullable;
    #[cfg(feature = "mmap")]
    use crate::ReadBinaryMmap;
    use crate::suffix_to_protein_index::{SuffixToProteinIndex, SuffixToProteinMapping};
    use super::{DenseSuffixToProtein, write_dense_mapping, read_dense_mapping};

    fn build_text() -> ProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        ProteinText::from_string(&text)
    }

    #[cfg(feature = "mmap")]
    fn write_to_tempfile(buf: &[u8]) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(buf).unwrap();
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn test_dense_build() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        let expected = DenseSuffixToProtein {
            mapping: vec![0, 0, 0, u32::NULL, 1, 1, u32::NULL, 2, 2, 2, u32::NULL]
        };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_dense() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_dense_roundtrip() {
        let text = build_text();
        let original = DenseSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_dense_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_dense_mapping(&mut cursor).unwrap();
        assert_eq!(original, restored);
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn test_mmap_dense_roundtrip() {
        let text = build_text();
        let original = DenseSuffixToProtein::new(&text);

        let mut buf = Vec::new();
        buf.push(0u8); // type byte
        write_dense_mapping(&original, &mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;

        for i in 0..text.len() as i64 {
            assert_eq!(
                original.suffix_to_protein(i),
                loaded.suffix_to_protein(i),
                "mismatch at suffix {}",
                i
            );
        }
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn test_search_mmap_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        buf.push(0u8);
        write_dense_mapping(&DenseSuffixToProtein::new(&text), &mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL); // SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL); // TERMINATION_CHARACTER
    }
}
