use std::error::Error;

use memmap2::Mmap;

use super::super::SuffixToProteinMappingBackend;

/// Mapping backed by a memory-mapped Dense binary file.
/// Format: [1 byte type=0x00] [8 bytes count (u64 LE)] [count × 4 bytes (u32 LE)]
pub struct MmapDenseSuffixToProtein {
    mmap: Mmap,
    data_offset: usize // 9 = 1 (type) + 8 (count)
}

impl SuffixToProteinMappingBackend for MmapDenseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let off = self.data_offset + suffix as usize * 4;
        u32::from_le_bytes(self.mmap[off..off + 4].try_into().unwrap())
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let off = self.data_offset + suffix as usize * 4;
        if off < self.mmap.len() {
            // Bounds checked above, so this indexes the mapping directly (Mmap derefs to [u8]).
            prefetch::prefetch_read(&self.mmap[off] as *const u8);
        }
    }

    fn touch_all_pages(&self) {
        let end = self.mmap.len();
        text_compression::mmap::touch_all_pages(&self.mmap, self.data_offset..end);
    }
}

pub(super) fn read_dense_mmap(mmap: Mmap) -> Result<MmapDenseSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 9 {
        return Err("Dense mapping file is truncated: missing count header".into());
    }
    let _count = u64::from_le_bytes(mmap[1..9].try_into()?) as usize;
    Ok(MmapDenseSuffixToProtein { mmap, data_offset: 9 })
}

#[cfg(test)]
mod tests {
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::{
        Nullable, ReadBinaryMmap, WriteBinary,
        suffix_to_protein_index::{
            SuffixToProteinMappingBackend, mmap::MmapBackedSuffixToProteinMapping, preloaded::DenseSuffixToProtein
        }
    };

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    fn write_to_tempfile(buf: &[u8]) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(buf).unwrap();
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn test_mmap_dense_roundtrip() {
        let text = build_text();
        let mut buf = Vec::new();
        DenseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        let original = DenseSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), loaded.suffix_to_protein(i), "mismatch at suffix {}", i);
        }
    }

    #[test]
    fn test_search_mmap_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        DenseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }
}
