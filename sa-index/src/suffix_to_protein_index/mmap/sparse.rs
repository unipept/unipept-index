use std::error::Error;

use memmap2::Mmap;

use crate::Nullable;
use super::super::SuffixToProteinMappingBackend;

/// Mapping backed by a memory-mapped Sparse binary file.
/// Format: [1 byte type=0x01] [8 bytes count (u64 LE)] [count × 8 bytes (i64 LE)]
pub struct MmapSparseSuffixToProtein {
    mmap: Mmap,
    data_offset: usize, // 9 = 1 (type) + 8 (count)
    count: usize,
}

impl SuffixToProteinMappingBackend for MmapSparseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let read_val = |i: usize| -> i64 {
            let off = self.data_offset + i * 8;
            i64::from_le_bytes(self.mmap[off..off + 8].try_into().unwrap())
        };

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

        if read_val(protein_index + 1) == suffix + 1 {
            return u32::NULL;
        }
        protein_index as u32
    }

    fn touch_all_pages(&self) {
        let end = self.data_offset + self.count * 8;
        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Sequential);

        for chunk in self.mmap[self.data_offset..end].chunks(4096) {
            std::hint::black_box(chunk[0]);
        }

        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Random);
    }
}

pub(super) fn read_sparse_mmap(mmap: Mmap) -> Result<MmapSparseSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 9 {
        return Err("Sparse mapping file is truncated: missing count header".into());
    }
    let count = u64::from_le_bytes(mmap[1..9].try_into()?) as usize;
    Ok(MmapSparseSuffixToProtein { mmap, count, data_offset: 9 })
}

#[cfg(test)]
mod tests {
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::InMemoryProteinText;

    use text_compression::ProteinTextBackend;

    use crate::{Nullable, ReadBinaryMmap, WriteBinary};
    use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
    use crate::suffix_to_protein_index::preloaded::SparseSuffixToProtein;
    use crate::suffix_to_protein_index::mmap::MmapBackedSuffixToProteinMapping;

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
    fn test_mmap_sparse_roundtrip() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        let original = SparseSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(
                original.suffix_to_protein(i),
                loaded.suffix_to_protein(i),
                "mismatch at suffix {}",
                i
            );
        }
    }

    #[test]
    fn test_search_mmap_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }
}
