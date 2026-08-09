use std::error::Error;

use memmap2::Mmap;

use super::super::SuffixToProteinMappingBackend;
use crate::Nullable;

/// Mapping backed by a memory-mapped BitVec binary file.
///
/// Format (type 0x02):
/// - 1 byte  type = 0x02
/// - 8 bytes bit_len (u64 LE)
/// - 8 bytes block_count (u64 LE)
/// - block_count × 8 bytes: raw u64 blocks (bit 0 = LSB of block 0)
/// - (block_count/8 + 1) × 16 bytes: superblock cells
///   each cell: [level1: u64 LE] [packed_level2: u64 LE]
///   packed_level2 bits (w-1)*9..(w-1)*9+9 hold cumulative count before word w (w=1..7)
pub struct MmapBitVecSuffixToProtein {
    mmap: Mmap,
    bit_len: u64,
    bits_offset: usize,
    counts_offset: usize
}

impl MmapBitVecSuffixToProtein {
    #[inline]
    fn get_bit(&self, pos: u64) -> bool {
        let block_idx = (pos / 64) as usize;
        let bit_idx = pos % 64;
        let off = self.bits_offset + block_idx * 8;
        let block = u64::from_le_bytes(self.mmap[off..off + 8].try_into().unwrap());
        (block >> bit_idx) & 1 == 1
    }

    #[inline]
    fn rank1(&self, position: u64) -> u64 {
        let bb_index = (position / 512) as usize;
        let word_index = (position / 64) as usize;
        let word_offset = word_index % 8;
        let bit_offset = position % 64;

        let cell_off = self.counts_offset + bb_index * 16;
        let level1 = u64::from_le_bytes(self.mmap[cell_off..cell_off + 8].try_into().unwrap());
        let packed = u64::from_le_bytes(self.mmap[cell_off + 8..cell_off + 16].try_into().unwrap());

        let level2 = if word_offset == 0 { 0u64 } else { (packed >> ((word_offset - 1) * 9)) & 0x1FF };

        let block_off = self.bits_offset + word_index * 8;
        let block = u64::from_le_bytes(self.mmap[block_off..block_off + 8].try_into().unwrap());
        let bit_count = (block << (63 - bit_offset)).count_ones() as u64;

        level1 + level2 + bit_count
    }
}

impl SuffixToProteinMappingBackend for MmapBitVecSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let pos = suffix as u64;
        if pos >= self.bit_len {
            return u32::NULL;
        }
        if self.get_bit(pos) {
            return u32::NULL;
        }
        self.rank1(pos).try_into().unwrap()
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let pos = suffix as usize;
        let bit_off = self.bits_offset + (pos / 64) * 8;
        let sb_off = self.counts_offset + (pos / 512) * 16;
        if bit_off < self.mmap.len() {
            // Bounds checked above, so this indexes the mapping directly (Mmap derefs to [u8]).
            prefetch::prefetch_read(&self.mmap[bit_off] as *const u8);
        }
        if sb_off + 16 <= self.mmap.len() {
            // Bounds checked above, so this indexes the mapping directly (Mmap derefs to [u8]).
            prefetch::prefetch_read(&self.mmap[sb_off] as *const u8);
        }
    }

    fn touch_all_pages(&self) {
        // The bits and counts regions are contiguous from `bits_offset` to the end of the file.
        let end = self.mmap.len();
        text_compression::mmap::touch_all_pages(&self.mmap, self.bits_offset..end);
    }
}

pub(super) fn read_bitvec_mmap(mmap: Mmap) -> Result<MmapBitVecSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 17 {
        return Err("Bitvec mapping file is truncated: missing header fields".into());
    }
    let bit_len = u64::from_le_bytes(mmap[1..9].try_into()?);
    let block_count = u64::from_le_bytes(mmap[9..17].try_into()?) as usize;
    let expected_size = 17 + block_count * 8 + (block_count / 8 + 1) * 16;
    if mmap.len() < expected_size {
        return Err(
            format!("Bitvec mapping file is truncated: expected {} bytes, got {}", expected_size, mmap.len()).into()
        );
    }
    let bits_offset = 17;
    let counts_offset = bits_offset + block_count * 8;
    Ok(MmapBitVecSuffixToProtein { mmap, bit_len, bits_offset, counts_offset })
}

#[cfg(test)]
mod tests {
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::{
        Nullable, ReadBinaryMmap, WriteBinary,
        suffix_to_protein_index::{
            SuffixToProteinMappingBackend, mmap::MmapBackedSuffixToProteinMapping, preloaded::BitVecSuffixToProtein
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
    fn test_mmap_bitvec_roundtrip() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        let original = BitVecSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), loaded.suffix_to_protein(i), "mismatch at suffix {}", i);
        }
    }

    #[test]
    fn test_mmap_bitvec_crosses_superblock_boundary() {
        let segment = "ACGKL-";
        let repeat = 100;
        let raw: String = segment.chars().cycle().take(repeat * segment.len()).collect();
        let mut raw = raw;
        let last = raw.len() - 1;
        raw.replace_range(last..=last, "$");

        let text = InMemoryProteinText::from_string(&raw);
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();

        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        let original = BitVecSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), loaded.suffix_to_protein(i), "mismatch at suffix {}", i);
        }
    }

    #[test]
    fn test_mmap_bitvec_random_equivalence() {
        use std::{
            collections::hash_map::DefaultHasher,
            hash::{Hash, Hasher}
        };

        let len = 2000usize;
        let mut raw: String = (0..len)
            .map(|i| {
                let mut h = DefaultHasher::new();
                i.hash(&mut h);
                if h.finish().is_multiple_of(8) { '-' } else { 'A' }
            })
            .collect();
        raw.pop();
        raw.push('$');

        let text = InMemoryProteinText::from_string(&raw);
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let mmap_idx = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        let original = BitVecSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), mmap_idx.suffix_to_protein(i), "mismatch at position {}", i);
        }
    }

    #[test]
    fn test_search_mmap_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }
}
