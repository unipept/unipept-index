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
    fn implied_text_len(&self) -> Option<usize> {
        Some(self.bit_len as usize)
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

    fn touch_all_pages(&self) -> u64 {
        // The bits and counts regions are contiguous from `bits_offset` to the end of the file.
        let end = self.mmap.len();
        text_compression::mmap::touch_all_pages(&self.mmap, self.bits_offset..end)
    }
}

/// Maps a bitvec mapping file. Unlike the other two readers this also checks the body, in the two
/// ways that between them make every in-range lookup addressable:
///
/// * the header's `block_count` fixes the size of both the bits and the superblock cells exactly,
///   so a file too short to hold them is rejected here;
/// * `bit_len` is what `suffix_to_protein` bounds a position against, so a `block_count` too small
///   to cover it would let an in-range position index past the bits region. That is the same check
///   `preloaded::bitvec::read_bitvec_mapping` makes, and must stay in step with it: without it a
///   crafted header loads cleanly and panics inside a request handler instead, which is what
///   [`ReadBinaryMmap`](binary_traits::ReadBinaryMmap)'s contract forbids.
pub(super) fn read_bitvec_mmap(mmap: Mmap) -> Result<MmapBitVecSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 17 {
        return Err("The bitvec mapping file is too small to contain the header".into());
    }
    let bit_len = u64::from_le_bytes(mmap[1..9].try_into()?);
    let block_count = u64::from_le_bytes(mmap[9..17].try_into()?) as usize;

    let needed = (bit_len as usize).div_ceil(64);
    if block_count < needed {
        return Err(format!(
            "Bitvec mapping declares {bit_len} bits but holds only {block_count} of the {needed} blocks that needs"
        )
        .into());
    }

    // Checked, because `block_count` is untrusted: an unchecked `* 8` wraps for a header near
    // `usize::MAX`, which would make `expected_size` small enough to pass the length check below
    // and leave `counts_offset` pointing anywhere.
    let too_many = || -> Box<dyn Error> { "The bitvec mapping header declares too many blocks".into() };
    let bits_bytes = block_count.checked_mul(8).ok_or_else(too_many)?;
    let counts_bytes = (block_count / 8 + 1).checked_mul(16).ok_or_else(too_many)?;
    let expected_size =
        17usize.checked_add(bits_bytes).and_then(|n| n.checked_add(counts_bytes)).ok_or_else(too_many)?;
    if mmap.len() < expected_size {
        return Err(format!(
            "The bitvec mapping file is too small to contain the mapping data: expected {} bytes, got {}",
            expected_size,
            mmap.len()
        )
        .into());
    }
    let bits_offset = 17;
    let counts_offset = bits_offset + bits_bytes;
    Ok(MmapBitVecSuffixToProtein { mmap, bit_len, bits_offset, counts_offset })
}

#[cfg(test)]
mod tests {
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::suffix_to_protein_index::{
        mmap::test_utils::write_and_map,
        preloaded::BitVecSuffixToProtein,
        test_utils::{assert_agree, many_proteins_text, sample_text}
    };

    /// The absolute answers are pinned by `mmap::tests::test_load_mmap_bitvec`; what this adds is
    /// that the two-level rank read out of the mapped file agrees with `Rank9` everywhere. The
    /// second text spans several superblocks, so it also exercises the level1 counts, which stay
    /// zero for anything shorter than 512 bits.
    #[test]
    fn test_mmap_bitvec_roundtrip() {
        for text in [sample_text(), many_proteins_text(300, 5)] {
            let (loaded, _tmp) = write_and_map(BitVecSuffixToProtein::new(&text));
            assert_agree(&BitVecSuffixToProtein::new(&text), &loaded, text.len());
        }
    }

    /// Proteins of regular length leave the separators at regular offsets, which a rank bug can
    /// happen to line up with. This scatters them irregularly across many superblocks instead.
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
        let (loaded, _tmp) = write_and_map(BitVecSuffixToProtein::new(&text));
        assert_agree(&BitVecSuffixToProtein::new(&text), &loaded, text.len());
    }
}
