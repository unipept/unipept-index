use std::io::{Read, Write};
use std::error::Error;

use memmap2::Mmap;
use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector, Rank9};
use text_compression::ProteinText;

use crate::Nullable;
use super::SuffixToProteinIndex;

/// Mapping that uses O(n) memory (1-2 bits per suffix) with n the size of the input text, with retrieval
/// of the protein in O(1)
#[derive(Debug)]
pub struct BitVecSuffixToProtein {
    rank: Rank9<BitVector<u64>>
}

/// Mapping backed by a memory-mapped BitVec binary file.
///
/// Format (type 0x02):
/// - 1 byte  type = 0x02
/// - 8 bytes bit_len (u64 LE)
/// - 8 bytes block_count (u64 LE)
/// - block_count × 8 bytes: raw u64 blocks (bit 0 = LSB of block 0)
/// - (block_count/8 + 1) × 16 bytes: superblock cells
///     each cell: [level1: u64 LE] [packed_level2: u64 LE]
///     packed_level2 bits (w-1)*9..(w-1)*9+9 hold cumulative count before word w (w=1..7)
pub struct MmapBitVecSuffixToProtein {
    pub(super) mmap: Mmap,
    pub(super) bit_len: u64,
    pub(super) bits_offset: usize,   // byte offset to first raw block
    pub(super) counts_offset: usize, // byte offset to first superblock cell
}

impl SuffixToProteinIndex for BitVecSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let suffix: u64 = suffix.try_into().unwrap();
        if self.rank.get_bit(suffix) {
            return u32::NULL;
        }
        self.rank.rank1(suffix).try_into().unwrap()
    }
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

    /// Count of 1-bits in [0..=position], O(1).
    #[inline]
    fn rank1(&self, position: u64) -> u64 {
        let bb_index = (position / 512) as usize;
        let word_index = (position / 64) as usize;
        let word_offset = word_index % 8;
        let bit_offset = position % 64;

        let cell_off = self.counts_offset + bb_index * 16;
        let level1 = u64::from_le_bytes(self.mmap[cell_off..cell_off + 8].try_into().unwrap());
        let packed = u64::from_le_bytes(self.mmap[cell_off + 8..cell_off + 16].try_into().unwrap());

        let level2 = if word_offset == 0 {
            0u64
        } else {
            (packed >> ((word_offset - 1) * 9)) & 0x1FF
        };

        // Count 1-bits in positions 0..=bit_offset of this block.
        // Shift left by (63 - bit_offset) so that bits 0..=bit_offset land in the
        // upper portion; count_ones() then gives the rank within the block.
        let block_off = self.bits_offset + word_index * 8;
        let block = u64::from_le_bytes(self.mmap[block_off..block_off + 8].try_into().unwrap());
        let bit_count = (block << (63 - bit_offset)).count_ones() as u64;

        level1 + level2 + bit_count
    }
}

impl SuffixToProteinIndex for MmapBitVecSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let pos = suffix as u64;
        if self.get_bit(pos) {
            return u32::NULL;
        }
        self.rank1(pos).try_into().unwrap()
    }
}

impl BitVecSuffixToProtein {
    /// Creates a new BitVecSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - the text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new BitVecSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
        let num_bits = text.len();

        // Create a BitVec (dynamic) first
        let mut bits = BitVector::with_capacity(num_bits as u64);

        // Set bits
        for c in text.iter() {
            bits.push_bit(c == SEPARATION_CHARACTER || c == TERMINATION_CHARACTER);
        }

        let rank = Rank9::new(bits);

        BitVecSuffixToProtein { rank }
    }
}

pub(super) fn write_bitvec_mapping<W: Write>(mapping: &BitVecSuffixToProtein, writer: &mut W) -> Result<(), Box<dyn Error>> {
    let bit_len = mapping.rank.bit_len();
    let block_count = mapping.rank.block_len();
    writer.write_all(&bit_len.to_le_bytes())?;
    writer.write_all(&(block_count as u64).to_le_bytes())?;

    // Write raw 64-bit blocks
    for i in 0..block_count {
        let block: u64 = mapping.rank.get_block(i);
        writer.write_all(&block.to_le_bytes())?;
    }

    // Write superblock array: block_count/8 + 1 entries of 16 bytes each.
    // Each entry: [level1: u64 LE] [packed_level2: u64 LE]
    //   level1        = cumulative 1-bit count before this 512-bit superblock
    //   packed_level2 = 7 × 9-bit counts (before words 1..=7 in the superblock)
    //                   bits (w-1)*9..(w-1)*9+9 hold the count before word w
    let sb_count = block_count / 8 + 1;
    let mut level1: u64 = 0;

    for sb in 0..sb_count {
        let word_start = sb * 8;
        let mut packed_level2: u64 = 0;
        let mut running: u64 = 0; // cumulative count within this superblock

        for w in 0..8usize {
            if w > 0 {
                // Store count-before-word-w at bits (w-1)*9 .. (w-1)*9+8
                packed_level2 |= (running & 0x1FF) << ((w - 1) * 9);
            }
            let word_idx = word_start + w;
            if word_idx < block_count {
                running += mapping.rank.get_block(word_idx).count_ones() as u64;
            }
        }

        writer.write_all(&level1.to_le_bytes())?;
        writer.write_all(&packed_level2.to_le_bytes())?;
        level1 += running;
    }

    Ok(())
}

pub(super) fn read_bitvec_mapping<R: Read>(reader: &mut R) -> Result<BitVecSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let bit_len = u64::from_le_bytes(buf8);
    reader.read_exact(&mut buf8)?;
    let block_count = u64::from_le_bytes(buf8) as usize;

    let mut bits = BitVector::with_capacity(bit_len);
    let mut bits_pushed: u64 = 0;

    for _ in 0..block_count {
        reader.read_exact(&mut buf8)?;
        let block = u64::from_le_bytes(buf8);
        let bits_in_block = std::cmp::min(64, bit_len - bits_pushed);
        for bit_pos in 0..bits_in_block {
            bits.push_bit((block >> bit_pos) & 1 == 1);
            bits_pushed += 1;
        }
    }

    // Read and discard the superblock array written by write_bitvec_mapping
    let sb_count = block_count / 8 + 1;
    let mut discard = [0u8; 16];
    for _ in 0..sb_count {
        reader.read_exact(&mut discard)?;
    }

    Ok(BitVecSuffixToProtein { rank: Rank9::new(bits) })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::ProteinText;

    use crate::{Nullable, ReadBinaryMmap};
    use crate::suffix_to_protein_index::{SuffixToProteinIndex, SuffixToProteinMapping};
    use super::{BitVecSuffixToProtein, write_bitvec_mapping, read_bitvec_mapping};

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
    fn test_search_bitvec() {
        let u8_text = &build_text();
        let index = BitVecSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_bitvec_roundtrip() {
        let text = build_text();
        let original = BitVecSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_bitvec_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_bitvec_mapping(&mut cursor).unwrap();
        // BitVecSuffixToProtein does not derive PartialEq, compare via suffix_to_protein
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), restored.suffix_to_protein(i));
        }
    }

    #[test]
    fn test_mmap_bitvec_roundtrip() {
        let text = build_text();
        let original = BitVecSuffixToProtein::new(&text);

        let mut buf = Vec::new();
        buf.push(2u8); // type byte
        write_bitvec_mapping(&original, &mut buf).unwrap();

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

    #[test]
    fn test_mmap_bitvec_crosses_superblock_boundary() {
        // Build a text long enough to span multiple 512-bit superblocks (>= 8 blocks = 512 chars)
        // We repeat a short pattern to get 600 characters
        let segment = "ACGKL-";
        let repeat = 100; // 600 chars total → 600 bits → 10 blocks → 2 superblocks
        let raw: String = segment.chars().cycle().take(repeat * segment.len()).collect();
        // Make a valid text: replace last char with '$'
        let mut raw = raw;
        let last = raw.len() - 1;
        raw.replace_range(last..=last, "$");

        let text = ProteinText::from_string(&raw);
        let original = BitVecSuffixToProtein::new(&text);

        let mut buf = Vec::new();
        buf.push(2u8);
        write_bitvec_mapping(&original, &mut buf).unwrap();

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

    #[test]
    fn test_search_mmap_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        buf.push(2u8);
        write_bitvec_mapping(&BitVecSuffixToProtein::new(&text), &mut buf).unwrap();
        let tmp = write_to_tempfile(&buf);
        let index = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;

        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL); // SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL); // TERMINATION_CHARACTER
    }
}
