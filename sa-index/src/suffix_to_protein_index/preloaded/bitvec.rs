use std::io::{Read, Write};
use std::error::Error;

use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector, Rank9};
use text_compression::ProteinTextBackend;

use crate::{Nullable, WriteBinary};
use super::super::SuffixToProteinMappingBackend;

/// Mapping that uses O(n) memory (1-2 bits per suffix) with n the size of the input text, with retrieval
/// of the protein in O(1)
#[derive(Debug)]
pub struct BitVecSuffixToProtein {
    rank: Rank9<BitVector<u64>>
}

impl SuffixToProteinMappingBackend for BitVecSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let suffix: u64 = suffix.try_into().unwrap();
        if self.rank.get_bit(suffix) {
            return u32::NULL;
        }
        self.rank.rank1(suffix).try_into().unwrap()
    }
}

impl BitVecSuffixToProtein {
    /// Creates a new BitVecSuffixToProtein mapping
    pub fn new<T: ProteinTextBackend>(text: &T) -> Self {
        Self::from_text_parts(text.len(), |i| text.get(i))
    }

    pub fn from_text_parts(text_len: usize, get_char: impl Fn(usize) -> u8) -> Self {
        let mut bits = BitVector::with_capacity(text_len as u64);
        for i in 0..text_len {
            let c = get_char(i);
            bits.push_bit(c == SEPARATION_CHARACTER || c == TERMINATION_CHARACTER);
        }
        BitVecSuffixToProtein { rank: Rank9::new(bits) }
    }
}

/// On-disk format for the BitVec mapping (type byte `0x02`).
///
/// ```text
/// [ type: u8 = 0x02 ]
/// [ bit_len: u64 LE ]        one bit per text position
/// [ block_count: u64 LE ]
/// [ blocks ]                 block_count * u64 LE, bit 0 = LSB of block 0
/// [ superblocks ]            (block_count / 8 + 1) cells of 16 bytes:
///                              [ level1: u64 LE ] [ packed_level2: u64 LE ]
/// ```
///
/// # The rank structure
///
/// A bit marks each text position that is *not* part of a protein (a separator or the
/// terminator). The protein index for a position is then the number of unset bits before it,
/// which is what makes `rank1` the lookup. Answering that in constant time needs precomputed
/// counts, stored in two levels:
///
/// * **level1** — the cumulative count before this superblock, i.e. before every one of its 8
///   words (512 bits). A full `u64`, since it can reach `bit_len`.
/// * **packed_level2** — seven 9-bit sub-counts, one per word after the first, each the
///   cumulative count within the superblock before that word. Nine bits suffice because a count
///   within a 512-bit superblock cannot exceed 512, and seven of them fit in a `u64` with a byte
///   to spare, so the whole cell is exactly 16 bytes and one cache line holds four of them.
///
/// Hence the constants below: `& 0x1FF` masks a 9-bit sub-count, `(w - 1) * 9` places it, and
/// the loop covers `w = 1..8` because word 0's sub-count is always zero and is not stored.
///
/// Only the mmap reader consumes this — the preloaded reader rebuilds `Rank9` from the raw bits
/// and skips the superblocks entirely — so the layout must be kept in step by hand with
/// `suffix_to_protein_index::mmap::bitvec`, which documents the same structure from the reading
/// side.
impl WriteBinary for BitVecSuffixToProtein {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        writer.write_all(&[2u8])?;
        let bit_len = self.rank.bit_len();
        let block_count = self.rank.block_len();
        writer.write_all(&bit_len.to_le_bytes())?;
        writer.write_all(&(block_count as u64).to_le_bytes())?;

        for i in 0..block_count {
            let block: u64 = self.rank.get_block(i);
            writer.write_all(&block.to_le_bytes())?;
        }

        let sb_count = block_count / 8 + 1;
        let mut level1: u64 = 0;

        for sb in 0..sb_count {
            let word_start = sb * 8;
            let mut packed_level2: u64 = 0;
            let mut running: u64 = 0;

            for w in 0..8usize {
                if w > 0 {
                    packed_level2 |= (running & 0x1FF) << ((w - 1) * 9);
                }
                let word_idx = word_start + w;
                if word_idx < block_count {
                    running += self.rank.get_block(word_idx).count_ones() as u64;
                }
            }

            writer.write_all(&level1.to_le_bytes())?;
            writer.write_all(&packed_level2.to_le_bytes())?;
            level1 += running;
        }

        Ok(())
    }
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

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::{Nullable, WriteBinary};
    use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
    use super::{BitVecSuffixToProtein, read_bitvec_mapping};

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    #[test]
    fn test_search_bitvec() {
        let u8_text = &build_text();
        let index = BitVecSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_bitvec_roundtrip() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let mut cursor = Cursor::new(&buf[1..]);
        let restored = read_bitvec_mapping(&mut cursor).unwrap();
        let reference = BitVecSuffixToProtein::new(&text);
        for i in 0..text.len() as i64 {
            assert_eq!(reference.suffix_to_protein(i), restored.suffix_to_protein(i));
        }
    }
}
