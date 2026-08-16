use std::{
    error::Error,
    io::{Read, Write}
};

use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector, Rank9};
use text_compression::ProteinTextBackend;

use super::super::SuffixToProteinMappingBackend;
use crate::{Nullable, WriteBinary};

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

    /// Closure-based constructor — works with any text type that exposes `len()` + `get()`.
    ///
    /// Sets one bit per separator and terminator; `Rank9` then builds its counts over those bits.
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
/// terminator). Since exactly one such byte closes each protein, the protein index for a position
/// is the number of set bits before it, which is what makes `rank1` the lookup. Answering that in
/// constant time needs precomputed counts, stored in two levels:
///
/// * **level1** — the cumulative count before this superblock, i.e. before every one of its 8
///   words (512 bits). A full `u64`, since it can reach `bit_len`.
/// * **packed_level2** — seven 9-bit sub-counts, one per word after the first, each the
///   cumulative count within the superblock before that word. Nine bits suffice because a count
///   within a 512-bit superblock cannot exceed 512, and seven of them occupy 63 of a `u64`'s
///   bits, so the whole cell is exactly 16 bytes and one cache line holds four of them.
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

/// Bytes read per iteration below. Kept at or above `binary_traits::load_owned`'s `BufReader`
/// capacity for the reason spelled out at `bitarray::binary`'s constant of the same name: a
/// destination smaller than that capacity makes `BufReader::read` decline to bypass its own buffer,
/// and every byte gets copied twice.
const READ_CHUNK_BYTES: usize = 4 << 20;

/// Reads the body of a bitvec mapping, after the type byte
/// [`InMemorySuffixToProteinMapping::read_binary`](super::InMemorySuffixToProteinMapping) consumed,
/// and rebuilds `Rank9` from the raw bits.
///
/// Reads whole words, and pushes whole words. It used to read one `u64` per `read_exact` and then
/// `push_bit` once **per bit** — a call per position in the protein text, which at UniProt scale is
/// the most expensive thing in a preloaded startup. `push_block` appends the same word in one
/// operation, so the file's own layout (LSB of block 0 first) is the layout that goes in.
pub(super) fn read_bitvec_mapping<R: Read>(reader: &mut R) -> Result<BitVecSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let bit_len = u64::from_le_bytes(buf8);
    reader.read_exact(&mut buf8)?;
    let block_count = u64::from_le_bytes(buf8) as usize;

    let mut bits = BitVector::with_capacity(bit_len);
    let mut buffer = vec![0u8; READ_CHUNK_BYTES];

    let mut words_left = block_count;
    while words_left > 0 {
        let words = words_left.min(buffer.len() / 8);
        let chunk = &mut buffer[..words * 8];
        reader.read_exact(chunk)?;
        for word in chunk.chunks_exact(8) {
            bits.push_block(u64::from_le_bytes(word.try_into().unwrap()));
        }
        words_left -= words;
    }

    // `push_block` appends 64 bits at a time, so the vector now holds `block_count * 64` bits while
    // the file declared `bit_len` — up to 63 more than there are text positions. `truncate` drops
    // them *and* zeroes the tail of the final word (`clear_extra_bits`), which is what makes this
    // byte-for-byte the vector the old per-bit loop produced. `Rank9`'s counts and
    // `suffix_to_protein`'s `rank1` both read `bit_len`, so an over-long vector is not cosmetic.
    bits.truncate(bit_len);

    // Read and discard the superblock array the `WriteBinary` impl above emitted: `Rank9::new`
    // recomputes those counts from the raw bits. Only the mmap reader consumes them. Skipped in the
    // same large chunks rather than 16 bytes at a time: at two bytes per word it is a quarter of the
    // body's size, so a per-cell loop here would undo most of what the loop above just bought.
    let mut skip_left = (block_count / 8 + 1) * 16;
    while skip_left > 0 {
        let bytes = skip_left.min(buffer.len());
        reader.read_exact(&mut buffer[..bytes])?;
        skip_left -= bytes;
    }

    Ok(BitVecSuffixToProtein { rank: Rank9::new(bits) })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use text_compression::ProteinTextBackend;

    use super::{BitVecSuffixToProtein, read_bitvec_mapping};
    use crate::suffix_to_protein_index::test_utils::{
        assert_agree, assert_sample_lookups, many_proteins_text, sample_text, to_binary
    };

    #[test]
    fn test_search_bitvec() {
        assert_sample_lookups(&BitVecSuffixToProtein::new(&sample_text()));
    }

    /// The reader rebuilds `Rank9` from the raw bits and skips the superblocks the writer emitted,
    /// so the second text is long enough to make that several cells rather than one — skip the
    /// wrong number of them and the bits that follow decode as garbage.
    ///
    /// The last two texts straddle the word boundary the reader's `truncate` exists for.
    /// `many_proteins_text(n, l)` is `n * (l + 1)` positions long, so 32x8 = 256 bits fills its
    /// final word exactly (nothing to truncate) while 33x8 = 264 leaves 8 bits of it in use (56 to
    /// drop, and to zero). Re-serialising is what makes that a byte-level check rather than a
    /// behavioural one: a stale bit above `bit_len` cannot change any `rank1` answer below it, so
    /// `assert_agree` alone would pass with a vector 63 bits too long.
    #[test]
    fn test_bitvec_roundtrip() {
        for text in [sample_text(), many_proteins_text(300, 5), many_proteins_text(32, 7), many_proteins_text(33, 7)] {
            let buf = to_binary(BitVecSuffixToProtein::new(&text));
            assert_eq!(buf[0], 2u8);
            let restored = read_bitvec_mapping(&mut Cursor::new(&buf[1..])).unwrap();
            assert_agree(&BitVecSuffixToProtein::new(&text), &restored, text.len());
            assert_eq!(to_binary(restored), buf, "the reloaded mapping is not the one that was written");
        }
    }
}
