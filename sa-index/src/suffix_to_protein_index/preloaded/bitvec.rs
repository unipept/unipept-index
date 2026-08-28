use std::{
    error::Error,
    io::{Read, Write}
};

use binary_traits::WriteBinary;
use protein_metadata::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use protein_text::ProteinTextBackend;

use super::super::SuffixToProteinMappingBackend;
use crate::Nullable;

/// One 16-byte rank cell, covering 512 bits (8 words). Laid out exactly as the file stores it, so
/// reading and writing are copies rather than conversions.
///
/// * `level1` — the cumulative count of set bits before this superblock.
/// * `packed_level2` — seven 9-bit sub-counts, one per word after the first, each the cumulative
///   count *within* the superblock before that word.
///
/// See the format documentation on [`WriteBinary`] below for why nine bits are enough.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Superblock {
    level1: u64,
    packed_level2: u64
}

/// Mapping that uses O(n) memory (1-2 bits per suffix) with n the size of the input text, with retrieval
/// of the protein in O(1)
///
/// Owns its rank structure — the raw bits and the two-level counts — rather than delegating to a
/// general-purpose succinct library. Two things follow from that, and both are the reason it does:
/// the counts the file already carries are read instead of recomputed at load, and
/// [`Self::prefetch_for_suffix`] can name the two addresses a lookup will touch. A library type
/// that hides its storage can do neither. The rank algorithm is deliberately the same one
/// [`crate::suffix_to_protein_index::mmap::bitvec`] runs against a mapping, so the two backends
/// agree by construction rather than by test.
#[derive(Debug)]
pub struct BitVecSuffixToProtein {
    /// Number of text positions; `blocks` may hold up to 63 bits more, which are always zero.
    bit_len: u64,
    /// One bit per text position, bit 0 = LSB of block 0.
    blocks: Vec<u64>,
    /// One cell per 512 bits, plus a trailing cell. See [`Superblock`].
    counts: Vec<Superblock>
}

impl BitVecSuffixToProtein {
    /// Whether `position` is marked — i.e. holds a separator or the terminator, and so belongs to
    /// no protein. The caller bounds `position` against `bit_len` first.
    #[inline]
    fn get_bit(&self, position: u64) -> bool {
        let block = self.blocks[(position / 64) as usize];
        (block >> (position % 64)) & 1 == 1
    }

    /// Set bits strictly before `position` — which *is* the protein index, since exactly one
    /// marked byte closes each protein.
    ///
    /// Constant time, from three reads and a popcount: the superblock cell holding `position`
    /// gives the count before its 512-bit block (`level1`) and before its word within that block
    /// (the 9-bit `level2` sub-count, zero for word 0, which is not stored); the data word itself
    /// supplies the rest. `block << (63 - bit_offset)` discards every bit at or above `position`,
    /// so bits past `bit_len` in the final word cannot affect an in-range answer.
    ///
    /// Deliberately the same arithmetic as `mmap::bitvec`'s `rank1`, which runs it against a
    /// mapping, so the two backends agree by construction rather than by test.
    #[inline]
    fn rank1(&self, position: u64) -> u64 {
        let word_index = (position / 64) as usize;
        let word_offset = word_index % 8;
        let bit_offset = position % 64;

        let cell = self.counts[(position / 512) as usize];
        let level2 = if word_offset == 0 { 0u64 } else { (cell.packed_level2 >> ((word_offset - 1) * 9)) & 0x1FF };

        let block = self.blocks[word_index];
        let bit_count = (block << (63 - bit_offset)).count_ones() as u64;

        cell.level1 + level2 + bit_count
    }
}

impl SuffixToProteinMappingBackend for BitVecSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let position = suffix as u64;
        if position >= self.bit_len {
            return u32::NULL;
        }
        if self.get_bit(position) {
            return u32::NULL;
        }
        self.rank1(position).try_into().unwrap()
    }

    #[inline]
    fn implied_text_len(&self) -> Option<usize> {
        Some(self.bit_len as usize)
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let position = suffix as usize;
        let word_index = position / 64;
        if word_index < self.blocks.len() {
            memory_hints::prefetch::prefetch_read(&self.blocks[word_index] as *const u64);
        }
        let cell_index = position / 512;
        if cell_index < self.counts.len() {
            memory_hints::prefetch::prefetch_read(&self.counts[cell_index] as *const Superblock);
        }
    }
}

impl BitVecSuffixToProtein {
    /// Creates a new BitVecSuffixToProtein mapping
    pub fn new<T: ProteinTextBackend>(text: &T) -> Self {
        Self::from_text_parts(text.len(), |i| text.get(i))
    }

    /// Closure-based constructor — works with any text type that exposes `len()` + `get()`.
    ///
    /// Sets one bit per separator and terminator, then builds the counts over those bits.
    pub fn from_text_parts(text_len: usize, get_char: impl Fn(usize) -> u8) -> Self {
        let mut blocks = vec![0u64; text_len.div_ceil(64)];
        for i in 0..text_len {
            let c = get_char(i);
            if c == SEPARATION_CHARACTER || c == TERMINATION_CHARACTER {
                blocks[i / 64] |= 1u64 << (i % 64);
            }
        }
        let counts = build_counts(&blocks);
        Self { bit_len: text_len as u64, blocks, counts }
    }
}

/// Builds the two-level counts over `blocks`, in the layout the file stores and
/// [`BitVecSuffixToProtein::rank1`] reads.
///
/// The trailing cell (hence `+ 1`) is not redundant: a position in the last, partially-filled
/// superblock still indexes a cell, and `block_count / 8` rounds that superblock away.
fn build_counts(blocks: &[u64]) -> Vec<Superblock> {
    let block_count = blocks.len();
    let mut counts = Vec::with_capacity(block_count / 8 + 1);
    let mut level1: u64 = 0;

    for superblock in 0..block_count / 8 + 1 {
        let word_start = superblock * 8;
        let mut packed_level2: u64 = 0;
        let mut running: u64 = 0;

        for word in 0..8usize {
            if word > 0 {
                packed_level2 |= (running & 0x1FF) << ((word - 1) * 9);
            }
            if let Some(block) = blocks.get(word_start + word) {
                running += block.count_ones() as u64;
            }
        }

        counts.push(Superblock { level1, packed_level2 });
        level1 += running;
    }

    counts
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
/// Hence the constants in `build_counts`: `& 0x1FF` masks a 9-bit sub-count, `(w - 1) * 9`
/// places it, and the loop covers `w = 1..8` because word 0's sub-count is always zero and is not
/// stored.
///
/// Both backends consume this layout — this one and
/// `suffix_to_protein_index::mmap::bitvec`, which documents the same structure from the reading
/// side — so the two must be kept in step by hand.
impl WriteBinary for BitVecSuffixToProtein {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        writer.write_all(&[2u8])?;
        writer.write_all(&self.bit_len.to_le_bytes())?;
        writer.write_all(&(self.blocks.len() as u64).to_le_bytes())?;

        for block in &self.blocks {
            writer.write_all(&block.to_le_bytes())?;
        }

        for cell in &self.counts {
            writer.write_all(&cell.level1.to_le_bytes())?;
            writer.write_all(&cell.packed_level2.to_le_bytes())?;
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
/// [`InMemorySuffixToProteinMapping::read_binary`](super::InMemorySuffixToProteinMapping) consumed.
///
/// Reads whole words and whole cells, in large chunks. It used to read one `u64` per `read_exact`
/// and then `push_bit` once **per bit** — a call per position in the protein text, which at UniProt
/// scale is the most expensive thing in a preloaded startup — and then throw the file's superblocks
/// away and recompute them. Both halves of the body now land as they are stored.
pub(super) fn read_bitvec_mapping<R: Read>(reader: &mut R) -> Result<BitVecSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let bit_len = u64::from_le_bytes(buf8);
    reader.read_exact(&mut buf8)?;
    let block_count = u64::from_le_bytes(buf8) as usize;

    let needed = (bit_len as usize).div_ceil(64);
    if block_count < needed {
        // The mmap reader makes the same check with the same message; keep the two in step.
        return Err(format!(
            "Bitvec mapping declares {bit_len} bits but holds only {block_count} of the {needed} blocks that needs"
        )
        .into());
    }

    let mut buffer = vec![0u8; READ_CHUNK_BYTES];

    let mut blocks = super::try_alloc_exact(block_count, "bitvec")?;
    // Before the loop below touches a page of it, which is the only point at which the advice does
    // anything — see `memory_hints::hugepages`, and `array::preloaded::original::load_original`, which
    // does the same for the same reason. At UniProt scale this is the ~8 GB half of a ~10 GB
    // structure, and it was the one large preloaded allocation in the index that went unadvised.
    // The fill below pushes exactly `block_count` entries and the `truncate` at the end cannot
    // reallocate, so the advice stays with the allocation it was issued for.
    memory_hints::hugepages::advise_capacity(&blocks);
    let mut words_left = block_count;
    while words_left > 0 {
        let words = words_left.min(buffer.len() / 8);
        let chunk = &mut buffer[..words * 8];
        reader.read_exact(chunk)?;
        for word in chunk.as_chunks::<8>().0 {
            blocks.push(u64::from_le_bytes(*word));
        }
        words_left -= words;
    }

    let sb_count = block_count / 8 + 1;
    let mut counts = super::try_alloc_exact(sb_count, "bitvec superblock")?;
    // The other ~2 GB of the same structure; same argument as for `blocks` above.
    memory_hints::hugepages::advise_capacity(&counts);
    let mut cells_left = sb_count;
    while cells_left > 0 {
        let cells = cells_left.min(buffer.len() / 16);
        let chunk = &mut buffer[..cells * 16];
        reader.read_exact(chunk)?;
        for cell in chunk.as_chunks::<16>().0 {
            counts.push(Superblock {
                level1: u64::from_le_bytes(cell[0..8].try_into().unwrap()),
                packed_level2: u64::from_le_bytes(cell[8..16].try_into().unwrap())
            });
        }
        cells_left -= cells;
    }

    // The file may carry whole words past `bit_len`, and bits past it inside the last word. Neither
    // can change an in-range answer — `rank1` masks off everything at or above the position it is
    // asked about — but dropping them is what makes a re-serialised mapping byte-identical to the
    // one that was written, which is the invariant the roundtrip test checks. `truncate` is a no-op
    // on a well-formed file.
    blocks.truncate(needed);
    let tail = bit_len % 64;
    if tail != 0 {
        if let Some(last) = blocks.last_mut() {
            *last &= (1u64 << tail) - 1;
        }
    }
    counts.truncate(blocks.len() / 8 + 1);

    Ok(BitVecSuffixToProtein { bit_len, blocks, counts })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use protein_text::ProteinTextBackend;

    use super::{BitVecSuffixToProtein, read_bitvec_mapping};
    use crate::suffix_to_protein_index::test_utils::{
        assert_agree, assert_prefetch_is_harmless, assert_sample_lookups, many_proteins_text, sample_text, to_binary
    };

    #[test]
    fn test_search_bitvec() {
        assert_sample_lookups(&BitVecSuffixToProtein::new(&sample_text()));
    }

    #[test]
    fn test_bitvec_prefetch_is_harmless() {
        assert_prefetch_is_harmless(&BitVecSuffixToProtein::new(&sample_text()));
    }

    /// The reader reads the superblock cells straight out of the file rather than recomputing them,
    /// so the second text is long enough to make that several cells rather than one — misread the
    /// cell boundaries and the counts that follow decode as garbage.
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
