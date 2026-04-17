use memmap2::Mmap;

/// Reads a u64 value in little-endian byte order from the given mmap at the given byte offset.
pub(super) fn read_u64_le(mmap: &Mmap, byte_offset: usize) -> u64 {
    let bytes: [u8; 8] = mmap[byte_offset..byte_offset + 8].try_into().unwrap();
    u64::from_le_bytes(bytes)
}

/// Streaming sequential iterator over a contiguous range of a compressed or uncompressed
/// mmap-backed suffix array.
///
/// Keeps `current_word` and `next_word` in local variables (register-allocated by the
/// compiler) so that a new mmap read only occurs when crossing a 64-bit block boundary —
/// roughly once per 1.6 entries for a 40-bit SA, vs 1–2 reads per entry with `get_mmap`.
pub(crate) struct MmapSaRangeIter<'a> {
    mmap: &'a Mmap,
    data_offset: usize,
    bits_per_value: usize,
    mask: u64,
    current_word: u64, // u64 block containing the next value to yield
    next_word: u64,    // u64 block after current_word (pre-loaded)
    block_idx: usize,  // index of current_word within the data section
    bit_off: usize,    // bit offset of next value within current_word (0..64)
    remaining: usize,  // entries left to yield
}

impl<'a> MmapSaRangeIter<'a> {
    pub(crate) fn new(
        mmap: &'a Mmap,
        data_offset: usize,
        bits_per_value: usize,
        start: usize,
        end: usize,
    ) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                mmap, data_offset, bits_per_value,
                mask: 0, current_word: 0, next_word: 0,
                block_idx: 0, bit_off: 0, remaining: 0,
            };
        }

        // (1u64 << 64) overflows; use u64::MAX for the 64-bit uncompressed case
        let mask = if bits_per_value == 64 { u64::MAX } else { (1u64 << bits_per_value) - 1 };

        let bit_pos   = start * bits_per_value;
        let block_idx = bit_pos / 64;
        let bit_off   = bit_pos % 64;

        let current_word = read_u64_le(mmap, data_offset + block_idx * 8);
        let next_off     = data_offset + (block_idx + 1) * 8;
        let next_word    = if next_off + 8 <= mmap.len() { read_u64_le(mmap, next_off) } else { 0 };

        Self { mmap, data_offset, bits_per_value, mask, current_word, next_word, block_idx, bit_off, remaining }
    }
}

impl Iterator for MmapSaRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 { return None; }
        self.remaining -= 1;

        let val = if self.bit_off + self.bits_per_value <= 64 {
            // Value fits entirely within current_word
            (self.current_word >> (64 - self.bit_off - self.bits_per_value)) & self.mask
        } else {
            // Value spans current_word and next_word
            let end_off = (self.bit_off + self.bits_per_value) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & self.mask
        };

        // Advance bit cursor; load next word from mmap only on block-boundary crossing
        self.bit_off += self.bits_per_value;
        if self.bit_off >= 64 {
            self.bit_off   -= 64;
            self.block_idx += 1;
            self.current_word = self.next_word;
            let next_off = self.data_offset + (self.block_idx + 1) * 8;
            self.next_word = if next_off + 8 <= self.mmap.len() {
                read_u64_le(self.mmap, next_off)
            } else { 0 };
        }

        Some(val as i64)
    }
}

/// Returns the suffix array value at the given index from a memory-mapped file.
pub(super) fn get_mmap(mmap: &Mmap, data_offset: usize, bits_per_value: usize, index: usize) -> i64 {
    if bits_per_value == 64 {
        let offset = data_offset + index * 8;
        let bytes: [u8; 8] = mmap[offset..offset + 8].try_into().unwrap();
        i64::from_le_bytes(bytes)
    } else {
        let mask: u64 = (1u64 << bits_per_value) - 1;
        let bit_offset = index * bits_per_value;
        let start_block = bit_offset / 64;
        let start_block_offset = bit_offset % 64;
        let block_byte_offset = data_offset + start_block * 8;
        let start_val = read_u64_le(mmap, block_byte_offset);
        if start_block_offset + bits_per_value <= 64 {
            ((start_val >> (64 - start_block_offset - bits_per_value)) & mask) as i64
        } else {
            let end_block_offset = (index + 1) * bits_per_value % 64;
            let end_val = read_u64_le(mmap, block_byte_offset + 8);
            let a = start_val << end_block_offset;
            let b = end_val >> (64 - end_block_offset);
            ((a | b) & mask) as i64
        }
    }
}
