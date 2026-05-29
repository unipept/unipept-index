use std::io::{BufRead, Result, Write};

use crate::binary::{self, Binary};

// ── BitArray<const BITS> ──────────────────────────────────────────────────────

/// A bit array whose bits-per-value is fixed at compile time.
pub struct BitArray<const BITS: usize> {
    data: Vec<u64>,
    len: usize,
}

impl<const BITS: usize> BitArray<BITS> {
    const MASK: u64 = u64::MAX >> (64 - BITS);

    pub fn with_capacity(capacity: usize) -> Self {
        let extra = if (capacity * BITS).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * BITS / 64 + extra],
            len: capacity,
        }
    }

    #[inline]
    pub fn get(&self, index: usize) -> u64 {
        let bit_offset = index * BITS;
        let start_block = bit_offset / 64;
        let start_bit = bit_offset % 64;
        if start_bit + BITS <= 64 {
            (self.data[start_block] >> (64 - start_bit - BITS)) & Self::MASK
        } else {
            let end_bit = (index + 1) * BITS % 64;
            ((self.data[start_block] << end_bit) | (self.data[start_block + 1] >> (64 - end_bit))) & Self::MASK
        }
    }

    pub fn set(&mut self, index: usize, value: u64) {
        let start_block = index * BITS / 64;
        let start_block_offset = index * BITS % 64;

        if start_block_offset + BITS <= 64 {
            self.data[start_block] &= !(Self::MASK << (64 - start_block_offset - BITS));
            self.data[start_block] |= value << (64 - start_block_offset - BITS);
            return;
        }

        let end_block = (index + 1) * BITS / 64;
        let end_block_offset = (index + 1) * BITS % 64;

        self.data[start_block] &= !(Self::MASK >> start_block_offset);
        self.data[start_block] |= value >> end_block_offset;

        self.data[end_block] &= !(Self::MASK << (64 - end_block_offset));
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    pub fn bits_per_value(&self) -> usize { BITS }
    pub fn len(&self) -> usize { self.len }
    pub fn is_empty(&self) -> bool { self.len == 0 }
    pub fn clear(&mut self) { self.data.iter_mut().for_each(|x| *x = 0); }

    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    pub fn iter_range(&self, start: usize, end: usize) -> BitArrayRangeIter<'_, BITS> {
        BitArrayRangeIter::new(&self.data, start, end)
    }
}

impl<const BITS: usize> Binary for BitArray<BITS> {
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()> {
        binary::write_words(&self.data, writer)
    }

    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()> {
        binary::read_words_into(&mut self.data, reader)
    }
}

// ── BitArrayRangeIter<const BITS> ─────────────────────────────────────────────

pub struct BitArrayRangeIter<'a, const BITS: usize> {
    data: &'a [u64],
    current_word: u64,
    next_word: u64,
    block_idx: usize,
    bit_off: usize,
    remaining: usize,
}

impl<'a, const BITS: usize> BitArrayRangeIter<'a, BITS> {
    const MASK: u64 = u64::MAX >> (64 - BITS);

    fn new(data: &'a [u64], start: usize, end: usize) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                data,
                current_word: 0, next_word: 0,
                block_idx: 0, bit_off: 0, remaining: 0,
            };
        }

        let bit_pos   = start * BITS;
        let block_idx = bit_pos / 64;
        let bit_off   = bit_pos % 64;

        let current_word = data[block_idx];
        let next_word    = if block_idx + 1 < data.len() { data[block_idx + 1] } else { 0 };

        Self { data, current_word, next_word, block_idx, bit_off, remaining }
    }
}

impl<const BITS: usize> Iterator for BitArrayRangeIter<'_, BITS> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 { return None; }
        self.remaining -= 1;

        let val = if self.bit_off + BITS <= 64 {
            (self.current_word >> (64 - self.bit_off - BITS)) & Self::MASK
        } else {
            let end_off = (self.bit_off + BITS) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & Self::MASK
        };

        self.bit_off += BITS;
        if self.bit_off >= 64 {
            self.bit_off   -= 64;
            self.block_idx += 1;
            self.current_word = self.next_word;
            self.next_word = if self.block_idx + 1 < self.data.len() {
                self.data[self.block_idx + 1]
            } else {
                0
            };
        }

        Some(val as i64)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<const BITS: usize> ExactSizeIterator for BitArrayRangeIter<'_, BITS> {}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_with_capacity() {
        let ba = BitArray::<40>::with_capacity(4);
        assert_eq!(ba.data, vec![0, 0, 0]);
        assert_eq!(ba.len, 4);
    }

    #[test]
    fn test_get() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(ba.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(ba.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(ba.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(ba.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_set() {
        let mut ba = BitArray::<40>::with_capacity(4);

        ba.set(0, 0b0001110011111010110001000111111100110010_u64);
        ba.set(1, 0b1100001001010010011000010100110111001001_u64);
        ba.set(2, 0b1111001101001101101101101011101001010001_u64);
        ba.set(3, 0b0000100010010001010001001110101110011100_u64);

        assert_eq!(ba.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_bits_per_value() {
        assert_eq!(BitArray::<40>::with_capacity(4).bits_per_value(), 40);
    }

    #[test]
    fn test_len_and_empty() {
        assert_eq!(BitArray::<40>::with_capacity(4).len(), 4);
        assert!(BitArray::<40>::with_capacity(0).is_empty());
        assert!(!BitArray::<40>::with_capacity(4).is_empty());
    }

    #[test]
    fn test_clear() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];
        ba.clear();
        assert_eq!(ba.data, vec![0, 0, 0]);
    }

    // ── Binary impl ───────────────────────────────────────────────────────────

    #[test]
    fn test_write_binary() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.set(0, 0x1234567890_u64);
        ba.set(1, 0xabcdef0123_u64);
        ba.set(2, 0x4567890abc_u64);
        ba.set(3, 0xdef0123456_u64);

        let mut buf = Vec::new();
        ba.write_binary(&mut buf).unwrap();

        assert_eq!(buf, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45,
            0x23, 0x01, 0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0,
        ]);
    }

    #[test]
    fn test_read_binary() {
        let buf = [
            0xef_u8, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67,
            0x45, 0x23, 0x01, 0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0,
        ];
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.read_binary(&buf[..]).unwrap();

        assert_eq!(ba.get(0), 0x1234567890);
        assert_eq!(ba.get(1), 0xabcdef0123);
        assert_eq!(ba.get(2), 0x4567890abc);
        assert_eq!(ba.get(3), 0xdef0123456);
    }

    // ── BitArrayRangeIter ─────────────────────────────────────────────────────

    fn collect<const BITS: usize>(ba: &BitArray<BITS>, start: usize, end: usize) -> Vec<i64> {
        ba.iter_range(start, end).collect()
    }

    fn expected<const BITS: usize>(ba: &BitArray<BITS>, start: usize, end: usize) -> Vec<i64> {
        (start..end).map(|i| ba.get(i) as i64).collect()
    }

    #[test]
    fn test_iter_range_empty() {
        let ba = BitArray::<32>::with_capacity(8);
        assert!(collect(&ba, 3, 3).is_empty());
        assert!(collect(&ba, 5, 3).is_empty());
    }

    #[test]
    fn test_iter_range_single_entry() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.set(2, 0xABCDEF1234_u64);
        assert_eq!(collect(&ba, 2, 3), vec![0xABCDEF1234_i64]);
    }

    #[test]
    fn test_iter_range_mid_block_start() {
        let values: Vec<u64> = (0..8).map(|i| i * 111 + 7).collect();
        let mut ba = BitArray::<32>::with_capacity(8);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 1, 6), expected(&ba, 1, 6));
    }

    #[test]
    fn test_iter_range_crosses_block_boundary() {
        let values: Vec<u64> = (0..16).map(|i| i as u64 * 0x100000001 + 3).collect();
        let mut ba = BitArray::<40>::with_capacity(16);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 0, 16), expected(&ba, 0, 16));
        assert_eq!(collect(&ba, 3, 13), expected(&ba, 3, 13));
    }

    #[test]
    fn test_iter_range_bits_per_value_64() {
        let values: Vec<u64> = (0..8).map(|i| i as u64 * 0xDEAD_BEEF + 1).collect();
        let mut ba = BitArray::<64>::with_capacity(8);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 0, 8), expected(&ba, 0, 8));
        assert_eq!(collect(&ba, 2, 6), expected(&ba, 2, 6));
    }

    #[test]
    fn test_iter_range_bits_per_value_1() {
        let mut ba = BitArray::<1>::with_capacity(128);
        for i in (0..128).step_by(3) { ba.set(i, 1); }
        assert_eq!(collect(&ba, 0, 128), expected(&ba, 0, 128));
        assert_eq!(collect(&ba, 60, 70), expected(&ba, 60, 70));
    }

    #[test]
    fn test_iter_range_exact_size() {
        let mut ba = BitArray::<40>::with_capacity(10);
        for i in 0..10 { ba.set(i, i as u64 * 99); }
        assert_eq!(ba.iter_range(2, 8).len(), 6);
    }
}
