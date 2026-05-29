//! This module contains the `BitArray` and `DynBitArray` structs and their associated methods.

mod binary;

use std::{
    cmp::max,
    io::{Result, Write}
};

/// Re-export the `Binary` trait.
pub use binary::Binary;

// ── DynBitArray ───────────────────────────────────────────────────────────────

/// A bit array whose bits-per-value is determined at runtime.
pub struct DynBitArray {
    data: Vec<u64>,
    mask: u64,
    len: usize,
    bits_per_value: usize
}

impl DynBitArray {
    pub fn with_capacity(capacity: usize, bits_per_value: usize) -> Self {
        let extra = if (capacity * bits_per_value).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * bits_per_value / 64 + extra],
            mask: if bits_per_value == 64 { u64::MAX } else { (1 << bits_per_value) - 1 },
            len: capacity,
            bits_per_value
        }
    }

    #[inline]
    pub fn get(&self, index: usize) -> u64 {
        let start_block = index * self.bits_per_value / 64;
        let start_block_offset = index * self.bits_per_value % 64;

        if start_block_offset + self.bits_per_value <= 64 {
            return (self.data[start_block] >> (64 - start_block_offset - self.bits_per_value)) & self.mask;
        }

        let end_block = (index + 1) * self.bits_per_value / 64;
        let end_block_offset = (index + 1) * self.bits_per_value % 64;

        let a = self.data[start_block] << end_block_offset;
        let b = self.data[end_block] >> (64 - end_block_offset);

        (a | b) & self.mask
    }

    pub fn set(&mut self, index: usize, value: u64) {
        let start_block = index * self.bits_per_value / 64;
        let start_block_offset = index * self.bits_per_value % 64;

        if start_block_offset + self.bits_per_value <= 64 {
            self.data[start_block] &= !(self.mask << (64 - start_block_offset - self.bits_per_value));
            self.data[start_block] |= value << (64 - start_block_offset - self.bits_per_value);
            return;
        }

        let end_block = (index + 1) * self.bits_per_value / 64;
        let end_block_offset = (index + 1) * self.bits_per_value % 64;

        self.data[start_block] &= !(self.mask >> start_block_offset);
        self.data[start_block] |= value >> end_block_offset;

        self.data[end_block] &= !(self.mask << (64 - end_block_offset));
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    pub fn bits_per_value(&self) -> usize {
        self.bits_per_value
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|x| *x = 0);
    }

    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    pub fn iter_range(&self, start: usize, end: usize) -> DynBitArrayRangeIter<'_> {
        DynBitArrayRangeIter::new(&self.data, self.bits_per_value, self.mask, start, end)
    }
}

// ── DynBitArrayRangeIter ──────────────────────────────────────────────────────

pub struct DynBitArrayRangeIter<'a> {
    data: &'a [u64],
    bits_per_value: usize,
    mask: u64,
    current_word: u64,
    next_word: u64,
    block_idx: usize,
    bit_off: usize,
    remaining: usize,
}

impl<'a> DynBitArrayRangeIter<'a> {
    fn new(data: &'a [u64], bits_per_value: usize, mask: u64, start: usize, end: usize) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                data, bits_per_value, mask,
                current_word: 0, next_word: 0,
                block_idx: 0, bit_off: 0, remaining: 0,
            };
        }

        let bit_pos   = start * bits_per_value;
        let block_idx = bit_pos / 64;
        let bit_off   = bit_pos % 64;

        let current_word = data[block_idx];
        let next_word    = if block_idx + 1 < data.len() { data[block_idx + 1] } else { 0 };

        Self { data, bits_per_value, mask, current_word, next_word, block_idx, bit_off, remaining }
    }
}

impl Iterator for DynBitArrayRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 { return None; }
        self.remaining -= 1;

        let val = if self.bit_off + self.bits_per_value <= 64 {
            (self.current_word >> (64 - self.bit_off - self.bits_per_value)) & self.mask
        } else {
            let end_off = (self.bit_off + self.bits_per_value) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & self.mask
        };

        self.bit_off += self.bits_per_value;
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

impl ExactSizeIterator for DynBitArrayRangeIter<'_> {}

// ── BitArray<const BITS> ──────────────────────────────────────────────────────

/// A bit array whose bits-per-value is fixed at compile time.
pub struct BitArray<const BITS: usize> {
    data: Vec<u64>,
    len: usize,
}

impl<const BITS: usize> BitArray<BITS> {
    pub fn with_capacity(capacity: usize) -> Self {
        let extra = if (capacity * BITS).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * BITS / 64 + extra],
            len: capacity,
        }
    }

    #[inline]
    pub fn get(&self, index: usize) -> u64 {
        let mask: u64 = u64::MAX >> (64 - BITS);
        let bit_offset = index * BITS;
        let start_block = bit_offset / 64;
        let start_bit = bit_offset % 64;
        if start_bit + BITS <= 64 {
            (self.data[start_block] >> (64 - start_bit - BITS)) & mask
        } else {
            let end_bit = (index + 1) * BITS % 64;
            ((self.data[start_block] << end_bit) | (self.data[start_block + 1] >> (64 - end_bit))) & mask
        }
    }

    pub fn set(&mut self, index: usize, value: u64) {
        let mask: u64 = u64::MAX >> (64 - BITS);
        let start_block = index * BITS / 64;
        let start_block_offset = index * BITS % 64;

        if start_block_offset + BITS <= 64 {
            self.data[start_block] &= !(mask << (64 - start_block_offset - BITS));
            self.data[start_block] |= value << (64 - start_block_offset - BITS);
            return;
        }

        let end_block = (index + 1) * BITS / 64;
        let end_block_offset = (index + 1) * BITS % 64;

        self.data[start_block] &= !(mask >> start_block_offset);
        self.data[start_block] |= value >> end_block_offset;

        self.data[end_block] &= !(mask << (64 - end_block_offset));
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    pub fn bits_per_value(&self) -> usize {
        BITS
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|x| *x = 0);
    }

    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    pub fn iter_range(&self, start: usize, end: usize) -> BitArrayRangeIter<'_, BITS> {
        BitArrayRangeIter::new(&self.data, start, end)
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

        let mask: u64 = u64::MAX >> (64 - BITS);
        let val = if self.bit_off + BITS <= 64 {
            (self.current_word >> (64 - self.bit_off - BITS)) & mask
        } else {
            let end_off = (self.bit_off + BITS) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & mask
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

// ── data_to_writer ────────────────────────────────────────────────────────────

/// Writes packed bit data to a writer in chunks, minimising peak memory.
pub fn data_to_writer(
    data: Vec<i64>,
    bits_per_value: usize,
    max_capacity: usize,
    writer: &mut impl Write
) -> Result<()> {
    let greates_common_divisor = gcd(bits_per_value, 64);
    let capacity = max(greates_common_divisor, max_capacity / greates_common_divisor * greates_common_divisor);

    if data.len() <= capacity {
        let mut bitarray = DynBitArray::with_capacity(data.len(), bits_per_value);
        for (i, &value) in data.iter().enumerate() {
            bitarray.set(i, value as u64);
        }
        bitarray.write_binary(writer)?;
        return Ok(());
    }

    let mut bitarray = DynBitArray::with_capacity(capacity, bits_per_value);
    let chunks = data.chunks_exact(capacity);
    let remainder = chunks.remainder();

    for chunk in chunks {
        for (i, &value) in chunk.iter().enumerate() {
            bitarray.set(i, value as u64);
        }
        bitarray.write_binary(writer)?;
        bitarray.clear();
    }

    bitarray = DynBitArray::with_capacity(remainder.len(), bits_per_value);
    for (i, &value) in remainder.iter().enumerate() {
        bitarray.set(i, value as u64);
    }
    bitarray.write_binary(writer)?;

    Ok(())
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        if b < a {
            std::mem::swap(&mut b, &mut a);
        }
        b %= a;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── DynBitArray tests ─────────────────────────────────────────────────────

    #[test]
    fn test_dynbitarray_with_capacity() {
        let ba = DynBitArray::with_capacity(4, 40);
        assert_eq!(ba.data, vec![0, 0, 0]);
        assert_eq!(ba.mask, 0xff_ffff_ffff);
        assert_eq!(ba.len, 4);
    }

    #[test]
    fn test_dynbitarray_get() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(ba.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(ba.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(ba.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(ba.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_dynbitarray_set() {
        let mut ba = DynBitArray::with_capacity(4, 40);

        ba.set(0, 0b0001110011111010110001000111111100110010_u64);
        ba.set(1, 0b1100001001010010011000010100110111001001_u64);
        ba.set(2, 0b1111001101001101101101101011101001010001_u64);
        ba.set(3, 0b0000100010010001010001001110101110011100_u64);

        assert_eq!(ba.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_dynbitarray_bits_per_value() {
        let ba = DynBitArray::with_capacity(4, 40);
        assert_eq!(ba.bits_per_value(), 40);
    }

    #[test]
    fn test_dynbitarray_len() {
        let ba = DynBitArray::with_capacity(4, 40);
        assert_eq!(ba.len(), 4);
    }

    #[test]
    fn test_dynbitarray_is_empty() {
        let ba = DynBitArray::with_capacity(0, 40);
        assert!(ba.is_empty());
    }

    #[test]
    fn test_dynbitarray_is_not_empty() {
        let ba = DynBitArray::with_capacity(4, 40);
        assert!(!ba.is_empty());
    }

    #[test]
    fn test_dynbitarray_clear() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];
        ba.clear();
        assert_eq!(ba.data, vec![0, 0, 0]);
    }

    // ── BitArray<const BITS> tests ────────────────────────────────────────────

    #[test]
    fn test_bitarray_with_capacity() {
        let ba = BitArray::<40>::with_capacity(4);
        assert_eq!(ba.data, vec![0, 0, 0]);
        assert_eq!(ba.len, 4);
    }

    #[test]
    fn test_bitarray_get() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(ba.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(ba.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(ba.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(ba.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_bitarray_set() {
        let mut ba = BitArray::<40>::with_capacity(4);

        ba.set(0, 0b0001110011111010110001000111111100110010_u64);
        ba.set(1, 0b1100001001010010011000010100110111001001_u64);
        ba.set(2, 0b1111001101001101101101101011101001010001_u64);
        ba.set(3, 0b0000100010010001010001001110101110011100_u64);

        assert_eq!(ba.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_bitarray_bits_per_value() {
        let ba = BitArray::<40>::with_capacity(4);
        assert_eq!(ba.bits_per_value(), 40);
    }

    #[test]
    fn test_bitarray_len_and_empty() {
        assert_eq!(BitArray::<40>::with_capacity(4).len(), 4);
        assert!(BitArray::<40>::with_capacity(0).is_empty());
        assert!(!BitArray::<40>::with_capacity(4).is_empty());
    }

    #[test]
    fn test_bitarray_clear() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];
        ba.clear();
        assert_eq!(ba.data, vec![0, 0, 0]);
    }

    // ── data_to_writer tests ──────────────────────────────────────────────────

    #[test]
    fn test_data_to_writer_no_chunks_needed() {
        let data = vec![0x1234567890, 0xabcdef0123, 0x4567890abc, 0xdef0123456];
        let mut writer = Vec::new();

        data_to_writer(data, 40, 2, &mut writer).unwrap();

        assert_eq!(writer, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01, 0x00, 0x00,
            0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
        ]);
    }

    #[test]
    fn test_data_to_writer_chunks_needed_no_remainder() {
        let data = vec![
            0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999,
            0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222,
            0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb,
            0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333, 0x44444444,
            0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd,
            0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666,
            0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff,
            0x00000000,
        ];
        let mut writer = Vec::new();

        data_to_writer(data, 32, 8, &mut writer).unwrap();

        assert_eq!(writer, vec![
            0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66,
            0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa,
            0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd,
            0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11,
            0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88,
            0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc,
            0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff,
            0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33,
            0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa,
            0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee,
            0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11,
            0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55,
            0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc,
            0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00,
            0xff, 0xff, 0xff, 0xff
        ]);
    }

    #[test]
    fn test_data_to_writer_chunks_needed_plus_remainder() {
        let data = vec![
            0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999,
            0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222,
            0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb,
            0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333,
        ];
        let mut writer = Vec::new();

        data_to_writer(data, 32, 8, &mut writer).unwrap();

        assert_eq!(writer, vec![
            0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66,
            0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa,
            0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd,
            0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11,
            0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88,
            0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc,
            0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff,
            0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x00, 0x00, 0x00, 0x00, 0x33, 0x33, 0x33, 0x33
        ]);
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(40, 64), 8);
        assert_eq!(gcd(64, 40), 8);
        assert_eq!(gcd(64, 64), 64);
        assert_eq!(gcd(32, 64), 32);
    }

    // ── DynBitArrayRangeIter tests ────────────────────────────────────────────

    fn collect_range_dyn(ba: &DynBitArray, start: usize, end: usize) -> Vec<i64> {
        ba.iter_range(start, end).collect()
    }

    fn expected_range_dyn(ba: &DynBitArray, start: usize, end: usize) -> Vec<i64> {
        (start..end).map(|i| ba.get(i) as i64).collect()
    }

    #[test]
    fn test_dyn_iter_range_empty() {
        let ba = DynBitArray::with_capacity(8, 32);
        assert!(collect_range_dyn(&ba, 3, 3).is_empty());
        assert!(collect_range_dyn(&ba, 5, 3).is_empty());
    }

    #[test]
    fn test_dyn_iter_range_single_entry() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.set(2, 0xABCDEF1234_u64);
        assert_eq!(collect_range_dyn(&ba, 2, 3), vec![0xABCDEF1234_i64]);
    }

    #[test]
    fn test_dyn_iter_range_mid_block_start() {
        let values: Vec<u64> = (0..8).map(|i| i * 111 + 7).collect();
        let mut ba = DynBitArray::with_capacity(8, 32);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range_dyn(&ba, 1, 6), expected_range_dyn(&ba, 1, 6));
    }

    #[test]
    fn test_dyn_iter_range_crosses_block_boundary() {
        let values: Vec<u64> = (0..16).map(|i| i as u64 * 0x100000001 + 3).collect();
        let mut ba = DynBitArray::with_capacity(16, 40);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range_dyn(&ba, 0, 16), expected_range_dyn(&ba, 0, 16));
        assert_eq!(collect_range_dyn(&ba, 3, 13), expected_range_dyn(&ba, 3, 13));
    }

    #[test]
    fn test_dyn_iter_range_bits_per_value_64() {
        let values: Vec<u64> = (0..8).map(|i| i as u64 * 0xDEAD_BEEF + 1).collect();
        let mut ba = DynBitArray::with_capacity(8, 64);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range_dyn(&ba, 0, 8), expected_range_dyn(&ba, 0, 8));
        assert_eq!(collect_range_dyn(&ba, 2, 6), expected_range_dyn(&ba, 2, 6));
    }

    #[test]
    fn test_dyn_iter_range_bits_per_value_1() {
        let mut ba = DynBitArray::with_capacity(128, 1);
        for i in (0..128).step_by(3) { ba.set(i, 1); }
        assert_eq!(collect_range_dyn(&ba, 0, 128), expected_range_dyn(&ba, 0, 128));
        assert_eq!(collect_range_dyn(&ba, 60, 70), expected_range_dyn(&ba, 60, 70));
    }

    #[test]
    fn test_dyn_iter_range_exact_size_iterator() {
        let mut ba = DynBitArray::with_capacity(10, 40);
        for i in 0..10 { ba.set(i, i as u64 * 99); }
        let iter = ba.iter_range(2, 8);
        assert_eq!(iter.len(), 6);
    }

    // ── BitArrayRangeIter<const BITS> tests ───────────────────────────────────

    fn collect_range<const BITS: usize>(ba: &BitArray<BITS>, start: usize, end: usize) -> Vec<i64> {
        ba.iter_range(start, end).collect()
    }

    fn expected_range<const BITS: usize>(ba: &BitArray<BITS>, start: usize, end: usize) -> Vec<i64> {
        (start..end).map(|i| ba.get(i) as i64).collect()
    }

    #[test]
    fn test_iter_range_empty() {
        let ba = BitArray::<32>::with_capacity(8);
        assert!(collect_range(&ba, 3, 3).is_empty());
        assert!(collect_range(&ba, 5, 3).is_empty());
    }

    #[test]
    fn test_iter_range_single_entry() {
        let mut ba = BitArray::<40>::with_capacity(4);
        ba.set(2, 0xABCDEF1234_u64);
        assert_eq!(collect_range(&ba, 2, 3), vec![0xABCDEF1234_i64]);
    }

    #[test]
    fn test_iter_range_mid_block_start() {
        let values: Vec<u64> = (0..8).map(|i| i * 111 + 7).collect();
        let mut ba = BitArray::<32>::with_capacity(8);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range(&ba, 1, 6), expected_range(&ba, 1, 6));
    }

    #[test]
    fn test_iter_range_crosses_block_boundary() {
        let values: Vec<u64> = (0..16).map(|i| i as u64 * 0x100000001 + 3).collect();
        let mut ba = BitArray::<40>::with_capacity(16);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range(&ba, 0, 16), expected_range(&ba, 0, 16));
        assert_eq!(collect_range(&ba, 3, 13), expected_range(&ba, 3, 13));
    }

    #[test]
    fn test_iter_range_bits_per_value_64() {
        let values: Vec<u64> = (0..8).map(|i| i as u64 * 0xDEAD_BEEF + 1).collect();
        let mut ba = BitArray::<64>::with_capacity(8);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect_range(&ba, 0, 8), expected_range(&ba, 0, 8));
        assert_eq!(collect_range(&ba, 2, 6), expected_range(&ba, 2, 6));
    }

    #[test]
    fn test_iter_range_bits_per_value_1() {
        let mut ba = BitArray::<1>::with_capacity(128);
        for i in (0..128).step_by(3) { ba.set(i, 1); }
        assert_eq!(collect_range(&ba, 0, 128), expected_range(&ba, 0, 128));
        assert_eq!(collect_range(&ba, 60, 70), expected_range(&ba, 60, 70));
    }

    #[test]
    fn test_iter_range_exact_size_iterator() {
        let mut ba = BitArray::<40>::with_capacity(10);
        for i in 0..10 { ba.set(i, i as u64 * 99); }
        let iter = ba.iter_range(2, 8);
        assert_eq!(iter.len(), 6);
    }
}
