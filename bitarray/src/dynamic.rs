use std::io::{BufRead, Result, Write};

use crate::binary::{self, Binary};

// ── DynBitArray ───────────────────────────────────────────────────────────────

/// A bit array whose bits-per-value is determined at runtime.
pub struct DynBitArray {
    data: Vec<u64>,
    mask: u64,
    len: usize,
    bits_per_value: usize,
}

impl DynBitArray {
    pub fn with_capacity(capacity: usize, bits_per_value: usize) -> Self {
        let extra = if (capacity * bits_per_value).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * bits_per_value / 64 + extra],
            mask: if bits_per_value == 64 { u64::MAX } else { (1 << bits_per_value) - 1 },
            len: capacity,
            bits_per_value,
        }
    }

    /// Best-effort request for transparent huge pages over the backing data (see
    /// `crate::advise_hugepages`). Env-gated; no-op off Linux or when unset.
    #[inline]
    pub fn advise_hugepages(&self) {
        crate::advise_hugepages(&self.data);
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

    pub fn bits_per_value(&self) -> usize { self.bits_per_value }
    pub fn len(&self) -> usize { self.len }
    pub fn is_empty(&self) -> bool { self.len == 0 }
    pub fn clear(&mut self) { self.data.iter_mut().for_each(|x| *x = 0); }

    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    pub fn iter_range(&self, start: usize, end: usize) -> DynBitArrayRangeIter<'_> {
        DynBitArrayRangeIter::new(&self.data, self.bits_per_value, self.mask, start, end)
    }
}

impl Binary for DynBitArray {
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()> {
        binary::write_words(&self.data, writer)
    }

    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()> {
        binary::read_words_into(&mut self.data, reader)
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

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_with_capacity() {
        let ba = DynBitArray::with_capacity(4, 40);
        assert_eq!(ba.data, vec![0, 0, 0]);
        assert_eq!(ba.mask, 0xff_ffff_ffff);
        assert_eq!(ba.len, 4);
    }

    #[test]
    fn test_get() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(ba.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(ba.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(ba.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(ba.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_set() {
        let mut ba = DynBitArray::with_capacity(4, 40);

        ba.set(0, 0b0001110011111010110001000111111100110010_u64);
        ba.set(1, 0b1100001001010010011000010100110111001001_u64);
        ba.set(2, 0b1111001101001101101101101011101001010001_u64);
        ba.set(3, 0b0000100010010001010001001110101110011100_u64);

        assert_eq!(ba.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_bits_per_value() {
        assert_eq!(DynBitArray::with_capacity(4, 40).bits_per_value(), 40);
    }

    #[test]
    fn test_len_and_empty() {
        assert_eq!(DynBitArray::with_capacity(4, 40).len(), 4);
        assert!(DynBitArray::with_capacity(0, 40).is_empty());
        assert!(!DynBitArray::with_capacity(4, 40).is_empty());
    }

    #[test]
    fn test_clear() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];
        ba.clear();
        assert_eq!(ba.data, vec![0, 0, 0]);
    }

    // ── Binary impl ───────────────────────────────────────────────────────────

    #[test]
    fn test_write_binary() {
        let mut ba = DynBitArray::with_capacity(4, 40);
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
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.read_binary(&buf[..]).unwrap();

        assert_eq!(ba.get(0), 0x1234567890);
        assert_eq!(ba.get(1), 0xabcdef0123);
        assert_eq!(ba.get(2), 0x4567890abc);
        assert_eq!(ba.get(3), 0xdef0123456);
    }

    // ── DynBitArrayRangeIter ──────────────────────────────────────────────────

    fn collect(ba: &DynBitArray, start: usize, end: usize) -> Vec<i64> {
        ba.iter_range(start, end).collect()
    }

    fn expected(ba: &DynBitArray, start: usize, end: usize) -> Vec<i64> {
        (start..end).map(|i| ba.get(i) as i64).collect()
    }

    #[test]
    fn test_iter_range_empty() {
        let ba = DynBitArray::with_capacity(8, 32);
        assert!(collect(&ba, 3, 3).is_empty());
        assert!(collect(&ba, 5, 3).is_empty());
    }

    #[test]
    fn test_iter_range_single_entry() {
        let mut ba = DynBitArray::with_capacity(4, 40);
        ba.set(2, 0xABCDEF1234_u64);
        assert_eq!(collect(&ba, 2, 3), vec![0xABCDEF1234_i64]);
    }

    #[test]
    fn test_iter_range_mid_block_start() {
        let values: Vec<u64> = (0..8).map(|i| i * 111 + 7).collect();
        let mut ba = DynBitArray::with_capacity(8, 32);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 1, 6), expected(&ba, 1, 6));
    }

    #[test]
    fn test_iter_range_crosses_block_boundary() {
        let values: Vec<u64> = (0..16).map(|i| i as u64 * 0x100000001 + 3).collect();
        let mut ba = DynBitArray::with_capacity(16, 40);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 0, 16), expected(&ba, 0, 16));
        assert_eq!(collect(&ba, 3, 13), expected(&ba, 3, 13));
    }

    #[test]
    fn test_iter_range_bits_per_value_64() {
        let values: Vec<u64> = (0..8).map(|i| i as u64 * 0xDEAD_BEEF + 1).collect();
        let mut ba = DynBitArray::with_capacity(8, 64);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }
        assert_eq!(collect(&ba, 0, 8), expected(&ba, 0, 8));
        assert_eq!(collect(&ba, 2, 6), expected(&ba, 2, 6));
    }

    #[test]
    fn test_iter_range_bits_per_value_1() {
        let mut ba = DynBitArray::with_capacity(128, 1);
        for i in (0..128).step_by(3) { ba.set(i, 1); }
        assert_eq!(collect(&ba, 0, 128), expected(&ba, 0, 128));
        assert_eq!(collect(&ba, 60, 70), expected(&ba, 60, 70));
    }

    #[test]
    fn test_iter_range_exact_size() {
        let mut ba = DynBitArray::with_capacity(10, 40);
        for i in 0..10 { ba.set(i, i as u64 * 99); }
        assert_eq!(ba.iter_range(2, 8).len(), 6);
    }
}
