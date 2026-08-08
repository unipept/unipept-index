//! Compile-time-width bit array — the fast half of the pair.
//!
//! Identical in behaviour to [`crate::dynamic::DynBitArray`], and asserted so by the shared
//! parity test. It exists separately because making the width a const generic turns every
//! quantity in [`BitArray::get`] into a compile-time constant:
//!
//! * `MASK` becomes an immediate instead of a field load,
//! * `index * BITS`, `64 - start_bit - BITS` and friends fold into a shift by a constant rather
//!   than a variable-shift instruction,
//! * and for widths that divide 64 the `start_bit + BITS <= 64` branch folds away entirely, so
//!   the straddling path is never even branched over.
//!
//! `get` is called once per residue compared during candidate validation — the innermost loop in
//! the whole index — so this matters. Prefer this type whenever the width is a property of the
//! data rather than of the file being read; see the crate docs for the choice.

use std::io::{BufRead, Result, Write};

use crate::binary::{self, Binary};

// ── BitArray<const BITS> ──────────────────────────────────────────────────────

/// A bit array whose bits-per-value is fixed at compile time.
///
/// Values are packed most-significant-bit first within each little-endian `u64`, and may straddle
/// a word boundary. `BITS` must be in `1..=64`; `BITS == 0` fails to compile at `MASK`.
pub struct BitArray<const BITS: usize> {
    data: Vec<u64>,
    len: usize,
}

impl<const BITS: usize> BitArray<BITS> {
    /// Low `BITS` bits set. Const-evaluated, so it costs no runtime work and no register.
    const MASK: u64 = u64::MAX >> (64 - BITS);

    /// Allocates room for `capacity` values, all zero.
    pub fn with_capacity(capacity: usize) -> Self {
        let extra = if (capacity * BITS).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * BITS / 64 + extra],
            len: capacity,
        }
    }

    /// Requests transparent huge pages over the backing data (no-op off Linux).
    /// See [`crate::hugepages`].
    #[inline]
    pub fn advise_hugepages(&self) {
        crate::hugepages::advise(&self.data);
    }

    /// Returns the value at `index`.
    ///
    /// `#[inline]` is load-bearing: every caller is in another crate and the workspace sets no
    /// `[profile.release]`, so without this the innermost loop of the index pays a cross-crate
    /// call per residue. See the crate docs on `prefetch` for the same argument at length.
    ///
    /// # Panics
    ///
    /// If `index` is out of bounds, via the underlying slice index.
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

    /// Writes `value` at `index`. Only the low `BITS` bits of `value` are stored.
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

    /// Bits stored per value, i.e. `BITS`. Present so callers can be generic over both
    /// implementations.
    pub fn bits_per_value(&self) -> usize { BITS }
    /// Number of values stored.
    pub fn len(&self) -> usize { self.len }
    /// Whether the array holds no values.
    pub fn is_empty(&self) -> bool { self.len == 0 }
    /// Zeroes every value without reallocating, so the buffer can be refilled.
    pub fn clear(&mut self) { self.data.iter_mut().for_each(|x| *x = 0); }

    /// Borrows the raw backing words in `start_slice..end_slice`.
    ///
    /// Exposed for prefetching: callers translate a value index to its word index and take the
    /// address of that word. Word indices are not value indices.
    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    /// Iterates values in `start..end` (half-open).
    ///
    /// This is the sequential fast path and the reason to prefer it over a `get` loop: the
    /// iterator carries the current and next words, so consecutive values that share a word cost
    /// no reload, and a straddling value already has its second word to hand. `get` re-derives
    /// and re-loads both on every call.
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

/// Constructor shim for the shared suite: the width is a const generic here.
#[cfg(test)]
macro_rules! new_bitarray {
    ($capacity:expr, $bits:literal) => { BitArray::<$bits>::with_capacity($capacity) };
}

#[cfg(test)]
bitarray_test_suite!(new_bitarray);

#[cfg(test)]
mod constant_only_tests {
    use super::*;

    /// The point of this type: `MASK` is a compile-time constant, usable in const context.
    #[test]
    fn mask_is_a_compile_time_constant() {
        const M40: u64 = BitArray::<40>::MASK;
        const M64: u64 = BitArray::<64>::MASK;
        const M1: u64 = BitArray::<1>::MASK;

        assert_eq!(M40, 0xff_ffff_ffff);
        assert_eq!(M64, u64::MAX);
        assert_eq!(M1, 1);
    }
}
