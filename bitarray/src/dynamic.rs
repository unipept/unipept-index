//! Runtime-width bit array — the flexible half of the pair.
//!
//! Identical in behaviour to [`crate::constant::BitArray`], and asserted so by the shared parity
//! test. Use this one when the width is not known until runtime, which in practice means it came
//! out of a file header: the compressed suffix array picks its width from the text length at
//! build time and records it in the file, so the reader cannot be generic over it.
//!
//! The cost of that flexibility is that `mask` is a field load rather than an immediate and every
//! shift amount is computed rather than folded — see [`crate::constant`] for what that buys. If
//! the width is a property of the data rather than of the file, prefer the const-generic type.

use std::io::{BufRead, Result, Write};

use crate::binary::{self, Binary};

/// A bit array whose bits-per-value is determined at runtime.
///
/// Values are packed most-significant-bit first within each little-endian `u64`, and may straddle
/// a word boundary — byte-for-byte the same layout as [`crate::constant::BitArray`].
pub struct DynBitArray {
    data: Vec<u64>,
    mask: u64,
    len: usize,
    bits_per_value: usize
}

impl DynBitArray {
    /// Allocates room for `capacity` values of `bits_per_value` bits each, all zero.
    ///
    /// `bits_per_value` is expected to be in `1..=64`.
    ///
    /// The huge-page advice is issued here, between reserving the allocation and zeroing it, and
    /// not by the caller once it has filled it: see [`memory_hints::hugepages`] for why that ordering is
    /// the whole point.
    pub fn with_capacity(capacity: usize, bits_per_value: usize) -> Self {
        Self::try_with_capacity(capacity, bits_per_value)
            .expect("DynBitArray::with_capacity: capacity * bits_per_value does not fit, or allocation failed")
    }

    /// Like [`Self::with_capacity`], but reports failure instead of aborting the process.
    ///
    /// Readers take the capacity straight out of an untrusted file header, where a corrupt count
    /// can ask for more memory than the machine has — or, once `capacity * bits_per_value`
    /// overflows, for an arbitrary amount unrelated to what the header said. `vec![0; n]` aborts on
    /// allocation failure, which turns a damaged index into a dead process; this returns `None` so
    /// the caller can turn it into a load error instead.
    ///
    /// Returns `None` if `capacity * bits_per_value` overflows or the allocation fails.
    pub fn try_with_capacity(capacity: usize, bits_per_value: usize) -> Option<Self> {
        let words = capacity.checked_mul(bits_per_value)?.div_ceil(64);
        let mut data: Vec<u64> = Vec::new();
        data.try_reserve_exact(words).ok()?;
        // Between the reservation and the `resize` that writes it, and not after: `resize` zeroes
        // every word, which faults the whole buffer in. Advice issued after that arrives at a
        // region that is already populated with 4 KB pages, where it buys nothing but khugepaged
        // eligibility. Here it still governs the faults `resize` is about to take. See
        // [`memory_hints::hugepages`]. `resize` cannot reallocate, since the capacity is already reserved,
        // so the advice stays with the allocation it was issued for.
        memory_hints::hugepages::advise_capacity(&data);
        data.resize(words, 0);
        Some(Self {
            data,
            mask: if bits_per_value == 64 { u64::MAX } else { (1 << bits_per_value) - 1 },
            len: capacity,
            bits_per_value
        })
    }

    /// Number of backing words the array currently holds.
    ///
    /// Exposed so a reader can check a body it has just decoded against the item count its header
    /// declared: [`Binary::read_binary`] refills the backing store with however many words the
    /// reader actually yielded, which for a truncated file is fewer than `len` implies. Nothing
    /// else relates the two, so without this check `len()` reports the declared count over a
    /// short body and every lookup past the real data panics.
    pub fn word_len(&self) -> usize {
        self.data.len()
    }

    /// Number of 64-bit words `len()` values at this width occupy — what [`Self::word_len`] must
    /// be at least, for the array to be readable through its whole declared length.
    pub fn required_words(&self) -> usize {
        (self.len * self.bits_per_value).div_ceil(64)
    }

    /// Returns the value at `index`.
    ///
    /// `#[inline]` for the same reason as [`BitArray::get`](crate::BitArray::get): every caller is
    /// in another crate and the workspace sets no `[profile.release]`, so there is no cross-crate
    /// LTO to fall back on. Unlike the const-generic type, the shift amounts and the mask are
    /// runtime values here, so inlining is all there is to win.
    ///
    /// # Panics
    ///
    /// If `index` is out of bounds, via the underlying slice index.
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

    /// Writes `value` at `index`. Only the low `bits_per_value` bits of `value` are stored.
    ///
    /// # Panics
    ///
    /// If `index` is out of bounds, via the underlying slice index.
    pub fn set(&mut self, index: usize, value: u64) {
        // Masked, because the write below ORs `value` into place: a value wider than the field
        // would otherwise spill its high bits into the *neighbouring* entry rather than being
        // rejected or truncated. Callers pass values derived from file headers and suffix arrays,
        // where a width mismatch is a real possibility, so the containment happens here.
        let value = value & self.mask;

        let start_block = index * self.bits_per_value / 64;
        let start_block_offset = index * self.bits_per_value % 64;

        if start_block_offset + self.bits_per_value <= 64 {
            self.data[start_block] &= !(self.mask << (64 - start_block_offset - self.bits_per_value));
            self.data[start_block] |= value << (64 - start_block_offset - self.bits_per_value);
            return;
        }

        let end_block = (index + 1) * self.bits_per_value / 64;
        let end_block_offset = (index + 1) * self.bits_per_value % 64;

        // The straddling half of the value occupies the low `64 - start_block_offset` bits of this
        // word, which is *not* what `self.mask >> start_block_offset` describes — for a field of 5
        // bits starting at offset 62 that expression is zero, so the old clear cleared nothing and
        // overwriting an index left the previous value's bits ORed underneath the new one.
        // `u64::MAX >> start_block_offset` is the correct set of bits, and cannot shift by 64
        // because a straddling field always starts past offset 0.
        self.data[start_block] &= !(u64::MAX >> start_block_offset);
        self.data[start_block] |= value >> end_block_offset;

        self.data[end_block] &= !(self.mask << (64 - end_block_offset));
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    /// Bits stored per value.
    pub fn bits_per_value(&self) -> usize {
        self.bits_per_value
    }
    /// Number of values stored.
    pub fn len(&self) -> usize {
        self.len
    }
    /// Whether the array holds no values.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    /// Zeroes every value without reallocating, so the buffer can be refilled.
    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|x| *x = 0);
    }

    /// Borrows the raw backing words in `start_slice..end_slice`.
    ///
    /// Exposed for prefetching: callers translate a value index to its word index and take the
    /// address of that word. Word indices are not value indices.
    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    /// Iterates values in `start..end` (half-open).
    ///
    /// The sequential fast path: the iterator carries the current and next words, so consecutive
    /// values that share a word cost no reload. See [`crate::constant::BitArray::iter_range`].
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

/// Iterator over a range of values, returned by [`DynBitArray::iter_range`].
///
/// The runtime-width counterpart of [`BitArrayRangeIter`](crate::BitArrayRangeIter), and carries
/// the same current/next word pair; the width and mask travel with it as fields rather than as
/// constants. Yields `i64` because the suffix array stores signed positions.
pub struct DynBitArrayRangeIter<'a> {
    data: &'a [u64],
    bits_per_value: usize,
    mask: u64,
    current_word: u64,
    next_word: u64,
    block_idx: usize,
    bit_off: usize,
    remaining: usize
}

impl<'a> DynBitArrayRangeIter<'a> {
    fn new(data: &'a [u64], bits_per_value: usize, mask: u64, start: usize, end: usize) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                data,
                bits_per_value,
                mask,
                current_word: 0,
                next_word: 0,
                block_idx: 0,
                bit_off: 0,
                remaining: 0
            };
        }

        let bit_pos = start * bits_per_value;
        let block_idx = bit_pos / 64;
        let bit_off = bit_pos % 64;

        let current_word = data[block_idx];
        let next_word = if block_idx + 1 < data.len() { data[block_idx + 1] } else { 0 };

        Self {
            data,
            bits_per_value,
            mask,
            current_word,
            next_word,
            block_idx,
            bit_off,
            remaining
        }
    }
}

impl Iterator for DynBitArrayRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;

        let val = if self.bit_off + self.bits_per_value <= 64 {
            (self.current_word >> (64 - self.bit_off - self.bits_per_value)) & self.mask
        } else {
            let end_off = (self.bit_off + self.bits_per_value) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & self.mask
        };

        self.bit_off += self.bits_per_value;
        if self.bit_off >= 64 {
            self.bit_off -= 64;
            self.block_idx += 1;
            self.current_word = self.next_word;
            self.next_word = if self.block_idx + 1 < self.data.len() { self.data[self.block_idx + 1] } else { 0 };
        }

        Some(val as i64)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl ExactSizeIterator for DynBitArrayRangeIter<'_> {}

/// Constructor shim for the shared suite: the width is a runtime argument here.
#[cfg(test)]
macro_rules! new_bitarray {
    ($capacity:expr, $bits:literal) => {
        DynBitArray::with_capacity($capacity, $bits)
    };
}

#[cfg(test)]
bitarray_test_suite!(new_bitarray);

#[cfg(test)]
mod dynamic_only_tests {
    use super::*;

    /// The runtime counterpart of `BitArray::MASK`, including the `1 << 64` overflow case that
    /// the const version sidesteps by construction.
    #[test]
    fn mask_is_derived_from_the_runtime_width() {
        assert_eq!(DynBitArray::with_capacity(4, 40).mask, 0xff_ffff_ffff);
        assert_eq!(DynBitArray::with_capacity(4, 64).mask, u64::MAX);
        assert_eq!(DynBitArray::with_capacity(4, 1).mask, 1);
    }
}
