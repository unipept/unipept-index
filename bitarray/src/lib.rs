//! This module contains the `BitArray` struct and its associated methods.

mod binary;

use std::{
    cmp::max,
    io::{Result, Write}
};

/// Re-export the `Binary` trait.
pub use binary::Binary;

/// A fixed-size bit array implementation.
pub struct BitArray {
    /// The underlying data storage for the bit array.
    data: Vec<u64>,
    /// The mask used to extract the relevant bits from each element in the data vector.
    mask: u64,
    /// The length of the bit array.
    len: usize,
    /// The number of bits in a single element of the data vector.
    bits_per_value: usize
}

impl BitArray {
    /// Creates a new `BitArray` with the specified capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` - The number of bits the `BitArray` can hold.
    /// * `bits_per_value` - The number of bits in a single value.
    ///
    /// # Returns
    ///
    /// A new `BitArray` with the specified capacity.
    pub fn with_capacity(capacity: usize, bits_per_value: usize) -> Self {
        let extra = if (capacity * bits_per_value).is_multiple_of(64) { 0 } else { 1 };
        Self {
            data: vec![0; capacity * bits_per_value / 64 + extra],
            mask: if bits_per_value == 64 { u64::MAX } else { (1 << bits_per_value) - 1 },
            len: capacity,
            bits_per_value
        }
    }

    /// Retrieves the value at the specified index in the `BitArray`.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the value to retrieve.
    ///
    /// # Returns
    ///
    /// The value at the specified index.
    pub fn get(&self, index: usize) -> u64 {
        let start_block = index * self.bits_per_value / 64;
        let start_block_offset = index * self.bits_per_value % 64;

        // If the value is contained within a single block
        if start_block_offset + self.bits_per_value <= 64 {
            // Shift the value to the right so that the relevant bits are in the least significant
            // position Then mask out the irrelevant bits
            return (self.data[start_block] >> (64 - start_block_offset - self.bits_per_value)) & self.mask;
        }

        let end_block = (index + 1) * self.bits_per_value / 64;
        let end_block_offset = (index + 1) * self.bits_per_value % 64;

        // Extract the relevant bits from the start block and shift them {end_block_offset} bits to
        // the left
        let a = self.data[start_block] << end_block_offset;

        // Extract the relevant bits from the end block and shift them to the least significant
        // position
        let b = self.data[end_block] >> (64 - end_block_offset);

        // Paste the two values together and mask out the irrelevant bits
        (a | b) & self.mask
    }

    /// Sets the value at the specified index in the `BitArray`.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the value to set.
    /// * `value` - The value to set at the specified index.
    pub fn set(&mut self, index: usize, value: u64) {
        let value: u64 = value;
        let start_block = index * self.bits_per_value / 64;
        let start_block_offset = index * self.bits_per_value % 64;

        // If the value is contained within a single block
        if start_block_offset + self.bits_per_value <= 64 {
            // Clear the relevant bits in the start block
            self.data[start_block] &= !(self.mask << (64 - start_block_offset - self.bits_per_value));
            // Set the relevant bits in the start block
            self.data[start_block] |= value << (64 - start_block_offset - self.bits_per_value);
            return;
        }

        let end_block = (index + 1) * self.bits_per_value / 64;
        let end_block_offset = (index + 1) * self.bits_per_value % 64;

        // Clear the relevant bits in the start block
        self.data[start_block] &= !(self.mask >> start_block_offset);
        // Set the relevant bits in the start block
        self.data[start_block] |= value >> end_block_offset;

        // Clear the relevant bits in the end block
        self.data[end_block] &= !(self.mask << (64 - end_block_offset));
        // Set the relevant bits in the end block
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    /// Returns the number of bits in a single value.
    ///
    /// # Returns
    ///
    /// The number of bits in a single value.
    pub fn bits_per_value(&self) -> usize {
        self.bits_per_value
    }

    /// Returns the length of the `BitArray`.
    ///
    /// # Returns
    ///
    /// The length of the `BitArray`.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Checks if the `BitArray` is empty.
    ///
    /// # Returns
    ///
    /// `true` if the `BitArray` is empty, `false` otherwise.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Clears the `BitArray`, setting all bits to 0.
    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|x| *x = 0);
    }

    pub fn get_data_slice(&self, start_slice: usize, end_slice: usize) -> &[u64] {
        &self.data[start_slice..end_slice]
    }

    /// Returns a streaming iterator over entries `[start, end)`.
    ///
    /// Keeps `current_word` and `next_word` in local variables so a new slice read only occurs
    /// when crossing a 64-bit block boundary — roughly once per `64 / bits_per_value` entries.
    pub fn iter_range(&self, start: usize, end: usize) -> BitArrayRangeIter<'_> {
        BitArrayRangeIter::new(&self.data, self.bits_per_value, self.mask, start, end)
    }
}

/// Streaming sequential iterator over a contiguous range of a `BitArray`.
///
/// Keeps `current_word` and `next_word` in local variables (register-allocated by the compiler)
/// so that a new slice read only occurs when crossing a 64-bit block boundary — roughly once per
/// `64 / bits_per_value` entries, vs one slice read per entry with `BitArray::get`.
pub struct BitArrayRangeIter<'a> {
    data: &'a [u64],
    bits_per_value: usize,
    mask: u64,
    current_word: u64, // u64 block containing the next value to yield
    next_word: u64,    // u64 block after current_word (pre-loaded)
    block_idx: usize,  // index of current_word within data
    bit_off: usize,    // bit offset of next value within current_word (0..64)
    remaining: usize,  // entries left to yield
}

impl<'a> BitArrayRangeIter<'a> {
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

impl Iterator for BitArrayRangeIter<'_> {
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

        // Advance bit cursor; load next word from data only on block-boundary crossing
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

impl ExactSizeIterator for BitArrayRangeIter<'_> {}

/// Writes the data to a writer in a binary format using a bit array. The data is written
/// in chunks of the specified capacity, so memory usage is minimized.
///
/// # Arguments
///
/// * `data` - The data to write.
/// * `bits_per_value` - The number of bits in a single value.
/// * `max_capacity` - The maximum amount of elements that may be stored in the bit array.
/// * `writer` - The writer to write the data to.
///
/// # Returns
///
/// A `Result` indicating whether the write operation was successful or not.
pub fn data_to_writer(
    data: Vec<i64>,
    bits_per_value: usize,
    max_capacity: usize,
    writer: &mut impl Write
) -> Result<()> {
    // Update the max capacity to be a multiple of the greatest common divisor of the bits per value
    // and 64. This is done to ensure that the bit array can store the data entirely
    let greates_common_divisor = gcd(bits_per_value, 64);
    let capacity = max(greates_common_divisor, max_capacity / greates_common_divisor * greates_common_divisor);

    // If amount of data is less than the max capacity, write the data to the writer in a single
    // chunk
    if data.len() <= capacity {
        let mut bitarray = BitArray::with_capacity(data.len(), bits_per_value);

        for (i, &value) in data.iter().enumerate() {
            bitarray.set(i, value as u64);
        }
        bitarray.write_binary(writer)?;

        return Ok(());
    }

    // Create a bit array that can store a single chunk of data
    let mut bitarray = BitArray::with_capacity(capacity, bits_per_value);

    // Write the data to the writer in chunks of the specified capacity
    let chunks = data.chunks_exact(capacity);

    // Store the remainder before looping over the chunks
    let remainder = chunks.remainder();

    for chunk in chunks {
        for (i, &value) in chunk.iter().enumerate() {
            bitarray.set(i, value as u64);
        }
        bitarray.write_binary(writer)?;
        bitarray.clear();
    }

    // Create a new bit array with the remainder capacity
    bitarray = BitArray::with_capacity(remainder.len(), bits_per_value);

    for (i, &value) in remainder.iter().enumerate() {
        bitarray.set(i, value as u64);
    }
    bitarray.write_binary(writer)?;

    Ok(())
}

/// Calculates the greatest common divisor of two numbers.
///
/// # Arguments
///
/// * `a` - The first number.
/// * `b` - The second number.
///
/// # Returns
///
/// The greatest common divisor of the two numbers.
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

    #[test]
    fn test_bitarray_with_capacity() {
        let bitarray = BitArray::with_capacity(4, 40);
        assert_eq!(bitarray.data, vec![0, 0, 0]);
        assert_eq!(bitarray.mask, 0xff_ffff_ffff);
        assert_eq!(bitarray.len, 4);
    }

    #[test]
    fn test_bitarray_get() {
        let mut bitarray = BitArray::with_capacity(4, 40);
        bitarray.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(bitarray.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(bitarray.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(bitarray.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(bitarray.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_bitarray_set() {
        let mut bitarray = BitArray::with_capacity(4, 40);

        bitarray.set(0, 0b0001110011111010110001000111111100110010_u64);
        bitarray.set(1, 0b1100001001010010011000010100110111001001_u64);
        bitarray.set(2, 0b1111001101001101101101101011101001010001_u64);
        bitarray.set(3, 0b0000100010010001010001001110101110011100_u64);

        assert_eq!(bitarray.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_bitarray_bits_per_value() {
        let bitarray = BitArray::with_capacity(4, 40);
        assert_eq!(bitarray.bits_per_value(), 40);
    }

    #[test]
    fn test_bitarray_len() {
        let bitarray = BitArray::with_capacity(4, 40);
        assert_eq!(bitarray.len(), 4);
    }

    #[test]
    fn test_bitarray_is_empty() {
        let bitarray = BitArray::with_capacity(0, 40);
        assert!(bitarray.is_empty());
    }

    #[test]
    fn test_bitarray_is_not_empty() {
        let bitarray = BitArray::with_capacity(4, 40);
        assert!(!bitarray.is_empty());
    }

    #[test]
    fn test_bitarray_clear() {
        let mut bitarray = BitArray::with_capacity(4, 40);
        bitarray.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        bitarray.clear();

        assert_eq!(bitarray.data, vec![0, 0, 0]);
    }

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

    // --- BitArrayRangeIter unit tests ---

    fn collect_range(ba: &BitArray, start: usize, end: usize) -> Vec<i64> {
        ba.iter_range(start, end).collect()
    }

    fn expected_range(ba: &BitArray, start: usize, end: usize) -> Vec<i64> {
        (start..end).map(|i| ba.get(i) as i64).collect()
    }

    #[test]
    fn test_iter_range_empty() {
        let ba = BitArray::with_capacity(8, 32);
        assert!(collect_range(&ba, 3, 3).is_empty());
        assert!(collect_range(&ba, 5, 3).is_empty()); // inverted range — must not panic
    }

    #[test]
    fn test_iter_range_single_entry() {
        let mut ba = BitArray::with_capacity(4, 40);
        ba.set(2, 0xABCDEF1234_u64);
        assert_eq!(collect_range(&ba, 2, 3), vec![0xABCDEF1234_i64]);
    }

    #[test]
    fn test_iter_range_mid_block_start() {
        // Start at an index whose bit offset is non-zero within the first u64 block.
        // With bits_per_value=32: block boundary every 2 entries. Start=1 → bit_off=32.
        let values: Vec<u64> = (0..8).map(|i| i * 111 + 7).collect();
        let mut ba = BitArray::with_capacity(8, 32);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }

        assert_eq!(collect_range(&ba, 1, 6), expected_range(&ba, 1, 6));
    }

    #[test]
    fn test_iter_range_crosses_block_boundary() {
        // bits_per_value=40: 64/gcd(40,64)=8 entries per cycle, 5 entries per u64 block approx.
        // Any range spanning >5 entries will cross a 64-bit boundary.
        let values: Vec<u64> = (0..16).map(|i| i as u64 * 0x100000001 + 3).collect();
        let mut ba = BitArray::with_capacity(16, 40);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }

        assert_eq!(collect_range(&ba, 0, 16), expected_range(&ba, 0, 16));
        assert_eq!(collect_range(&ba, 3, 13), expected_range(&ba, 3, 13));
    }

    #[test]
    fn test_iter_range_bits_per_value_64() {
        // Each entry occupies exactly one u64 block — no boundary crossing needed.
        let values: Vec<u64> = (0..8).map(|i| i as u64 * 0xDEAD_BEEF + 1).collect();
        let mut ba = BitArray::with_capacity(8, 64);
        for (i, &v) in values.iter().enumerate() { ba.set(i, v); }

        assert_eq!(collect_range(&ba, 0, 8), expected_range(&ba, 0, 8));
        assert_eq!(collect_range(&ba, 2, 6), expected_range(&ba, 2, 6));
    }

    #[test]
    fn test_iter_range_bits_per_value_1() {
        // 64 entries per block; exercises many iterations before a boundary crossing.
        let mut ba = BitArray::with_capacity(128, 1);
        for i in (0..128).step_by(3) { ba.set(i, 1); }

        assert_eq!(collect_range(&ba, 0, 128), expected_range(&ba, 0, 128));
        assert_eq!(collect_range(&ba, 60, 70), expected_range(&ba, 60, 70)); // crosses block boundary at 64
    }

    #[test]
    fn test_iter_range_exact_size_iterator() {
        let mut ba = BitArray::with_capacity(10, 40);
        for i in 0..10 { ba.set(i, i as u64 * 99); }
        let iter = ba.iter_range(2, 8);
        assert_eq!(iter.len(), 6);
    }
}
