//! This module contains the `BitArray` struct and its associated methods.

mod binary;

use std::{io::{Result, Write}, ops::Index};

/// Re-export the `Binary` trait.
pub use binary::Binary;

/// A fixed-size bit array implementation.
pub struct BitArray<const B: usize> {
    /// The underlying data storage for the bit array.
    data: Vec<u64>,
    /// The mask used to extract the relevant bits from each element in the data vector.
    mask: u64,
    /// The length of the bit array.
    len:  usize
}

impl<const B: usize> BitArray<B> {
    /// Creates a new `BitArray` with the specified capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` - The number of bits the `BitArray` can hold.
    ///
    /// # Returns
    ///
    /// A new `BitArray` with the specified capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: vec![0; capacity * B / 64 + 1],
            mask: (1 << B) - 1,
            len:  capacity
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
        let start_block = index * B / 64;
        let start_block_offset = index * B % 64;

        // If the value is contained within a single block
        if start_block_offset + B <= 64 {
            // Shift the value to the right so that the relevant bits are in the least significant
            // position Then mask out the irrelevant bits
            return self.data[start_block] >> (64 - start_block_offset - B) & self.mask;
        }

        let end_block = (index + 1) * B / 64;
        let end_block_offset = (index + 1) * B % 64;

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
        let start_block = index * B / 64;
        let start_block_offset = index * B % 64;

        // If the value is contained within a single block
        if start_block_offset + B <= 64 {
            // Clear the relevant bits in the start block
            self.data[start_block] &= !(self.mask << (64 - start_block_offset - B));
            // Set the relevant bits in the start block
            self.data[start_block] |= value << (64 - start_block_offset - B);
            return;
        }

        let end_block = (index + 1) * B / 64;
        let end_block_offset = (index + 1) * B % 64;

        // Clear the relevant bits in the start block
        self.data[start_block] &= !(self.mask >> start_block_offset);
        // Set the relevant bits in the start block
        self.data[start_block] |= value >> end_block_offset;

        // Clear the relevant bits in the end block
        self.data[end_block] &= !(self.mask << (64 - end_block_offset));
        // Set the relevant bits in the end block
        self.data[end_block] |= value << (64 - end_block_offset);
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
}

/// Writes the data to a writer in a binary format using a bit array. This function is helpfull
/// when writing large amounts of data to a writer in chunks. The data is written in chunks of the
/// specified capacity, so memory usage is minimized.
///
/// # Arguments
///
/// * `data` - The data to write.
/// * `writer` - The writer to write the data to.
/// * `max_capacity` - The maximum amount of elements that may be stored in the bit array.
///
/// # Returns
///
/// A `Result` indicating whether the write operation was successful or not.
pub fn data_to_writer<const B: usize>(
    data: Vec<i64>, 
    writer: &mut impl Write,
    max_capacity: usize
) -> Result<()> {
    // Calculate the capacity of the bit array so the data buffer can be stored entirely
    // This makes the process of writing partial data to the writer easier as bounds checking is not needed
    let capacity = max_capacity / (B * 64) * B * 64;

    // Create a bit array that can store a single chunk of data
    let mut bitarray = BitArray::<B>::with_capacity(capacity);

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
    bitarray = BitArray::<B>::with_capacity(remainder.len());

    for (i, &value) in remainder.iter().enumerate() {
        bitarray.set(i, value as u64);
    }
    bitarray.write_binary(writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitarray_with_capacity() {
        let bitarray = BitArray::<40>::with_capacity(4);
        assert_eq!(bitarray.data, vec![0, 0, 0]);
        assert_eq!(bitarray.mask, 0xff_ffff_ffff);
        assert_eq!(bitarray.len, 4);
    }

    #[test]
    fn test_bitarray_get() {
        let mut bitarray = BitArray::<40>::with_capacity(4);
        bitarray.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

        assert_eq!(bitarray.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(bitarray.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(bitarray.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(bitarray.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_bitarray_set() {
        let mut bitarray = BitArray::<40>::with_capacity(4);

        bitarray.set(0, 0b0001110011111010110001000111111100110010);
        bitarray.set(1, 0b1100001001010010011000010100110111001001);
        bitarray.set(2, 0b1111001101001101101101101011101001010001);
        bitarray.set(3, 0b0000100010010001010001001110101110011100);

        assert_eq!(bitarray.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
    }

    #[test]
    fn test_bitarray_len() {
        let bitarray = BitArray::<40>::with_capacity(4);
        assert_eq!(bitarray.len(), 4);
    }

    #[test]
    fn test_bitarray_is_empty() {
        let bitarray = BitArray::<40>::with_capacity(0);
        assert!(bitarray.is_empty());
    }

    #[test]
    fn test_bitarray_is_not_empty() {
        let bitarray = BitArray::<40>::with_capacity(4);
        assert!(!bitarray.is_empty());
    }

    #[test]
    fn test_data_to_writer() {
        let data = vec![0x1234567890, 0xabcdef0123, 0x4567890abc, 0xdef0123456];
        let mut writer = Vec::new();

        data_to_writer::<40>(data, &mut writer, 2).unwrap();

        assert_eq!(writer, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12,
            0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01,
            0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
        ]);
    }
}
