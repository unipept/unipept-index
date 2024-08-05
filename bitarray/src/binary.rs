//! This module provides utilities for reading and writing the bitarray as binary.

use std::io::{BufRead, Read, Result, Write};

use crate::BitArray;

/// The `Binary` trait provides methods for reading and writing a struct as binary.
pub trait Binary {
    /// Writes the struct as binary to the given writer.
    ///
    /// # Arguments
    ///
    /// * `writer` - The writer to write the binary data to.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the write operation is successful, or an `Err` if an error occurs.
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()>;

    /// Reads binary data into a struct from the given reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - The reader to read the binary data from.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read operation is successful, or an `Err` if an error occurs.
    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()>;
}

/// Implementation of the `Binary` trait for the `BitArray` struct.
impl Binary for BitArray {
    /// Writes the binary representation of the `BitArray` to the given writer.
    ///
    /// # Arguments
    ///
    /// * `writer` - The writer to which the binary data will be written.
    ///
    /// # Errors
    ///
    /// Returns an error if there was a problem writing to the writer.
    fn write_binary<W: Write>(&self, writer: &mut W) -> Result<()> {
        for value in self.data.iter() {
            writer.write_all(&value.to_le_bytes())?;
        }

        Ok(())
    }

    /// Reads the binary representation of the `BitArray` from the given reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - The reader from which the binary data will be read.
    ///
    /// # Errors
    ///
    /// Returns an error if there was a problem reading from the reader.
    fn read_binary<R: BufRead>(&mut self, mut reader: R) -> Result<()> {
        self.data.clear();

        let mut buffer = vec![0; 8 * 1024];

        loop {
            let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer)?;
            for buffer_slice in buffer[..bytes_read].chunks_exact(8) {
                self.data.push(u64::from_le_bytes(buffer_slice.try_into().unwrap()));
            }

            if finished {
                break;
            }
        }

        Ok(())
    }
}

/// Fills the buffer with data read from the input.
///
/// # Arguments
///
/// * `input` - The input source to read data from.
/// * `buffer` - The buffer to fill with data.
///
/// # Returns
///
/// Returns a tuple `(finished, bytes_read)` where `finished` indicates whether the end of the input
/// is reached, and `bytes_read` is the number of bytes read into the buffer.
fn fill_buffer<T: Read>(input: &mut T, buffer: &mut Vec<u8>) -> std::io::Result<(bool, usize)> {
    // Store the buffer size in advance, because rust will complain
    // about the buffer being borrowed mutably while it's borrowed
    let buffer_size = buffer.len();

    let mut writable_buffer_space = buffer.as_mut();

    loop {
        match input.read(writable_buffer_space) {
            // No bytes written, which means we've completely filled the buffer
            // or we've reached the end of the file
            Ok(0) => {
                return Ok((!writable_buffer_space.is_empty(), buffer_size - writable_buffer_space.len()));
            }

            // We've read {bytes_read} bytes
            Ok(bytes_read) => {
                // Shrink the writable buffer slice
                writable_buffer_space = writable_buffer_space[bytes_read..].as_mut();
            }

            // An error occurred while reading
            Err(e) => {
                return Err(e);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    pub struct ErrorInput;

    impl Read for ErrorInput {
        fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::Error::new(std::io::ErrorKind::Other, "read error"))
        }
    }

    #[test]
    fn test_fill_buffer() {
        let input_str = "a".repeat(8_000);
        let mut input = input_str.as_bytes();

        let mut buffer = vec![0; 800];

        loop {
            let (finished, bytes_read) = fill_buffer(&mut input, &mut buffer).unwrap();

            if finished {
                assert!(bytes_read < 800);
                break;
            } else {
                assert_eq!(bytes_read, 800);
            }
        }
    }

    #[test]
    fn test_fill_buffer_read_error() {
        let mut input = ErrorInput;
        let mut buffer = vec![0; 800];

        assert!(fill_buffer(&mut input, &mut buffer).is_err());
    }

    #[test]
    fn test_write_binary() {
        let mut bitarray = BitArray::with_capacity(4, 40);
        bitarray.set(0, 0x1234567890);
        bitarray.set(1, 0xabcdef0123);
        bitarray.set(2, 0x4567890abc);
        bitarray.set(3, 0xdef0123456);

        let mut buffer = Vec::new();
        bitarray.write_binary(&mut buffer).unwrap();

        assert_eq!(buffer, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01, 0x00, 0x00,
            0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
        ]);
    }

    #[test]
    fn test_read_binary() {
        let buffer = vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01, 0x00, 0x00,
            0x00, 0x00, 0x56, 0x34, 0x12, 0xf0,
        ];

        let mut bitarray = BitArray::with_capacity(4, 40);
        bitarray.read_binary(&buffer[..]).unwrap();

        assert_eq!(bitarray.get(0), 0x1234567890);
        assert_eq!(bitarray.get(1), 0xabcdef0123);
        assert_eq!(bitarray.get(2), 0x4567890abc);
        assert_eq!(bitarray.get(3), 0xdef0123456);
    }
}
