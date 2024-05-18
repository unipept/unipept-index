use std::{
    error::Error,
    io::{
        BufRead,
        Read,
        Write
    }
};

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
    fn write_binary<W: Write>(&self, writer: &mut W) -> std::io::Result<()>;

    /// Reads binary data into a struct from the given reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - The reader to read the binary data from.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read operation is successful, or an `Err` if an error occurs.
    fn read_binary<R: BufRead>(&mut self, reader: R) -> std::io::Result<()>;
}

/// Implements the `Binary` trait for `Vec<i64>`.
impl Binary for Vec<i64> {
    /// Writes the elements of the vector to a binary file.
    ///
    /// # Arguments
    ///
    /// * `writer` - The writer to which the binary data will be written.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the write operation is successful, or an `std::io::Error` otherwise.
    fn write_binary<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        for value in self {
            writer.write_all(&value.to_le_bytes())?;
        }

        Ok(())
    }

    /// Reads binary data from a reader and populates the vector with the read values.
    ///
    /// # Arguments
    ///
    /// * `reader` - The reader from which the binary data will be read.
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` if the read operation is successful, or an `std::io::Error` otherwise.
    fn read_binary<R: BufRead>(&mut self, mut reader: R) -> std::io::Result<()> {
        self.clear();

        let mut buffer = vec![0; 8 * 1024];

        loop {
            let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer)?;
            for buffer_slice in buffer[.. bytes_read].chunks_exact(8) {
                self.push(i64::from_le_bytes(buffer_slice.try_into().unwrap()));
            }

            if finished {
                break;
            }
        }

        Ok(())
    }
}

/// Writes the suffix array to a binary file.
///
/// # Arguments
///
/// * `sa` - The suffix array to dump.
/// * `sparseness_factor` - The sparseness factor to write to the file.
/// * `writer` - The writer to write the binary data to.
///
/// # Returns
///
/// Returns `Ok(())` if the write operation is successful, or an `Err` if an error occurs.
pub fn dump_suffix_array(
    sa: &Vec<i64>,
    sparseness_factor: u8,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    // Write the required bits to the writer
    // 01000000 indicates that the suffix array is not compressed
    writer
        .write(&[64_u8])
        .map_err(|_| "Could not write the required bits to the writer")?;

    // Write the sparseness factor to the writer
    writer
        .write(&[sparseness_factor])
        .map_err(|_| "Could not write the sparseness factor to the writer")?;

    // Write the size of the suffix array to the writer
    let sa_len = sa.len();
    writer
        .write(&(sa_len).to_le_bytes())
        .map_err(|_| "Could not write the size of the suffix array to the writer")?;

    // Write the suffix array to the writer
    sa.write_binary(writer)
        .map_err(|_| "Could not write the suffix array to the writer")?;

    Ok(())
}

/// Loads the suffix array from the file with the given `filename`
///
/// # Arguments
/// * `filename` - The filename of the file where the suffix array is stored
///
/// # Returns
///
/// Returns the sample rate of the suffix array, together with the suffix array
///
/// # Errors
///
/// Returns any error from opening the file or reading the file
pub fn load_suffix_array(reader: &mut impl BufRead) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    // Read the sample rate from the binary file (1 byte)
    let mut sample_rate_buffer = [0_u8; 1];
    reader
        .read_exact(&mut sample_rate_buffer)
        .map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    // Read the size of the suffix array from the binary file (8 bytes)
    let mut size_buffer = [0_u8; 8];
    reader
        .read_exact(&mut size_buffer)
        .map_err(|_| "Could not read the size of the suffix array from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    let mut sa = Vec::with_capacity(size);
    sa.read_binary(reader)
        .map_err(|_| "Could not read the suffix array from the binary file")?;

    Ok((sample_rate, sa))
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
                return Ok((
                    !writable_buffer_space.is_empty(),
                    buffer_size - writable_buffer_space.len()
                ));
            }

            // We've read {bytes_read} bytes
            Ok(bytes_read) => {
                // Shrink the writable buffer slice
                writable_buffer_space = writable_buffer_space[bytes_read ..].as_mut();
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

    pub struct FailingWriter {
        /// The number of times the write function can be called before it fails.
        pub valid_write_count: usize
    }

    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> Result<usize, std::io::Error> {
            if self.valid_write_count == 0 {
                return Err(std::io::Error::new(std::io::ErrorKind::Other, "Write failed"));
            }

            self.valid_write_count -= 1;
            Ok(1)
        }

        fn flush(&mut self) -> Result<(), std::io::Error> {
            Ok(())
        }
    }

    pub struct FailingReader {
        /// The number of times the read function can be called before it fails.
        pub valid_read_count: usize
    }

    impl Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.valid_read_count == 0 {
                return Err(std::io::Error::new(std::io::ErrorKind::Other, "Read failed"));
            }

            self.valid_read_count -= 1;
            Ok(buf.len())
        }
    }

    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
            Ok(&[])
        }

        fn consume(&mut self, _: usize) {}
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
        let mut input = FailingReader { valid_read_count: 0 };
        let mut buffer = vec![0; 800];

        assert!(fill_buffer(&mut input, &mut buffer).is_err());
    }

    #[test]
    fn test_write_binary() {
        let mut buffer = Vec::new();
        let values = vec![1, 2, 3, 4, 5];

        values.write_binary(&mut buffer).unwrap();

        assert_eq!(
            buffer,
            vec![
                1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0
            ]
        );
    }

    #[test]
    fn test_read_binary() {
        let buffer = vec![
            1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0,
            0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
        ];

        let mut values = Vec::new();
        values.read_binary(buffer.as_slice()).unwrap();

        assert_eq!(values, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_dump_suffix_array() {
        let mut buffer = Vec::new();
        let sa = vec![1, 2, 3, 4, 5];

        dump_suffix_array(&sa, 1, &mut buffer).unwrap();

        assert_eq!(
            buffer,
            vec![
                // required bits
                64, // Sparseness factor
                1,  // Size of the suffix array
                5, 0, 0, 0, 0, 0, 0, 0, // Suffix array
                1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0
            ]
        );
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_suffix_array_fail_required_bits() {
        let mut writer = FailingWriter {
            valid_write_count: 0
        };

        dump_suffix_array(&vec![], 1, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the sparseness factor to the writer")]
    fn test_dump_suffix_array_fail_sparseness_factor() {
        let mut writer = FailingWriter {
            valid_write_count: 1
        };

        dump_suffix_array(&vec![], 1, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the suffix array to the writer")]
    fn test_dump_suffix_array_fail_size() {
        let mut writer = FailingWriter {
            valid_write_count: 2
        };

        dump_suffix_array(&vec![], 1, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the suffix array to the writer")]
    fn test_dump_suffix_array_fail_suffix_array() {
        let mut writer = FailingWriter {
            valid_write_count: 3
        };

        dump_suffix_array(&vec![ 1 ], 1, &mut writer).unwrap();
    }

    #[test]
    fn test_load_suffix_array() {
        let buffer = vec![
            // Sample rate
            1, // Size of the suffix array
            5, 0, 0, 0, 0, 0, 0, 0, // Suffix array
            1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0,
            0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
        ];

        let mut reader = buffer.as_slice();
        let (sample_rate, sa) = load_suffix_array(&mut reader).unwrap();

        assert_eq!(sample_rate, 1);
        assert_eq!(sa, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    #[should_panic(expected = "Could not read the sample rate from the binary file")]
    fn test_load_suffix_array_fail_sample_rate() {
        let mut reader = FailingReader {
            valid_read_count: 0
        };

        load_suffix_array(&mut reader).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the suffix array from the binary file")]
    fn test_load_suffix_array_fail_size() {
        let mut reader = FailingReader {
            valid_read_count: 1
        };

        load_suffix_array(&mut reader).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the suffix array from the binary file")]
    fn test_load_suffix_array_fail_suffix_array() {
        let mut reader = FailingReader {
            valid_read_count: 2
        };

        load_suffix_array(&mut reader).unwrap();
    }
}
