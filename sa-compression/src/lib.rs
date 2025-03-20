use std::{
    error::Error,
    io::{BufRead, Write}
};

use bitarray::{Binary, BitArray, data_to_writer};
use sa_index::SuffixArray;

/// Writes the compressed suffix array to a writer.
///
/// # Arguments
///
/// * `sa` - The suffix array to be compressed.
/// * `sparseness_factor` - The sparseness factor used for compression.
/// * `bits_per_value` - The number of bits used to represent each value in the compressed array.
/// * `writer` - The writer to which the compressed array will be written.
///
/// # Errors
///
/// Returns an error if writing to the writer fails.
pub fn dump_compressed_suffix_array(
    sa: Vec<i64>,
    sparseness_factor: u8,
    bits_per_value: usize,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    // Write the flags to the writer
    // 00000001 indicates that the suffix array is compressed
    writer
        .write(&[bits_per_value as u8])
        .map_err(|_| "Could not write the required bits to the writer")?;

    // Write the sparseness factor to the writer
    writer
        .write(&[sparseness_factor])
        .map_err(|_| "Could not write the sparseness factor to the writer")?;

    // Write the size of the suffix array to the writer
    writer
        .write(&(sa.len() as u64).to_le_bytes())
        .map_err(|_| "Could not write the size of the suffix array to the writer")?;

    // Compress the suffix array and write it to the writer
    data_to_writer(sa, bits_per_value, 8 * 1024, writer)
        .map_err(|_| "Could not write the compressed suffix array to the writer")?;

    Ok(())
}

/// Load the compressed suffix array from a reader.
///
/// # Arguments
///
/// * `reader` - The reader from which the compressed array will be read.
/// * `bits_per_value` - The number of bits used to represent each value in the compressed array.
///
/// # Errors
///
/// Returns an error if reading from the reader fails.
pub fn load_compressed_suffix_array(
    reader: &mut impl BufRead,
    bits_per_value: usize
) -> Result<SuffixArray, Box<dyn Error>> {
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

    // Read the compressed suffix array from the binary file
    let mut compressed_suffix_array = BitArray::with_capacity(size, bits_per_value);
    compressed_suffix_array
        .read_binary(reader)
        .map_err(|_| "Could not read the compressed suffix array from the binary file")?;

    Ok(SuffixArray::Compressed(compressed_suffix_array, sample_rate))
}

#[cfg(test)]
mod tests {
    use std::io::Read;

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
    fn test_dump_compressed_suffix_array() {
        let sa = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        let mut writer = vec![];
        dump_compressed_suffix_array(sa, 1, 8, &mut writer).unwrap();

        assert_eq!(writer, vec![
            // bits per value
            8, // sparseness factor
            1, // size of the suffix array
            10, 0, 0, 0, 0, 0, 0, 0, // compressed suffix array
            8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9
        ]);
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_compressed_suffix_array_fail_required_bits() {
        let mut writer = FailingWriter { valid_write_count: 0 };

        dump_compressed_suffix_array(vec![], 1, 8, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the sparseness factor to the writer")]
    fn test_dump_compressed_suffix_array_fail_sparseness_factor() {
        let mut writer = FailingWriter { valid_write_count: 1 };

        dump_compressed_suffix_array(vec![], 1, 8, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the suffix array to the writer")]
    fn test_dump_compressed_suffix_array_fail_size() {
        let mut writer = FailingWriter { valid_write_count: 2 };

        dump_compressed_suffix_array(vec![], 1, 8, &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the compressed suffix array to the writer")]
    fn test_dump_compressed_suffix_array_fail_compressed_suffix_array() {
        let mut writer = FailingWriter { valid_write_count: 3 };

        dump_compressed_suffix_array(vec![1], 1, 8, &mut writer).unwrap();
    }

    #[test]
    fn test_load_compressed_suffix_array() {
        let data = vec![
            // sparseness factor
            1, // size of the suffix array
            10, 0, 0, 0, 0, 0, 0, 0, // compressed suffix array
            8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9,
        ];

        let mut reader = std::io::BufReader::new(&data[..]);
        let compressed_suffix_array = load_compressed_suffix_array(&mut reader, 8).unwrap();

        assert_eq!(compressed_suffix_array.sample_rate(), 1);
        for i in 0..10 {
            assert_eq!(compressed_suffix_array.get(i), i as i64 + 1);
        }
    }

    #[test]
    #[should_panic(expected = "Could not read the sample rate from the binary file")]
    fn test_load_compressed_suffix_array_fail_sample_rate() {
        let mut reader = FailingReader { valid_read_count: 0 };

        load_compressed_suffix_array(&mut reader, 8).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the suffix array from the binary file")]
    fn test_load_compressed_suffix_array_fail_size() {
        let mut reader = FailingReader { valid_read_count: 1 };

        load_compressed_suffix_array(&mut reader, 8).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the compressed suffix array from the binary file")]
    fn test_load_compressed_suffix_array_fail_compressed_suffix_array() {
        let mut reader = FailingReader { valid_read_count: 2 };

        load_compressed_suffix_array(&mut reader, 8).unwrap();
    }

    #[test]
    fn test_failing_writer() {
        let mut writer = FailingWriter { valid_write_count: 0 };
        assert!(writer.flush().is_ok());
        assert!(writer.write(&[0]).is_err());
    }

    #[test]
    fn test_failing_reader() {
        let mut reader = FailingReader { valid_read_count: 0 };
        let right_buffer: [u8; 0] = [];
        assert_eq!(reader.fill_buf().unwrap(), &right_buffer);
        assert_eq!(reader.consume(0), ());
        let mut buffer = [0_u8; 1];
        assert!(reader.read(&mut buffer).is_err());
    }
}
