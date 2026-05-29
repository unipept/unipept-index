use std::{
    error::Error,
    io::{BufRead, Write}
};

use bitarray::{Binary, DynBitArray, data_to_writer};
use text_compression::WriteBinary;

/// Owned, in-memory compressed suffix array backend.
pub struct CompressedSA(pub DynBitArray, pub u8);

impl super::SuffixArrayBackend for CompressedSA {
    type RangeIter<'a> = bitarray::DynBitArrayRangeIter<'a>;

    fn len(&self) -> usize { self.0.len() }
    fn bits_per_value(&self) -> usize { self.0.bits_per_value() }
    fn sample_rate(&self) -> u8 { self.1 }
    #[inline]
    fn get(&self, index: usize) -> i64 { self.0.get(index) as i64 }

    fn iter_range(&self, start: usize, end: usize) -> bitarray::DynBitArrayRangeIter<'_> {
        self.0.iter_range(start, end)
    }

    #[inline]
    fn prefetch_sa_index(&self, index: usize) {
        if index < self.0.len() {
            let word_idx = (index * self.0.bits_per_value()) / 64;
            let ptr: *const u64 = self.0.get_data_slice(word_idx, word_idx + 1).as_ptr();
            prefetch::prefetch_read(ptr);
        }
    }
}

impl WriteBinary for CompressedSA {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let CompressedSA(bit_array, sample_rate) = self;
        writer.write_all(&[bit_array.bits_per_value() as u8])
            .map_err(|_| "Could not write the required bits")?;
        writer.write_all(&[sample_rate])
            .map_err(|_| "Could not write the sparseness factor")?;
        writer.write_all(&(bit_array.len() as u64).to_le_bytes())
            .map_err(|_| "Could not write the size")?;
        bit_array.write_binary(writer)
            .map_err(|_| "Could not write the compressed suffix array")?;
        Ok(())
    }
}

// ── I/O helpers ──────────────────────────────────────────────────────────────

/// Writes the compressed suffix array to a writer.
pub fn dump_compressed_suffix_array(
    sa: Vec<i64>,
    sparseness_factor: u8,
    bits_per_value: usize,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    writer.write(&[bits_per_value as u8]).map_err(|_| "Could not write the required bits to the writer")?;
    writer.write(&[sparseness_factor]).map_err(|_| "Could not write the sparseness factor to the writer")?;
    writer.write(&(sa.len() as u64).to_le_bytes()).map_err(|_| "Could not write the size of the suffix array to the writer")?;
    data_to_writer(sa, bits_per_value, 8 * 1024, writer).map_err(|_| "Could not write the compressed suffix array to the writer")?;
    Ok(())
}

/// Load the compressed suffix array, reading the sample_rate + size header first.
/// Returns a `CompressedSA`.
pub fn load_compressed_suffix_array(
    reader: &mut impl BufRead,
    bits_per_value: usize
) -> Result<CompressedSA, Box<dyn Error>> {
    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    let mut size_buffer = [0_u8; 8];
    reader.read_exact(&mut size_buffer).map_err(|_| "Could not read the size of the suffix array from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    let sa = load_compressed(reader, bits_per_value, size)?;
    Ok(CompressedSA(sa, sample_rate))
}

/// Inner helper: load the compressed BitArray body (no header).
pub(super) fn load_compressed(
    reader: &mut impl BufRead,
    bits_per_value: usize,
    size: usize
) -> Result<DynBitArray, Box<dyn Error>> {
    let mut compressed_suffix_array = DynBitArray::with_capacity(size, bits_per_value);
    compressed_suffix_array.read_binary(reader).map_err(|_| "Could not read the compressed suffix array from the binary file")?;
    Ok(compressed_suffix_array)
}

#[cfg(test)]
mod tests {
    use std::io::Read;
    use super::*;
    use super::super::SuffixArrayBackend;

    pub struct FailingWriter { pub valid_write_count: usize }

    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> Result<usize, std::io::Error> {
            if self.valid_write_count == 0 { return Err(std::io::Error::other("Write failed")); }
            self.valid_write_count -= 1;
            Ok(1)
        }
        fn flush(&mut self) -> Result<(), std::io::Error> { Ok(()) }
    }

    pub struct FailingReader { pub valid_read_count: usize }

    impl Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.valid_read_count == 0 { return Err(std::io::Error::other("Read failed")); }
            self.valid_read_count -= 1;
            Ok(buf.len())
        }
    }

    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> { Ok(&[]) }
        fn consume(&mut self, _: usize) {}
    }

    #[test]
    fn test_dump_compressed_suffix_array() {
        let sa = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut writer = vec![];
        dump_compressed_suffix_array(sa, 1, 8, &mut writer).unwrap();
        assert_eq!(writer, vec![
            8, 1, 10, 0, 0, 0, 0, 0, 0, 0,
            8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9
        ]);
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_compressed_suffix_array_fail_required_bits() {
        dump_compressed_suffix_array(vec![], 1, 8, &mut FailingWriter { valid_write_count: 0 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the sparseness factor to the writer")]
    fn test_dump_compressed_suffix_array_fail_sparseness_factor() {
        dump_compressed_suffix_array(vec![], 1, 8, &mut FailingWriter { valid_write_count: 1 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the suffix array to the writer")]
    fn test_dump_compressed_suffix_array_fail_size() {
        dump_compressed_suffix_array(vec![], 1, 8, &mut FailingWriter { valid_write_count: 2 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the compressed suffix array to the writer")]
    fn test_dump_compressed_suffix_array_fail_compressed_suffix_array() {
        dump_compressed_suffix_array(vec![1], 1, 8, &mut FailingWriter { valid_write_count: 3 }).unwrap();
    }

    #[test]
    fn test_load_compressed_suffix_array() {
        let data = [
            1,
            10, 0, 0, 0, 0, 0, 0, 0,
            8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9
        ];
        let mut reader = std::io::BufReader::new(&data[..]);
        let sa = load_compressed_suffix_array(&mut reader, 8).unwrap();
        assert_eq!(sa.sample_rate(), 1);
        for i in 0..10 {
            assert_eq!(sa.get(i), i as i64 + 1);
        }
    }

    #[test]
    #[should_panic(expected = "Could not read the sample rate from the binary file")]
    fn test_load_compressed_suffix_array_fail_sample_rate() {
        load_compressed_suffix_array(&mut FailingReader { valid_read_count: 0 }, 8).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the suffix array from the binary file")]
    fn test_load_compressed_suffix_array_fail_size() {
        load_compressed_suffix_array(&mut FailingReader { valid_read_count: 1 }, 8).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the compressed suffix array from the binary file")]
    fn test_load_compressed_suffix_array_fail_compressed_suffix_array() {
        load_compressed_suffix_array(&mut FailingReader { valid_read_count: 2 }, 8).unwrap();
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
        assert_eq!(reader.fill_buf().unwrap(), &[] as &[u8]);
        let mut buffer = [0_u8; 1];
        assert!(reader.read(&mut buffer).is_err());
    }
}
