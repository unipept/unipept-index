use std::{
    error::Error,
    io::{BufRead, Read, Write}
};
use text_compression::WriteBinary;

/// Owned, in-memory original (uncompressed) suffix array backend.
pub struct OriginalSA(pub Vec<i64>, pub u8);

/// Wrapper around `slice::Iter` that yields `i64` values directly.
pub struct OriginalRangeIter<'a>(pub std::slice::Iter<'a, i64>);

impl Iterator for OriginalRangeIter<'_> {
    type Item = i64;
    #[inline]
    fn next(&mut self) -> Option<i64> { self.0.next().copied() }
    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) { let n = self.0.len(); (n, Some(n)) }
}

impl ExactSizeIterator for OriginalRangeIter<'_> {}

impl super::SuffixArrayBackend for OriginalSA {
    type RangeIter<'a> = OriginalRangeIter<'a>;

    fn len(&self) -> usize { self.0.len() }
    fn bits_per_value(&self) -> usize { 64 }
    fn sample_rate(&self) -> u8 { self.1 }
    #[inline]
    fn get(&self, index: usize) -> i64 { self.0[index] }

    fn iter_range(&self, start: usize, end: usize) -> OriginalRangeIter<'_> {
        OriginalRangeIter(self.0.get(start..end).unwrap_or(&[]).iter())
    }

    #[inline]
    fn prefetch_sa_index(&self, index: usize) {
        if index < self.0.len() {
            let ptr: *const i64 = &self.0[index];
            prefetch::prefetch_read(ptr);
        }
    }
}

impl WriteBinary for OriginalSA {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        dump_suffix_array(self.0, self.1, writer)
    }
}

// ── I/O helpers ──────────────────────────────────────────────────────────────

fn fill_buffer<T: Read>(input: &mut T, buffer: &mut Vec<u8>) -> std::io::Result<(bool, usize)> {
    let buffer_size = buffer.len();
    let mut writable_buffer_space = buffer.as_mut();

    loop {
        match input.read(writable_buffer_space) {
            Ok(0) => {
                return Ok((!writable_buffer_space.is_empty(), buffer_size - writable_buffer_space.len()));
            }
            Ok(bytes_read) => {
                writable_buffer_space = writable_buffer_space[bytes_read..].as_mut();
            }
            Err(e) => {
                return Err(e);
            }
        }
    }
}

fn write_vec_i64(vec: Vec<i64>, writer: &mut impl Write) -> std::io::Result<()> {
    for value in vec {
        writer.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

fn read_vec_i64(vec: &mut Vec<i64>, mut reader: impl BufRead) -> std::io::Result<()> {
    vec.clear();
    let mut buffer = vec![0; 8 * 1024];

    loop {
        let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer)?;
        for buffer_slice in buffer[..bytes_read].chunks_exact(8) {
            vec.push(i64::from_le_bytes(buffer_slice.try_into().unwrap()));
        }

        if finished {
            break;
        }
    }

    Ok(())
}

/// Writes the suffix array to a binary file.
pub fn dump_suffix_array(sa: Vec<i64>, sparseness_factor: u8, writer: &mut impl Write) -> Result<(), Box<dyn Error>> {
    writer.write(&[64_u8]).map_err(|_| "Could not write the required bits to the writer")?;
    writer.write(&[sparseness_factor]).map_err(|_| "Could not write the sparseness factor to the writer")?;
    writer.write(&(sa.len()).to_le_bytes()).map_err(|_| "Could not write the size of the suffix array to the writer")?;
    write_vec_i64(sa, writer).map_err(|_| "Could not write the suffix array to the writer")?;
    Ok(())
}

/// Inner helper: load the original (uncompressed) suffix array body after the header is already read.
pub(super) fn load_original(reader: &mut impl BufRead, _sample_rate: u8, size: usize) -> Result<Vec<i64>, Box<dyn Error>> {
    let mut sa = Vec::with_capacity(size);
    read_vec_i64(&mut sa, reader).map_err(|_| "Could not read the suffix array from the binary file")?;
    Ok(sa)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_fill_buffer() {
        let input_str = "a".repeat(8_000);
        let mut input = input_str.as_bytes();
        let mut buffer = vec![0; 800];
        loop {
            let (finished, bytes_read) = fill_buffer(&mut input, &mut buffer).unwrap();
            if finished { assert!(bytes_read < 800); break; }
            else { assert_eq!(bytes_read, 800); }
        }
    }

    #[test]
    fn test_fill_buffer_read_error() {
        let mut input = FailingReader { valid_read_count: 0 };
        let mut buffer = vec![0; 800];
        assert!(fill_buffer(&mut input, &mut buffer).is_err());
    }

    #[test]
    fn test_dump_suffix_array() {
        let mut buffer = Vec::new();
        dump_suffix_array(vec![1, 2, 3, 4, 5], 1, &mut buffer).unwrap();
        assert_eq!(buffer, vec![
            64, 1, 5, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
            3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0,
            5, 0, 0, 0, 0, 0, 0, 0,
        ]);
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_suffix_array_fail_required_bits() {
        dump_suffix_array(vec![], 1, &mut FailingWriter { valid_write_count: 0 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the sparseness factor to the writer")]
    fn test_dump_suffix_array_fail_sparseness_factor() {
        dump_suffix_array(vec![], 1, &mut FailingWriter { valid_write_count: 1 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the suffix array to the writer")]
    fn test_dump_suffix_array_fail_size() {
        dump_suffix_array(vec![], 1, &mut FailingWriter { valid_write_count: 2 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the suffix array to the writer")]
    fn test_dump_suffix_array_fail_suffix_array() {
        dump_suffix_array(vec![1], 1, &mut FailingWriter { valid_write_count: 3 }).unwrap();
    }
}
