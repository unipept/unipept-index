use std::{
    error::Error,
    io::{BufRead, Read, Write}
};

use binary_traits::WriteBinary;

use super::super::{SuffixArrayBackend, write_sa_header};

/// Suffix array in owned memory, one `i64` per entry.
///
/// Field 0 is the entries; field 1 is the sparseness factor the array was built with.
pub struct OriginalSA(pub Vec<i64>, pub u8);

/// Wrapper around `slice::Iter` that yields `i64` values directly.
pub struct OriginalRangeIter<'a>(pub std::slice::Iter<'a, i64>);

impl Iterator for OriginalRangeIter<'_> {
    type Item = i64;
    #[inline]
    fn next(&mut self) -> Option<i64> {
        self.0.next().copied()
    }
    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.0.len();
        (n, Some(n))
    }
}

impl ExactSizeIterator for OriginalRangeIter<'_> {}

impl SuffixArrayBackend for OriginalSA {
    type RangeIter<'a> = OriginalRangeIter<'a>;

    fn len(&self) -> usize {
        self.0.len()
    }
    fn bits_per_value(&self) -> usize {
        64
    }
    fn sample_rate(&self) -> u8 {
        self.1
    }
    #[inline]
    fn get(&self, index: usize) -> i64 {
        self.0[index]
    }

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

/// Reads until `buffer` is full or the input is exhausted, since one `read` may return less than
/// asked for. Returns `(input_exhausted, bytes_read)`.
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

/// Writes entries as little-endian `i64`s, the uncompressed body format.
fn write_vec_i64(vec: Vec<i64>, writer: &mut impl Write) -> std::io::Result<()> {
    for value in vec {
        writer.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

/// Bytes pulled out of the reader per iteration of the refill loop below.
///
/// Must stay at or above `binary_traits::load_owned`'s `BufReader` capacity, for the reason spelled
/// out at `bitarray::binary`'s constant of the same name: below it, `BufReader::read` declines to
/// bypass its own buffer and every byte of the body is copied twice. The two constants are separate
/// only because the two `fill_buffer` helpers are; they belong in step.
const READ_CHUNK_BYTES: usize = 4 << 20;

/// Reads little-endian `i64`s into `vec` until the reader is exhausted, in [`READ_CHUNK_BYTES`]
/// chunks so a UniProt-scale body does not need a second copy of itself in memory. Callers bound
/// the reader to the body they want; a trailing partial entry is dropped.
fn read_vec_i64(vec: &mut Vec<i64>, mut reader: impl BufRead) -> std::io::Result<()> {
    vec.clear();
    let mut buffer = vec![0; READ_CHUNK_BYTES];

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

/// Writes an uncompressed suffix array: the shared header at 64 bits per value, then the entries
/// as plain little-endian `i64`s. See the [module docs](super) for the layout.
pub fn dump_suffix_array(sa: Vec<i64>, sparseness_factor: u8, writer: &mut impl Write) -> Result<(), Box<dyn Error>> {
    write_sa_header(64, sparseness_factor, sa.len(), writer)?;
    write_vec_i64(sa, writer).map_err(|_| "Could not write the suffix array to the writer")?;
    Ok(())
}

/// Loads an uncompressed suffix array body, the header having already been read.
///
/// Reads exactly the `size` entries the header declared: a body holding fewer is an error, and one
/// holding more — a header that disagrees with its file — stops at `size` rather than returning
/// entries nothing accounted for.
pub(super) fn load_original(reader: &mut impl BufRead, size: usize) -> Result<Vec<i64>, Box<dyn Error>> {
    // `size` is eight bytes straight out of the file header, so it is related to something real
    // *before* it becomes an allocation: `checked_mul` rejects a count whose body would overflow,
    // and `try_reserve_exact` reports an impossible allocation instead of aborting the process the
    // way `Vec::with_capacity` does. `read_binary_mmap` gets there by never allocating at all.
    let body_bytes = size.checked_mul(8).ok_or("The SA header declares too many items")? as u64;
    let mut sa: Vec<i64> = Vec::new();
    sa.try_reserve_exact(size)
        .map_err(|_| "The SA header declares more entries than can be allocated")?;
    // Before `read_vec_i64` touches a page of it: the advice only shapes the faults that populate
    // the buffer. `read_vec_i64` clears and refills within this capacity, so it never reallocates.
    bitarray::hugepages::advise_capacity(&sa);
    read_vec_i64(&mut sa, reader.take(body_bytes))
        .map_err(|_| "Could not read the suffix array from the binary file")?;
    if sa.len() != size {
        return Err(format!("The SA header declares {} entries but the file holds {}", size, sa.len()).into());
    }
    Ok(sa)
}

#[cfg(test)]
mod tests {
    use super::{
        super::test_utils::{FailingReader, FailingWriter},
        *
    };

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
    fn test_dump_suffix_array() {
        let mut buffer = Vec::new();
        dump_suffix_array(vec![1, 2, 3, 4, 5], 1, &mut buffer).unwrap();
        assert_eq!(buffer, vec![
            64, 1, 5, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4,
            0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
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
        // 1 call for the width, 1 for the sparseness factor, 8 for the count field, since the
        // writer claims one byte per call.
        dump_suffix_array(vec![1], 1, &mut FailingWriter { valid_write_count: 10 }).unwrap();
    }
}
