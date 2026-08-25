use std::{
    error::Error,
    io::{BufRead, Write}
};

use binary_traits::WriteBinary;
use bitarray::{Binary, DynBitArray, data_to_writer};

use super::super::{SuffixArrayBackend, write_sa_header};

/// Suffix array in owned memory, bit-packed at the minimum width the text length requires.
///
/// Field 0 is the packed entries; field 1 is the sparseness factor the array was built with.
pub struct CompressedSA(pub DynBitArray, pub u8);

impl SuffixArrayBackend for CompressedSA {
    type RangeIter<'a> = bitarray::DynBitArrayRangeIter<'a>;

    fn len(&self) -> usize {
        self.0.len()
    }
    fn bits_per_value(&self) -> usize {
        self.0.bits_per_value()
    }
    fn sample_rate(&self) -> u8 {
        self.1
    }
    #[inline]
    fn get(&self, index: usize) -> i64 {
        self.0.get(index) as i64
    }

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
    /// Goes through `write_sa_header` rather than emitting the ten bytes itself, so that this and
    /// [`dump_compressed_suffix_array`] — the two routes to a compressed file — cannot drift apart
    /// in layout or in what they report when a write fails.
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let CompressedSA(bit_array, sample_rate) = self;
        write_sa_header(bit_array.bits_per_value(), sample_rate, bit_array.len(), writer)?;
        bit_array
            .write_binary(writer)
            .map_err(|_| "Could not write the compressed suffix array to the writer")?;
        Ok(())
    }
}

// ── I/O helpers ──────────────────────────────────────────────────────────────

/// Writes a bit-packed suffix array: the shared header at `bits_per_value` bits per value — see
/// the [module docs](super) for the layout — then the body packed by `bitarray`. The body goes out
/// in chunks so peak memory stays bounded at index-build scale.
///
/// # Errors
///
/// Returns an error if any value does not fit in `bits_per_value`. That check is the difference
/// between a build that fails and one that silently produces a wrong index: `bitarray::set` now
/// masks an over-wide value to its field, so without this the offending entry would be written
/// truncated — a plausible-looking but wrong suffix position — and the caller would never learn.
/// (Before the mask it was worse still: the excess bits spilled into the *preceding* entry,
/// corrupting a value that was already correct.)
///
/// `sa-builder` derives `bits_per_value` from the text length, so it cannot trip this; the check is
/// for every other caller of a `pub` function.
pub fn dump_compressed_suffix_array(
    sa: Vec<i64>,
    sparseness_factor: u8,
    bits_per_value: usize,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    super::super::check_bits_per_value(bits_per_value)?;
    let max_representable: i64 = if bits_per_value >= 63 { i64::MAX } else { (1i64 << bits_per_value) - 1 };
    if let Some((index, &value)) = sa.iter().enumerate().find(|&(_, &v)| v < 0 || v > max_representable) {
        return Err(format!(
            "suffix array entry {index} is {value}, which does not fit in {bits_per_value} bits \
             (representable range 0..={max_representable})"
        )
        .into());
    }

    write_sa_header(bits_per_value, sparseness_factor, sa.len(), writer)?;
    data_to_writer(sa, bits_per_value, 8 * 1024, writer)
        .map_err(|_| "Could not write the compressed suffix array to the writer")?;
    Ok(())
}

/// Loads a compressed suffix array whose width has already been read, taking the sparseness factor
/// and item count that follow it before the body.
pub fn load_compressed_suffix_array(
    reader: &mut impl BufRead,
    bits_per_value: usize
) -> Result<CompressedSA, Box<dyn Error>> {
    let mut sample_rate_buffer = [0_u8; 1];
    reader
        .read_exact(&mut sample_rate_buffer)
        .map_err(|_| "Could not read the sample rate from the binary file")?;
    let sample_rate = sample_rate_buffer[0];

    let mut size_buffer = [0_u8; 8];
    reader
        .read_exact(&mut size_buffer)
        .map_err(|_| "Could not read the size of the suffix array from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    let sa = load_compressed(reader, bits_per_value, size)?;
    Ok(CompressedSA(sa, sample_rate))
}

/// Loads a compressed body, the header having already been read. Unlike its uncompressed
/// counterpart the entry count is a capacity rather than a bound: `DynBitArray::read_binary`
/// consumes the reader to its end.
pub(super) fn load_compressed(
    reader: &mut impl BufRead,
    bits_per_value: usize,
    size: usize
) -> Result<DynBitArray, Box<dyn Error>> {
    // `size` is eight bytes straight out of the file header, so it is allocated through the
    // fallible constructor: a corrupt count that overflows `size * bits_per_value`, or that simply
    // asks for more memory than exists, becomes a load error rather than an aborted process.
    // `read_binary_mmap` reaches the same conclusion by never allocating at all.
    let mut compressed_suffix_array = DynBitArray::try_with_capacity(size, bits_per_value)
        .ok_or("The SA header declares more entries than can be allocated")?;
    compressed_suffix_array
        .read_binary(reader)
        .map_err(|_| "Could not read the compressed suffix array from the binary file")?;

    // `read_binary` refills the backing store with however many words the reader yielded, which
    // says nothing about the count the header declared. Without this, a truncated `sa.bin` loads
    // cleanly, `len()` reports the declared count, and the first binary-search probe past the real
    // data panics inside a request handler. `load_original` and `read_binary_mmap` both check the
    // equivalent; this was the one sibling that did not.
    if compressed_suffix_array.word_len() < compressed_suffix_array.required_words() {
        return Err(format!(
            "The SA header declares {} entries at {} bits ({} words) but the file holds {} words",
            size,
            bits_per_value,
            compressed_suffix_array.required_words(),
            compressed_suffix_array.word_len()
        )
        .into());
    }
    // The huge-page advice is not issued here: `DynBitArray::with_capacity` already did it, on the
    // untouched allocation, which is the only point at which it does anything. See
    // `bitarray::hugepages`.
    Ok(compressed_suffix_array)
}

#[cfg(test)]
mod tests {

    /// A value too wide for the declared width must fail the write, not be packed truncated.
    ///
    /// `bitarray::set` masks the value to its field, so an over-wide entry would otherwise be
    /// written as a different, plausible-looking suffix position and the caller would never learn.
    #[test]
    fn a_value_wider_than_the_declared_width_is_rejected() {
        // 2^10 - 1 is the widest value 10 bits can hold.
        let mut buf = Vec::new();
        assert!(
            dump_compressed_suffix_array(vec![0, 5, 1023], 1, 10, &mut buf).is_ok(),
            "values inside the width must be written"
        );

        let mut buf = Vec::new();
        let err =
            dump_compressed_suffix_array(vec![0, 5, 1024], 1, 10, &mut buf).expect_err("1024 does not fit in 10 bits");
        let msg = err.to_string();
        assert!(msg.contains("entry 2"), "the error should name the offending entry: {msg}");
        assert!(msg.contains("10 bits"), "the error should name the width: {msg}");

        // A negative entry cannot be represented either.
        let mut buf = Vec::new();
        assert!(dump_compressed_suffix_array(vec![0, -1, 3], 1, 10, &mut buf).is_err());
    }
    use super::{
        super::test_utils::{FailingReader, FailingWriter},
        *
    };

    #[test]
    fn test_dump_compressed_suffix_array() {
        let sa = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut writer = vec![];
        dump_compressed_suffix_array(sa, 1, 8, &mut writer).unwrap();
        assert_eq!(writer, vec![8, 1, 10, 0, 0, 0, 0, 0, 0, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9]);
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
        // 1 call for the width, 1 for the sparseness factor, 8 for the count field, since the
        // writer claims one byte per call.
        dump_compressed_suffix_array(vec![1], 1, 8, &mut FailingWriter { valid_write_count: 10 }).unwrap();
    }

    #[test]
    fn test_load_compressed_suffix_array() {
        let data = [1, 10, 0, 0, 0, 0, 0, 0, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 10, 9];
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
}
