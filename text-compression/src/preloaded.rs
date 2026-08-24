//! In-memory protein text: the whole text decompressed into owned RAM.
//!
//! This module owns the `WriteBinary` implementation, so `sa-builder` uses it to produce the file
//! that *either* backend later reads — which is why the builder never names a backend at all.
//!
//! See [`crate::mmap`] for the counterpart that decodes straight out of a mapping.
use std::{
    collections::HashMap,
    error::Error,
    io::{BufRead, Read, Write}
};

use binary_traits::{ReadBinary, WriteBinary};
use bitarray::{Binary, BitArray};

use crate::{BIT5_TO_CHAR, ProteinTextBackend, bit_array_byte_size};

// ── InMemoryProteinText ───────────────────────────────────────────────────────

/// The protein text held in owned memory, packed at 5 bits per residue.
///
/// Carries the ASCII → 5-bit table alongside the packed data so that `set` can encode; decoding
/// needs only [`crate::BIT5_TO_CHAR`].
pub struct InMemoryProteinText {
    pub(crate) bit_array: BitArray<5>,
    pub(crate) char_to_5bit: HashMap<u8, u8>
}

impl InMemoryProteinText {
    pub(crate) fn create_char_to_5bit_hashmap() -> HashMap<u8, u8> {
        BIT5_TO_CHAR.iter().enumerate().map(|(i, &c)| (c, i as u8)).collect()
    }

    /// Encodes `input_string`, one residue per character.
    ///
    /// # Panics
    ///
    /// If any character is outside the alphabet in [`crate::BIT5_TO_CHAR`].
    pub fn from_string(input_string: &str) -> Self {
        let char_to_5bit = Self::create_char_to_5bit_hashmap();
        let mut bit_array = BitArray::<5>::with_capacity(input_string.len());
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 =
                *char_to_5bit.get(&(c as u8)).unwrap_or_else(|| panic!("Input character '{}' not in alphabet", c));
            bit_array.set(i, char_5bit as u64);
        }
        Self { bit_array, char_to_5bit }
    }

    /// Encodes `input_vec`, one residue per byte. Same alphabet constraint as
    /// [`Self::from_string`].
    ///
    /// # Panics
    ///
    /// If any byte is outside the alphabet in [`crate::BIT5_TO_CHAR`].
    pub fn from_vec(input_vec: &[u8]) -> Self {
        let char_to_5bit = Self::create_char_to_5bit_hashmap();
        let mut bit_array = BitArray::<5>::with_capacity(input_vec.len());
        for (i, e) in input_vec.iter().enumerate() {
            let char_5bit: u8 =
                *char_to_5bit.get(e).unwrap_or_else(|| panic!("Input character '{}' not in alphabet", e));
            bit_array.set(i, char_5bit as u64);
        }
        Self { bit_array, char_to_5bit }
    }

    /// Wraps an already-packed bit array, as produced by `read_binary`.
    pub fn new(bit_array: BitArray<5>) -> Self {
        Self { bit_array, char_to_5bit: Self::create_char_to_5bit_hashmap() }
    }

    /// Allocates room for `capacity` residues, all decoding to the first alphabet entry.
    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(BitArray::<5>::with_capacity(capacity))
    }

    /// Encodes and stores the ASCII residue `value` at `index`.
    ///
    /// # Panics
    ///
    /// If `value` is outside the alphabet, or `index` is out of bounds.
    pub fn set(&mut self, index: usize, value: u8) {
        let char_5bit: u8 = *self
            .char_to_5bit
            .get(&value)
            .unwrap_or_else(|| panic!("Input character '{}' not in alphabet", value));
        self.bit_array.set(index, char_5bit as u64);
    }

    /// Zeroes the text without reallocating.
    pub fn clear(&mut self) {
        self.bit_array.clear();
    }
}

impl ProteinTextBackend for InMemoryProteinText {
    #[inline]
    fn get(&self, index: usize) -> u8 {
        BIT5_TO_CHAR[self.bit_array.get(index) as usize]
    }

    #[inline]
    fn len(&self) -> usize {
        self.bit_array.len()
    }

    /// Prefetches the backing word holding `index`.
    ///
    /// Note the index conversion: `get_data_slice` is indexed by `u64` *word*, not by residue,
    /// so the residue index has to be scaled by the 5-bit width first. Out-of-range indices are
    /// skipped rather than clamped, per the trait contract.
    #[inline]
    fn prefetch_at(&self, index: usize) {
        if index < self.bit_array.len() {
            let word_idx = (index * 5) / 64;
            let ptr: *const u64 = self.bit_array.get_data_slice(word_idx, word_idx + 1).as_ptr();
            prefetch::prefetch_read(ptr);
        }
    }
}

/// On-disk format for the protein text — written here, read by both backends.
///
/// ```text
/// [ text_length: u64 little-endian ][ packed residues: ceil(len*5/64) * 8 bytes ]
/// ```
///
/// The payload is exactly `BitArray<5>`'s backing words, so it is packed
/// most-significant-bit-first within each little-endian `u64` and values may straddle a word
/// boundary. `crate::bit_array_byte_size` computes the payload length and is what the mmap
/// readers bounds-check against.
impl WriteBinary for InMemoryProteinText {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let text_length = self.bit_array.len() as u64;
        writer.write_all(&text_length.to_le_bytes())?;
        self.bit_array.write_binary(writer)?;
        Ok(())
    }
}

impl ReadBinary for InMemoryProteinText {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8).map_err(|_| "Could not read text_length from binary file")?;
        let text_length = u64::from_le_bytes(buf8) as usize;

        let n_bytes = bit_array_byte_size(text_length);
        let mut bit_array = BitArray::<5>::with_capacity(text_length);
        let mut limited = <&mut R as Read>::take(reader, n_bytes as u64);
        bit_array.read_binary(&mut limited).map_err(|_| "Could not parse BitArray data from binary file")?;
        // The huge-page advice is not issued here: `BitArray::with_capacity` already did it, on the
        // untouched allocation, which is the only point at which it does anything. See
        // `bitarray::hugepages`.

        Ok(Self::new(bit_array))
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::io::{BufRead, Read};

    use super::*;

    pub struct FailingWriter {
        pub valid_write_count: usize
    }
    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> Result<usize, std::io::Error> {
            if self.valid_write_count == 0 {
                return Err(std::io::Error::other("Write failed"));
            }
            self.valid_write_count -= 1;
            Ok(1)
        }
        fn flush(&mut self) -> Result<(), std::io::Error> {
            Ok(())
        }
    }

    pub struct FailingReader {
        pub valid_read_count: usize
    }
    impl Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.valid_read_count == 0 {
                return Err(std::io::Error::other("Read failed"));
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
    fn test_u8_5bit_conversion() {
        let char_to_5bit = InMemoryProteinText::create_char_to_5bit_hashmap();
        for (i, &c) in BIT5_TO_CHAR.iter().enumerate() {
            let idx = *char_to_5bit.get(&c).unwrap() as usize;
            assert_eq!(i, idx);
            assert_eq!(BIT5_TO_CHAR[idx], c);
        }
    }

    #[test]
    fn test_build_from_string() {
        let text = InMemoryProteinText::from_string("ACACA-CAC$");
        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_from_vec() {
        let vec = vec![b'A', b'C', b'A', b'C', b'A', b'-', b'C', b'A', b'C', b'$'];
        let text = InMemoryProteinText::from_vec(&vec);
        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_from_bitarray() {
        let input_string = "ACACA-CAC$";
        let char_to_5bit = InMemoryProteinText::create_char_to_5bit_hashmap();
        let mut bit_array = BitArray::<5>::with_capacity(input_string.len());
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(&(c as u8)).unwrap();
            bit_array.set(i, char_5bit as u64);
        }
        let text = InMemoryProteinText::new(bit_array);
        for (i, c) in input_string.chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_with_capacity() {
        let input_string = "ACACA-CAC$";
        let mut text = InMemoryProteinText::with_capacity(input_string.len());
        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            text.set(i, c as u8);
        }
        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_text_slice() {
        let text = InMemoryProteinText::from_string("ACACA-CAC$");
        let text_slice = text.slice(1, 5);
        for (i, c) in "ACACA-CAC$"[1..5].chars().enumerate() {
            assert_eq!(c as u8, text_slice.get(i));
        }
    }

    #[test]
    fn test_equals_slice() {
        let text = InMemoryProteinText::from_string("ACICA-CAC$");
        let text_slice = text.slice(1, 5);
        assert!(text_slice.equals_slice(&b"CICA"[..], false));
        assert!(!text_slice.equals_slice(&b"CCCA"[..], false));
        assert!(text_slice.equals_slice(&b"CLCA"[..], true));
    }

    #[test]
    fn test_check_il_locations() {
        let text = InMemoryProteinText::from_string("ACILA-CAC$");
        let text_slice = text.slice(1, 5);
        let il_locations = [1, 2];
        assert!(text_slice.check_il_locations(0, &il_locations, &b"CILA"[..]));
        assert!(!text_slice.check_il_locations(0, &il_locations, &b"CICA"[..]));
    }

    #[test]
    fn test_write_and_read_binary() {
        let input = "ACACA-CAC$";
        let text = InMemoryProteinText::from_string(input);
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();
        let mut reader = std::io::BufReader::new(buf.as_slice());
        let loaded = InMemoryProteinText::read_binary(&mut reader).unwrap();
        for (i, c) in input.chars().enumerate() {
            assert_eq!(loaded.get(i), c as u8);
        }
        assert_eq!(loaded.len(), input.len());
    }

    /// The preloaded half of the hint contract; the mapped half is the twin of this test.
    #[test]
    fn prefetch_hints_are_harmless() {
        let input = "ACACA-CAC$MLPGLALLLL$";
        let text = InMemoryProteinText::from_string(input);
        crate::test_utils::assert_prefetch_is_harmless(&text, input);
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
