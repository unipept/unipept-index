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
    /// The input is assumed to be ASCII, and nothing here enforces it: the allocation is sized
    /// from `str::len` (bytes) while the loop walks `chars()` (code points), and each `char` is
    /// narrowed with `as u8`. Non-ASCII input therefore does *not* reliably hit the panic below —
    /// `'ł'` (U+0142) narrows to `0x42`, a perfectly valid `'B'`, and is encoded silently. Use
    /// [`Self::from_vec`] for anything not known to be ASCII.
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
        // Not redundant with `bitarray`'s own bounds check: that one indexes the backing *word*
        // slice, so it lets an index in the final word's zero padding through. See the trait.
        debug_assert!(index < self.bit_array.len(), "residue index {index} is past the end of the text");
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

/// Reads the format documented on the [`WriteBinary`] impl above, validating the header against
/// what the file actually holds.
///
/// [`crate::MmapBackedProteinText`]'s `read_binary_mmap` is the sibling that must reject exactly
/// the same files; `tests::both_backends_reject_the_same_damaged_text_files` pins the two
/// together.
impl ReadBinary for InMemoryProteinText {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8).map_err(|_| "Could not read text_length from binary file")?;
        let text_length = u64::from_le_bytes(buf8) as usize;

        let n_bytes =
            bit_array_byte_size(text_length).ok_or("The protein text header declares an implausible text length")?;

        // Fallibly, because `text_length` is eight bytes straight out of the file: a corrupt header
        // asking for more memory than exists becomes a load error rather than an aborted process.
        let mut bit_array = BitArray::<5>::try_with_capacity(text_length)
            .ok_or("The protein text header declares more residues than can be allocated")?;
        let mut limited = <&mut R as Read>::take(reader, n_bytes as u64);
        bit_array.read_binary(&mut limited).map_err(|_| "Could not parse BitArray data from binary file")?;
        // The huge-page advice is not issued here: `BitArray::with_capacity` already did it, on the
        // untouched allocation, which is the only point at which it does anything. See
        // `bitarray::hugepages`.

        // `read_binary` refills the backing store with however many words the reader supplied, which
        // says nothing about the length the header declared. Without this, a truncated or
        // over-declared `proteins.bin` loaded cleanly, `len()` reported the declared length, and the
        // first `get` past the real data panicked inside `bitarray` — on a live query rather than at
        // load. `read_binary_mmap` rejects both cases; this is the sibling that did not.
        if bit_array.word_len() < bit_array.required_words() {
            return Err(format!(
                "The protein text header declares {} residues ({} words) but the file holds {} words",
                text_length,
                bit_array.required_words(),
                bit_array.word_len()
            )
            .into());
        }

        Ok(Self::new(bit_array))
    }
}

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

    /// The preloaded reader must refuse the same damaged files the mmap reader refuses.
    ///
    /// `mmap::tests` has pinned truncation and over-declared headers for a while; this side had no
    /// equivalent, and no check either. A short or over-declared `proteins.bin` loaded cleanly,
    /// `len()` reported the declared length, and the first `get` past the real data panicked inside
    /// `bitarray` — on a live query rather than at load, which is exactly what the
    /// `ReadBinaryMmap` contract exists to prevent. Both backends now answer the same way.
    #[test]
    fn both_backends_reject_the_same_damaged_text_files() {
        use binary_traits::ReadBinaryMmap;

        let text = InMemoryProteinText::from_string("ACACA-CAC$MLPGLALLLL$");
        let text_len = text.len();
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();

        let map_bytes = |bytes: &[u8]| {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            std::io::Write::write_all(&mut tmp, bytes).unwrap();
            std::io::Write::flush(&mut tmp).unwrap();
            (crate::MmapBackedProteinText::read_binary_mmap(tmp.path()).is_ok(), tmp)
        };

        // Every truncation: both backends must reject every prefix.
        for cut in 0..buf.len() {
            let preloaded_ok = InMemoryProteinText::read_binary(&mut &buf[..cut]).is_ok();
            let (mmap_ok, _tmp) = map_bytes(&buf[..cut]);
            assert!(!preloaded_ok, "preloaded accepted a {cut}-byte prefix of {}", buf.len());
            assert_eq!(preloaded_ok, mmap_ok, "backends disagree on a {cut}-byte prefix");
        }

        // A header claiming far more text than the body holds.
        let mut overlong = buf.clone();
        overlong[0..8].copy_from_slice(&1_000_000_u64.to_le_bytes());
        let preloaded_ok = InMemoryProteinText::read_binary(&mut overlong.as_slice()).is_ok();
        let (mmap_ok, _tmp) = map_bytes(&overlong);
        assert!(!preloaded_ok, "preloaded accepted an over-declared header");
        assert_eq!(preloaded_ok, mmap_ok, "backends disagree on an over-declared header");

        // A header whose length overflows the byte-size computation must be an error on both,
        // not a panic inside the check itself.
        let mut absurd = buf.clone();
        absurd[0..8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(InMemoryProteinText::read_binary(&mut absurd.as_slice()).is_err());
        let (mmap_ok, _tmp) = map_bytes(&absurd);
        assert!(!mmap_ok, "mmap accepted a text length of u64::MAX");

        // The intact file loads on both and reads back, so the rejections above are specific.
        let loaded = InMemoryProteinText::read_binary(&mut buf.as_slice()).expect("intact file must load");
        assert_eq!(loaded.len(), text_len);
        let (mmap_ok, _tmp) = map_bytes(&buf);
        assert!(mmap_ok, "mmap rejected an intact file");
    }
}
