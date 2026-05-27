// Non-mmap build: in-memory protein text backed by a BitArray.
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, Read, Write};

use bitarray::{Binary, BitArray, data_to_writer};

use binary_traits::{ReadBinary, WriteBinary};
use crate::{bit_array_byte_size, BIT5_TO_CHAR, ProteinTextBackend};

// ── InMemoryProteinText ───────────────────────────────────────────────────────

pub struct InMemoryProteinText {
    pub(crate) bit_array: BitArray,
    pub(crate) char_to_5bit: HashMap<u8, u8>,
}

impl InMemoryProteinText {
    pub(crate) fn create_char_to_5bit_hashmap() -> HashMap<u8, u8> {
        BIT5_TO_CHAR.iter().enumerate().map(|(i, &c)| (c, i as u8)).collect()
    }

    pub fn from_string(input_string: &str) -> Self {
        let char_to_5bit = Self::create_char_to_5bit_hashmap();
        let mut bit_array = BitArray::with_capacity(input_string.len(), 5);
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(&(c as u8))
                .unwrap_or_else(|| panic!("Input character '{}' not in alphabet", c));
            bit_array.set(i, char_5bit as u64);
        }
        Self { bit_array, char_to_5bit }
    }

    pub fn from_vec(input_vec: &[u8]) -> Self {
        let char_to_5bit = Self::create_char_to_5bit_hashmap();
        let mut bit_array = BitArray::with_capacity(input_vec.len(), 5);
        for (i, e) in input_vec.iter().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(e)
                .unwrap_or_else(|| panic!("Input character '{}' not in alphabet", e));
            bit_array.set(i, char_5bit as u64);
        }
        Self { bit_array, char_to_5bit }
    }

    pub fn new(bit_array: BitArray) -> Self {
        Self { bit_array, char_to_5bit: Self::create_char_to_5bit_hashmap() }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(BitArray::with_capacity(capacity, 5))
    }

    pub fn set(&mut self, index: usize, value: u8) {
        let char_5bit: u8 = *self.char_to_5bit.get(&value)
            .unwrap_or_else(|| panic!("Input character '{}' not in alphabet", value));
        self.bit_array.set(index, char_5bit as u64);
    }

    pub fn clear(&mut self) { self.bit_array.clear(); }
}

impl ProteinTextBackend for InMemoryProteinText {
    #[inline]
    fn get(&self, index: usize) -> u8 {
        BIT5_TO_CHAR[self.bit_array.get(index) as usize]
    }

    fn len(&self) -> usize { self.bit_array.len() }

    #[inline]
    fn prefetch_at(&self, index: usize) {
        if index < self.bit_array.len() {
            let word_idx = (index * 5) / 64;
            let ptr: *const u64 = self.bit_array.get_data_slice(word_idx, word_idx + 1).as_ptr();
            prefetch::prefetch_read(ptr);
        }
    }
}

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
        let mut bit_array = BitArray::with_capacity(text_length, 5);
        let mut limited = <&mut R as Read>::take(reader, n_bytes as u64);
        bit_array
            .read_binary(&mut limited)
            .map_err(|_| "Could not parse BitArray data from binary file")?;

        Ok(Self::new(bit_array))
    }
}

// ── I/O helpers ──────────────────────────────────────────────────────────────

pub fn dump_compressed_text(text: Vec<u8>, writer: &mut impl Write) -> Result<(), Box<dyn Error>> {
    let bits_per_value = 5;
    writer.write(&[bits_per_value as u8]).map_err(|_| "Could not write the required bits to the writer")?;
    writer.write(&(text.len() as u64).to_le_bytes()).map_err(|_| "Could not write the size of the text to the writer")?;
    let text_writer: Vec<i64> = text.iter().map(|item| <i64>::from(*item)).collect();
    data_to_writer(text_writer, bits_per_value, 8 * 1024, writer)
        .map_err(|_| "Could not write the compressed text to the writer")?;
    Ok(())
}

pub fn load_compressed_text(reader: &mut impl BufRead) -> Result<InMemoryProteinText, Box<dyn Error>> {
    let bits_per_value: usize = 5;
    let mut size_buffer = [0_u8; 8];
    reader.read_exact(&mut size_buffer).map_err(|_| "Could not read the size of the text from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;
    let mut compressed_text = BitArray::with_capacity(size, bits_per_value);
    compressed_text.read_binary(reader).map_err(|_| "Could not read the compressed text from the binary file")?;
    Ok(InMemoryProteinText::new(compressed_text))
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::io::{Read, BufRead};
    use super::*;

    pub struct FailingWriter { pub valid_write_count: usize }
    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> Result<usize, std::io::Error> {
            if self.valid_write_count == 0 { return Err(std::io::Error::other("Write failed")); }
            self.valid_write_count -= 1; Ok(1)
        }
        fn flush(&mut self) -> Result<(), std::io::Error> { Ok(()) }
    }

    pub struct FailingReader { pub valid_read_count: usize }
    impl Read for FailingReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.valid_read_count == 0 { return Err(std::io::Error::other("Read failed")); }
            self.valid_read_count -= 1; Ok(buf.len())
        }
    }
    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> { Ok(&[]) }
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
        let mut bit_array = BitArray::with_capacity(input_string.len(), 5);
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
        for (i, c) in "ACACA-CAC$".chars().enumerate() { text.set(i, c as u8); }
        for (i, c) in "ACACA-CAC$".chars().enumerate() { assert_eq!(c as u8, text.get(i)); }
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
        assert!(text_slice.equals_slice(&[b'C', b'I', b'C', b'A'], false));
        assert!(!text_slice.equals_slice(&[b'C', b'C', b'C', b'A'], false));
        assert!(text_slice.equals_slice(&[b'C', b'L', b'C', b'A'], true));
    }

    #[test]
    fn test_check_il_locations() {
        let text = InMemoryProteinText::from_string("ACILA-CAC$");
        let text_slice = text.slice(1, 5);
        let il_locations = [1, 2];
        assert!(text_slice.check_il_locations(0, &il_locations, &[b'C', b'I', b'L', b'A']));
        assert!(!text_slice.check_il_locations(0, &il_locations, &[b'C', b'I', b'C', b'A']));
    }

    #[test]
    fn test_dump_compressed_text() {
        let text: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut writer = vec![];
        dump_compressed_text(text, &mut writer).unwrap();
        assert_eq!(writer, vec![5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 128, 74, 232, 152, 66, 134, 8]);
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_compressed_text_fail_required_bits() {
        dump_compressed_text(vec![], &mut FailingWriter { valid_write_count: 0 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the text to the writer")]
    fn test_dump_compressed_text_fail_size() {
        dump_compressed_text(vec![], &mut FailingWriter { valid_write_count: 1 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the compressed text to the writer")]
    fn test_dump_compressed_text_fail_compressed_text() {
        dump_compressed_text(vec![1], &mut FailingWriter { valid_write_count: 3 }).unwrap();
    }

    #[test]
    fn test_load_compressed_text() {
        let data = [10, 0, 0, 0, 0, 0, 0, 0, 0, 128, 74, 232, 152, 66, 134, 8];
        let mut reader = std::io::BufReader::new(&data[..]);
        let compressed_text = load_compressed_text(&mut reader).unwrap();
        for (i, c) in "BCDEFGHIKL".chars().enumerate() {
            assert_eq!(compressed_text.get(i), c as u8);
        }
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the text from the binary file")]
    fn test_load_compressed_text_fail_size() {
        load_compressed_text(&mut FailingReader { valid_read_count: 0 }).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the compressed text from the binary file")]
    fn test_load_compressed_text_fail_compressed_text() {
        load_compressed_text(&mut FailingReader { valid_read_count: 2 }).unwrap();
    }

    #[test]
    fn test_write_and_read_binary() {
        let input = "ACACA-CAC$";
        let text = InMemoryProteinText::from_string(input);
        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();
        let mut reader = std::io::BufReader::new(buf.as_slice());
        let loaded = InMemoryProteinText::read_binary(&mut reader).unwrap();
        for (i, c) in input.chars().enumerate() { assert_eq!(loaded.get(i), c as u8); }
        assert_eq!(loaded.len(), input.len());
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
