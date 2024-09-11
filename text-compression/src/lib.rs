use std::{
    error::Error,
    io::{BufRead, Write}
};
use std::collections::HashMap;

use bitarray::{data_to_writer, Binary, BitArray};

pub struct ProteinText {
    bit_array: BitArray,
    char_to_5bit: HashMap<u8, u8>,
    bit5_to_char: Vec<u8>,
}

impl ProteinText {

    fn create_char_to_5bit_hashmap() -> HashMap<u8, u8> {
        let mut hashmap = HashMap::<u8, u8>::new();
        for (i, c) in "ACDEFGHIKLMNPQRSTVWY-$".chars().enumerate() {
            hashmap.insert(c as u8, i as u8);
        }

        hashmap
    }

    fn create_bit5_to_char() -> Vec<u8> {
        let mut vec = Vec::<u8>::new();
        for c in "ACDEFGHIKLMNPQRSTVWY-$".chars() {
            vec.push(c as u8);
        }
        vec
    }
    
    pub fn from_string(input_string: &str) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        let mut bit_array = BitArray::with_capacity(input_string.len(), 5);
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(&(c as u8)).expect("Input character not in alphabet");
            bit_array.set(i, char_5bit as u64);
        }

        Self { bit_array, char_to_5bit, bit5_to_char }
    }

    pub fn from_vec(input_vec: &Vec<u8>) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        let mut bit_array = BitArray::with_capacity(input_vec.len(), 5);
        for (i, e) in input_vec.iter().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(e).expect("Input character not in alphabet");
            bit_array.set(i, char_5bit as u64);
        }

        Self { bit_array, char_to_5bit, bit5_to_char }
    }

    pub fn new(bit_array: BitArray) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();
        Self { bit_array, char_to_5bit, bit5_to_char }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(BitArray::with_capacity(capacity, 5))
    }

    pub fn get(&self, index: usize) -> u8 {
        let char_5bit = self.bit_array.get(index) as usize;
        self.bit5_to_char[char_5bit]
    }

    pub fn set(&mut self, index: usize, value: u8) {
        let char_5bit: u8 = *self.char_to_5bit.get(&value).expect("Input character not in alphabet");
        self.bit_array.set(index, char_5bit as u64);
    }

    pub fn len(&self) -> usize {
        self.bit_array.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bit_array.len() == 0
    }

    /// Clears the `BitArray`, setting all bits to 0.
    pub fn clear(&mut self) {
        self.bit_array.clear()
    }

    pub fn iter(&self) -> ProteinTextIterator {
        ProteinTextIterator {protein_text: self, index: 0, }
    }

}

pub struct ProteinTextSlice<'a> {
    text: &'a ProteinText,
    start: usize, // included
    end: usize,   // excluded
}

impl<'a> ProteinTextSlice<'a> {

    pub fn new(text: &'a ProteinText, start: usize, end: usize) -> ProteinTextSlice {
        Self {text, start, end }
    }

    pub fn get(&self, index: usize) -> u8 {
        self.text.get(self.start + index)
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    #[inline]
    pub fn equals_slice(&self, other: &[u8], equate_il: bool) -> bool {
        if equate_il {
            other.iter().zip(self.iter()).all(|(&search_character, text_character)| {
                search_character == text_character
                    || (search_character == b'I' && text_character == b'L')
                    || (search_character == b'L' && text_character == b'I')
            })
        } else {
            other.iter().zip(self.iter()).all(|(&search_character, text_character)| search_character == text_character)
        }
    }

    pub fn check_il_locations(
        &self,
        skip: usize,
        il_locations: &[usize],
        search_string: &[u8],
    ) -> bool {
        for &il_location in il_locations {
            let index = il_location - skip;
            if search_string[index] != self.get(index) {
                return false;
            }
        }
        true
    }

    pub fn iter(&self) -> ProteinTextSliceIterator {
        ProteinTextSliceIterator {text_slice: self, index: 0, }
    }
}

pub struct ProteinTextIterator<'a> {
    protein_text: &'a ProteinText,
    index: usize,
}

pub struct ProteinTextSliceIterator<'a> {
    text_slice: &'a ProteinTextSlice<'a>,
    index: usize,
}

impl<'a> Iterator for ProteinTextSliceIterator<'a> {

    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.text_slice.len() {
            return None;
        }

        self.index += 1;
        Some(self.text_slice.get(self.index - 1))
    }
}

impl<'a> Iterator for ProteinTextIterator<'a> {

    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.protein_text.len() {
            return None;
        }

        self.index += 1;
        Some(self.protein_text.get(self.index - 1))
    }
}

/// Writes the compressed text to a writer.
///
/// # Arguments
///
/// * `text` - The text to be compressed.
/// * `writer` - The writer to which the compressed text will be written.
///
/// # Errors
///
/// Returns an error if writing to the writer fails.
pub fn dump_compressed_text(
    text: Vec<u8>,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    let bits_per_value = 5;

    // Write the flags to the writer
    // 00000001 indicates that the text is compressed
    writer
        .write(&[bits_per_value as u8])
        .map_err(|_| "Could not write the required bits to the writer")?;

    // Write the size of the text to the writer
    writer
        .write(&(text.len() as u64).to_le_bytes())
        .map_err(|_| "Could not write the size of the text to the writer")?;

    // Compress the text and write it to the writer
    let text_writer: Vec<i64> = text.iter().map(|item| <i64>::from(*item)).collect();
    data_to_writer(text_writer, bits_per_value, 8 * 1024, writer)
        .map_err(|_| "Could not write the compressed text to the writer")?;

    Ok(())
}

/// Load the compressed text from a reader.
///
/// # Arguments
///
/// * `reader` - The reader from which the compressed text will be read.
///
/// # Errors
///
/// Returns an error if reading from the reader fails.
pub fn load_compressed_text(
    reader: &mut impl BufRead
) -> Result<ProteinText, Box<dyn Error>> {
    let bits_per_value: usize = 5;
    // Read the size of the text from the binary file (8 bytes)
    let mut size_buffer = [0_u8; 8];
    reader
        .read_exact(&mut size_buffer)
        .map_err(|_| "Could not read the size of the text from the binary file")?;
    let size = u64::from_le_bytes(size_buffer) as usize;

    // Read the compressed text from the binary file
    let mut compressed_text = BitArray::with_capacity(size, bits_per_value);
    compressed_text
        .read_binary(reader)
        .map_err(|_| "Could not read the compressed text from the binary file")?;

    Ok(ProteinText::new(compressed_text))
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
    fn test_dump_compressed_text() {
        let text: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        let mut writer = vec![];
        dump_compressed_text(text, &mut writer).unwrap();

        assert_eq!(writer, vec![
            // bits per value
            5, // size of the text
            10, 0, 0, 0, 0, 0, 0, 0, // compressed text
            0, 128, 74, 232, 152, 66, 134, 8
        ]);
    }

    #[test]
    #[should_panic(expected = "Could not write the required bits to the writer")]
    fn test_dump_compressed_text_fail_required_bits() {
        let mut writer = FailingWriter { valid_write_count: 0 };

        dump_compressed_text(vec![], &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the size of the text to the writer")]
    fn test_dump_compressed_text_fail_size() {
        let mut writer = FailingWriter { valid_write_count: 1 };

        dump_compressed_text(vec![], &mut writer).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not write the compressed text to the writer")]
    fn test_dump_compressed_text_fail_compressed_text() {
        let mut writer = FailingWriter { valid_write_count: 3 };

        dump_compressed_text(vec![1], &mut writer).unwrap();
    }

    #[test]
    fn test_load_compressed_text() {
        let data = vec![
             // size of the text
            10, 0, 0, 0, 0, 0, 0, 0, // compressed text
            0, 128, 74, 232, 152, 66, 134, 8
        ];

        let mut reader = std::io::BufReader::new(&data[..]);
        let compressed_text = load_compressed_text(&mut reader).unwrap();

        for (i, c) in "CDEFGHIKLM".chars().enumerate() {
            assert_eq!(compressed_text.get(i), c as u8);
        }
    }

    #[test]
    #[should_panic(expected = "Could not read the size of the text from the binary file")]
    fn test_load_compressed_text_fail_size() {
        let mut reader = FailingReader { valid_read_count: 0 };

        load_compressed_text(&mut reader).unwrap();
    }

    #[test]
    #[should_panic(expected = "Could not read the compressed text from the binary file")]
    fn test_load_compressed_text_fail_compressed_text() {
        let mut reader = FailingReader { valid_read_count: 2 };

        load_compressed_text(&mut reader).unwrap();
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
