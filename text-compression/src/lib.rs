use std::{
    error::Error,
    io::{BufRead, Write}
};
use std::collections::HashMap;

use bitarray::{data_to_writer, Binary, BitArray};

/// Structure representing the proteins, stored in a bit array using 5 bits per amino acid.
pub struct ProteinText {
    /// Bit array holding the sequence of amino acids
    bit_array: BitArray,
    /// Hashmap storing the mapping between the character as `u8` and a 5 bit number.
    char_to_5bit: HashMap<u8, u8>,
    /// Vector storing the mapping between the 5 bit number and the character as `u8`.
    bit5_to_char: Vec<u8>,
}

impl ProteinText {

    /// Creates the hashmap storing the mappings between the characters as `u8` and 5 bit numbers.
    ///
    /// # Returns
    ///
    /// Returns the hashmap
    fn create_char_to_5bit_hashmap() -> HashMap<u8, u8> {
        let mut hashmap = HashMap::<u8, u8>::new();
        for (i, c) in "ACDEFGHIKLMNPQRSTVWY-$".chars().enumerate() {
            hashmap.insert(c as u8, i as u8);
        }

        hashmap
    }

    /// Creates the vector storing the mappings between the 5 bit numbers and the characters as `u8`.
    ///
    /// # Returns
    ///
    /// Returns the vector
    fn create_bit5_to_char() -> Vec<u8> {
        let mut vec = Vec::<u8>::new();
        for c in "ACDEFGHIKLMNPQRSTVWY-$".chars() {
            vec.push(c as u8);
        }
        vec
    }
    
    /// Creates the compressed text from a string.
    /// 
    /// # Arguments
    /// * `input_string` - The text (proteins) in string format
    ///
    /// # Returns
    ///
    /// An instance of `ProteinText`
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

    /// Creates the compressed text from a vector.
    /// 
    /// # Arguments
    /// * `input_vec` - The text (proteins) in a vector with elements of type `u8` representing the amino acids.
    ///
    /// # Returns
    ///
    /// An instance of `ProteinText`
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

    /// Creates the compressed text from a bit array.
    /// 
    /// # Arguments
    /// * `bit_array` - The text (proteins) in a bit array using 5 bits for each amino acid.
    ///
    /// # Returns
    ///
    /// An instance of `ProteinText`
    pub fn new(bit_array: BitArray) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();
        Self { bit_array, char_to_5bit, bit5_to_char }
    }

    /// Creates an instance of `ProteinText` with a given capacity.
    /// 
    /// # Arguments
    /// * `capacity` - The amount of characters in the text.
    ///
    /// # Returns
    ///
    /// An instance of `ProteinText`
    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(BitArray::with_capacity(capacity, 5))
    }

    /// Search the character at a given position in the compressed text.
    /// 
    /// # Arguments
    /// * `index` - The index of the character to search.
    ///
    /// # Returns
    ///
    /// the character at position `index` as `u8`.
    pub fn get(&self, index: usize) -> u8 {
        let char_5bit = self.bit_array.get(index) as usize;
        self.bit5_to_char[char_5bit]
    }

    /// Set the character at a given index.
    /// 
    /// # Arguments
    /// * `index` - The index of the character to change.
    /// * `value` - The character to fill in as `u8`.
    pub fn set(&mut self, index: usize, value: u8) {
        let char_5bit: u8 = *self.char_to_5bit.get(&value).expect("Input character not in alphabet");
        self.bit_array.set(index, char_5bit as u64);
    }

    /// Queries the length of the text.
    ///
    /// # Returns
    /// 
    /// the length of the text
    pub fn len(&self) -> usize {
        self.bit_array.len()
    }

    /// Check if the text is empty (length 0).
    ///
    /// # Returns
    /// 
    /// true if the the text has length 0, false otherwise.
    pub fn is_empty(&self) -> bool {
        self.bit_array.len() == 0
    }

    /// Clears the `BitArray`, setting all bits to 0.
    pub fn clear(&mut self) {
        self.bit_array.clear()
    }

    /// Get an iterator over the characters of the text.
    ///
    /// # Returns
    /// 
    /// A `ProteinTextIterator`, which can iterate over the characters of the text.
    pub fn iter(&self) -> ProteinTextIterator {
        ProteinTextIterator {protein_text: self, index: 0, }
    }

    /// Get a slice of the text
    ///
    /// # Returns
    /// 
    /// An `ProteinTextSlice` representing a slice of the text.
    pub fn slice(&self, start: usize, end:usize) -> ProteinTextSlice {
        ProteinTextSlice::new(self, start, end)
    }

}

/// Structure representing a slice of a `ProteinText`.
pub struct ProteinTextSlice<'a> {
    /// The `Proteintext` of whih to take a slice.
    text: &'a ProteinText,
    /// The start of the slice.
    start: usize, // included
    /// The end of the slice.
    end: usize,   // excluded
}

impl<'a> ProteinTextSlice<'a> {

    /// Creates an instance of `ProteintextSlice`, given the text and boundaries.
    /// 
    /// # Arguments
    /// * `text` - The `Proteintext` representing the text of proteins with 5 bits per amino acid.
    /// * `start` - The start of the slice.
    /// * `end` - The end of the slice.
    ///
    /// # Returns
    ///
    /// An instance of `ProteinTextSlice`
    pub fn new(text: &'a ProteinText, start: usize, end: usize) -> ProteinTextSlice {
        Self {text, start, end }
    }

    /// Get a character (amino acid) in the slice.
    /// 
    /// # Arguments
    /// * `index` - The index in the slice of the character to get.
    ///
    /// # Returns
    ///
    /// The character as `u8`.
    pub fn get(&self, index: usize) -> u8 {
        self.text.get(self.start + index)
    }

    /// Get the length of the slice.
    ///
    /// # Returns
    ///
    /// The length of the slice.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Checks if the slice and a given array of `u8` are equal.
    /// I and L can be equated.
    /// 
    /// # Arguments
    /// * `other` - the array of `u8` to compare the slice with.
    /// * `equate_il` - true if I and L need to be equated, false otherwise.
    ///
    /// # Returns
    ///
    /// True if the slice is equal to the given array, false otherwise.
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

    /// Check if the slice and a given array of `u8` are equal on the I and L positions.
    /// 
    /// # Arguments
    /// * `skip` - The amount of positions this slice skipped, this has an influence on the I and L positions.
    /// * `il_locations` - The positions where I and L occur.
    /// * `search_string` -  An array of `u8` to compare the slice with.
    ///
    /// # Returns
    ///
    /// True if the slice and `search_string` have the same contents on the I and L positions, false otherwise.
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

    /// Get an iterator over the slice.
    ///
    /// # Returns
    ///
    /// An iterator over the slice.
    pub fn iter(&self) -> ProteinTextSliceIterator {
        ProteinTextSliceIterator {text_slice: self, index: 0, }
    }
}

/// Structure representing an iterator over a `ProteinText` instance, iterating the characters of the text.
pub struct ProteinTextIterator<'a> {
    protein_text: &'a ProteinText,
    index: usize,
}

/// Structure representing an iterator over a `ProteintextSlice` instance, iterating the characters of the slice.
pub struct ProteinTextSliceIterator<'a> {
    text_slice: &'a ProteinTextSlice<'a>,
    index: usize,
}

impl<'a> Iterator for ProteinTextSliceIterator<'a> {

    type Item = u8;
    
    /// Get the next character in the `ProteinTextSlice`.
    /// 
    /// # Returns
    /// 
    /// The next character in the slice.
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
    
    /// Get the next character in the `ProteinText`.
    /// 
    /// # Returns
    /// 
    /// The next character in the text.
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
    fn test_u8_5bit_conversion() {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        for c in "ACDEFGHIKLMNPQRSTVWY-$".chars() {
            let char_5bit = char_to_5bit.get(&(c as u8)).unwrap();
            assert_eq!(c as u8, bit5_to_char[*char_5bit as usize]);
        }
    }

    #[test]
    fn test_build_from_string() {
        let text = ProteinText::from_string("ACACA-CAC$");

        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_from_vec() {
        let vec = vec![b'A', b'C', b'A', b'C', b'A', b'-', b'C', b'A', b'C', b'$'];
        let text = ProteinText::from_vec(&vec);

        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_from_bitarray() {
        let input_string = "ACACA-CAC$";
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();

        let mut bit_array = BitArray::with_capacity(input_string.len(), 5);
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 = *char_to_5bit.get(&(c as u8)).expect("Input character not in alphabet");
            bit_array.set(i, char_5bit as u64);
        }

        let text = ProteinText::new(bit_array);

        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_build_with_capacity() {
        let input_string = "ACACA-CAC$";

        let mut text = ProteinText::with_capacity(input_string.len());
        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            text.set(i, c as u8);
        }

        for (i, c) in "ACACA-CAC$".chars().enumerate() {
            assert_eq!(c as u8, text.get(i));
        }
    }

    #[test]
    fn test_text_slice() {
        let input_string = "ACACA-CAC$";
        let start = 1;
        let end  = 5;
        let text = ProteinText::from_string(&input_string);
        let text_slice = text.slice(start, end);

        for (i, c) in input_string[start..end].chars().enumerate() {
            assert_eq!(c as u8, text_slice.get(i));
        }
    }

    #[test]
    fn test_equals_slice() {
        let input_string = "ACICA-CAC$";
        let text = ProteinText::from_string(&input_string);
        let text_slice = text.slice(1, 5);
        let eq_slice_true = [b'C', b'I', b'C', b'A'];
        let eq_slice_false = [b'C', b'C', b'C', b'A'];
        let eq_slice_il_true = [b'C', b'L', b'C', b'A'];

        assert!(text_slice.equals_slice(&eq_slice_true, false));
        assert!(! text_slice.equals_slice(&eq_slice_false, false));
        assert!(text_slice.equals_slice(&eq_slice_il_true, true));
    }

    #[test]
    fn test_check_il_locations() {
        let input_string = "ACILA-CAC$";
        let text = ProteinText::from_string(&input_string);
        let text_slice = text.slice(1, 5);
        let il_locations = [1, 2];
        let il_true = [b'C', b'I', b'L', b'A'];
        let il_false = [b'C', b'I', b'C', b'A'];

        assert!(text_slice.check_il_locations(0, &il_locations, &il_true));
        assert!(! text_slice.check_il_locations(0, &il_locations, &il_false));
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
