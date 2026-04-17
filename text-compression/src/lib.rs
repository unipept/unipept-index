use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, Write},
    path::Path,
    sync::Arc
};
use std::io::Read;
use bitarray::{Binary, BitArray, data_to_writer};
use memmap2::Mmap;

pub mod traits;
pub use traits::{WriteBinary, ReadBinary, ReadBinaryMmap};

/// The 5-bit-to-char lookup table for mmap-backed ProteinText.
const BIT5_TO_CHAR: &[u8; 27] = b"ABCDEFGHIKLMNOPQRSTUVWXYZ-$";

/// Returns the number of bytes the BitArray data occupies for a given text length at 5 bits/value.
pub fn bit_array_byte_size(text_length: usize) -> usize {
    let extra = if (text_length * 5).is_multiple_of(64) { 0 } else { 1 };
    (text_length * 5 / 64 + extra) * 8
}

/// Structure representing the proteins, stored in a bit array using 5 bits per amino acid.
pub enum ProteinText {
    /// In-memory representation using a BitArray.
    InMemory {
        bit_array: BitArray,
        char_to_5bit: HashMap<u8, u8>,
        bit5_to_char: Vec<u8>,
    },
    /// Memory-mapped representation backed by a file.
    MmapBacked {
        mmap: Arc<Mmap>,
        data_offset: usize,
        len: usize,
    },
}

impl ProteinText {
    /// Creates the hashmap storing the mappings between the characters as `u8` and 5 bit numbers.
    fn create_char_to_5bit_hashmap() -> HashMap<u8, u8> {
        let mut hashmap = HashMap::<u8, u8>::new();
        for (i, c) in "ABCDEFGHIKLMNOPQRSTUVWXYZ-$".chars().enumerate() {
            hashmap.insert(c as u8, i as u8);
        }
        hashmap
    }

    /// Creates the vector storing the mappings between the 5 bit numbers and the characters as `u8`.
    fn create_bit5_to_char() -> Vec<u8> {
        let mut vec = Vec::<u8>::new();
        for c in "ABCDEFGHIKLMNOPQRSTUVWXYZ-$".chars() {
            vec.push(c as u8);
        }
        vec
    }

    /// Creates the compressed text from a string.
    pub fn from_string(input_string: &str) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        let mut bit_array = BitArray::with_capacity(input_string.len(), 5);
        for (i, c) in input_string.chars().enumerate() {
            let char_5bit: u8 =
                *char_to_5bit.get(&(c as u8)).unwrap_or_else(|| panic!("Input character '{}' not in alphabet", c));
            bit_array.set(i, char_5bit as u64);
        }

        ProteinText::InMemory { bit_array, char_to_5bit, bit5_to_char }
    }

    /// Creates the compressed text from a vector.
    pub fn from_vec(input_vec: &[u8]) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        let mut bit_array = BitArray::with_capacity(input_vec.len(), 5);
        for (i, e) in input_vec.iter().enumerate() {
            let char_5bit: u8 =
                *char_to_5bit.get(e).unwrap_or_else(|| panic!("Input character '{}' not in alphabet", e));
            bit_array.set(i, char_5bit as u64);
        }

        ProteinText::InMemory { bit_array, char_to_5bit, bit5_to_char }
    }

    /// Creates the compressed text from a bit array.
    pub fn new(bit_array: BitArray) -> ProteinText {
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();
        ProteinText::InMemory { bit_array, char_to_5bit, bit5_to_char }
    }

    /// Creates an instance of `ProteinText` with a given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self::new(BitArray::with_capacity(capacity, 5))
    }

    /// Creates a memory-mapped ProteinText backed by an existing mmap.
    pub fn from_mmap(mmap: Arc<Mmap>, data_offset: usize, len: usize) -> Self {
        ProteinText::MmapBacked { mmap, data_offset, len }
    }

    /// Extract a 5-bit value from the mmap and convert to a character.
    fn get_mmap(mmap: &Mmap, data_offset: usize, index: usize) -> u8 {
        const BITS: usize = 5;
        const MASK: u64 = (1u64 << BITS) - 1;

        let bit_offset = index * BITS;
        let start_block = bit_offset / 64;
        let start_bit = bit_offset % 64;
        let byte_off = data_offset + start_block * 8;

        let lo = u64::from_le_bytes(mmap[byte_off..byte_off + 8].try_into().unwrap());

        let raw = if start_bit + BITS <= 64 {
            (lo >> (64 - start_bit - BITS)) & MASK
        } else {
            let end_bit = (index + 1) * BITS % 64;
            let hi = u64::from_le_bytes(mmap[byte_off + 8..byte_off + 16].try_into().unwrap());
            ((lo << end_bit) | (hi >> (64 - end_bit))) & MASK
        };

        BIT5_TO_CHAR[raw as usize]
    }

    /// Search the character at a given position in the compressed text.
    pub fn get(&self, index: usize) -> u8 {
        match self {
            ProteinText::InMemory { bit_array, bit5_to_char, .. } => {
                let char_5bit = bit_array.get(index) as usize;
                bit5_to_char[char_5bit]
            }
            ProteinText::MmapBacked { mmap, data_offset, .. } => {
                Self::get_mmap(mmap, *data_offset, index)
            }
        }
    }

    /// Set the character at a given index. Only valid for InMemory variant.
    pub fn set(&mut self, index: usize, value: u8) {
        match self {
            ProteinText::InMemory { bit_array, char_to_5bit, .. } => {
                let char_5bit: u8 = *char_to_5bit
                    .get(&value)
                    .unwrap_or_else(|| panic!("Input character '{}' not in alphabet", value));
                bit_array.set(index, char_5bit as u64);
            }
            ProteinText::MmapBacked { .. } => {
                panic!("set() is not supported on MmapBacked ProteinText");
            }
        }
    }

    /// Queries the length of the text.
    pub fn len(&self) -> usize {
        match self {
            ProteinText::InMemory { bit_array, .. } => bit_array.len(),
            ProteinText::MmapBacked { len, .. } => *len,
        }
    }

    /// Check if the text is empty (length 0).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Clears the `BitArray`, setting all bits to 0. Only valid for InMemory variant.
    pub fn clear(&mut self) {
        match self {
            ProteinText::InMemory { bit_array, .. } => bit_array.clear(),
            ProteinText::MmapBacked { .. } => {
                panic!("clear() is not supported on MmapBacked ProteinText");
            }
        }
    }

    /// Get an iterator over the characters of the text.
    pub fn iter(&self) -> ProteinTextIterator<'_> {
        ProteinTextIterator { protein_text: self, index: 0 }
    }

    /// Get a slice of the text.
    pub fn slice(&self, start: usize, end: usize) -> ProteinTextSlice<'_> {
        ProteinTextSlice::new(self, start, end)
    }

}

impl WriteBinary for ProteinText {
    /// Writes this `ProteinText` to a writer in the binary proteins format.
    ///
    /// Format:
    /// - 8 bytes: text_length (u64 le)
    /// - N bytes: BitArray data where N = ceil(text_length * 5 / 64) * 8
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        match self {
            ProteinText::InMemory { bit_array, .. } => {
                let text_length = bit_array.len() as u64;
                writer.write_all(&text_length.to_le_bytes())?;
                bit_array.write_binary(writer)?;
            }
            ProteinText::MmapBacked { mmap, data_offset, len } => {
                let text_length = len as u64;
                writer.write_all(&text_length.to_le_bytes())?;
                let n_bytes = bit_array_byte_size(len);
                writer.write_all(&mmap[data_offset..data_offset + n_bytes])?;
            }
        }
        Ok(())
    }
}

impl ReadBinary for ProteinText {
    /// Reads a `ProteinText` from a reader in the binary proteins format.
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8).map_err(|_| "Could not read text_length from binary file")?;
        let text_length = u64::from_le_bytes(buf8) as usize;

        let n_bytes = bit_array_byte_size(text_length);
        let mut bit_array = BitArray::with_capacity(text_length, 5);

        let mut limited = reader.take(n_bytes as u64);
        bit_array
            .read_binary(&mut limited)
            .map_err(|_| "Could not parse BitArray data from binary file")?;

        Ok(ProteinText::new(bit_array))
    }
}

impl ReadBinaryMmap for ProteinText {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let f = File::open(path)?;
        let mmap = Arc::new(unsafe { Mmap::map(&f)? });

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        // Ensure the file is large enough to contain the 8-byte header.
        if mmap.len() < 8 {
            return Err("File is too small to contain ProteinText header (8 bytes required)".into());
        }

        let text_length = u64::from_le_bytes(mmap[0..8].try_into()
            .map_err(|_| "Failed to parse ProteinText header")?) as usize;

        // Ensure the file is large enough to contain the BitArray data for the declared text length.
        if mmap.len() < 8 + bit_array_byte_size(text_length) {
            return Err("File is too small to contain ProteinText BitArray data for declared length".into());
        }

        Ok(ProteinText::from_mmap(mmap, 8, text_length))
    }
}

/// Structure representing a slice of a `ProteinText`.
pub struct ProteinTextSlice<'a> {
    /// The `Proteintext` of which to take a slice.
    text: &'a ProteinText,
    /// The start of the slice.
    start: usize, // included
    /// The end of the slice.
    end: usize // excluded
}

impl<'a> ProteinTextSlice<'a> {
    /// Creates an instance of `ProteintextSlice`, given the text and boundaries.
    pub fn new(text: &'a ProteinText, start: usize, end: usize) -> ProteinTextSlice<'a> {
        Self { text, start, end }
    }

    /// Get a character (amino acid) in the slice.
    pub fn get(&self, index: usize) -> u8 {
        self.text.get(self.start + index)
    }

    /// Get the length of the slice.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Checks if the slice and a given array of `u8` are equal.
    /// I and L can be equated.
    #[inline]
    pub fn equals_slice(&self, other: &[u8], equate_il: bool) -> bool {
        if equate_il {
            other.iter().zip(self.iter()).all(|(&search_character, text_character)| {
                search_character == text_character
                    || (search_character == b'I' && text_character == b'L')
                    || (search_character == b'L' && text_character == b'I')
            })
        } else {
            other
                .iter()
                .zip(self.iter())
                .all(|(&search_character, text_character)| search_character == text_character)
        }
    }

    /// Check if the slice and a given array of `u8` are equal on the I and L positions.
    pub fn check_il_locations(&self, skip: usize, il_locations: &[usize], search_string: &[u8]) -> bool {
        for &il_location in il_locations {
            let index = il_location - skip;
            if search_string[index] != self.get(index) {
                return false;
            }
        }
        true
    }

    /// Get an iterator over the slice.
    pub fn iter(&self) -> ProteinTextSliceIterator<'_> {
        ProteinTextSliceIterator { text_slice: self, index: 0 }
    }
}

/// Structure representing an iterator over a `ProteinText` instance.
pub struct ProteinTextIterator<'a> {
    protein_text: &'a ProteinText,
    index: usize
}

/// Structure representing an iterator over a `ProteintextSlice` instance.
pub struct ProteinTextSliceIterator<'a> {
    text_slice: &'a ProteinTextSlice<'a>,
    index: usize
}

impl Iterator for ProteinTextSliceIterator<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.text_slice.len() {
            return None;
        }

        self.index += 1;
        Some(self.text_slice.get(self.index - 1))
    }
}

impl Iterator for ProteinTextIterator<'_> {
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
pub fn dump_compressed_text(text: Vec<u8>, writer: &mut impl Write) -> Result<(), Box<dyn Error>> {
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
pub fn load_compressed_text(reader: &mut impl BufRead) -> Result<ProteinText, Box<dyn Error>> {
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
        /// The number of times the read function can be called before it fails.
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
        let char_to_5bit = ProteinText::create_char_to_5bit_hashmap();
        let bit5_to_char = ProteinText::create_bit5_to_char();

        for c in "ABCDEFGHIKLMNOPQRSTUVWXYZ-$".chars() {
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
            let char_5bit: u8 =
                *char_to_5bit.get(&(c as u8)).unwrap_or_else(|| panic!("Input character '{}' not in alphabet", c));
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
        let end = 5;
        let text = ProteinText::from_string(input_string);
        let text_slice = text.slice(start, end);

        for (i, c) in input_string[start..end].chars().enumerate() {
            assert_eq!(c as u8, text_slice.get(i));
        }
    }

    #[test]
    fn test_equals_slice() {
        let input_string = "ACICA-CAC$";
        let text = ProteinText::from_string(input_string);
        let text_slice = text.slice(1, 5);
        let eq_slice_true = [b'C', b'I', b'C', b'A'];
        let eq_slice_false = [b'C', b'C', b'C', b'A'];
        let eq_slice_il_true = [b'C', b'L', b'C', b'A'];

        assert!(text_slice.equals_slice(&eq_slice_true, false));
        assert!(!text_slice.equals_slice(&eq_slice_false, false));
        assert!(text_slice.equals_slice(&eq_slice_il_true, true));
    }

    #[test]
    fn test_check_il_locations() {
        let input_string = "ACILA-CAC$";
        let text = ProteinText::from_string(input_string);
        let text_slice = text.slice(1, 5);
        let il_locations = [1, 2];
        let il_true = [b'C', b'I', b'L', b'A'];
        let il_false = [b'C', b'I', b'C', b'A'];

        assert!(text_slice.check_il_locations(0, &il_locations, &il_true));
        assert!(!text_slice.check_il_locations(0, &il_locations, &il_false));
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
        let data = [10, 0, 0, 0, 0, 0, 0, 0, // compressed text
            0, 128, 74, 232, 152, 66, 134, 8];

        let mut reader = std::io::BufReader::new(&data[..]);
        let compressed_text = load_compressed_text(&mut reader).unwrap();

        for (i, c) in "BCDEFGHIKL".chars().enumerate() {
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
        let mut buffer = [0_u8; 1];
        assert!(reader.read(&mut buffer).is_err());
    }

    #[test]
    fn test_write_and_read_binary() {
        let input = "ACACA-CAC$";
        let text = ProteinText::from_string(input);

        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();

        let mut reader = std::io::BufReader::new(buf.as_slice());
        let loaded = ProteinText::read_binary(&mut reader).unwrap();

        for (i, c) in input.chars().enumerate() {
            assert_eq!(loaded.get(i), c as u8);
        }
        assert_eq!(loaded.len(), input.len());
    }

    #[test]
    fn test_mmap_roundtrip() {
        use std::fs::File;
        use memmap2::Mmap;

        let input = "ACACA-CAC$MLPGLALLLL$";
        let text = ProteinText::from_string(input);

        // Write to a temp file
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        text.write_binary(tmp.as_file_mut()).unwrap();
        tmp.as_file_mut().flush().unwrap();

        // Mmap the file
        let f = File::open(tmp.path()).unwrap();
        let mmap = Arc::new(unsafe { Mmap::map(&f).unwrap() });

        let text_length = u64::from_le_bytes(mmap[0..8].try_into().unwrap()) as usize;
        let mmap_text = ProteinText::from_mmap(Arc::clone(&mmap), 8, text_length);

        assert_eq!(mmap_text.len(), input.len());
        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8, "mismatch at index {}", i);
        }
    }

    #[test]
    fn test_mmap_block_boundary() {
        // 13 characters: 13*5=65 bits, crosses a u64 boundary at index 12 (bit 60..65)
        let input = "ABCDEFGHIKLMN";
        let text = ProteinText::from_string(input);

        let mut buf: Vec<u8> = Vec::new();
        text.write_binary(&mut buf).unwrap();

        let mmap = Arc::new(unsafe {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            tmp.write_all(&buf).unwrap();
            tmp.flush().unwrap();
            let f = std::fs::File::open(tmp.path()).unwrap();
            Mmap::map(&f).unwrap()
        });

        let text_length = u64::from_le_bytes(mmap[0..8].try_into().unwrap()) as usize;
        let mmap_text = ProteinText::from_mmap(Arc::clone(&mmap), 8, text_length);

        for (i, c) in input.chars().enumerate() {
            assert_eq!(mmap_text.get(i), c as u8, "boundary mismatch at index {}", i);
        }
    }
}
