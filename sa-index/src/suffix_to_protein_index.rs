use std::io::{Read, Write};
use std::error::Error;

use clap::ValueEnum;
use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinText;
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector, Rank9};
use crate::Nullable;

/// Enum used to define the commandline arguments and choose which index style is used
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SuffixToProteinMappingStyle {
    Dense,
    Sparse,
    BitVec
}

/// Trait implemented by the SuffixToProtein mappings
pub trait SuffixToProteinIndex: Send + Sync {
    /// Returns the index of the protein in the protein list for the given suffix
    ///
    /// # Arguments
    /// * `suffix` - The suffix of which we want to know of which protein it is a part
    ///
    /// # Returns
    ///
    /// Returns the index of the protein in the proteins list of which the suffix is a part
    fn suffix_to_protein(&self, suffix: i64) -> u32;
}

/// Mapping that uses O(n) memory with n the size of the input text, but retrieval of the protein is
/// in O(1)
#[derive(Debug, PartialEq)]
pub struct DenseSuffixToProtein {
    // UniProtKB does not have more that u32::MAX proteins, so a larger type is not needed
    mapping: Vec<u32>
}

/// Mapping that uses O(m) memory with m the number of proteins, but retrieval of the protein is
/// O(log m)
#[derive(Debug, PartialEq)]
pub struct SparseSuffixToProtein {
    mapping: Vec<i64>
}

/// Mapping that uses O(n) memory (1-2 bits per suffix) with n the size of the input text, with retrieval
/// of the protein in O(1)
#[derive(Debug)]
pub struct BitVecSuffixToProtein {
    rank: Rank9<BitVector<u64>>
}

impl SuffixToProteinIndex for DenseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        self.mapping[suffix as usize]
    }
}

impl SuffixToProteinIndex for SparseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let protein_index = self.mapping.binary_search(&suffix).unwrap_or_else(|index| index - 1);
        // if the next value in the mapping is 1 larger than the current suffix, that means that the
        // current suffix starts with a SEPARATION_CHARACTER or TERMINATION_CHARACTER
        // this means it does not belong to a protein
        if self.mapping[protein_index + 1] == suffix + 1 {
            return u32::NULL;
        }
        protein_index as u32
    }
}

impl SuffixToProteinIndex for BitVecSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let suffix: u64 = suffix.try_into().unwrap();
        if self.rank.get_bit(suffix) {
            return u32::NULL;
        }
        self.rank.rank1(suffix).try_into().unwrap()
    }
}

impl DenseSuffixToProtein {
    /// Creates a new DenseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new DenseSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
        let mut current_protein_index: u32 = 0;
        let mut suffix_index_to_protein: Vec<u32> = vec![];
        for char in text.iter() {
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                current_protein_index += 1;
                suffix_index_to_protein.push(u32::NULL);
            } else {
                assert_ne!(current_protein_index, u32::NULL);
                suffix_index_to_protein.push(current_protein_index);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        DenseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

impl SparseSuffixToProtein {
    /// Creates a new SparseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new SparseSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
        let mut suffix_index_to_protein: Vec<i64> = vec![0];
        for (index, char) in text.iter().enumerate() {
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                suffix_index_to_protein.push(index as i64 + 1);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        SparseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

impl BitVecSuffixToProtein {
    /// Creates a new BitVecSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - the text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new BitVecSuffixToProtein build over the provided text
    pub fn new(text: &ProteinText) -> Self {
        let num_bits = text.len();

        // Create a BitVec (dynamic) first
        let mut bits = BitVector::with_capacity(num_bits as u64);

        // Set bits
        for c in text.iter() {
            bits.push_bit(c == SEPARATION_CHARACTER || c == TERMINATION_CHARACTER);
        }

        let rank = Rank9::new(bits);

        BitVecSuffixToProtein { rank }
    }
}

fn write_dense_mapping<W: Write>(mapping: &DenseSuffixToProtein, writer: &mut W) -> Result<(), Box<dyn Error>> {
    let count = mapping.mapping.len() as u64;
    writer.write_all(&count.to_le_bytes())?;
    for &val in &mapping.mapping {
        writer.write_all(&val.to_le_bytes())?;
    }
    Ok(())
}

fn read_dense_mapping<R: Read>(reader: &mut R) -> Result<DenseSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let count = u64::from_le_bytes(buf8) as usize;
    let mut mapping = Vec::with_capacity(count);
    for _ in 0..count {
        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        mapping.push(u32::from_le_bytes(buf4));
    }
    Ok(DenseSuffixToProtein { mapping })
}

fn write_sparse_mapping<W: Write>(mapping: &SparseSuffixToProtein, writer: &mut W) -> Result<(), Box<dyn Error>> {
    let count = mapping.mapping.len() as u64;
    writer.write_all(&count.to_le_bytes())?;
    for &val in &mapping.mapping {
        writer.write_all(&val.to_le_bytes())?;
    }
    Ok(())
}

fn read_sparse_mapping<R: Read>(reader: &mut R) -> Result<SparseSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let count = u64::from_le_bytes(buf8) as usize;
    let mut mapping = Vec::with_capacity(count);
    for _ in 0..count {
        reader.read_exact(&mut buf8)?;
        mapping.push(i64::from_le_bytes(buf8));
    }
    Ok(SparseSuffixToProtein { mapping })
}

fn write_bitvec_mapping<W: Write>(mapping: &BitVecSuffixToProtein, writer: &mut W) -> Result<(), Box<dyn Error>> {
    let bit_len = mapping.rank.bit_len();
    let block_count = mapping.rank.block_len() as u64;
    writer.write_all(&bit_len.to_le_bytes())?;
    writer.write_all(&block_count.to_le_bytes())?;
    for i in 0..block_count as usize {
        let block: u64 = mapping.rank.get_block(i);
        writer.write_all(&block.to_le_bytes())?;
    }
    Ok(())
}

fn read_bitvec_mapping<R: Read>(reader: &mut R) -> Result<BitVecSuffixToProtein, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let bit_len = u64::from_le_bytes(buf8);
    reader.read_exact(&mut buf8)?;
    let block_count = u64::from_le_bytes(buf8);

    let mut bits = BitVector::with_capacity(bit_len);
    let mut bits_pushed: u64 = 0;

    for _ in 0..block_count {
        reader.read_exact(&mut buf8)?;
        let block = u64::from_le_bytes(buf8);
        let bits_in_block = std::cmp::min(64, bit_len - bits_pushed);
        for bit_pos in 0..bits_in_block {
            bits.push_bit((block >> bit_pos) & 1 == 1);
            bits_pushed += 1;
        }
    }

    Ok(BitVecSuffixToProtein { rank: Rank9::new(bits) })
}

/// Constructs the appropriate mapping from text and writes type byte + data to the writer.
///
/// # Arguments
/// * `style` - The style of mapping to construct
/// * `text` - The protein text over which to build the mapping
/// * `writer` - The writer to write the mapping to
///
/// # Returns
///
/// Returns Ok(()) on success, or an error if writing fails
pub fn dump_mapping<W: Write>(
    style: &SuffixToProteinMappingStyle,
    text: &ProteinText,
    writer: &mut W
) -> Result<(), Box<dyn Error>> {
    match style {
        SuffixToProteinMappingStyle::Dense => {
            writer.write_all(&[0u8])?;
            write_dense_mapping(&DenseSuffixToProtein::new(text), writer)
        }
        SuffixToProteinMappingStyle::Sparse => {
            writer.write_all(&[1u8])?;
            write_sparse_mapping(&SparseSuffixToProtein::new(text), writer)
        }
        SuffixToProteinMappingStyle::BitVec => {
            writer.write_all(&[2u8])?;
            write_bitvec_mapping(&BitVecSuffixToProtein::new(text), writer)
        }
    }
}

/// Reads a type byte from the reader and reconstructs the appropriate mapping.
///
/// # Arguments
/// * `reader` - The reader to read the mapping from
///
/// # Returns
///
/// Returns a boxed SuffixToProteinIndex on success, or an error if reading fails or the type byte
/// is unknown
pub fn load_mapping<R: Read>(reader: &mut R) -> Result<Box<dyn SuffixToProteinIndex>, Box<dyn Error>> {
    let mut type_buf = [0u8; 1];
    reader.read_exact(&mut type_buf)?;
    match type_buf[0] {
        0 => Ok(Box::new(read_dense_mapping(reader)?)),
        1 => Ok(Box::new(read_sparse_mapping(reader)?)),
        2 => Ok(Box::new(read_bitvec_mapping(reader)?)),
        t => Err(format!("Unknown mapping type byte: {}", t).into()),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use clap::ValueEnum;
    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::ProteinText;

    use crate::{
        Nullable,
        suffix_to_protein_index::{
            DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein, SuffixToProteinIndex,
            SuffixToProteinMappingStyle, dump_mapping, load_mapping,
            write_dense_mapping, read_dense_mapping,
            write_sparse_mapping, read_sparse_mapping,
            write_bitvec_mapping, read_bitvec_mapping,
        }
    };

    fn build_text() -> ProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        ProteinText::from_string(&text)
    }

    #[test]
    fn test_suffix_to_protein_mapping_style() {
        assert_eq!(SuffixToProteinMappingStyle::Dense, SuffixToProteinMappingStyle::from_str("dense", false).unwrap());
        assert_eq!(
            SuffixToProteinMappingStyle::Sparse,
            SuffixToProteinMappingStyle::from_str("sparse", false).unwrap()
        );
        assert_eq!(
            SuffixToProteinMappingStyle::BitVec,
            SuffixToProteinMappingStyle::from_str("bit-vec", false).unwrap()
        );
    }

    #[test]
    fn test_dense_build() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        let expected = DenseSuffixToProtein {
            mapping: vec![0, 0, 0, u32::NULL, 1, 1, u32::NULL, 2, 2, 2, u32::NULL]
        };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_sparse_build() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        let expected = SparseSuffixToProtein { mapping: vec![0, 4, 7, 11] };
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_dense() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_search_sparse() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_search_bitvec() {
        let u8_text = &build_text();
        let index = BitVecSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_dense_roundtrip() {
        let text = build_text();
        let original = DenseSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_dense_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_dense_mapping(&mut cursor).unwrap();
        assert_eq!(original, restored);
    }

    #[test]
    fn test_sparse_roundtrip() {
        let text = build_text();
        let original = SparseSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_sparse_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_sparse_mapping(&mut cursor).unwrap();
        assert_eq!(original, restored);
    }

    #[test]
    fn test_bitvec_roundtrip() {
        let text = build_text();
        let original = BitVecSuffixToProtein::new(&text);
        let mut buf = Vec::new();
        write_bitvec_mapping(&original, &mut buf).unwrap();
        let mut cursor = Cursor::new(buf);
        let restored = read_bitvec_mapping(&mut cursor).unwrap();
        // BitVecSuffixToProtein does not derive PartialEq, compare via suffix_to_protein
        for i in 0..text.len() as i64 {
            assert_eq!(original.suffix_to_protein(i), restored.suffix_to_protein(i));
        }
    }

    #[test]
    fn test_dump_and_load_mapping_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::Dense, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let mut cursor = Cursor::new(buf);
        let loaded = load_mapping(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mapping_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::Sparse, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let mut cursor = Cursor::new(buf);
        let loaded = load_mapping(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(7), 2);
        assert_eq!(loaded.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mapping_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::BitVec, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let mut cursor = Cursor::new(buf);
        let loaded = load_mapping(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_load_mapping_unknown_type() {
        let buf = vec![99u8];
        let mut cursor = Cursor::new(buf);
        let result = load_mapping(&mut cursor);
        assert!(result.is_err());
    }
}
