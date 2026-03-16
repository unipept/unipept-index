use std::io::Write;
use std::error::Error;
use std::path::Path;
use std::fs::File;

use clap::ValueEnum;
use memmap2::Mmap;
use text_compression::ProteinText;
use crate::{ReadBinary, ReadBinaryMmap};

pub mod dense;
pub mod sparse;
pub mod bitvec;

pub use dense::DenseSuffixToProtein;
pub use sparse::SparseSuffixToProtein;
pub use bitvec::BitVecSuffixToProtein;

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
            dense::write_dense_mapping(&DenseSuffixToProtein::new(text), writer)
        }
        SuffixToProteinMappingStyle::Sparse => {
            writer.write_all(&[1u8])?;
            sparse::write_sparse_mapping(&SparseSuffixToProtein::new(text), writer)
        }
        SuffixToProteinMappingStyle::BitVec => {
            writer.write_all(&[2u8])?;
            bitvec::write_bitvec_mapping(&BitVecSuffixToProtein::new(text), writer)
        }
    }
}

/// A newtype wrapping a boxed `SuffixToProteinIndex` to enable trait-based loading.
pub struct SuffixToProteinMapping(pub Box<dyn SuffixToProteinIndex>);

impl ReadBinary for SuffixToProteinMapping {
    fn read_binary<R: std::io::BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut type_buf = [0u8; 1];
        reader.read_exact(&mut type_buf)?;
        let index: Box<dyn SuffixToProteinIndex> = match type_buf[0] {
            0 => Box::new(dense::read_dense_mapping(reader)?),
            1 => Box::new(sparse::read_sparse_mapping(reader)?),
            2 => Box::new(bitvec::read_bitvec_mapping(reader)?),
            t => return Err(format!("Unknown mapping type byte: {}", t).into()),
        };
        Ok(SuffixToProteinMapping(index))
    }
}

impl ReadBinaryMmap for SuffixToProteinMapping {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        let index: Box<dyn SuffixToProteinIndex> = match mmap[0] {
            0 => Box::new(dense::MmapDenseSuffixToProtein { mmap, data_offset: 9 }),
            1 => {
                let count = u64::from_le_bytes(mmap[1..9].try_into().unwrap()) as usize;
                Box::new(sparse::MmapSparseSuffixToProtein { mmap, data_offset: 9, count })
            }
            2 => {
                let bit_len = u64::from_le_bytes(mmap[1..9].try_into().unwrap());
                let block_count = u64::from_le_bytes(mmap[9..17].try_into().unwrap()) as usize;
                let bits_offset = 17;
                let counts_offset = bits_offset + block_count * 8;
                Box::new(bitvec::MmapBitVecSuffixToProtein { mmap, bit_len, bits_offset, counts_offset })
            }
            t => return Err(format!("Unknown mapping type byte: {}", t).into()),
        };
        Ok(SuffixToProteinMapping(index))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use clap::ValueEnum;
    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::ProteinText;

    use crate::{Nullable, ReadBinary, ReadBinaryMmap};
    use crate::suffix_to_protein_index::{SuffixToProteinMapping, SuffixToProteinMappingStyle, dump_mapping};

    fn build_text() -> ProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        ProteinText::from_string(&text)
    }

    fn write_to_tempfile(buf: &[u8]) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(buf).unwrap();
        tmp.flush().unwrap();
        tmp
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
    fn test_dump_and_load_mapping_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::Dense, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let mut cursor = Cursor::new(buf);
        let loaded = SuffixToProteinMapping::read_binary(&mut cursor).unwrap().0;
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
        let loaded = SuffixToProteinMapping::read_binary(&mut cursor).unwrap().0;
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
        let loaded = SuffixToProteinMapping::read_binary(&mut cursor).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_load_mapping_unknown_type() {
        let buf = vec![99u8];
        let mut cursor = Cursor::new(buf);
        let result = SuffixToProteinMapping::read_binary(&mut cursor);
        assert!(result.is_err());
    }

    #[test]
    fn test_mmap_unknown_type() {
        let buf = vec![99u8];
        let tmp = write_to_tempfile(&buf);
        let result = SuffixToProteinMapping::read_binary_mmap(tmp.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_dump_and_load_mmap_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::Dense, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mmap_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::Sparse, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(7), 2);
        assert_eq!(loaded.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mmap_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        dump_mapping(&SuffixToProteinMappingStyle::BitVec, &text, &mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }
}
