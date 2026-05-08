use std::error::Error;

use clap::ValueEnum;
use crate::ReadBinary;
#[cfg(feature = "mmap")]
use std::{path::Path, fs::File};
#[cfg(feature = "mmap")]
use memmap2::Mmap;
#[cfg(feature = "mmap")]
use crate::ReadBinaryMmap;

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

    /// Non-blocking hardware prefetch hint for the mmap data that
    /// `suffix_to_protein(suffix)` will access.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn prefetch_for_suffix(&self, _suffix: i64) {}

    /// Reads at least one byte from every OS page in the mmap backing this mapping,
    /// ensuring all pages are resident in the page cache.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn touch_all_pages(&self) {}
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

#[cfg(feature = "mmap")]
impl ReadBinaryMmap for SuffixToProteinMapping {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        if mmap.is_empty() {
            return Err("Mapping file is too small to contain a type header byte".into());
        }

        let index: Box<dyn SuffixToProteinIndex> = match mmap[0] {
            0 => Box::new(dense::read_dense_mmap(mmap)?),
            1 => Box::new(sparse::read_sparse_mmap(mmap)?),
            2 => Box::new(bitvec::read_bitvec_mmap(mmap)?),
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
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::{Nullable, ReadBinary, WriteBinary};
    #[cfg(feature = "mmap")]
    use crate::ReadBinaryMmap;
    use crate::suffix_to_protein_index::{SuffixToProteinMapping, SuffixToProteinMappingStyle, DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein};

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    #[cfg(feature = "mmap")]
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
        DenseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
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
        SparseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
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
        BitVecSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
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

    #[cfg(feature = "mmap")]
    #[test]
    fn test_mmap_unknown_type() {
        let buf = vec![99u8];
        let tmp = write_to_tempfile(&buf);
        let result = SuffixToProteinMapping::read_binary_mmap(tmp.path());
        assert!(result.is_err());
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn test_dump_and_load_mmap_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        DenseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn test_dump_and_load_mmap_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(7), 2);
        assert_eq!(loaded.suffix_to_protein(10), u32::NULL);
    }

    #[cfg(feature = "mmap")]
    #[test]
    fn test_dump_and_load_mmap_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = SuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap().0;
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }
}
