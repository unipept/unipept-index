use std::error::Error;

use crate::ReadBinary;
use super::SuffixToProteinMappingBackend;

pub mod dense;
pub mod sparse;
pub mod bitvec;

pub use dense::DenseSuffixToProtein;
pub use sparse::SparseSuffixToProtein;
pub use bitvec::BitVecSuffixToProtein;

pub enum InMemorySuffixToProteinMapping {
    Dense(DenseSuffixToProtein),
    Sparse(SparseSuffixToProtein),
    BitVec(BitVecSuffixToProtein),
}

impl SuffixToProteinMappingBackend for InMemorySuffixToProteinMapping {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        match self {
            Self::Dense(m) => m.suffix_to_protein(suffix),
            Self::Sparse(m) => m.suffix_to_protein(suffix),
            Self::BitVec(m) => m.suffix_to_protein(suffix),
        }
    }

    fn prefetch_for_suffix(&self, suffix: i64) {
        match self {
            Self::Dense(m) => m.prefetch_for_suffix(suffix),
            Self::Sparse(m) => m.prefetch_for_suffix(suffix),
            Self::BitVec(m) => m.prefetch_for_suffix(suffix),
        }
    }
}

impl ReadBinary for InMemorySuffixToProteinMapping {
    fn read_binary<R: std::io::BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut type_buf = [0u8; 1];
        reader.read_exact(&mut type_buf)?;
        match type_buf[0] {
            0 => Ok(InMemorySuffixToProteinMapping::Dense(dense::read_dense_mapping(reader)?)),
            1 => Ok(InMemorySuffixToProteinMapping::Sparse(sparse::read_sparse_mapping(reader)?)),
            2 => Ok(InMemorySuffixToProteinMapping::BitVec(bitvec::read_bitvec_mapping(reader)?)),
            t => Err(format!("Unknown mapping type byte: {}", t).into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::{InMemoryProteinText, ProteinTextBackend};

    use crate::{Nullable, ReadBinary, WriteBinary};
    use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
    use crate::suffix_to_protein_index::preloaded::{
        InMemorySuffixToProteinMapping, DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein,
    };

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    #[test]
    fn test_dump_and_load_mapping_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        DenseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let mut cursor = std::io::Cursor::new(buf);
        let loaded = InMemorySuffixToProteinMapping::read_binary(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mapping_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let mut cursor = std::io::Cursor::new(buf);
        let loaded = InMemorySuffixToProteinMapping::read_binary(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(7), 2);
        assert_eq!(loaded.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_dump_and_load_mapping_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let mut cursor = std::io::Cursor::new(buf);
        let loaded = InMemorySuffixToProteinMapping::read_binary(&mut cursor).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_load_mapping_unknown_type() {
        let buf = vec![99u8];
        let mut cursor = std::io::Cursor::new(buf);
        let result = InMemorySuffixToProteinMapping::read_binary(&mut cursor);
        assert!(result.is_err());
    }
}
