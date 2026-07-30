use std::error::Error;
use std::path::Path;

use memmap2::MmapOptions;

use crate::ReadBinaryMmap;
use super::SuffixToProteinMappingBackend;

pub mod dense;
pub mod sparse;
pub mod bitvec;

pub use dense::MmapDenseSuffixToProtein;
pub use sparse::MmapSparseSuffixToProtein;
pub use bitvec::MmapBitVecSuffixToProtein;

pub enum MmapBackedSuffixToProteinMapping {
    Dense(MmapDenseSuffixToProtein),
    Sparse(MmapSparseSuffixToProtein),
    BitVec(MmapBitVecSuffixToProtein),
}

impl SuffixToProteinMappingBackend for MmapBackedSuffixToProteinMapping {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        match self {
            Self::Dense(m) => m.suffix_to_protein(suffix),
            Self::Sparse(m) => m.suffix_to_protein(suffix),
            Self::BitVec(m) => m.suffix_to_protein(suffix),
        }
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        match self {
            Self::Dense(m) => m.prefetch_for_suffix(suffix),
            Self::Sparse(m) => m.prefetch_for_suffix(suffix),
            Self::BitVec(m) => m.prefetch_for_suffix(suffix),
        }
    }

    fn touch_all_pages(&self) {
        match self {
            Self::Dense(m) => m.touch_all_pages(),
            Self::Sparse(m) => m.touch_all_pages(),
            Self::BitVec(m) => m.touch_all_pages(),
        }
    }
}

impl ReadBinaryMmap for MmapBackedSuffixToProteinMapping {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        // Experiment: env-gated transparent huge pages (see sa-index/array/mmap.rs).
        #[cfg(target_os = "linux")]
        if std::env::var_os("SA_MADV_HUGEPAGE").is_some() {
            let _ = mmap.advise(memmap2::Advice::HugePage);
        }

        if mmap.is_empty() {
            return Err("Mapping file is empty".into());
        }
        match mmap[0] {
            0 => Ok(Self::Dense(dense::read_dense_mmap(mmap)?)),
            1 => Ok(Self::Sparse(sparse::read_sparse_mmap(mmap)?)),
            2 => Ok(Self::BitVec(bitvec::read_bitvec_mmap(mmap)?)),
            t => Err(format!("Unknown mapping type byte: {}", t).into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write as IoWrite;

    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use text_compression::InMemoryProteinText;

    use crate::{Nullable, ReadBinaryMmap, WriteBinary};
    use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
    use crate::suffix_to_protein_index::preloaded::{DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein};
    use super::MmapBackedSuffixToProteinMapping;

    fn build_text() -> InMemoryProteinText {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        InMemoryProteinText::from_string(&text)
    }

    fn write_to_tempfile(buf: &[u8]) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(buf).unwrap();
        tmp.flush().unwrap();
        tmp
    }

    #[test]
    fn test_load_mmap_dense() {
        let text = build_text();
        let mut buf = Vec::new();
        DenseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 0u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_load_mmap_sparse() {
        let text = build_text();
        let mut buf = Vec::new();
        SparseSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 1u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();
        assert_eq!(loaded.suffix_to_protein(7), 2);
        assert_eq!(loaded.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_load_mmap_bitvec() {
        let text = build_text();
        let mut buf = Vec::new();
        BitVecSuffixToProtein::new(&text).write_binary(&mut buf).unwrap();
        assert_eq!(buf[0], 2u8);
        let tmp = write_to_tempfile(&buf);
        let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();
        assert_eq!(loaded.suffix_to_protein(5), 1);
        assert_eq!(loaded.suffix_to_protein(3), u32::NULL);
    }

    #[test]
    fn test_load_mmap_unknown_type() {
        let buf = vec![99u8, 0, 0, 0, 0, 0, 0, 0, 0];
        let tmp = write_to_tempfile(&buf);
        let result = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path());
        assert!(result.is_err());
    }
}
