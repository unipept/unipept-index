//! Runtime dispatch over the three mmap-backed mapping representations.
//!
//! Each answers out of the mapping itself rather than an owned structure, so a lookup is a page
//! the kernel may still have to fault in. The layouts they read are the ones the matching
//! `preloaded` types write, and are documented on both sides.
//!
//! This entire module is mmap-only.

use std::{error::Error, path::Path};

use memmap2::MmapOptions;

use crate::ReadBinaryMmap;

pub mod bitvec;
pub mod dense;
pub mod sparse;
#[cfg(test)]
pub(super) mod test_utils;

pub use bitvec::MmapBitVecSuffixToProtein;
pub use dense::MmapDenseSuffixToProtein;
pub use sparse::MmapSparseSuffixToProtein;

/// Wraps whichever of the three mappings a file holds, picked at load time from its type byte.
pub enum MmapBackedSuffixToProteinMapping {
    Dense(MmapDenseSuffixToProtein),
    Sparse(MmapSparseSuffixToProtein),
    BitVec(MmapBitVecSuffixToProtein)
}

delegate_suffix_to_protein_mapping!(MmapBackedSuffixToProteinMapping {
    fn suffix_to_protein(&self, suffix: i64) -> u32;
    fn prefetch_for_suffix(&self, suffix: i64);
    fn touch_all_pages(&self);
});

impl ReadBinaryMmap for MmapBackedSuffixToProteinMapping {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = std::fs::File::open(path)?;
        // SAFETY: see the note in `text_compression::mmap` — an index file is written once by
        // sa-builder and is read-only for the lifetime of the process, so the mapping cannot be
        // truncated or written underneath us.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        if mmap.is_empty() {
            return Err("Mapping file is empty".into());
        }
        match mmap[0] {
            0 => Ok(Self::Dense(dense::read_dense_mmap(mmap)?)),
            1 => Ok(Self::Sparse(sparse::read_sparse_mmap(mmap)?)),
            2 => Ok(Self::BitVec(bitvec::read_bitvec_mmap(mmap)?)),
            t => Err(format!("Unknown mapping type byte: {}", t).into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        MmapBackedSuffixToProteinMapping,
        test_utils::{assert_hints_are_harmless, assert_load_mmap, write_to_tempfile}
    };
    use crate::{
        ReadBinaryMmap,
        suffix_to_protein_index::{
            preloaded::{BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein},
            test_utils::sample_text
        }
    };

    #[test]
    fn test_load_mmap_dense() {
        assert_load_mmap(DenseSuffixToProtein::new(&sample_text()), 0u8);
        assert_hints_are_harmless(DenseSuffixToProtein::new(&sample_text()));
    }

    #[test]
    fn test_load_mmap_sparse() {
        assert_load_mmap(SparseSuffixToProtein::new(&sample_text()), 1u8);
        assert_hints_are_harmless(SparseSuffixToProtein::new(&sample_text()));
    }

    #[test]
    fn test_load_mmap_bitvec() {
        assert_load_mmap(BitVecSuffixToProtein::new(&sample_text()), 2u8);
        assert_hints_are_harmless(BitVecSuffixToProtein::new(&sample_text()));
    }

    /// A file that is not an index — or one whose header promises more than it holds — must come
    /// back as an error. These readers index the mapping by offsets taken from the header, so
    /// anything they fail to reject here becomes a panic on the first lookup instead.
    #[test]
    fn test_load_mmap_rejects_malformed_files() {
        let cases: [(&str, Vec<u8>); 6] = [
            ("empty file", vec![]),
            ("unknown type byte", vec![99u8, 0, 0, 0, 0, 0, 0, 0, 0]),
            ("dense without a count header", vec![0u8, 0, 0]),
            ("sparse without a count header", vec![1u8, 0, 0]),
            ("bitvec without a full header", vec![2u8, 0, 0]),
            (
                "bitvec header promising blocks the file does not hold",
                [vec![2u8], 6400u64.to_le_bytes().to_vec(), 100u64.to_le_bytes().to_vec()].concat()
            )
        ];

        for (case, buf) in cases {
            let tmp = write_to_tempfile(&buf);
            let result = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path());
            assert!(result.is_err(), "{} was accepted", case);
        }
    }
}
