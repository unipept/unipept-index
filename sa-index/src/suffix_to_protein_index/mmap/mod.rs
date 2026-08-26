//! Runtime dispatch over the three mmap-backed mapping representations.
//!
//! Each answers out of the mapping itself rather than an owned structure, so a lookup is a page
//! the kernel may still have to fault in. The layouts they read are the ones the matching
//! `preloaded` types write, and are documented on both sides.

use std::{error::Error, path::Path};

use binary_traits::{LoadIndex, ReadBinaryMmap};
use memmap2::Mmap;

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
    fn implied_text_len(&self) -> Option<usize>;
    fn touch_all_pages(&self) -> u64;
});

impl ReadBinaryMmap for MmapBackedSuffixToProteinMapping {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = std::fs::File::open(path)?;
        // SAFETY: see the note in `protein_text::mmap` — an index file is written once by
        // sa-builder and is read-only for the lifetime of the process, so the mapping cannot be
        // truncated or written underneath us.
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.is_empty() {
            return Err("Mapping file is empty".into());
        }

        // As every other mmap loader in the workspace does. A lookup lands on a position the suffix
        // array chose, so the access order is one the kernel cannot predict, and the default
        // readahead drags in neighbouring pages that will not be used. This loader was the only one
        // missing the advice; the omission is invisible to any benchmark that warms first, because
        // `touch_all_pages` sets `Random` on the way out.
        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        match mmap[0] {
            0 => Ok(Self::Dense(dense::read_dense_mmap(mmap)?)),
            1 => Ok(Self::Sparse(sparse::read_sparse_mmap(mmap)?)),
            2 => Ok(Self::BitVec(bitvec::read_bitvec_mmap(mmap)?)),
            t => Err(format!("Unknown mapping type byte: {}", t).into())
        }
    }
}

impl LoadIndex for MmapBackedSuffixToProteinMapping {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        Self::read_binary_mmap(path)
    }
}

#[cfg(test)]
mod tests {
    use binary_traits::ReadBinaryMmap;

    use super::{
        MmapBackedSuffixToProteinMapping,
        test_utils::{assert_hints_are_harmless, assert_load_mmap, write_to_tempfile}
    };
    use crate::suffix_to_protein_index::{
        preloaded::{BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein},
        test_utils::sample_text
    };

    /// A truncated mapping file must be refused by *all three* representations, on both backends.
    ///
    /// Dense and sparse used to validate only their 9-byte header, so a file whose body was short
    /// loaded cleanly and panicked on the first lookup — inside a request handler, since this is a
    /// server startup path. Bitvec always checked. The contract in `binary_traits` requires all
    /// three to behave the same way, and now they do.
    #[test]
    fn every_mmap_representation_rejects_a_truncated_body() {
        use binary_traits::ReadBinary;

        use crate::suffix_to_protein_index::{InMemorySuffixToProteinMapping, test_utils::to_binary};

        let text = sample_text();
        let cases: [(&str, Vec<u8>); 3] = [
            ("dense", to_binary(DenseSuffixToProtein::new(&text))),
            ("sparse", to_binary(SparseSuffixToProtein::new(&text))),
            ("bitvec", to_binary(BitVecSuffixToProtein::new(&text)))
        ];

        for (name, full) in cases {
            // Header intact, body cut in half — unambiguously short for every layout.
            let cut = 9 + (full.len() - 9) / 2;
            assert!(cut < full.len(), "{name}: fixture too small to truncate meaningfully");

            let tmp = write_to_tempfile(&full[..cut]);
            assert!(
                MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).is_err(),
                "{name}: the mmap reader accepted a truncated body"
            );

            // The preloaded sibling has always rejected these; assert it stays that way, so the
            // two backends are pinned to the same answer rather than merely both erroring today.
            assert!(
                InMemorySuffixToProteinMapping::read_binary(&mut &full[..cut]).is_err(),
                "{name}: the preloaded reader accepted a truncated body"
            );

            // And the intact file still loads on both, so the assertions above reject for the
            // right reason.
            let tmp_full = write_to_tempfile(&full);
            assert!(
                MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp_full.path()).is_ok(),
                "{name}: the mmap reader rejected an intact file"
            );
            assert!(
                InMemorySuffixToProteinMapping::read_binary(&mut full.as_slice()).is_ok(),
                "{name}: the preloaded reader rejected an intact file"
            );
        }
    }

    /// A corrupt count must be a load error on both backends, not a process abort.
    ///
    /// All three preloaded readers used to size a `Vec` straight from the header, so an
    /// implausible count reached `handle_alloc_error` and `abort()`ed — un-catchable, and the
    /// opposite of the `Err` the mmap readers return for the same nine bytes.
    #[test]
    fn an_implausible_count_is_an_error_not_an_abort() {
        use binary_traits::ReadBinary;

        use crate::suffix_to_protein_index::InMemorySuffixToProteinMapping;

        for tag in [0u8, 1u8] {
            let mut bytes = vec![tag];
            bytes.extend_from_slice(&(1u64 << 60).to_le_bytes());

            assert!(
                InMemorySuffixToProteinMapping::read_binary(&mut bytes.as_slice()).is_err(),
                "tag {tag}: preloaded accepted a count of 2^60"
            );

            let tmp = write_to_tempfile(&bytes);
            assert!(
                MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).is_err(),
                "tag {tag}: mmap accepted a count of 2^60"
            );
        }
    }

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
        let cases: [(&str, Vec<u8>); 7] = [
            ("empty file", vec![]),
            ("unknown type byte", vec![99u8, 0, 0, 0, 0, 0, 0, 0, 0]),
            ("dense without a count header", vec![0u8, 0, 0]),
            ("sparse without a count header", vec![1u8, 0, 0]),
            ("bitvec without a full header", vec![2u8, 0, 0]),
            (
                "bitvec header promising blocks the file does not hold",
                [vec![2u8], 6400u64.to_le_bytes().to_vec(), 100u64.to_le_bytes().to_vec()].concat()
            ),
            (
                // A well-formed 41-byte file whose `bit_len` needs 157 blocks but whose
                // `block_count` is 1. Every length in it is honest, so the size check passes; only
                // relating the two fields catches it. Before that check this loaded cleanly and
                // `suffix_to_protein(1000)` panicked on `mmap[137..]` of a 41-byte mapping.
                "bitvec bit_len needing more blocks than the header declares",
                [vec![2u8], 10_000u64.to_le_bytes().to_vec(), 1u64.to_le_bytes().to_vec(), vec![0u8; 8], vec![0u8; 16]]
                    .concat()
            )
        ];

        for (case, buf) in cases {
            let tmp = write_to_tempfile(&buf);
            let result = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path());
            assert!(result.is_err(), "{} was accepted", case);
        }
    }
}
