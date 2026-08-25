//! Runtime dispatch over the three owned-memory mapping representations.
//!
//! These types own the `WriteBinary` implementations that produce the files *either* backend
//! reads, which is why `sa-builder` names only them. The tag-sniffing `read_binary` below is the
//! single place that decides which of the three a file holds.

use std::{error::Error, path::Path};

use binary_traits::{LoadIndex, ReadBinary};

pub mod bitvec;
pub mod dense;
pub mod sparse;
#[cfg(test)]
pub(super) mod test_utils;

pub use bitvec::BitVecSuffixToProtein;
pub use dense::DenseSuffixToProtein;
pub use sparse::SparseSuffixToProtein;

/// Allocates a `Vec` with exactly `capacity`, reporting failure instead of aborting the process.
///
/// Every reader below takes its capacity from a count in the file header. `Vec::with_capacity`
/// reaches `handle_alloc_error` on an implausible one, which `abort()`s — un-catchable, and a very
/// different outcome from the `Err` the mmap readers return for the very same bytes. A crafted
/// 9-byte header declaring `2^60` entries used to kill the process here and load-error there.
///
/// The capacity stays exact on purpose: `bitarray::hugepages::advise_capacity` is issued against
/// these allocations, and a rounded-up or capped capacity would detach the advice from the
/// allocation it describes.
pub(super) fn try_alloc_exact<T>(capacity: usize, what: &str) -> Result<Vec<T>, Box<dyn Error>> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(capacity)
        .map_err(|_| format!("The {what} mapping header declares {capacity} entries, which cannot be allocated"))?;
    Ok(buffer)
}

/// Wraps whichever of the three mappings a file holds, picked at load time from its type byte.
pub enum InMemorySuffixToProteinMapping {
    Dense(DenseSuffixToProtein),
    Sparse(SparseSuffixToProtein),
    BitVec(BitVecSuffixToProtein)
}

// `touch_all_pages` is left to the trait default: there are no pages to fault in here.
delegate_suffix_to_protein_mapping!(InMemorySuffixToProteinMapping {
    fn suffix_to_protein(&self, suffix: i64) -> u32;
    fn prefetch_for_suffix(&self, suffix: i64);
    fn implied_text_len(&self) -> Option<usize>;
});

impl ReadBinary for InMemorySuffixToProteinMapping {
    fn read_binary<R: std::io::BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut type_buf = [0u8; 1];
        reader.read_exact(&mut type_buf)?;
        match type_buf[0] {
            0 => Ok(InMemorySuffixToProteinMapping::Dense(dense::read_dense_mapping(reader)?)),
            1 => Ok(InMemorySuffixToProteinMapping::Sparse(sparse::read_sparse_mapping(reader)?)),
            2 => Ok(InMemorySuffixToProteinMapping::BitVec(bitvec::read_bitvec_mapping(reader)?)),
            t => Err(format!("Unknown mapping type byte: {}", t).into())
        }
    }
}

impl LoadIndex for InMemorySuffixToProteinMapping {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        binary_traits::load_owned(path)
    }
}

#[cfg(test)]
mod tests {
    use binary_traits::ReadBinary;
    use text_compression::ProteinTextBackend;

    use crate::suffix_to_protein_index::{
        preloaded::{
            BitVecSuffixToProtein, DenseSuffixToProtein, InMemorySuffixToProteinMapping, SparseSuffixToProtein,
            test_utils::assert_dump_and_load
        },
        test_utils::sample_text
    };

    // The mappings below are built through `from_text_parts`, the constructor `sa-builder` uses.

    #[test]
    fn test_dump_and_load_mapping_dense() {
        let text = sample_text();
        assert_dump_and_load(DenseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)), 0u8);
    }

    #[test]
    fn test_dump_and_load_mapping_sparse() {
        let text = sample_text();
        assert_dump_and_load(SparseSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)), 1u8);
    }

    #[test]
    fn test_dump_and_load_mapping_bitvec() {
        let text = sample_text();
        assert_dump_and_load(BitVecSuffixToProtein::from_text_parts(text.len(), |i| text.get(i)), 2u8);
    }

    /// A file that is not an index — or one cut short of what its header promises — must come back
    /// as an error rather than a mapping that answers nonsense.
    #[test]
    fn test_load_mapping_rejects_malformed_input() {
        let cases: [(&str, Vec<u8>); 6] = [
            ("empty input", vec![]),
            ("unknown type byte", vec![99u8]),
            ("dense without a count header", vec![0u8, 0, 0]),
            ("sparse promising four values it does not hold", vec![1u8, 4, 0, 0, 0, 0, 0, 0, 0]),
            ("bitvec without a full header", vec![2u8, 8, 0, 0, 0, 0, 0, 0, 0]),
            (
                // The same crafted header `mmap::tests` rejects, so both readers stay pinned to the
                // one check that relates `bit_len` to `block_count`.
                "bitvec bit_len needing more blocks than the header declares",
                [vec![2u8], 10_000u64.to_le_bytes().to_vec(), 1u64.to_le_bytes().to_vec(), vec![0u8; 8], vec![0u8; 16]]
                    .concat()
            )
        ];

        for (case, buf) in cases {
            let result = InMemorySuffixToProteinMapping::read_binary(&mut std::io::Cursor::new(buf));
            assert!(result.is_err(), "{} was accepted", case);
        }
    }
}
