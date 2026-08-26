use std::error::Error;

use memmap2::Mmap;

use super::super::SuffixToProteinMappingBackend;
use crate::Nullable;

/// Mapping backed by a memory-mapped Sparse binary file.
/// Format: [1 byte type=0x01] [8 bytes count (u64 LE)] [count × 8 bytes (i64 LE)]
pub struct MmapSparseSuffixToProtein {
    mmap: Mmap,
    data_offset: usize, // 9 = 1 (type) + 8 (count)
    count: usize
}

impl SuffixToProteinMappingBackend for MmapSparseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let read_val = |i: usize| -> i64 {
            let off = self.data_offset + i * 8;
            i64::from_le_bytes(self.mmap[off..off + 8].try_into().unwrap())
        };

        let mut lo = 0usize;
        let mut hi = self.count;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if read_val(mid) <= suffix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let protein_index = lo - 1;

        if read_val(protein_index + 1) == suffix + 1 {
            return u32::NULL;
        }
        protein_index as u32
    }

    fn touch_all_pages(&self) -> u64 {
        let end = self.data_offset + self.count * 8;
        memory_hints::warmup::touch_all_pages(&self.mmap, self.data_offset..end)
    }
}

/// Maps a sparse mapping file, validating the header *and* the body it declares.
///
/// The count is checked against the mapping's actual length before the struct is built, so every
/// lookup below can index the mapping without re-checking. It used to be validated only as far as
/// the 9-byte header, which meant a file whose body was short — a build killed part-way, a partial
/// copy, a full disk — loaded cleanly and then panicked on the first lookup, inside a request
/// handler. That is what [`ReadBinaryMmap`](binary_traits::ReadBinaryMmap)'s contract forbids,
/// and what the bitvec reader next door already did correctly.
///
/// This also bounds `touch_all_pages`, whose sweep is `data_offset .. data_offset + count * 8`:
/// with the count validated, that range is guaranteed to lie inside the mapping, so the sweep no
/// longer needs a clamp the dense sibling has and this one lacked.
pub(super) fn read_sparse_mmap(mmap: Mmap) -> Result<MmapSparseSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 9 {
        return Err("The sparse mapping file is too small to contain the count header".into());
    }
    let count = u64::from_le_bytes(mmap[1..9].try_into()?) as usize;

    // Checked, because `count` is untrusted: an unchecked `* 8` wraps for a header near
    // `usize::MAX` and would make `expected` small enough to pass the comparison below.
    let expected = count
        .checked_mul(8)
        .and_then(|n| n.checked_add(9))
        .ok_or("The sparse mapping header declares too many entries")?;
    if mmap.len() < expected {
        return Err(format!(
            "The sparse mapping file is too small to contain the mapping data: expected {expected} bytes, got {}",
            mmap.len()
        )
        .into());
    }

    Ok(MmapSparseSuffixToProtein { mmap, count, data_offset: 9 })
}

#[cfg(test)]
mod tests {
    use protein_text::ProteinTextBackend;

    use crate::suffix_to_protein_index::{
        mmap::test_utils::write_and_map,
        preloaded::SparseSuffixToProtein,
        test_utils::{assert_agree, many_proteins_text, sample_text}
    };

    /// The absolute answers are pinned by `mmap::tests::test_load_mmap_sparse`; what this adds is
    /// that the hand-written binary search over the mapped file agrees with `Vec::binary_search`
    /// in the preloaded mapping, at every position and over enough proteins to recurse a few
    /// levels deep.
    #[test]
    fn test_mmap_sparse_roundtrip() {
        for text in [sample_text(), many_proteins_text(300, 5)] {
            let (loaded, _tmp) = write_and_map(SparseSuffixToProtein::new(&text));
            assert_agree(&SparseSuffixToProtein::new(&text), &loaded, text.len());
        }
    }
}
