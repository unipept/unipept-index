use std::error::Error;

use memmap2::Mmap;

use super::super::SuffixToProteinMappingBackend;

/// Mapping backed by a memory-mapped Dense binary file.
/// Format: [1 byte type=0x00] [8 bytes count (u64 LE)] [count × 4 bytes (u32 LE)]
pub struct MmapDenseSuffixToProtein {
    mmap: Mmap,
    data_offset: usize, // 9 = 1 (type) + 8 (count)
    /// Entries the header declares. A lookup does not need it — it addresses its entry directly —
    /// but [`Self::touch_all_pages`] does, so that the sweep covers this structure's own entries
    /// rather than everything to the end of the file, the way the sparse one already did.
    count: usize
}

impl SuffixToProteinMappingBackend for MmapDenseSuffixToProtein {
    #[inline]
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let off = self.data_offset + suffix as usize * 4;
        u32::from_le_bytes(self.mmap[off..off + 4].try_into().unwrap())
    }

    #[inline]
    fn implied_text_len(&self) -> Option<usize> {
        Some(self.count)
    }

    #[inline]
    fn prefetch_for_suffix(&self, suffix: i64) {
        let off = self.data_offset + suffix as usize * 4;
        if off < self.mmap.len() {
            // Bounds checked above, so this indexes the mapping directly (Mmap derefs to [u8]).
            prefetch::prefetch_read(&self.mmap[off] as *const u8);
        }
    }

    fn touch_all_pages(&self) -> u64 {
        // This structure's own entries, the way the sparse backend bounds its sweep. Clamped to the
        // mapping, because the header is untrusted and `read_dense_mmap` does not check the body
        // against it.
        let end = (self.data_offset + self.count * 4).min(self.mmap.len());
        text_compression::mmap::touch_all_pages(&self.mmap, self.data_offset..end)
    }
}

/// Maps a dense mapping file, validating the header *and* the body it declares.
///
/// The count is checked against the mapping's actual length before the struct is built, so every
/// lookup below can index the mapping without re-checking. It used to be validated only as far as
/// the 9-byte header, which meant a file whose body was short — a build killed part-way, a partial
/// copy, a full disk — loaded cleanly and then panicked on the first lookup, inside a request
/// handler. That is what [`ReadBinaryMmap`](text_compression::ReadBinaryMmap)'s contract forbids,
/// and what the bitvec reader next door already did correctly.
pub(super) fn read_dense_mmap(mmap: Mmap) -> Result<MmapDenseSuffixToProtein, Box<dyn Error>> {
    if mmap.len() < 9 {
        return Err("The dense mapping file is too small to contain the count header".into());
    }
    let count = u64::from_le_bytes(mmap[1..9].try_into()?) as usize;

    // Checked, because `count` is untrusted: an unchecked `* 4` wraps for a header near
    // `usize::MAX` and would make `expected` small enough to pass the comparison below.
    let expected = count
        .checked_mul(4)
        .and_then(|n| n.checked_add(9))
        .ok_or("The dense mapping header declares too many entries")?;
    if mmap.len() < expected {
        return Err(format!(
            "The dense mapping file is too small to contain the mapping data: expected {expected} bytes, got {}",
            mmap.len()
        )
        .into());
    }

    Ok(MmapDenseSuffixToProtein { mmap, data_offset: 9, count })
}

#[cfg(test)]
mod tests {
    use text_compression::ProteinTextBackend;

    use crate::suffix_to_protein_index::{
        mmap::test_utils::write_and_map,
        preloaded::DenseSuffixToProtein,
        test_utils::{assert_agree, sample_text}
    };

    /// The absolute answers are pinned by `mmap::tests::test_load_mmap_dense`; what this adds is
    /// that reading them out of the mapped file agrees with the preloaded mapping everywhere.
    #[test]
    fn test_mmap_dense_roundtrip() {
        let text = sample_text();
        let (loaded, _tmp) = write_and_map(DenseSuffixToProtein::new(&text));
        assert_agree(&DenseSuffixToProtein::new(&text), &loaded, text.len());
    }
}
