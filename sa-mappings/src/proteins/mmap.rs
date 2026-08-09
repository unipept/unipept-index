//! Protein metadata borrowed from a memory mapping.
//!
//! The counterpart to [`super::preloaded`]: same file, but the accession strings and the encoded
//! annotations are returned as slices into the mapping instead of owned `String`s and `Vec`s.
//! That is what keeps a multi-gigabyte protein table servable in bounded RSS.
//!
//! The file itself is written by `preloaded`'s `WriteBinary`; see there for the format.

use std::{error::Error, path::Path, sync::Arc};

use memmap2::Mmap;
use text_compression::{ReadBinaryMmap, ProteinText, bit_array_byte_size};

use super::{ProteinRef, ProteinsBackend};

// ── MmapBackedProteins ────────────────────────────────────────────────────────

/// Protein table borrowed from a memory mapping.
///
/// The mapping is shared with [`Self::text`], which borrows the text section of the same file.
pub struct MmapBackedProteins {
    /// The mapping of `proteins.bin`, shared with `text`.
    pub mmap: Arc<Mmap>,
    /// The concatenated protein text, borrowing the same mapping.
    pub text: ProteinText,
    /// Number of entries in the fixed-size table.
    pub protein_count: usize,
    pub(crate) fixed_table_offset: usize,
    pub(crate) uid_data_offset: usize,
    pub(crate) fa_data_offset: usize,
}

/// Field layout of one fixed-size entry in the protein table.
///
/// Must stay in lockstep with the writer in [`super::preloaded`], which emits these fields in
/// this order; nothing but this comment ties the two together.
mod entry_offsets {
    /// NCBI taxon id.
    pub const TAXON_ID:   std::ops::Range<usize> = 0..4;
    /// Byte offset of this protein's accession within the UID blob.
    pub const UID_OFFSET: std::ops::Range<usize> = 4..8;
    /// Length of the accession in bytes.
    pub const UID_LEN:    std::ops::Range<usize> = 8..10;
    /// Byte offset of this protein's encoded annotations within the FA blob.
    pub const FA_OFFSET:  std::ops::Range<usize> = 10..14;
    /// Length of the encoded annotations in bytes.
    pub const FA_LEN:     std::ops::Range<usize> = 14..16;
    /// Total size of one entry. Entries are fixed-size so that `get` is O(1).
    pub const ENTRY_SIZE: usize = 16;
}

impl MmapBackedProteins {
    /// Borrows the fixed-size table entry for `index`, or `None` if it falls outside the mapping.
    ///
    /// Both `get` and `prefetch_strings` used to compute this offset independently, in two
    /// different indexing styles; keeping it in one place is what stops them drifting apart.
    #[inline]
    fn entry(&self, index: usize) -> Option<&[u8]> {
        let off = self.fixed_table_offset + index * entry_offsets::ENTRY_SIZE;
        self.mmap.get(off..off + entry_offsets::ENTRY_SIZE)
    }
}

impl ProteinsBackend for MmapBackedProteins {
    type Text = ProteinText;

    #[inline]
    fn text(&self) -> &ProteinText { &self.text }
    fn len(&self) -> usize { self.protein_count }

    fn touch_all_pages(&self) {
        // The whole file: the text, the entry table and both string blobs are all needed.
        let end = self.mmap.len();
        text_compression::mmap::touch_all_pages(&self.mmap, 0..end);
    }

    /// Prefetches the fixed-size table entry for `index`.
    #[inline]
    fn prefetch(&self, index: usize) {
        if let Some(entry) = self.entry(index) {
            prefetch::prefetch_read(entry.as_ptr());
        }
    }

    /// Prefetches the accession and annotation bytes for `index`.
    ///
    /// Separate from [`Self::prefetch`] because they are two extra dependent loads: the entry has
    /// to arrive before its offsets can be followed. Retrieval issues this one batch ahead so the
    /// string data is in flight while the current batch is being decoded.
    #[inline]
    fn prefetch_strings(&self, index: usize) {
        let Some(entry) = self.entry(index) else { return };
        use entry_offsets as eo;
        let uid_off = u32::from_le_bytes(entry[eo::UID_OFFSET].try_into().unwrap()) as usize;
        let fa_off  = u32::from_le_bytes(entry[eo::FA_OFFSET ].try_into().unwrap()) as usize;
        let uid_ptr = self.uid_data_offset + uid_off;
        let fa_ptr  = self.fa_data_offset  + fa_off;
        if uid_ptr < self.mmap.len() { prefetch::prefetch_read(&self.mmap[uid_ptr] as *const u8); }
        if fa_ptr  < self.mmap.len() { prefetch::prefetch_read(&self.mmap[fa_ptr]  as *const u8); }
    }

    /// Decodes the entry at `index` into slices borrowed from the mapping.
    ///
    /// # Panics
    ///
    /// If `index` is out of range, or the file's offsets point outside the mapping. The bounds
    /// are only `debug_assert`ed, so a corrupt file panics in release rather than erroring —
    /// tracked as a known issue.
    #[inline]
    fn get(&self, index: usize) -> ProteinRef<'_> {
        use entry_offsets as eo;
        let entry_off = self.fixed_table_offset + index * eo::ENTRY_SIZE;
        debug_assert!(entry_off + eo::ENTRY_SIZE <= self.mmap.len(), "protein index {index} out of range");
        let entry = &self.mmap[entry_off..entry_off + eo::ENTRY_SIZE];

        let taxon_id   = u32::from_le_bytes(entry[eo::TAXON_ID  ].try_into().unwrap());
        let uid_offset = u32::from_le_bytes(entry[eo::UID_OFFSET].try_into().unwrap()) as usize;
        let uid_len    = u16::from_le_bytes(entry[eo::UID_LEN   ].try_into().unwrap()) as usize;
        let fa_offset  = u32::from_le_bytes(entry[eo::FA_OFFSET ].try_into().unwrap()) as usize;
        let fa_len     = u16::from_le_bytes(entry[eo::FA_LEN    ].try_into().unwrap()) as usize;

        let uid_start = self.uid_data_offset + uid_offset;
        let fa_start  = self.fa_data_offset  + fa_offset;

        ProteinRef {
            uniprot_id: std::str::from_utf8(&self.mmap[uid_start..uid_start + uid_len])
                .expect("invalid UTF-8 in protein UID data"),
            taxon_id,
            functional_annotations: &self.mmap[fa_start..fa_start + fa_len],
        }
    }
}

impl ReadBinaryMmap for MmapBackedProteins {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        use std::fs::File;
        let f = File::open(path)?;
        // SAFETY: see the note in `text_compression::mmap` — an index file is written once by
        // sa-builder and is read-only for the lifetime of the process, so the mapping cannot be
        // truncated or written underneath us.
        let mmap = Arc::new(unsafe { Mmap::map(&f)? });

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        let mmap_len = mmap.len();
        if mmap_len < 8 { return Err("proteins file too short to contain text header".into()); }

        let text_length = u64::from_le_bytes(mmap[0..8].try_into()?) as usize;
        let text_data_offset: usize = 8;
        let bit_array_bytes = bit_array_byte_size(text_length);

        let meta_offset = text_data_offset.checked_add(bit_array_bytes)
            .ok_or_else(|| "overflow while computing metadata offset".to_string())?;
        let meta_end = meta_offset.checked_add(24)
            .ok_or_else(|| "overflow while computing metadata end offset".to_string())?;
        if meta_end > mmap_len { return Err("proteins file too short to contain metadata section".into()); }

        let text = ProteinText::from_mmap(Arc::clone(&mmap), text_data_offset, text_length);

        let protein_count = u64::from_le_bytes(mmap[meta_offset..meta_offset + 8].try_into()?) as usize;
        let uid_bytes_total = u64::from_le_bytes(mmap[meta_offset + 8..meta_offset + 16].try_into()?) as usize;

        let fixed_table_offset = meta_offset.checked_add(24)
            .ok_or_else(|| "overflow while computing fixed table offset".to_string())?;
        let protein_entry_bytes = protein_count.checked_mul(16)
            .ok_or_else(|| "overflow while computing protein table size".to_string())?;
        let uid_data_offset = fixed_table_offset.checked_add(protein_entry_bytes)
            .ok_or_else(|| "overflow while computing uid data offset".to_string())?;
        let fa_data_offset = uid_data_offset.checked_add(uid_bytes_total)
            .ok_or_else(|| "overflow while computing fa data offset".to_string())?;

        if uid_data_offset > mmap_len || fa_data_offset > mmap_len {
            return Err("proteins file truncated: data section offsets exceed file length".into());
        }

        Ok(Self { mmap, text, protein_count, fixed_table_offset, uid_data_offset, fa_data_offset })
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::{fs::File, path::PathBuf};
    use tempdir::TempDir;
    use text_compression::{ReadBinaryMmap, WriteBinary};
    use super::MmapBackedProteins;
    use crate::proteins::{ProteinsBackend as _, preloaded::InMemoryProteins};
    use crate::proteins::test_fixtures::{TEST_PROTEINS, write_database_file};
    use text_compression::{ProteinTextBackend as _, bit_array_byte_size};

    fn write_binary_to_tempfile(tmp_dir: &TempDir) -> PathBuf {
        let db = write_database_file(tmp_dir, &TEST_PROTEINS[..3]);
        let original = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();
        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        bin_path
    }

    #[test]
    fn test_mmap_roundtrip_len() {
        let tmp_dir = TempDir::new("test_mmap_len").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = MmapBackedProteins::read_binary_mmap(&bin_path).unwrap();
        assert_eq!(mmap.len(), 3);
    }

    #[test]
    fn test_mmap_roundtrip_taxon() {
        let tmp_dir = TempDir::new("test_mmap_taxon").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = MmapBackedProteins::read_binary_mmap(&bin_path).unwrap();
        for (i, &taxon) in [1u32, 2, 6].iter().enumerate() {
            assert_eq!(mmap.get(i).taxon_id, taxon);
        }
    }

    #[test]
    fn test_mmap_roundtrip_uniprot_id() {
        let tmp_dir = TempDir::new("test_mmap_uid").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = MmapBackedProteins::read_binary_mmap(&bin_path).unwrap();
        assert_eq!(mmap.get(0).uniprot_id, "P12345");
        assert_eq!(mmap.get(1).uniprot_id, "P54321");
        assert_eq!(mmap.get(2).uniprot_id, "P67890");
    }

    #[test]
    fn test_mmap_roundtrip_functional_annotations() {
        let tmp_dir = TempDir::new("test_mmap_fa").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = MmapBackedProteins::read_binary_mmap(&bin_path).unwrap();
        for i in 0..mmap.len() {
            assert_eq!(mmap.get(i).get_functional_annotations(), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }

    /// Every backend agrees on every field. This is the check that would catch a drift between
    /// the writer's field order in `preloaded` and the reader's `entry_offsets` here, which are
    /// two independent statements of the same 16-byte layout.
    #[test]
    fn matches_the_preloaded_backend_field_for_field() {
        let tmp_dir = TempDir::new("test_mmap_parity").unwrap();
        let db = write_database_file(&tmp_dir, &TEST_PROTEINS);
        let preloaded = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();

        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap()
            .write_binary(&mut bin_file).unwrap();
        let mapped = MmapBackedProteins::read_binary_mmap(&bin_path).unwrap();

        assert_eq!(mapped.len(), preloaded.len());
        for i in 0..preloaded.len() {
            let (a, b) = (preloaded.get(i), mapped.get(i));
            assert_eq!(a.uniprot_id, b.uniprot_id, "uniprot_id differs at {i}");
            assert_eq!(a.taxon_id, b.taxon_id, "taxon_id differs at {i}");
            assert_eq!(a.functional_annotations, b.functional_annotations, "annotations differ at {i}");
        }

        assert_eq!(mapped.text().len(), preloaded.text().len());
        for i in 0..preloaded.text().len() {
            assert_eq!(mapped.text().get(i), preloaded.text().get(i), "text differs at {i}");
        }
    }

    /// `read_binary_mmap` parses an untrusted header. Truncation anywhere up to the end of the
    /// UID blob must error rather than panic or read out of bounds; none of its length checks had
    /// coverage before.
    ///
    /// The sweep deliberately stops at `fa_data_offset`. Beyond that the reader validates
    /// nothing: it parses `fa_bytes_total` from the header but never checks
    /// `fa_data_offset + fa_bytes_total <= mmap.len()`, so a file truncated inside the FA blob
    /// loads happily and only fails later, as an out-of-bounds slice inside `get` — in a request
    /// handler. That is a known issue, reported rather than fixed in this pass; when it is fixed,
    /// extend this sweep to `full.len()` and it should stay green.
    #[test]
    fn truncated_files_error_rather_than_panicking() {
        let tmp_dir = TempDir::new("test_mmap_truncated").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let full = std::fs::read(&bin_path).unwrap();

        // Recompute the section boundaries the way the reader does, rather than hard-coding them.
        let text_length = u64::from_le_bytes(full[0..8].try_into().unwrap()) as usize;
        let meta_offset = 8 + bit_array_byte_size(text_length);
        let protein_count =
            u64::from_le_bytes(full[meta_offset..meta_offset + 8].try_into().unwrap()) as usize;
        let uid_bytes_total =
            u64::from_le_bytes(full[meta_offset + 8..meta_offset + 16].try_into().unwrap()) as usize;
        let fa_data_offset = meta_offset + 24 + protein_count * 16 + uid_bytes_total;
        assert!(fa_data_offset < full.len(), "fixture should have a non-empty FA blob");

        for cut in 0..fa_data_offset {
            let path = tmp_dir.path().join(format!("truncated_{cut}.bin"));
            std::fs::write(&path, &full[..cut]).unwrap();
            assert!(
                MmapBackedProteins::read_binary_mmap(&path).is_err(),
                "{cut} of {} bytes should not load", full.len()
            );
        }
    }

    /// A header declaring more proteins than the file holds must be rejected.
    #[test]
    fn overlong_protein_count_is_rejected() {
        let tmp_dir = TempDir::new("test_mmap_overlong").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mut bytes = std::fs::read(&bin_path).unwrap();

        // The protein count sits just past the text section; locate it the same way the reader
        // does rather than hard-coding an offset.
        let text_length = u64::from_le_bytes(bytes[0..8].try_into().unwrap()) as usize;
        let meta_offset = 8 + bit_array_byte_size(text_length);
        bytes[meta_offset..meta_offset + 8].copy_from_slice(&1_000_000_u64.to_le_bytes());

        let path = tmp_dir.path().join("overlong.bin");
        std::fs::write(&path, &bytes).unwrap();
        assert!(MmapBackedProteins::read_binary_mmap(&path).is_err());
    }
}
