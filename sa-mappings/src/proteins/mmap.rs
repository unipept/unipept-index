// This entire file is mmap-only.
use std::{error::Error, path::Path, sync::Arc};

use memmap2::Mmap;
use text_compression::{ReadBinaryMmap, ProteinText, bit_array_byte_size};

use super::{ProteinRef, ProteinsBackend};

// ── MmapBackedProteins ────────────────────────────────────────────────────────

pub struct MmapBackedProteins {
    pub mmap: Arc<Mmap>,
    pub text: ProteinText,
    pub protein_count: usize,
    pub(crate) fixed_table_offset: usize,
    pub(crate) uid_data_offset: usize,
    pub(crate) fa_data_offset: usize,
}

mod entry_offsets {
    pub const TAXON_ID:   std::ops::Range<usize> = 0..4;
    pub const UID_OFFSET: std::ops::Range<usize> = 4..8;
    pub const UID_LEN:    std::ops::Range<usize> = 8..10;
    pub const FA_OFFSET:  std::ops::Range<usize> = 10..14;
    pub const FA_LEN:     std::ops::Range<usize> = 14..16;
    pub const ENTRY_SIZE: usize = 16;
}

impl ProteinsBackend for MmapBackedProteins {
    type Text = ProteinText;

    #[inline]
    fn text(&self) -> &ProteinText { &self.text }
    fn len(&self) -> usize { self.protein_count }

    fn touch_all_pages(&self) {
        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Sequential);
        for chunk in self.mmap.chunks(4096) { std::hint::black_box(chunk[0]); }
        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Random);
    }

    #[inline]
    fn prefetch(&self, index: usize) {
        let off = self.fixed_table_offset + index * entry_offsets::ENTRY_SIZE;
        if off + entry_offsets::ENTRY_SIZE <= self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[off] as *const u8);
        }
    }

    #[inline]
    fn prefetch_strings(&self, index: usize) {
        use entry_offsets as eo;
        let entry_off = self.fixed_table_offset + index * eo::ENTRY_SIZE;
        if entry_off + eo::ENTRY_SIZE > self.mmap.len() { return; }
        let uid_off = u32::from_le_bytes(self.mmap[entry_off + eo::UID_OFFSET.start..entry_off + eo::UID_OFFSET.end].try_into().unwrap()) as usize;
        let fa_off  = u32::from_le_bytes(self.mmap[entry_off + eo::FA_OFFSET.start ..entry_off + eo::FA_OFFSET.end ].try_into().unwrap()) as usize;
        let uid_ptr = self.uid_data_offset + uid_off;
        let fa_ptr  = self.fa_data_offset  + fa_off;
        if uid_ptr < self.mmap.len() { prefetch::prefetch_read(&self.mmap[uid_ptr] as *const u8); }
        if fa_ptr  < self.mmap.len() { prefetch::prefetch_read(&self.mmap[fa_ptr]  as *const u8); }
    }

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
    use std::{fs::File, io::Write, path::PathBuf};
    use tempdir::TempDir;
    use text_compression::{ReadBinaryMmap, WriteBinary};
    use super::MmapBackedProteins;
    use crate::proteins::{ProteinsBackend as _, preloaded::InMemoryProteins};

    fn create_database_file(tmp_dir: &TempDir) -> PathBuf {
        let path = tmp_dir.path().join("database.tsv");
        let mut f = File::create(&path).unwrap();
        f.write_all("P12345\t1\tMLPGLALLLLAAWTARALEV\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        f.write_all("P54321\t2\tPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        f.write_all("P67890\t6\tKWDSDPSGTKTCIDT\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        path
    }

    fn write_binary_to_tempfile(tmp_dir: &TempDir) -> PathBuf {
        let db = create_database_file(tmp_dir);
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
}
