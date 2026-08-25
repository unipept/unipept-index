//! Protein metadata borrowed from a memory mapping.
//!
//! The counterpart to [`super::preloaded`]: same file, but the accession strings and the encoded
//! annotations are returned as slices into the mapping instead of owned `String`s and `Vec`s.
//! That is what keeps a multi-gigabyte protein table servable in bounded RSS.
//!
//! The file itself is written by `preloaded`'s `WriteBinary`; see there for the format.
//!
//! This module also owns every loader that needs the mapping — including the one producing an
//! *owned* metadata table alongside a mapped text, which lives here rather than in `preloaded`
//! because it is the mapping that has to be opened and kept alive. The three of them share
//! `layout` for the header and `preloaded`'s `read_metadata_section` for the entries, so the
//! storage combinations cannot drift apart in how they read the same file.

use std::{error::Error, path::Path, sync::Arc};

use memmap2::Mmap;
use text_compression::{
    InMemoryProteinText, LoadIndex, MmapBackedProteinText, ProteinTextBackend, ReadBinary, ReadBinaryMmap,
    bit_array_byte_size
};

use super::{ProteinRef, ProteinsBackend, preloaded::read_metadata_section};

// ── MmapBackedProteins ────────────────────────────────────────────────────────

/// Protein table borrowed from a memory mapping.
///
/// `T` is the text backend. At [`MmapBackedProteinText`] the text borrows the text section of the
/// very same mapping; at [`InMemoryProteinText`] it is copied into owned RAM while the much larger
/// metadata table stays mapped. The text is the hottest structure in the index and the metadata
/// the biggest, so that second pairing is the point of the parameter.
pub struct MmapBackedProteins<T> {
    /// The mapping of `proteins.bin`, shared with `text` when the text is mapped too.
    pub mmap: Arc<Mmap>,
    /// The concatenated protein text.
    pub text: T,
    /// Number of entries in the fixed-size table.
    pub protein_count: usize,
    pub(crate) fixed_table_offset: usize,
    pub(crate) uid_data_offset: usize,
    pub(crate) fa_data_offset: usize,
    /// First byte [`ProteinsBackend::touch_all_pages`] warms.
    ///
    /// `0` when the text borrows the mapping, so the sweep covers the text section too; the start
    /// of the metadata section when the text was copied into owned memory, since those text pages
    /// are then never read again and faulting them in would be pure startup cost.
    pub(crate) warm_from: usize
}

/// Field layout of one fixed-size entry in the protein table.
///
/// Must stay in lockstep with the writer in [`super::preloaded`], which emits these fields in
/// this order. The two are independent statements of the same 16-byte layout; what catches a
/// drift between them is the `matches_the_preloaded_backend_field_for_field` test below, which
/// writes a file with one backend and reads it back with the other.
mod entry_offsets {
    /// NCBI taxon id.
    pub const TAXON_ID: std::ops::Range<usize> = 0..4;
    /// Byte offset of this protein's accession within the UID blob.
    pub const UID_OFFSET: std::ops::Range<usize> = 4..8;
    /// Length of the accession in bytes.
    pub const UID_LEN: std::ops::Range<usize> = 8..10;
    /// Byte offset of this protein's encoded annotations within the FA blob.
    pub const FA_OFFSET: std::ops::Range<usize> = 10..14;
    /// Length of the encoded annotations in bytes.
    pub const FA_LEN: std::ops::Range<usize> = 14..16;
    /// Total size of one entry. Entries are fixed-size so that `get` is O(1).
    pub const ENTRY_SIZE: usize = 16;
}

impl<T> MmapBackedProteins<T> {
    /// Borrows the fixed-size table entry for `index`, or `None` if it falls outside the mapping.
    ///
    /// Shared by [`Self::prefetch`] and [`Self::prefetch_strings`], which need the `Option`: a
    /// hint is speculative and must skip an index it cannot address rather than fault on it.
    /// [`Self::get`] computes the same offset itself instead of going through here, because it
    /// wants the panicking slice — an unaddressable entry there means a corrupt file, which is
    /// not something to pass over silently.
    #[inline]
    fn entry(&self, index: usize) -> Option<&[u8]> {
        let off = self.fixed_table_offset + index * entry_offsets::ENTRY_SIZE;
        self.mmap.get(off..off + entry_offsets::ENTRY_SIZE)
    }
}

impl<T: ProteinTextBackend + Send + Sync> ProteinsBackend for MmapBackedProteins<T> {
    type Text = T;

    #[inline]
    fn text(&self) -> &T {
        &self.text
    }
    fn len(&self) -> usize {
        self.protein_count
    }

    fn touch_all_pages(&self) -> u64 {
        // Everything from `warm_from` on: the entry table and both string blobs always, and the
        // text section too when the text is mapped rather than owned. See the field's doc.
        let end = self.mmap.len();
        text_compression::mmap::touch_all_pages(&self.mmap, self.warm_from..end)
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
    /// Separate from [`Self::prefetch`] because they are two extra *dependent* loads: the entry
    /// has to arrive before its offsets can be followed.
    ///
    /// That dependency is why nothing calls this today. `Searcher::retrieve_proteins` in
    /// `sa-index/src/sa_searcher/retrieval.rs` deliberately issues only [`Self::prefetch`]: at its
    /// look-ahead distance the entry has usually not landed yet when this would read it, so the
    /// call stalls on the very load it is meant to hide. Its doc comment records the reasoning.
    /// A caller with a longer look-ahead could use this; the current one cannot.
    #[inline]
    fn prefetch_strings(&self, index: usize) {
        let Some(entry) = self.entry(index) else { return };
        use entry_offsets as eo;
        let uid_off = u32::from_le_bytes(entry[eo::UID_OFFSET].try_into().unwrap()) as usize;
        let fa_off = u32::from_le_bytes(entry[eo::FA_OFFSET].try_into().unwrap()) as usize;
        let uid_ptr = self.uid_data_offset + uid_off;
        let fa_ptr = self.fa_data_offset + fa_off;
        if uid_ptr < self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[uid_ptr] as *const u8);
        }
        if fa_ptr < self.mmap.len() {
            prefetch::prefetch_read(&self.mmap[fa_ptr] as *const u8);
        }
    }

    /// Decodes the entry at `index` into slices borrowed from the mapping.
    ///
    /// # Panics
    ///
    /// Three ways, all of them on a *deliberately edited* file rather than on bad input:
    ///
    /// * `index` is out of range. Only `debug_assert`ed, so a release build panics on the slice
    ///   below instead of erroring.
    /// * the entry's `uid_offset` / `fa_offset` point outside their blob.
    /// * the accession bytes are not valid UTF-8.
    ///
    /// # Why this is not checked here
    ///
    /// A *truncated* file no longer reaches this point: `layout` bounds all four sections —
    /// including the annotation blob, whose length header it used to skip — so the realistic
    /// damage (a build killed part-way, a partial copy, a full disk) is now a load error rather
    /// than a panic in a request handler. What remains is an entry whose offsets were edited to
    /// point elsewhere while the section headers stayed consistent, which truncation cannot
    /// produce.
    ///
    /// Closing that last gap is a deliberate non-goal, not an oversight. Validating every entry at
    /// load is O(protein_count) and would fault in the whole entry table, which is precisely the
    /// lazy startup this backend exists to provide; checking on each lookup would put two bounds
    /// checks on the retrieval hot path. The preloaded sibling *does* check per entry, because it
    /// materialises both blobs anyway and the check is free there — so the two backends differ
    /// here on purpose.
    ///
    /// See `truncation_inside_the_annotation_blob_is_rejected_by_both_backends` for what the load
    /// path now catches.
    #[inline]
    fn get(&self, index: usize) -> ProteinRef<'_> {
        use entry_offsets as eo;
        let entry_off = self.fixed_table_offset + index * eo::ENTRY_SIZE;
        debug_assert!(entry_off + eo::ENTRY_SIZE <= self.mmap.len(), "protein index {index} out of range");
        let entry = &self.mmap[entry_off..entry_off + eo::ENTRY_SIZE];

        let taxon_id = u32::from_le_bytes(entry[eo::TAXON_ID].try_into().unwrap());
        let uid_offset = u32::from_le_bytes(entry[eo::UID_OFFSET].try_into().unwrap()) as usize;
        let uid_len = u16::from_le_bytes(entry[eo::UID_LEN].try_into().unwrap()) as usize;
        let fa_offset = u32::from_le_bytes(entry[eo::FA_OFFSET].try_into().unwrap()) as usize;
        let fa_len = u16::from_le_bytes(entry[eo::FA_LEN].try_into().unwrap()) as usize;

        let uid_start = self.uid_data_offset + uid_offset;
        let fa_start = self.fa_data_offset + fa_offset;

        ProteinRef {
            uniprot_id: std::str::from_utf8(&self.mmap[uid_start..uid_start + uid_len])
                .expect("invalid UTF-8 in protein UID data"),
            taxon_id,
            functional_annotations: &self.mmap[fa_start..fa_start + fa_len]
        }
    }
}

// ── Loading ───────────────────────────────────────────────────────────────────

/// Where each section of `proteins.bin` starts, once the header has been validated.
///
/// Every loader below needs these offsets, so they are computed once here rather than three times.
struct Layout {
    /// Residues in the text, i.e. its length before 5-bit packing.
    text_length: usize,
    /// Byte offset of the packed text data. Constant, but named so the loaders do not repeat it.
    text_data_offset: usize,
    /// Byte offset of the `protein_count` header, i.e. the start of the metadata section.
    meta_offset: usize,
    protein_count: usize,
    fixed_table_offset: usize,
    uid_data_offset: usize,
    fa_data_offset: usize
}

/// Parses and bounds-checks the header of a mapped `proteins.bin`.
///
/// Split out of the loaders so that all three read the same header the same way — a divergence
/// here would mean one storage combination validating less than another. It deliberately does
/// *all* the checking before any loader builds a text or copies metadata, which is what keeps the
/// guarantee that a truncated file errors rather than panicking.
///
/// The checking still stops at `fa_data_offset`; see `truncated_files_error_rather_than_panicking`
/// for what is not validated beyond it.
fn layout(mmap: &Mmap) -> Result<Layout, Box<dyn Error>> {
    let mmap_len = mmap.len();
    if mmap_len < 8 {
        return Err("The proteins file is too small to contain the text header".into());
    }

    let text_length = u64::from_le_bytes(mmap[0..8].try_into()?) as usize;
    let text_data_offset: usize = 8;
    let bit_array_bytes = bit_array_byte_size(text_length)
        .ok_or_else(|| "The proteins header declares an implausible text length".to_string())?;

    let meta_offset = text_data_offset
        .checked_add(bit_array_bytes)
        .ok_or_else(|| "overflow while computing metadata offset".to_string())?;
    let meta_end = meta_offset
        .checked_add(24)
        .ok_or_else(|| "overflow while computing metadata end offset".to_string())?;
    if meta_end > mmap_len {
        return Err("The proteins file is too small to contain the metadata header".into());
    }

    let protein_count = u64::from_le_bytes(mmap[meta_offset..meta_offset + 8].try_into()?) as usize;
    let uid_bytes_total = u64::from_le_bytes(mmap[meta_offset + 8..meta_offset + 16].try_into()?) as usize;
    // The third header field. It used to be skipped entirely, which left the annotation blob as
    // the one section never bounded against the file: a `proteins.bin` truncated anywhere inside
    // it mapped cleanly and panicked later, in a query for one of the last proteins.
    let fa_bytes_total = u64::from_le_bytes(mmap[meta_offset + 16..meta_offset + 24].try_into()?) as usize;

    let fixed_table_offset = meta_offset
        .checked_add(24)
        .ok_or_else(|| "overflow while computing fixed table offset".to_string())?;
    let protein_entry_bytes = protein_count
        .checked_mul(16)
        .ok_or_else(|| "overflow while computing protein table size".to_string())?;
    let uid_data_offset = fixed_table_offset
        .checked_add(protein_entry_bytes)
        .ok_or_else(|| "overflow while computing uid data offset".to_string())?;
    let fa_data_offset = uid_data_offset
        .checked_add(uid_bytes_total)
        .ok_or_else(|| "overflow while computing fa data offset".to_string())?;

    let fa_data_end = fa_data_offset
        .checked_add(fa_bytes_total)
        .ok_or_else(|| "overflow while computing fa data end offset".to_string())?;

    if uid_data_offset > mmap_len || fa_data_offset > mmap_len || fa_data_end > mmap_len {
        return Err("The proteins file is too small to contain the data sections its header declares".into());
    }

    Ok(Layout {
        text_length,
        text_data_offset,
        meta_offset,
        protein_count,
        fixed_table_offset,
        uid_data_offset,
        fa_data_offset
    })
}

/// Maps `path` read-only and advises random access, the way every loader here wants it.
fn map_file(path: &Path) -> Result<Arc<Mmap>, Box<dyn Error>> {
    use std::fs::File;
    let f = File::open(path)?;
    // SAFETY: see the note in `text_compression::mmap` — an index file is written once by
    // sa-builder and is read-only for the lifetime of the process, so the mapping cannot be
    // truncated or written underneath us.
    let mmap = Arc::new(unsafe { Mmap::map(&f)? });

    #[cfg(unix)]
    mmap.advise(memmap2::Advice::Random)?;

    Ok(mmap)
}

/// Metadata and text both borrowed from the mapping — the smallest resident footprint.
impl ReadBinaryMmap for MmapBackedProteins<MmapBackedProteinText> {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let mmap = map_file(path)?;
        let l = layout(&mmap)?;
        let text = MmapBackedProteinText::from_mmap(Arc::clone(&mmap), l.text_data_offset, l.text_length);

        Ok(Self {
            mmap,
            text,
            protein_count: l.protein_count,
            fixed_table_offset: l.fixed_table_offset,
            uid_data_offset: l.uid_data_offset,
            fa_data_offset: l.fa_data_offset,
            // The text borrows the mapping, so its pages are worth warming.
            warm_from: 0
        })
    }
}

/// Metadata mapped, text copied into owned RAM.
///
/// The text is read through the ordinary [`ReadBinary`] path with the mapping itself as the
/// source — `&[u8]` is a `BufRead`, and the text section starts at byte 0 — so the 5-bit unpacking
/// and the huge-page advice come from the one implementation in `text_compression::preloaded`
/// rather than a second copy here.
impl ReadBinaryMmap for MmapBackedProteins<InMemoryProteinText> {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let mmap = map_file(path)?;
        let l = layout(&mmap)?;
        let text = InMemoryProteinText::read_binary(&mut &mmap[..])?;

        Ok(Self {
            mmap,
            text,
            protein_count: l.protein_count,
            fixed_table_offset: l.fixed_table_offset,
            uid_data_offset: l.uid_data_offset,
            fa_data_offset: l.fa_data_offset,
            // The text now lives in owned memory; warming its pages would fault in ~190 MB at
            // UniProt scale that nothing reads again.
            warm_from: l.meta_offset
        })
    }
}

/// Metadata copied into owned RAM, text mapped.
///
/// This is a `ReadBinaryMmap` impl on the *preloaded* struct, which reads oddly until you notice
/// that its text still lives in the mapping: something has to hold the file open, and the metadata
/// is what gets copied out. `read_metadata_section` is the same parser the fully-preloaded reader
/// uses, applied to the mapping instead of a file handle.
impl ReadBinaryMmap for super::InMemoryProteins<MmapBackedProteinText> {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let mmap = map_file(path)?;
        let l = layout(&mmap)?;
        let text = MmapBackedProteinText::from_mmap(Arc::clone(&mmap), l.text_data_offset, l.text_length);

        // Stream the metadata section in before parsing it, rather than letting the parser demand
        // page it. `map_file` sets `MADV_RANDOM` for the query workload this mapping exists to
        // serve, and `read_metadata_section` then walks the section with one `read_exact` per
        // 16-byte entry — so with readahead disabled it faults 8 GB in through ~200 M sixteen-byte
        // reads. Cold, that measured 38.8 MB/s against the ~530 MB/s the same parser reaches over
        // the same bytes from a `BufReader`, and it cost this configuration 218.9 s of a 536.7 s
        // startup (`startup`, b530143049). It went unnoticed for as long as it did because the
        // suite ran this arm behind the mapped one, whose page sweep left the section warm.
        //
        // `touch_all_pages` is the right helper and not just a convenient one: it brackets the walk
        // with `MADV_SEQUENTIAL` and restores `MADV_RANDOM` afterwards, which is exactly the
        // temporary reversal wanted for one bulk copy out of a mapping tuned for random access.
        let _ = text_compression::mmap::touch_all_pages(&mmap, l.meta_offset..mmap.len());
        let proteins = read_metadata_section(&mut &mmap[l.meta_offset..])?;

        Ok(Self::new(text, proteins))
    }
}

// ── LoadIndex ─────────────────────────────────────────────────────────────────

// Three of the four pairings load through the mapping, and this is where that fact is recorded.
// `proteins.bin` holds the text and the metadata in one file, so it has to be *mapped* whenever
// either section is mapped — not merely when the metadata is. The odd one out is the last impl
// below, whose metadata is owned; the fourth pairing, with nothing mapped, takes the owned route
// in `super::preloaded`. This used to be a `#[cfg]` predicate spelled out at every loader.

impl LoadIndex for MmapBackedProteins<MmapBackedProteinText> {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        Self::read_binary_mmap(path)
    }
}

impl LoadIndex for MmapBackedProteins<InMemoryProteinText> {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        Self::read_binary_mmap(path)
    }
}

impl LoadIndex for super::InMemoryProteins<MmapBackedProteinText> {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        Self::read_binary_mmap(path)
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {

    /// The property every `Mmap::map` in the workspace depends on: replacing an index file by
    /// *rename* leaves a live mapping intact.
    ///
    /// `sa-builder` writes each section to a temporary sibling and renames it over the destination,
    /// so a rebuild never truncates a path a running server may have mapped. A rename only unlinks
    /// the old name; the inode survives for as long as anything still holds it open, which is what
    /// makes the mapping safe to keep reading.
    ///
    /// The contrast was measured rather than assumed: truncating the same file *in place* under a
    /// live mapping and then reading it dies with `signal: 10, SIGBUS: access to undefined memory`.
    /// That is the hazard the SAFETY comment describes, and it is exactly one `cp` away — see the
    /// note there.
    #[test]
    fn replacing_an_index_by_rename_keeps_a_live_mapping_valid() {
        use std::io::Write;

        let dir = TempDir::new("rename_probe").unwrap();
        let live = dir.path().join("index.bin");
        std::fs::write(&live, vec![0xAAu8; 4096]).unwrap();

        let f = File::open(&live).unwrap();
        let mapped = unsafe { memmap2::Mmap::map(&f) }.unwrap();
        assert_eq!(mapped[0], 0xAA);

        // What the builder now does: write a sibling, rename it over the destination.
        let tmp = dir.path().join("index.bin.tmp");
        let mut t = File::create(&tmp).unwrap();
        t.write_all(&vec![0xBBu8; 8192]).unwrap();
        t.sync_all().unwrap();
        drop(t);
        std::fs::rename(&tmp, &live).unwrap();

        // The mapping must still see the *old* bytes, at the old length, with no fault.
        assert_eq!(mapped.len(), 4096, "mapping length changed under us");
        let sum: u64 = mapped.iter().map(|&b| b as u64).sum();
        assert_eq!(sum, 0xAA * 4096, "mapping content changed under us");
    }

    use std::{fs::File, path::PathBuf};

    use tempdir::TempDir;
    use text_compression::{ProteinTextBackend as _, ReadBinaryMmap, WriteBinary, bit_array_byte_size};

    use super::{InMemoryProteinText, MmapBackedProteinText, MmapBackedProteins};
    use crate::proteins::{
        ProteinsBackend,
        preloaded::InMemoryProteins,
        test_fixtures::{TEST_PROTEINS, write_database_file}
    };

    /// The four ways `proteins.bin` can be loaded, named once so the tests read as prose.
    ///
    /// Default type parameters are not applied in expression position, and there are three
    /// `ReadBinaryMmap` impls to choose between, so every call site has to say which it means.
    type Mapped = MmapBackedProteins<MmapBackedProteinText>;
    type MappedMetaOwnedText = MmapBackedProteins<InMemoryProteinText>;
    type OwnedMetaMappedText = InMemoryProteins<MmapBackedProteinText>;

    fn write_binary_to_tempfile(tmp_dir: &TempDir) -> PathBuf {
        let db = write_database_file(tmp_dir, &TEST_PROTEINS[..3]);
        let original = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();
        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        bin_path
    }

    /// Both hint methods, walked over every entry and well past the last one.
    ///
    /// Neither had a test. `prefetch` is on the retrieval hot path; `prefetch_strings` is called by
    /// nothing at all today — see its doc for why — which is exactly what makes it worth pinning:
    /// it follows two `u32` offsets read out of the mapping and indexes with them, so a mistake
    /// there would sit unnoticed until someone with a longer look-ahead started calling it. Both
    /// must tolerate an index no entry exists for, and neither may disturb the `get` that follows.
    #[test]
    fn prefetch_hints_are_harmless() {
        let tmp_dir = TempDir::new("test_mmap_prefetch").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let proteins = <Mapped as ReadBinaryMmap>::read_binary_mmap(&bin_path).unwrap();

        for index in 0..proteins.len() * 4 + 16 {
            proteins.prefetch(index);
            proteins.prefetch_strings(index);
        }

        assert_eq!(proteins.len(), 3);
        for (index, (uid, taxon, _, _)) in TEST_PROTEINS[..3].iter().enumerate() {
            let protein = proteins.get(index);
            assert_eq!(protein.uniprot_id, *uid, "uid {index} differs after the hints");
            assert_eq!(protein.taxon_id, *taxon, "taxon {index} differs after the hints");
        }
    }

    #[test]
    fn test_mmap_roundtrip_len() {
        let tmp_dir = TempDir::new("test_mmap_len").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = Mapped::read_binary_mmap(&bin_path).unwrap();
        assert_eq!(mmap.len(), 3);
    }

    #[test]
    fn test_mmap_roundtrip_taxon() {
        let tmp_dir = TempDir::new("test_mmap_taxon").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = Mapped::read_binary_mmap(&bin_path).unwrap();
        for (i, &taxon) in [1u32, 2, 6].iter().enumerate() {
            assert_eq!(mmap.get(i).taxon_id, taxon);
        }
    }

    #[test]
    fn test_mmap_roundtrip_uniprot_id() {
        let tmp_dir = TempDir::new("test_mmap_uid").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = Mapped::read_binary_mmap(&bin_path).unwrap();
        assert_eq!(mmap.get(0).uniprot_id, "P12345");
        assert_eq!(mmap.get(1).uniprot_id, "P54321");
        assert_eq!(mmap.get(2).uniprot_id, "P67890");
    }

    #[test]
    fn test_mmap_roundtrip_functional_annotations() {
        let tmp_dir = TempDir::new("test_mmap_fa").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let mmap = Mapped::read_binary_mmap(&bin_path).unwrap();
        for i in 0..mmap.len() {
            assert_eq!(mmap.get(i).get_functional_annotations(), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }

    /// Asserts `actual` reports exactly what `expected` does, field by field and residue by
    /// residue. `what` names the combination so a failure says which one broke.
    fn assert_agrees_with<A: ProteinsBackend, E: ProteinsBackend>(what: &str, actual: &A, expected: &E) {
        assert_eq!(actual.len(), expected.len(), "{what}: protein count differs");
        for i in 0..expected.len() {
            let (a, b) = (expected.get(i), actual.get(i));
            assert_eq!(a.uniprot_id, b.uniprot_id, "{what}: uniprot_id differs at {i}");
            assert_eq!(a.taxon_id, b.taxon_id, "{what}: taxon_id differs at {i}");
            assert_eq!(a.functional_annotations, b.functional_annotations, "{what}: annotations differ at {i}");
        }

        assert_eq!(actual.text().len(), expected.text().len(), "{what}: text length differs");
        for i in 0..expected.text().len() {
            assert_eq!(actual.text().get(i), expected.text().get(i), "{what}: text differs at {i}");
        }
    }

    /// Every combination agrees on every field. Two independent things ride on this:
    ///
    /// * the writer's field order in `preloaded` against the reader's `entry_offsets` here, which
    ///   are two separate statements of the same 16-byte layout;
    /// * the four text × metadata pairings against each other. They share `layout` and
    ///   `read_metadata_section`, but each assembles the pieces itself, and a wrong offset in one
    ///   would otherwise surface only as wrong answers from a differently-built server.
    #[test]
    fn matches_the_preloaded_backend_field_for_field() {
        let tmp_dir = TempDir::new("test_mmap_parity").unwrap();
        let db = write_database_file(&tmp_dir, &TEST_PROTEINS);
        let preloaded = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();

        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap().write_binary(&mut bin_file).unwrap();

        assert_agrees_with("mapped metadata, mapped text", &Mapped::read_binary_mmap(&bin_path).unwrap(), &preloaded);
        assert_agrees_with(
            "mapped metadata, owned text",
            &MappedMetaOwnedText::read_binary_mmap(&bin_path).unwrap(),
            &preloaded
        );
        assert_agrees_with(
            "owned metadata, mapped text",
            &OwnedMetaMappedText::read_binary_mmap(&bin_path).unwrap(),
            &preloaded
        );
    }

    /// The warm range must skip the text section exactly when the text is not in the mapping.
    ///
    /// Getting this wrong is silent: the pages fault in, nothing ever reads them again, and the
    /// only symptoms are startup time and page-cache pressure — ~190 MB of it at UniProt scale.
    /// That would quietly cancel the benefit `preloaded-text` exists to deliver.
    #[test]
    fn warm_range_skips_the_text_only_when_the_text_is_owned() {
        let tmp_dir = TempDir::new("test_mmap_warm_from").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);

        assert_eq!(Mapped::read_binary_mmap(&bin_path).unwrap().warm_from, 0, "mapped text: warm the whole file");

        let owned_text = MappedMetaOwnedText::read_binary_mmap(&bin_path).unwrap();
        let meta_offset = 8 + bit_array_byte_size(owned_text.text.len()).unwrap();
        assert!(meta_offset > 8, "fixture should have a non-empty text section");
        assert_eq!(owned_text.warm_from, meta_offset, "owned text: warming must start at the metadata section");
    }

    /// Owned metadata over a mapped text is the one pairing where nothing else can reach the text
    /// pages: `MmapBackedProteins` sweeps its own mapping and would cover them, but this pairing
    /// does not have that mapping. It inherited the trait's no-op sweep until it was given one,
    /// which left the whole text section unwarmed in the `mmap + preloaded-proteins` build.
    ///
    /// Residency is not observable from inside the process without `mincore`, so what this pins is
    /// the bound: the sweep slices the mapping by a range it computes from the header, so a wrong
    /// end offset is a panic here rather than a silently short sweep.
    #[test]
    fn owned_metadata_warms_its_mapped_text() {
        let tmp_dir = TempDir::new("test_mmap_warm_owned_meta").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);

        let proteins = OwnedMetaMappedText::read_binary_mmap(&bin_path).unwrap();
        proteins.touch_all_pages();

        assert_eq!(proteins.len(), 3, "the sweep must not disturb what follows it");
        assert_eq!(proteins.text().len(), proteins.text.len());
    }

    /// A named loader, reduced to "did it reject this file?" so the malformed-header tests can
    /// hold all three in one list despite their different return types.
    type RejectsFile = (&'static str, fn(&std::path::Path) -> bool);

    /// Every loader that maps the file, so a malformed-header test covers all of them.
    ///
    /// They share `layout`, which is exactly why this matters: the shared header check has to run
    /// *before* each loader touches the sections it owns, or one of them reads out of bounds on a
    /// file the others reject.
    fn mmap_loaders() -> Vec<RejectsFile> {
        vec![
            ("mapped metadata, mapped text", |p| Mapped::read_binary_mmap(p).is_err()),
            ("mapped metadata, owned text", |p| MappedMetaOwnedText::read_binary_mmap(p).is_err()),
            ("owned metadata, mapped text", |p| OwnedMetaMappedText::read_binary_mmap(p).is_err()),
        ]
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
        let meta_offset = 8 + bit_array_byte_size(text_length).unwrap();
        let protein_count = u64::from_le_bytes(full[meta_offset..meta_offset + 8].try_into().unwrap()) as usize;
        let uid_bytes_total = u64::from_le_bytes(full[meta_offset + 8..meta_offset + 16].try_into().unwrap()) as usize;
        let fa_data_offset = meta_offset + 24 + protein_count * 16 + uid_bytes_total;
        assert!(fa_data_offset < full.len(), "fixture should have a non-empty FA blob");

        for cut in 0..fa_data_offset {
            let path = tmp_dir.path().join(format!("truncated_{cut}.bin"));
            std::fs::write(&path, &full[..cut]).unwrap();
            for (what, load_is_err) in mmap_loaders() {
                assert!(load_is_err(&path), "{what}: {cut} of {} bytes should not load", full.len());
            }
        }
    }

    /// A file truncated inside the annotation blob must be refused by *both* backends.
    ///
    /// `layout` bounded the entry table and the UID blob but never parsed `fa_bytes_total`, so the
    /// annotation section was the one part of the file no reader checked. A `proteins.bin` cut
    /// anywhere inside it mapped cleanly, served the early proteins, and panicked on one of the
    /// last — inside a request handler, long after startup. The owned loader rejected the same
    /// file with "failed to fill whole buffer", so the two disagreed.
    #[test]
    fn truncation_inside_the_annotation_blob_is_rejected_by_both_backends() {
        use text_compression::ReadBinary;

        let tmp_dir = TempDir::new("test_fa_truncation").unwrap();
        let bin_path = write_binary_to_tempfile(&tmp_dir);
        let buf = std::fs::read(&bin_path).unwrap();

        // The annotation blob is the last section, so the final bytes lie inside it.
        for cut in [buf.len() - 1, buf.len() - 2] {
            let short = tmp_dir.path().join(format!("short{cut}.bin"));
            std::fs::write(&short, &buf[..cut]).unwrap();

            assert!(
                Mapped::read_binary_mmap(&short).is_err(),
                "mmap accepted a file truncated to {cut} of {} bytes",
                buf.len()
            );
            assert!(
                InMemoryProteins::<InMemoryProteinText>::read_binary(&mut &buf[..cut]).is_err(),
                "preloaded accepted a file truncated to {cut} of {} bytes",
                buf.len()
            );
        }

        // The intact file still loads on both, so the rejections above are specific.
        assert!(Mapped::read_binary_mmap(&bin_path).is_ok());
        assert!(InMemoryProteins::<InMemoryProteinText>::read_binary(&mut buf.as_slice()).is_ok());
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
        let meta_offset = 8 + bit_array_byte_size(text_length).unwrap();
        bytes[meta_offset..meta_offset + 8].copy_from_slice(&1_000_000_u64.to_le_bytes());

        let path = tmp_dir.path().join("overlong.bin");
        std::fs::write(&path, &bytes).unwrap();
        for (what, load_is_err) in mmap_loaders() {
            assert!(load_is_err(&path), "{what}: an overlong protein count should be rejected");
        }
    }
}
