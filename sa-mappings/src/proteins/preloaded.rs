//! Protein metadata held in owned memory, plus the TSV loader that builds an index.
//!
//! This module owns the `WriteBinary` implementation, so `sa-builder` uses it to produce the
//! `proteins.bin` that either backend later reads — which is why the builder never mentions a
//! backend at all.
//!
//! The struct is generic over its text backend, so "owned metadata" does not imply "owned text" —
//! see [`super`] for the two axes. The reader here handles the both-owned case; the pairings that
//! involve the mapping live in [`super::mmap`], which has the mapping to hand.

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
    str::from_utf8
};

use bytelines::ByteLines;
use fa_compression::algorithm1::encode;
use text_compression::{InMemoryProteinText, ProteinTextBackend};
pub use text_compression::{LoadIndex, ReadBinary, WriteBinary};

use super::{Protein, ProteinRef, ProteinsBackend, SEPARATION_CHARACTER, TERMINATION_CHARACTER};

// ── InMemoryProteins ──────────────────────────────────────────────────────────

/// All protein metadata held in owned memory, plus the concatenated text.
///
/// `T` is the text backend. At [`InMemoryProteinText`] metadata and text are both owned; it may
/// instead be instantiated at `MmapBackedProteinText` to keep the (much larger) metadata in RAM
/// while the text stays mapped. The two axes are independent.
pub struct InMemoryProteins<T> {
    /// The concatenated protein text the suffix array is built over.
    pub text: T,
    /// Metadata per protein, in the same order as the runs in `text`.
    pub proteins: Vec<Protein>
}

impl<T> InMemoryProteins<T> {
    /// Pairs an already-built text with the protein table describing it. The two must agree:
    /// `proteins[i]` describes the i-th `-`-separated run in `text`.
    pub fn new(text: T, proteins: Vec<Protein>) -> Self {
        Self { text, proteins }
    }
}

impl InMemoryProteins<InMemoryProteinText> {
    // ── TSV loader ────────────────────────────────────────────────────────────

    /// Builds an index from a UniProt TSV: `uniprot_id\ttaxon_id\tsequence\tannotations`.
    ///
    /// Sequences are upper-cased and concatenated with [`SEPARATION_CHARACTER`] between them and
    /// [`TERMINATION_CHARACTER`] at the end, which is the layout the suffix array is built over.
    ///
    /// # Errors
    ///
    /// Malformed UTF-8 or an unparseable taxon id. Note that a mid-file I/O error currently ends
    /// the loop and returns the proteins read so far rather than failing, and a row with fewer
    /// than four fields panics — both tracked as known issues.
    pub fn load_from_tsv(file: &str) -> Result<Self, Box<dyn Error>> {
        let mut input_string = String::new();
        let mut proteins: Vec<Protein> = Vec::new();
        let file = File::open(file)?;
        let mut lines = ByteLines::new(BufReader::new(file));
        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');
            let uniprot_id = from_utf8(fields.next().unwrap())?;
            let taxon_id = from_utf8(fields.next().unwrap())?.parse()?;
            let sequence = from_utf8(fields.next().unwrap())?;
            let functional_annotations: Vec<u8> = encode(from_utf8(fields.next().unwrap())?);
            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
            proteins.push(Protein {
                uniprot_id: uniprot_id.to_string(),
                taxon_id,
                functional_annotations
            });
        }
        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());
        proteins.shrink_to_fit();
        let text = InMemoryProteinText::from_string(&input_string);
        Ok(Self { text, proteins })
    }
}

impl<T: ProteinTextBackend + Send + Sync> ProteinsBackend for InMemoryProteins<T> {
    type Text = T;

    #[inline]
    fn text(&self) -> &T {
        &self.text
    }
    fn len(&self) -> usize {
        self.proteins.len()
    }

    #[inline]
    fn get(&self, index: usize) -> ProteinRef<'_> {
        let p = &self.proteins[index];
        ProteinRef {
            uniprot_id: &p.uniprot_id,
            taxon_id: p.taxon_id,
            functional_annotations: &p.functional_annotations
        }
    }

    /// Warms the text, which may still be mapped even though the metadata is owned — that is
    /// exactly the `mmap + preloaded-proteins` build. The `Vec<Protein>` itself is anonymous
    /// memory and has nothing to fault, so it contributes nothing; without this delegation the
    /// whole text section went unwarmed in that build, since no other structure can reach it.
    fn touch_all_pages(&self) -> u64 {
        self.text.touch_all_pages()
    }

    /// Prefetches the `Protein` struct at `index`. Only the struct header is hinted; the
    /// `String` and `Vec` it owns live elsewhere and are not reachable from one cache line.
    #[inline]
    fn prefetch(&self, index: usize) {
        if index < self.proteins.len() {
            prefetch::prefetch_read(&self.proteins[index] as *const Protein);
        }
    }
}

/// Largest byte offset the fixed-size entry table can address, since offsets are stored as `u32`.
const MAX_BLOB_BYTES: u64 = u32::MAX as u64;

/// Rejects a database whose accession or annotation blob would overflow the `u32` entry offsets.
///
/// Kept separate from `write_binary` and pure in its arguments so it can be tested at the limit
/// without materialising four gigabytes of input.
fn check_blob_sizes(uid_bytes_total: u64, fa_bytes_total: u64) -> Result<(), String> {
    for (what, total, flag) in
        [("accession (UID)", uid_bytes_total, "uid_offset"), ("annotation (FA)", fa_bytes_total, "fa_offset")]
    {
        if total > MAX_BLOB_BYTES {
            return Err(format!(
                "the {what} blob is {total} bytes, which exceeds the {MAX_BLOB_BYTES} byte limit \
                 imposed by the u32 `{flag}` field in the proteins.bin entry table. Writing it \
                 would silently wrap the offsets and corrupt every entry past the wrap. Raising \
                 this limit means widening the offsets to u64 and rebuilding existing indexes."
            ));
        }
    }
    Ok(())
}

/// Rejects any single protein whose accession or annotations would overflow the `u16` lengths.
fn check_entry_lengths(proteins: &[Protein]) -> Result<(), String> {
    const MAX: usize = u16::MAX as usize;
    for protein in proteins {
        for (what, len) in
            [("accession", protein.uniprot_id.len()), ("annotations", protein.functional_annotations.len())]
        {
            if len > MAX {
                return Err(format!(
                    "protein {}: {what} is {len} bytes, which exceeds the {MAX} byte limit \
                     imposed by the u16 length field in the proteins.bin entry table.",
                    protein.uniprot_id
                ));
            }
        }
    }
    Ok(())
}

/// On-disk format for `proteins.bin` — written here, read by both backends.
///
/// ```text
/// [ protein text        ]  see text_compression::preloaded's WriteBinary
/// [ protein_count: u64  ]
/// [ uid_bytes_total: u64]
/// [ fa_bytes_total: u64 ]
/// [ fixed table         ]  protein_count entries of 16 bytes:
///                            taxon_id:   u32   bytes  0..4
///                            uid_offset: u32   bytes  4..8
///                            uid_len:    u16   bytes  8..10
///                            fa_offset:  u32   bytes 10..14
///                            fa_len:     u16   bytes 14..16
/// [ UID blob            ]  uid_bytes_total bytes, all accessions concatenated
/// [ FA blob             ]  fa_bytes_total bytes, all encoded annotations concatenated
/// ```
///
/// All integers little-endian. Entries are fixed-size so the mmap reader can index the table
/// directly; the variable-length strings live in the two trailing blobs, addressed by offset.
///
/// The reader's copy of this layout is `mmap::entry_offsets` — keep the two in step.
///
/// # Size limits, and why they are checked
///
/// `uid_offset` and `fa_offset` are `u32`, so neither blob can exceed 4 GiB, and `uid_len` /
/// `fa_len` are `u16`, so no single entry can exceed 64 KiB. Exceeding either used to be silent:
/// the offset accumulators wrap in release builds and the length casts truncate, producing a
/// file whose entries alias each other. Nothing detected it — the index simply returned the
/// wrong accession and annotations for every protein past the wrap.
///
/// Both are now validated up front, before a single byte is written, so an over-large database
/// fails the build loudly instead and leaves no partial file behind.
///
/// For scale: measured over a 573,911-protein UniProt release, accessions run 6.04 bytes and
/// encoded annotations 41.58 bytes per protein. The accession blob therefore has room for ~711M
/// proteins; the annotation blob for ~103M. Raising the annotation ceiling means widening the
/// offsets to `u64` and rebuilding every index, so it is deliberately not done pre-emptively —
/// the check exists so that the need announces itself.
impl<T: WriteBinary> WriteBinary for InMemoryProteins<T> {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        let protein_count = self.proteins.len() as u64;
        let uid_bytes_total: u64 = self.proteins.iter().map(|p| p.uniprot_id.len() as u64).sum();
        let fa_bytes_total: u64 = self.proteins.iter().map(|p| p.functional_annotations.len() as u64).sum();

        // Validate before writing anything, including the text section, so a rejected database
        // leaves no partial file. Passing these two checks is what makes the `u32` accumulators
        // and `u16` casts below provably safe.
        check_blob_sizes(uid_bytes_total, fa_bytes_total)?;
        check_entry_lengths(&self.proteins)?;

        WriteBinary::write_binary(self.text, writer)?;
        writer.write_all(&protein_count.to_le_bytes())?;
        writer.write_all(&uid_bytes_total.to_le_bytes())?;
        writer.write_all(&fa_bytes_total.to_le_bytes())?;
        let mut uid_offset: u32 = 0;
        let mut fa_offset: u32 = 0;
        for protein in &self.proteins {
            // Both casts are checked above.
            let uid_len = protein.uniprot_id.len() as u16;
            let fa_len = protein.functional_annotations.len() as u16;
            writer.write_all(&protein.taxon_id.to_le_bytes())?;
            writer.write_all(&uid_offset.to_le_bytes())?;
            writer.write_all(&uid_len.to_le_bytes())?;
            writer.write_all(&fa_offset.to_le_bytes())?;
            writer.write_all(&fa_len.to_le_bytes())?;
            uid_offset += uid_len as u32;
            fa_offset += fa_len as u32;
        }
        for protein in &self.proteins {
            writer.write_all(protein.uniprot_id.as_bytes())?;
        }
        for protein in self.proteins {
            writer.write_all(&protein.functional_annotations)?;
        }
        Ok(())
    }
}

/// Reads the metadata section — everything after the text — into owned `Protein`s.
///
/// Split out of [`ReadBinary::read_binary`] because the metadata half of `proteins.bin` is read
/// the same way regardless of what the *text* half is doing. The mmap backend reads it straight
/// out of a mapping (`&[u8]` is a `BufRead`) when the text is mapped but the metadata is not, so
/// both storage choices share this one parser rather than each keeping a copy that could drift.
///
/// `reader` must be positioned at the start of the metadata section: immediately after the text,
/// at the `protein_count` header. See the format documented on `impl WriteBinary` above.
pub(super) fn read_metadata_section<R: BufRead>(reader: &mut R) -> Result<Vec<Protein>, Box<dyn Error>> {
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8)?;
    let protein_count = u64::from_le_bytes(buf8) as usize;
    reader.read_exact(&mut buf8)?;
    let uid_bytes_total = u64::from_le_bytes(buf8) as usize;
    reader.read_exact(&mut buf8)?;
    let fa_bytes_total = u64::from_le_bytes(buf8) as usize;
    let mut table = vec![[0u8; 16]; protein_count];
    for entry in table.iter_mut() {
        reader.read_exact(entry)?;
    }
    let mut uid_data = vec![0u8; uid_bytes_total];
    reader.read_exact(&mut uid_data)?;
    let mut fa_data = vec![0u8; fa_bytes_total];
    reader.read_exact(&mut fa_data)?;
    let mut proteins = Vec::with_capacity(protein_count);
    for entry in &table {
        let taxon_id = u32::from_le_bytes(entry[0..4].try_into()?);
        let uid_offset = u32::from_le_bytes(entry[4..8].try_into()?) as usize;
        let uid_len = u16::from_le_bytes(entry[8..10].try_into()?) as usize;
        let fa_offset = u32::from_le_bytes(entry[10..14].try_into()?) as usize;
        let fa_len = u16::from_le_bytes(entry[14..16].try_into()?) as usize;
        let uniprot_id = String::from_utf8(uid_data[uid_offset..uid_offset + uid_len].to_vec())?;
        let functional_annotations = fa_data[fa_offset..fa_offset + fa_len].to_vec();
        proteins.push(Protein { uniprot_id, taxon_id, functional_annotations });
    }
    Ok(proteins)
}

/// Reads the whole file into owned memory.
///
/// The bound admits only [`InMemoryProteinText`] in practice: the mmap text implements
/// [`ReadBinaryMmap`](text_compression::ReadBinaryMmap) instead, since it needs a path to map
/// rather than a stream to consume. The combination "mapped text, owned metadata" therefore lives
/// in the mmap module, which has the mapping to hand.
impl<T: ReadBinary> ReadBinary for InMemoryProteins<T> {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let text = T::read_binary(reader)?;
        let proteins = read_metadata_section(reader)?;
        Ok(Self { text, proteins })
    }
}

/// The only pairing with nothing mapped, and so the only one loaded through the owned route.
///
/// The other three are in the mmap module: this impl is deliberately on the concrete
/// `InMemoryProteins<InMemoryProteinText>` rather than on a generic `T`, because
/// `InMemoryProteins<MmapBackedProteinText>` must map the file even though its metadata is owned.
impl LoadIndex for InMemoryProteins<InMemoryProteinText> {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        text_compression::load_owned(path)
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod size_limit_tests {
    use text_compression::InMemoryProteinText;

    use super::*;

    // ── check_blob_sizes ──────────────────────────────────────────────────────
    //
    // Pure in its arguments precisely so the boundary can be tested without building a 4 GB
    // database. The real write path cannot reach these values in any practical test.

    #[test]
    fn blob_exactly_at_the_limit_is_accepted() {
        assert!(check_blob_sizes(MAX_BLOB_BYTES, MAX_BLOB_BYTES).is_ok());
    }

    #[test]
    fn oversized_uid_blob_is_rejected_and_named() {
        let err = check_blob_sizes(MAX_BLOB_BYTES + 1, 0).unwrap_err();
        assert!(err.contains("accession"), "should name the offending blob: {err}");
        assert!(err.contains("uid_offset"), "should name the field: {err}");
        assert!(!err.contains("annotation"), "should not blame the other blob: {err}");
    }

    #[test]
    fn oversized_fa_blob_is_rejected_and_named() {
        let err = check_blob_sizes(0, MAX_BLOB_BYTES + 1).unwrap_err();
        assert!(err.contains("annotation"), "should name the offending blob: {err}");
        assert!(err.contains("fa_offset"), "should name the field: {err}");
    }

    /// The blobs are checked independently — a healthy one must not mask an over-large one.
    #[test]
    fn each_blob_is_checked_independently() {
        assert!(check_blob_sizes(MAX_BLOB_BYTES, MAX_BLOB_BYTES + 1).is_err());
        assert!(check_blob_sizes(MAX_BLOB_BYTES + 1, MAX_BLOB_BYTES).is_err());
        assert!(check_blob_sizes(0, 0).is_ok());
    }

    // ── check_entry_lengths ───────────────────────────────────────────────────

    fn protein_with(uid: &str, fa_len: usize) -> Protein {
        Protein {
            uniprot_id: uid.to_string(),
            taxon_id: 1,
            functional_annotations: vec![0u8; fa_len]
        }
    }

    #[test]
    fn annotation_at_the_u16_limit_is_accepted() {
        assert!(check_entry_lengths(&[protein_with("P12345", u16::MAX as usize)]).is_ok());
    }

    #[test]
    fn oversized_annotation_is_rejected_and_names_the_protein() {
        let err = check_entry_lengths(&[protein_with("P12345", u16::MAX as usize + 1)]).unwrap_err();
        assert!(err.contains("P12345"), "should name the protein: {err}");
        assert!(err.contains("annotations"), "should name the field: {err}");
    }

    #[test]
    fn oversized_accession_is_rejected() {
        let err = check_entry_lengths(&[protein_with(&"A".repeat(u16::MAX as usize + 1), 0)]).unwrap_err();
        assert!(err.contains("accession"), "should name the field: {err}");
    }

    // ── end to end ────────────────────────────────────────────────────────────

    /// A rejected database must leave the writer untouched: validation runs before the text
    /// section, so a failed build cannot produce a half-written index that later loads.
    #[test]
    fn rejection_writes_nothing_at_all() {
        let proteins = InMemoryProteins::new(InMemoryProteinText::from_string("MA-$"), vec![protein_with(
            "P12345",
            u16::MAX as usize + 1
        )]);

        let mut buf: Vec<u8> = Vec::new();
        let err = proteins.write_binary(&mut buf).unwrap_err().to_string();

        assert!(err.contains("P12345"), "{err}");
        assert!(buf.is_empty(), "a rejected database wrote {} bytes", buf.len());
    }

    /// The guards are inert for anything realistic — this is the case that must keep working.
    #[test]
    fn ordinary_database_is_unaffected() {
        let proteins = InMemoryProteins::new(InMemoryProteinText::from_string("MA-CD$"), vec![
            protein_with("P12345", 40),
            protein_with("Q6GZX4", 120),
        ]);

        let mut buf: Vec<u8> = Vec::new();
        proteins.write_binary(&mut buf).unwrap();
        assert!(!buf.is_empty());
    }
}

#[cfg(test)]
mod tests {
    use tempdir::TempDir;
    use text_compression::ProteinTextBackend as _;

    use super::*;
    use crate::proteins::test_fixtures::{TEST_PROTEINS, write_database_file};

    #[test]
    fn test_new_protein() {
        let protein = Protein {
            uniprot_id: "P12345".to_string(),
            taxon_id: 1,
            functional_annotations: vec![0xD1, 0x11]
        };
        assert_eq!(protein.uniprot_id, "P12345");
        assert_eq!(protein.taxon_id, 1);
    }

    #[test]
    fn test_new_proteins() {
        let text = InMemoryProteinText::from_string("MLPGLALLLLAAWTARALEV-PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG");
        let proteins = InMemoryProteins::new(text, vec![
            Protein {
                uniprot_id: "P12345".to_string(),
                taxon_id: 1,
                functional_annotations: vec![0xD1, 0x11]
            },
            Protein {
                uniprot_id: "P54321".to_string(),
                taxon_id: 2,
                functional_annotations: vec![0xD1, 0x11]
            },
        ]);
        assert_eq!(proteins.len(), 2);
        assert_eq!(proteins.get(0).uniprot_id, "P12345");
        assert_eq!(proteins.get(1).taxon_id, 2);
    }

    #[test]
    fn test_get_taxon() {
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();
        let proteins =
            InMemoryProteins::load_from_tsv(write_database_file(&tmp_dir, &TEST_PROTEINS).to_str().unwrap()).unwrap();
        for (i, &taxon) in [1u32, 2, 6, 17].iter().enumerate() {
            assert_eq!(proteins.get(i).taxon_id, taxon);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();
        let proteins =
            InMemoryProteins::load_from_tsv(write_database_file(&tmp_dir, &TEST_PROTEINS).to_str().unwrap()).unwrap();
        for i in 0..proteins.len() {
            assert_eq!(proteins.get(i).get_functional_annotations(), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }

    #[test]
    fn test_write_and_read_binary_buffered() {
        use std::io::BufReader;
        let tmp_dir = TempDir::new("test_binary_roundtrip").unwrap();
        let db = write_database_file(&tmp_dir, &TEST_PROTEINS);
        let original = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();
        let original_save = InMemoryProteins::load_from_tsv(db.to_str().unwrap()).unwrap();
        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        drop(bin_file);
        // The text type has to be named: `read_binary` is generic over it, and inference will not
        // pick the one type implementing `ReadBinary` on its own.
        let loaded =
            InMemoryProteins::<InMemoryProteinText>::read_binary(&mut BufReader::new(File::open(&bin_path).unwrap()))
                .unwrap();
        assert_eq!(loaded.len(), original_save.len());
        for i in 0..original_save.len() {
            assert_eq!(original_save.get(i).uniprot_id, loaded.get(i).uniprot_id);
        }
        assert_eq!(loaded.text().len(), original_save.text().len());
    }
}
