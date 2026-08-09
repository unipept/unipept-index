//! Protein metadata held in owned memory, plus the TSV loader that builds an index.
//!
//! Compiled in *both* configurations: this module owns the `WriteBinary` implementation, so
//! `sa-builder` uses it to produce the `proteins.bin` that the mmap backend later reads. Only the
//! reading half is configuration-specific.

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    str::from_utf8
};

use bytelines::ByteLines;
use fa_compression::algorithm1::encode;
use text_compression::InMemoryProteinText;
pub use text_compression::{ReadBinary, WriteBinary};

use super::{Protein, ProteinRef, ProteinsBackend, SEPARATION_CHARACTER, TERMINATION_CHARACTER};

// ── InMemoryProteins ──────────────────────────────────────────────────────────

/// All protein metadata plus the concatenated text, held in owned memory.
pub struct InMemoryProteins {
    /// The concatenated protein text the suffix array is built over.
    pub text: InMemoryProteinText,
    /// Metadata per protein, in the same order as the runs in `text`.
    pub proteins: Vec<Protein>
}

impl InMemoryProteins {
    /// Pairs an already-built text with the protein table describing it. The two must agree:
    /// `proteins[i]` describes the i-th `-`-separated run in `text`.
    pub fn new(text: InMemoryProteinText, proteins: Vec<Protein>) -> Self {
        Self { text, proteins }
    }

    // ── TSV loaders (non-mmap only) ───────────────────────────────────────────

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

    /// Builds only the concatenated text from a TSV, skipping the protein metadata.
    pub fn text_from_tsv(database_file: &str) -> Result<InMemoryProteinText, Box<dyn Error>> {
        Ok(InMemoryProteinText::from_string(&Self::read_sequences_from_tsv(database_file)?))
    }

    fn read_sequences_from_tsv(database_file: &str) -> Result<String, Box<dyn Error>> {
        let mut input_string = String::new();
        let file = File::open(database_file)?;
        let mut lines = ByteLines::new(BufReader::new(file));
        while let Some(Ok(line)) = lines.next() {
            let sequence = from_utf8(line.split(|b| *b == b'\t').nth(2).unwrap())?;
            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
        }
        Ok(input_string)
    }
}

impl ProteinsBackend for InMemoryProteins {
    type Text = InMemoryProteinText;

    #[inline]
    fn text(&self) -> &InMemoryProteinText {
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

    /// Prefetches the `Protein` struct at `index`. Only the struct header is hinted; the
    /// `String` and `Vec` it owns live elsewhere and are not reachable from one cache line.
    #[inline]
    fn prefetch(&self, index: usize) {
        if index < self.proteins.len() {
            prefetch::prefetch_read(&self.proteins[index] as *const Protein);
        }
    }
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
/// # Known limit
///
/// `uid_offset` and `fa_offset` are `u32`, so the format cannot address blobs beyond 4 GiB, and
/// the accumulators below wrap silently in release once they pass it. At full UniProt scale the
/// UID blob is within an order of magnitude of that ceiling. Tracked as a known issue; fixing it
/// means widening the fields and rebuilding every index.
impl WriteBinary for InMemoryProteins {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        WriteBinary::write_binary(self.text, writer)?;
        let protein_count = self.proteins.len() as u64;
        let uid_bytes_total: u64 = self.proteins.iter().map(|p| p.uniprot_id.len() as u64).sum();
        let fa_bytes_total: u64 = self.proteins.iter().map(|p| p.functional_annotations.len() as u64).sum();
        writer.write_all(&protein_count.to_le_bytes())?;
        writer.write_all(&uid_bytes_total.to_le_bytes())?;
        writer.write_all(&fa_bytes_total.to_le_bytes())?;
        let mut uid_offset: u32 = 0;
        let mut fa_offset: u32 = 0;
        for protein in &self.proteins {
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

impl ReadBinary for InMemoryProteins {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let text = InMemoryProteinText::read_binary(reader)?;
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
        Ok(Self { text, proteins })
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use tempdir::TempDir;
    use text_compression::ProteinTextBackend as _;

    use super::{
        super::test_fixtures::{TEST_PROTEINS, write_database_file},
        *
    };

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
    fn test_get_concatenated_proteins() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();
        let text =
            InMemoryProteins::text_from_tsv(write_database_file(&tmp_dir, &TEST_PROTEINS).to_str().unwrap()).unwrap();
        assert_eq!(text.get(4), b'L');
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
        let loaded = InMemoryProteins::read_binary(&mut BufReader::new(File::open(&bin_path).unwrap())).unwrap();
        assert_eq!(loaded.len(), original_save.len());
        for i in 0..original_save.len() {
            assert_eq!(original_save.get(i).uniprot_id, loaded.get(i).uniprot_id);
        }
        assert_eq!(loaded.text().len(), original_save.text().len());
    }
}
