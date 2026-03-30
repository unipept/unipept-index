//! This module contains the `Protein` and `Proteins` structs, which are used to represent proteins
//! and collections of proteins, respectively.

use std::{error::Error, fs::File, io::{BufRead, BufReader, Write}, path::Path, str::from_utf8, sync::Arc};

use bytelines::ByteLines;
use fa_compression::algorithm1::{decode, encode};
use memmap2::Mmap;
pub use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};
use text_compression::{ProteinText, bit_array_byte_size};

/// The separation character used in the input string
pub static SEPARATION_CHARACTER: u8 = b'-';

/// The termination character used in the input string
/// This character should be smaller than the separation character
pub static TERMINATION_CHARACTER: u8 = b'$';

/// A struct that represents a protein and its linked information
pub struct Protein {
    /// The id of the protein
    pub uniprot_id: String,

    /// the taxon id of the protein
    pub taxon_id: u32,

    /// The encoded functional annotations of the protein
    pub functional_annotations: Vec<u8>
}

impl Protein {
    /// Returns the decoded functional annotations of the protein
    pub fn get_functional_annotations(&self) -> String {
        decode(&self.functional_annotations)
    }
}

/// A zero-copy view into a single protein's metadata, borrowing from an in-memory `Vec<Protein>`.
#[derive(Clone, Copy)]
pub struct ProteinRef<'a> {
    /// The UniProt accession ID of the protein.
    pub uniprot_id: &'a str,
    /// The taxon ID of the protein.
    pub taxon_id: u32,
    /// The encoded functional annotations of the protein.
    pub functional_annotations: &'a [u8]
}

impl<'a> ProteinRef<'a> {
    /// Returns the decoded functional annotations of the protein
    pub fn get_functional_annotations(&self) -> String {
        decode(self.functional_annotations)
    }
}

/// A collection of proteins, either fully in memory or backed by a memory-mapped file.
pub enum Proteins {
    /// All data loaded into memory.
    InMemory {
        text: ProteinText,
        proteins: Vec<Protein>,
    },
    /// Data accessed via memory-mapped file.
    MmapBacked {
        mmap: Arc<Mmap>,
        text: ProteinText,
        protein_count: usize,
        fixed_table_offset: usize,
        uid_data_offset: usize,
        fa_data_offset: usize,
    },
}

impl Proteins {
    /// Creates a new in-memory `Proteins` collection.
    pub fn new(text: ProteinText, proteins: Vec<Protein>) -> Self {
        Proteins::InMemory { text, proteins }
    }

    /// Returns a reference to the underlying `ProteinText`.
    pub fn text(&self) -> &ProteinText {
        match self {
            Proteins::InMemory { text, .. } => text,
            Proteins::MmapBacked { text, .. } => text,
        }
    }

    /// Returns the number of proteins in the collection.
    pub fn len(&self) -> usize {
        match self {
            Proteins::InMemory { proteins, .. } => proteins.len(),
            Proteins::MmapBacked { protein_count, .. } => *protein_count,
        }
    }

    /// Returns true if there are no proteins.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns a zero-copy view of the protein at `index`.
    pub fn get(&self, index: usize) -> ProteinRef<'_> {
        match self {
            Proteins::InMemory { proteins, .. } => {
                let p = &proteins[index];
                ProteinRef {
                    uniprot_id: &p.uniprot_id,
                    taxon_id: p.taxon_id,
                    functional_annotations: &p.functional_annotations
                }
            }
            Proteins::MmapBacked { mmap, fixed_table_offset, uid_data_offset, fa_data_offset, .. } => {
                let entry_off = fixed_table_offset + index * 16;
                let entry = &mmap[entry_off..entry_off + 16];

                let taxon_id = u32::from_le_bytes(entry[0..4].try_into().unwrap());
                let uid_offset = u32::from_le_bytes(entry[4..8].try_into().unwrap()) as usize;
                let uid_len = u16::from_le_bytes(entry[8..10].try_into().unwrap()) as usize;
                let fa_offset = u32::from_le_bytes(entry[10..14].try_into().unwrap()) as usize;
                let fa_len = u16::from_le_bytes(entry[14..16].try_into().unwrap()) as usize;

                ProteinRef {
                    uniprot_id: std::str::from_utf8(
                        &mmap[uid_data_offset + uid_offset..uid_data_offset + uid_offset + uid_len]
                    ).unwrap(),
                    taxon_id,
                    functional_annotations: &mmap[fa_data_offset + fa_offset..fa_data_offset + fa_offset + fa_len],
                }
            }
        }
    }

    /// Creates a new `Proteins` struct from a database file.
    pub fn load_from_tsv(file: &str) -> Result<Self, Box<dyn Error>> {
        let mut input_string: String = String::new();
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

        let text = ProteinText::from_string(&input_string);
        Ok(Proteins::InMemory { text, proteins })
    }

    /// Creates a `ProteinText` which represents all the proteins concatenated from the database file.
    pub fn text_from_tsv(database_file: &str) -> Result<ProteinText, Box<dyn Error>> {
        let input_string = Self::read_sequences_from_tsv(database_file)?;
        let text = ProteinText::from_string(&input_string);
        Ok(text)
    }

    /// Creates a `Vec<u8>` which represents all the proteins concatenated from the database file.
    pub fn raw_text_from_tsv(database_file: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let mut input_string = Self::read_sequences_from_tsv(database_file)?;

        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());

        input_string.shrink_to_fit();
        Ok(input_string.into_bytes())
    }

    /// Reads concatenated sequences from a TSV database file.
    fn read_sequences_from_tsv(database_file: &str) -> Result<String, Box<dyn Error>> {
        let mut input_string: String = String::new();

        let file = File::open(database_file)?;

        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            let sequence = from_utf8(fields.nth(2).unwrap())?;

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
        }

        Ok(input_string)
    }

}

impl WriteBinary for Proteins {
    /// Writes this `Proteins` to a writer in the binary proteins format.
    ///
    /// Format:
    /// - ProteinText section (text_length u64 + BitArray bytes)
    /// - protein_count (u64 le)
    /// - uid_bytes_total (u64 le)
    /// - fa_bytes_total (u64 le)
    /// - protein_count × 16-byte fixed entries: taxon_id u32, uid_offset u32, uid_len u16,
    ///   fa_offset u32, fa_len u16
    /// - uid bytes (concatenated uniprot_ids)
    /// - fa bytes (concatenated functional_annotations)
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        match self {
            Proteins::InMemory { text, proteins } => {
                WriteBinary::write_binary(text, writer)?;

                let protein_count = proteins.len() as u64;
                let uid_bytes_total: u64 =
                    proteins.iter().map(|p| p.uniprot_id.len() as u64).sum();
                let fa_bytes_total: u64 =
                    proteins.iter().map(|p| p.functional_annotations.len() as u64).sum();

                writer.write_all(&protein_count.to_le_bytes())?;
                writer.write_all(&uid_bytes_total.to_le_bytes())?;
                writer.write_all(&fa_bytes_total.to_le_bytes())?;

                let mut uid_offset: u32 = 0;
                let mut fa_offset: u32 = 0;
                for protein in &proteins {
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

                for protein in &proteins {
                    writer.write_all(protein.uniprot_id.as_bytes())?;
                }

                for protein in proteins {
                    writer.write_all(&protein.functional_annotations)?;
                }

                Ok(())
            }
            Proteins::MmapBacked { .. } => {
                panic!("write_binary() is not supported on MmapBacked Proteins");
            }
        }
    }
}

impl ReadBinary for Proteins {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let text = ProteinText::read_binary(reader)?;

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

        Ok(Proteins::InMemory { text, proteins })
    }
}

impl ReadBinaryMmap for Proteins {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let f = File::open(path)?;
        let mmap = Arc::new(unsafe { Mmap::map(&f)? });

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        let mmap_len = mmap.len();
        if mmap_len < 8 {
            return Err("proteins file too short to contain text header".into());
        }

        // Parse text section header
        let text_length = u64::from_le_bytes(mmap[0..8].try_into()?) as usize;
        let text_data_offset: usize = 8;
        let bit_array_bytes = bit_array_byte_size(text_length);

        // Compute metadata offset with overflow and bounds checks
        let meta_offset = text_data_offset
            .checked_add(bit_array_bytes)
            .ok_or_else(|| "overflow while computing metadata offset".to_string())?;
        let meta_end = meta_offset
            .checked_add(24)
            .ok_or_else(|| "overflow while computing metadata end offset".to_string())?;
        if meta_end > mmap_len {
            return Err("proteins file too short to contain metadata section".into());
        }

        // Build ProteinText::MmapBacked
        let text = ProteinText::from_mmap(Arc::clone(&mmap), text_data_offset, text_length);

        // Parse metadata section
        let protein_count = u64::from_le_bytes(mmap[meta_offset..meta_offset + 8].try_into()?) as usize;
        let uid_bytes_total = u64::from_le_bytes(mmap[meta_offset + 8..meta_offset + 16].try_into()?, ) as usize;

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

        if uid_data_offset > mmap_len || fa_data_offset > mmap_len {
            return Err("proteins file truncated: data section offsets exceed file length".into());
        }

        Ok(Proteins::MmapBacked {
            mmap,
            text,
            protein_count,
            fixed_table_offset,
            uid_data_offset,
            fa_data_offset,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Write, path::PathBuf};

    use tempdir::TempDir;

    use super::*;

    fn create_database_file(tmp_dir: &TempDir) -> PathBuf {
        let database_file = tmp_dir.path().join("database.tsv");
        let mut file = File::create(&database_file).unwrap();

        file.write_all("P12345\t1\tMLPGLALLLLAAWTARALEV\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes())
            .unwrap();
        file.write_all("P54321\t2\tPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes())
            .unwrap();
        file.write_all("P67890\t6\tKWDSDPSGTKTCIDT\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes())
            .unwrap();
        file.write_all(
            "P13579\t17\tKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n"
                .as_bytes()
        )
        .unwrap();

        database_file
    }

    #[test]
    fn test_new_protein() {
        let protein = Protein {
            uniprot_id: "P12345".to_string(),
            taxon_id: 1,
            functional_annotations: vec![0xD1, 0x11]
        };

        assert_eq!(protein.uniprot_id, "P12345");
        assert_eq!(protein.taxon_id, 1);
        assert_eq!(protein.functional_annotations, vec![0xD1, 0x11]);
    }

    #[test]
    fn test_new_proteins() {
        let input_string = "MLPGLALLLLAAWTARALEV-PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG";
        let text = ProteinText::from_string(input_string);
        let proteins = Proteins::new(text, vec![
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
        assert_eq!(proteins.get(0).taxon_id, 1);
        assert_eq!(proteins.get(0).functional_annotations, &[0xD1u8, 0x11u8][..]);
        assert_eq!(proteins.get(1).uniprot_id, "P54321");
        assert_eq!(proteins.get(1).taxon_id, 2);
        assert_eq!(proteins.get(1).functional_annotations, &[0xD1u8, 0x11u8][..]);
    }

    #[test]
    fn test_get_taxon() {
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();

        let taxa = [1u32, 2, 6, 17];
        for (i, &taxon) in taxa.iter().enumerate() {
            assert_eq!(proteins.get(i).taxon_id, taxon);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();

        for i in 0..proteins.len() {
            assert_eq!(
                proteins.get(i).get_functional_annotations(),
                "GO:0009279;IPR:IPR016364;IPR:IPR008816"
            );
        }
    }

    #[test]
    fn test_get_concatenated_proteins() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins =
            Proteins::text_from_tsv(database_file.to_str().unwrap())
                .unwrap();

        let expected = b'L';
        assert_eq!(proteins.get(4), expected);
    }

    #[test]
    fn test_write_and_read_binary_buffered() {
        let tmp_dir = TempDir::new("test_binary_roundtrip").unwrap();
        let database_file = create_database_file(&tmp_dir);

        let original = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();
        let original_save = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();

        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        drop(bin_file);

        let f = File::open(&bin_path).unwrap();
        let mut reader = BufReader::new(f);
        let loaded = Proteins::read_binary(&mut reader).unwrap();

        assert_eq!(loaded.len(), original_save.len());
        for i in 0..original_save.len() {
            let orig = original_save.get(i);
            let load = loaded.get(i);
            assert_eq!(orig.uniprot_id, load.uniprot_id, "uniprot_id mismatch at {}", i);
            assert_eq!(orig.taxon_id, load.taxon_id, "taxon_id mismatch at {}", i);
            assert_eq!(
                orig.functional_annotations,
                load.functional_annotations,
                "fa mismatch at {}",
                i
            );
        }
        assert_eq!(loaded.text().len(), original_save.text().len());
        for i in 0..original_save.text().len() {
            assert_eq!(
                loaded.text().get(i),
                original_save.text().get(i),
                "text mismatch at {}",
                i
            );
        }
    }

    #[test]
    fn test_load_from_binary_mmap() {
        let tmp_dir = TempDir::new("test_mmap_roundtrip").unwrap();
        let database_file = create_database_file(&tmp_dir);

        let original = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();
        let original_save = Proteins::load_from_tsv(database_file.to_str().unwrap()).unwrap();

        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        drop(bin_file);

        let mmap_loaded = Proteins::read_binary_mmap(bin_path.as_path()).unwrap();

        assert_eq!(mmap_loaded.len(), original_save.len());
        for i in 0..original_save.len() {
            let orig = original_save.get(i);
            let mmap = mmap_loaded.get(i);
            assert_eq!(orig.uniprot_id, mmap.uniprot_id, "uniprot_id mismatch at {}", i);
            assert_eq!(orig.taxon_id, mmap.taxon_id, "taxon_id mismatch at {}", i);
            assert_eq!(
                orig.functional_annotations,
                mmap.functional_annotations,
                "fa mismatch at {}",
                i
            );
        }
        assert_eq!(mmap_loaded.text().len(), original_save.text().len());
        for i in 0..original_save.text().len() {
            assert_eq!(
                mmap_loaded.text().get(i),
                original_save.text().get(i),
                "text mismatch at {}",
                i
            );
        }
    }
}
