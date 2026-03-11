//! This module contains the `Protein` and `Proteins` structs, which are used to represent proteins
//! and collections of proteins, respectively.

use std::{error::Error, fs::File, io::{BufReader, Read, Write}, ops::Index, path::Path, str::from_utf8};

use bytelines::ByteLines;
use fa_compression::algorithm1::{decode, encode};
use memmap2::Mmap;
use text_compression::ProteinText;

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

/// A struct that represents a collection of proteins
pub struct Proteins {
    /// The input string containing all proteins
    pub text: ProteinText,

    /// The proteins in the input string
    pub proteins: Vec<Protein>
}

impl Protein {
    /// Returns the decoded functional annotations of the protein
    pub fn get_functional_annotations(&self) -> String {
        decode(&self.functional_annotations)
    }
}

impl Proteins {
    /// Creates a new `Proteins` struct from a database file and a `TaxonAggregator`
    ///
    /// # Arguments
    /// * `file` - The path to the database file
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `Proteins` struct
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the database file
    pub fn try_from_database_file(file: &str) -> Result<Self, Box<dyn Error>> {
        let mut input_string: String = String::new();
        let mut proteins: Vec<Protein> = Vec::new();

        let file = File::open(file)?;

        // Read the lines as bytes, since the input string is not guaranteed to be utf8
        // because of the encoded functional annotations
        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            // uniprot_id, taxon_id and sequence should always contain valid utf8
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
        Ok(Self { text, proteins })
    }

    /// Creates a `ProteinText` which represents all the proteins concatenated from the database file
    ///
    /// # Arguments
    /// * `file` - The path to the database file
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `ProteinText`
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the database file
    pub fn try_from_database_file_without_annotations(database_file: &str) -> Result<ProteinText, Box<dyn Error>> {
        let mut input_string: String = String::new();

        let file = File::open(database_file)?;

        // Read the lines as bytes, since the input string is not guaranteed to be utf8
        // because of the encoded functional annotations
        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            // only get the taxon id and sequence from each line, we don't need the other parts
            let sequence = from_utf8(fields.nth(2).unwrap())?;

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
        }

        let text = ProteinText::from_string(&input_string);

        Ok(text)
    }

    /// Creates a `vec<u8>` which represents all the proteins concatenated from the database file
    ///
    /// # Arguments
    /// * `file` - The path to the database file
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `Vec<u8>`
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the database file
    pub fn try_from_database_file_uncompressed(database_file: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let mut input_string: String = String::new();

        let file = File::open(database_file)?;

        // Read the lines as bytes, since the input string is not guaranteed to be utf8
        // because of the encoded functional annotations
        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            // only get the taxon id and sequence from each line, we don't need the other parts
            let sequence = from_utf8(fields.nth(2).unwrap())?;

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
        }

        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());

        input_string.shrink_to_fit();
        Ok(input_string.into_bytes())
    }
}

impl Proteins {
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
    pub fn write_binary(&self, writer: &mut impl Write) -> Result<(), Box<dyn Error>> {
        self.text.write_binary(writer)?;

        let protein_count = self.proteins.len() as u64;
        let uid_bytes_total: u64 = self.proteins.iter().map(|p| p.uniprot_id.len() as u64).sum();
        let fa_bytes_total: u64 = self.proteins.iter().map(|p| p.functional_annotations.len() as u64).sum();

        writer.write_all(&protein_count.to_le_bytes())?;
        writer.write_all(&uid_bytes_total.to_le_bytes())?;
        writer.write_all(&fa_bytes_total.to_le_bytes())?;

        // Write fixed table
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

        // Write uid data
        for protein in &self.proteins {
            writer.write_all(protein.uniprot_id.as_bytes())?;
        }

        // Write fa data
        for protein in &self.proteins {
            writer.write_all(&protein.functional_annotations)?;
        }

        Ok(())
    }

    /// Loads a `Proteins` from a binary file produced by `write_binary`.
    ///
    /// If `mmap` is true, the ProteinText is backed by a memory-mapped file (zero-copy);
    /// otherwise it is read into heap memory.
    pub fn try_from_binary_file(file: &str, mmap: bool) -> Result<Self, Box<dyn Error>> {
        if mmap {
            Self::load_binary_mmap(file)
        } else {
            Self::load_binary_buffered(file)
        }
    }

    fn load_binary_buffered(file: &str) -> Result<Self, Box<dyn Error>> {
        let f = File::open(file)?;
        let mut reader = BufReader::new(f);

        let text = ProteinText::read_binary(&mut reader)?;
        let proteins = Self::read_protein_metadata(&mut reader)?;

        Ok(Self { text, proteins })
    }

    fn load_binary_mmap(file: &str) -> Result<Self, Box<dyn Error>> {
        let f = File::open(Path::new(file))?;
        let mmap = unsafe { Mmap::map(&f)? };

        let (text, meta_offset) = ProteinText::from_mmap(mmap)?;

        // Read protein metadata directly from the mmap bytes
        let mmap_ref = match &text {
            ProteinText::MmapBacked { mmap, .. } => mmap,
            _ => unreachable!()
        };
        let proteins = Self::read_protein_metadata_from_bytes(&mmap_ref[meta_offset..])?;

        Ok(Self { text, proteins })
    }

    /// Reads the protein metadata section from a reader.
    fn read_protein_metadata(reader: &mut impl Read) -> Result<Vec<Protein>, Box<dyn Error>> {
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf8)?;
        let protein_count = u64::from_le_bytes(buf8) as usize;

        reader.read_exact(&mut buf8)?;
        let uid_bytes_total = u64::from_le_bytes(buf8) as usize;

        reader.read_exact(&mut buf8)?;
        let fa_bytes_total = u64::from_le_bytes(buf8) as usize;

        // Read fixed table
        let mut table = vec![[0u8; 16]; protein_count];
        for entry in table.iter_mut() {
            reader.read_exact(entry)?;
        }

        // Read uid and fa byte pools
        let mut uid_data = vec![0u8; uid_bytes_total];
        reader.read_exact(&mut uid_data)?;

        let mut fa_data = vec![0u8; fa_bytes_total];
        reader.read_exact(&mut fa_data)?;

        Self::reconstruct_proteins(&table, &uid_data, &fa_data)
    }

    /// Reads the protein metadata section from a byte slice (e.g. from an mmap).
    fn read_protein_metadata_from_bytes(bytes: &[u8]) -> Result<Vec<Protein>, Box<dyn Error>> {
        let mut cursor = std::io::Cursor::new(bytes);
        Self::read_protein_metadata(&mut cursor)
    }

    fn reconstruct_proteins(
        table: &[[u8; 16]],
        uid_data: &[u8],
        fa_data: &[u8]
    ) -> Result<Vec<Protein>, Box<dyn Error>> {
        let mut proteins = Vec::with_capacity(table.len());

        for entry in table {
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
}

impl Index<usize> for Proteins {
    type Output = Protein;

    fn index(&self, index: usize) -> &Self::Output {
        &self.proteins[index]
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
        let proteins = Proteins {
            text,
            proteins: vec![
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
            ]
        };

        assert_eq!(proteins.proteins.len(), 2);
        assert_eq!(proteins[0].uniprot_id, "P12345");
        assert_eq!(proteins[0].taxon_id, 1);
        assert_eq!(proteins[0].functional_annotations, vec![0xD1, 0x11]);
        assert_eq!(proteins[1].uniprot_id, "P54321");
        assert_eq!(proteins[1].taxon_id, 2);
        assert_eq!(proteins[1].functional_annotations, vec![0xD1, 0x11]);
    }

    #[test]
    fn test_get_taxon() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins = Proteins::try_from_database_file(database_file.to_str().unwrap()).unwrap();

        let taxa = [1, 2, 6, 17];
        for (i, protein) in proteins.proteins.iter().enumerate() {
            assert_eq!(protein.taxon_id, taxa[i]);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins = Proteins::try_from_database_file(database_file.to_str().unwrap()).unwrap();

        for protein in proteins.proteins.iter() {
            assert_eq!(protein.get_functional_annotations(), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }

    #[test]
    fn test_get_concatenated_proteins() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);

        let proteins = Proteins::try_from_database_file_without_annotations(database_file.to_str().unwrap()).unwrap();

        let expected = b'L';
        assert_eq!(proteins.get(4), expected);
    }

    #[test]
    fn test_write_and_read_binary_buffered() {
        let tmp_dir = TempDir::new("test_binary_roundtrip").unwrap();
        let database_file = create_database_file(&tmp_dir);

        let original = Proteins::try_from_database_file(database_file.to_str().unwrap()).unwrap();

        // Write to binary file
        let bin_path = tmp_dir.path().join("proteins.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        drop(bin_file);

        // Load back (non-mmap)
        let loaded = Proteins::try_from_binary_file(bin_path.to_str().unwrap(), false).unwrap();

        assert_eq!(loaded.proteins.len(), original.proteins.len());
        for (i, (orig, load)) in original.proteins.iter().zip(loaded.proteins.iter()).enumerate() {
            assert_eq!(orig.uniprot_id, load.uniprot_id, "uniprot_id mismatch at {}", i);
            assert_eq!(orig.taxon_id, load.taxon_id, "taxon_id mismatch at {}", i);
            assert_eq!(orig.functional_annotations, load.functional_annotations, "fa mismatch at {}", i);
        }
        assert_eq!(loaded.text.len(), original.text.len());
        for i in 0..original.text.len() {
            assert_eq!(loaded.text.get(i), original.text.get(i), "text mismatch at {}", i);
        }
    }

    #[test]
    fn test_write_and_read_binary_mmap() {
        let tmp_dir = TempDir::new("test_binary_mmap_roundtrip").unwrap();
        let database_file = create_database_file(&tmp_dir);

        let original = Proteins::try_from_database_file(database_file.to_str().unwrap()).unwrap();

        // Write to binary file
        let bin_path = tmp_dir.path().join("proteins_mmap.bin");
        let mut bin_file = File::create(&bin_path).unwrap();
        original.write_binary(&mut bin_file).unwrap();
        drop(bin_file);

        // Load back (mmap)
        let loaded = Proteins::try_from_binary_file(bin_path.to_str().unwrap(), true).unwrap();

        assert_eq!(loaded.proteins.len(), original.proteins.len());
        for (i, (orig, load)) in original.proteins.iter().zip(loaded.proteins.iter()).enumerate() {
            assert_eq!(orig.uniprot_id, load.uniprot_id, "uniprot_id mismatch at {}", i);
            assert_eq!(orig.taxon_id, load.taxon_id, "taxon_id mismatch at {}", i);
            assert_eq!(orig.functional_annotations, load.functional_annotations, "fa mismatch at {}", i);
        }
        assert_eq!(loaded.text.len(), original.text.len());
        for i in 0..original.text.len() {
            assert_eq!(loaded.text.get(i), original.text.get(i), "text mismatch at {}", i);
        }
    }
}
