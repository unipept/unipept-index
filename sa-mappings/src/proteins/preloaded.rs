use std::{error::Error, fs::File};
use std::io::{BufRead, BufReader, Write};
use std::str::from_utf8;

use bytelines::ByteLines;
use fa_compression::algorithm1::encode;
pub use text_compression::{WriteBinary, ReadBinary};
use text_compression::InMemoryProteinText;

use super::{Protein, ProteinRef, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use super::ProteinsBackend;

// ── InMemoryProteins ──────────────────────────────────────────────────────────

pub struct InMemoryProteins {
    pub text: InMemoryProteinText,
    pub proteins: Vec<Protein>,
}

impl InMemoryProteins {
    pub fn new(text: InMemoryProteinText, proteins: Vec<Protein>) -> Self {
        Self { text, proteins }
    }

    // ── TSV loaders (non-mmap only) ───────────────────────────────────────────

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
            proteins.push(Protein { uniprot_id: uniprot_id.to_string(), taxon_id, functional_annotations });
        }
        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());
        proteins.shrink_to_fit();
        let text = InMemoryProteinText::from_string(&input_string);
        Ok(Self { text, proteins })
    }

    pub fn text_from_tsv(database_file: &str) -> Result<InMemoryProteinText, Box<dyn Error>> {
        Ok(InMemoryProteinText::from_string(&Self::read_sequences_from_tsv(database_file)?))
    }

    pub fn raw_text_from_tsv(database_file: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let mut s = Self::read_sequences_from_tsv(database_file)?;
        s.pop();
        s.push(TERMINATION_CHARACTER.into());
        s.shrink_to_fit();
        Ok(s.into_bytes())
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
    fn text(&self) -> &InMemoryProteinText { &self.text }
    fn len(&self) -> usize { self.proteins.len() }

    fn get(&self, index: usize) -> ProteinRef<'_> {
        let p = &self.proteins[index];
        ProteinRef { uniprot_id: &p.uniprot_id, taxon_id: p.taxon_id, functional_annotations: &p.functional_annotations }
    }

    #[inline]
    fn prefetch(&self, index: usize) {
        if index < self.proteins.len() {
            prefetch::prefetch_read(&self.proteins[index] as *const Protein);
        }
    }
}

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
        for protein in &self.proteins { writer.write_all(protein.uniprot_id.as_bytes())?; }
        for protein in self.proteins { writer.write_all(&protein.functional_annotations)?; }
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
        for entry in table.iter_mut() { reader.read_exact(entry)?; }
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
    use std::{fs::File, io::Write, path::PathBuf};
    use tempdir::TempDir;
    use super::*;
    use text_compression::ProteinTextBackend as _;

    fn create_database_file(tmp_dir: &TempDir) -> PathBuf {
        let path = tmp_dir.path().join("database.tsv");
        let mut f = File::create(&path).unwrap();
        f.write_all("P12345\t1\tMLPGLALLLLAAWTARALEV\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        f.write_all("P54321\t2\tPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        f.write_all("P67890\t6\tKWDSDPSGTKTCIDT\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        f.write_all("P13579\t17\tKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH\tGO:0009279;IPR:IPR016364;IPR:IPR008816\n".as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_new_protein() {
        let protein = Protein { uniprot_id: "P12345".to_string(), taxon_id: 1, functional_annotations: vec![0xD1, 0x11] };
        assert_eq!(protein.uniprot_id, "P12345");
        assert_eq!(protein.taxon_id, 1);
    }

    #[test]
    fn test_new_proteins() {
        let text = InMemoryProteinText::from_string("MLPGLALLLLAAWTARALEV-PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG");
        let proteins = InMemoryProteins::new(text, vec![
            Protein { uniprot_id: "P12345".to_string(), taxon_id: 1, functional_annotations: vec![0xD1, 0x11] },
            Protein { uniprot_id: "P54321".to_string(), taxon_id: 2, functional_annotations: vec![0xD1, 0x11] },
        ]);
        assert_eq!(proteins.len(), 2);
        assert_eq!(proteins.get(0).uniprot_id, "P12345");
        assert_eq!(proteins.get(1).taxon_id, 2);
    }

    #[test]
    fn test_get_taxon() {
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();
        let proteins = InMemoryProteins::load_from_tsv(create_database_file(&tmp_dir).to_str().unwrap()).unwrap();
        for (i, &taxon) in [1u32, 2, 6, 17].iter().enumerate() {
            assert_eq!(proteins.get(i).taxon_id, taxon);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();
        let proteins = InMemoryProteins::load_from_tsv(create_database_file(&tmp_dir).to_str().unwrap()).unwrap();
        for i in 0..proteins.len() {
            assert_eq!(proteins.get(i).get_functional_annotations(), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }

    #[test]
    fn test_get_concatenated_proteins() {
        let tmp_dir = TempDir::new("test_get_fa").unwrap();
        let text = InMemoryProteins::text_from_tsv(create_database_file(&tmp_dir).to_str().unwrap()).unwrap();
        assert_eq!(text.get(4), b'L');
    }

    #[test]
    fn test_write_and_read_binary_buffered() {
        use std::io::BufReader;
        let tmp_dir = TempDir::new("test_binary_roundtrip").unwrap();
        let db = create_database_file(&tmp_dir);
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
