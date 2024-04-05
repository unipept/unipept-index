use std::{error::Error, fs::File, io::{BufRead, BufReader}, ops::Index};

use memchr::memchr_iter;
use umgap::taxon::TaxonId;

use crate::{taxonomy::TaxonAggregator, DatabaseFormatError};

pub static SEPARATION_CHARACTER: u8 = b'-';
pub static TERMINATION_CHARACTER: u8 = b'$';

#[derive(Debug)]
pub struct Protein {
    /// The id of the protein
    pub uniprot_id: String,

    /// start position and length of the protein in the input string
    pub sequence: (usize, u32),

    /// the taxon id of the protein
    pub taxon_id: TaxonId,

    // /// The encoded functional annotations of the protein
    functional_annotations: Vec<u8>,
}

#[derive(Debug)]
pub struct Proteins {
    /// The input string containing all proteins
    input_string: Vec<u8>,

    /// The proteins in the input string
    proteins: Vec<Protein>,
}

impl Proteins {
    pub fn try_from_database_file(file: &str, taxon_aggregator: &TaxonAggregator) -> Result<Self, Box<dyn Error>> {
        let mut input_string: String = String::new();
        let mut proteins: Vec<Protein> = Vec::new();

        let file = File::open(file)?;

        let mut start_index = 0;

        let mut reader = BufReader::new(file);

        let mut buffer = Vec::new();
        println!("{:?}", reader.read_until(b'\n', &mut buffer));

        println!("{:?}", buffer);

        for line in reader.lines().into_iter().map_while(Result::ok) {
            println!("{:?}", line);
            let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            let [uniprot_id, taxon_id, sequence, fa]: [String; 4] = fields.try_into().map_err(DatabaseFormatError::new)?;
            println!("{:?}", taxon_id);
            let taxon_id = taxon_id.parse::<TaxonId>()?;

            if !taxon_aggregator.taxon_exists(taxon_id) {
                continue;
            }

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());

            proteins.push(Protein {
                uniprot_id,
                sequence: (start_index, sequence.len() as u32),
                taxon_id,
                functional_annotations: fa.as_bytes().to_vec(),
            });

            start_index += sequence.len() + 1;
        }

        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());

        Ok(Self { input_string: input_string.into_bytes(), proteins })
    }

    pub fn get_sequence(&self, protein: &Protein) -> &str {
        let (start, length) = protein.sequence;
        let end = start + length as usize;

        std::str::from_utf8(&self.input_string[start..end]).unwrap() // should never fail since the input string will always be utf8
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
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;

    use fa_compression::decode;
    use tempdir::TempDir;

    use crate::taxonomy::AggregationMethod;

    use super::*;

    fn create_database_file(tmp_dir: &TempDir) -> PathBuf {
        let database_file = tmp_dir.path().join("database.tsv");
        let mut file = File::create(&database_file).unwrap();

        file.write("P12345\t1\tMLPGLALLLLAAWTARALEV\t".as_bytes()).unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27]).unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P54321\t2\tPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG\t".as_bytes()).unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27]).unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P67890\t6\tKWDSDPSGTKTCIDT\t".as_bytes()).unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27]).unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P13579\t17\tKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH\t".as_bytes()).unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27]).unwrap();
        file.write("\n".as_bytes()).unwrap();

        database_file
    }

    fn create_taxonomy_file(tmp_dir: &TempDir) -> PathBuf {
        let taxonomy_file = tmp_dir.path().join("taxonomy.tsv");
        let mut file = File::create(&taxonomy_file).unwrap();

        writeln!(file, "1\troot\tno rank\t1\t\x01").unwrap();
        writeln!(file, "2\tBacteria\tsuperkingdom\t1\t\x01").unwrap();
        writeln!(file, "6\tAzorhizobium\tgenus\t1\t\x01").unwrap();
        writeln!(file, "7\tAzorhizobium caulinodans\tspecies\t6\t\x01").unwrap();
        writeln!(file, "9\tBuchnera aphidicola\tspecies\t6\t\x01").unwrap();
        writeln!(file, "10\tCellvibrio\tgenus\t6\t\x01").unwrap();
        writeln!(file, "11\tCellulomonas gilvus\tspecies\t10\t\x01").unwrap();
        writeln!(file, "13\tDictyoglomus\tgenus\t11\t\x01").unwrap();
        writeln!(file, "14\tDictyoglomus thermophilum\tspecies\t10\t\x01").unwrap();
        writeln!(file, "16\tMethylophilus\tgenus\t14\t\x01").unwrap();
        writeln!(file, "17\tMethylophilus methylotrophus\tspecies\t16\t\x01").unwrap();
        writeln!(file, "18\tPelobacter\tgenus\t17\t\x01").unwrap();
        writeln!(file, "19\tSyntrophotalea carbinolica\tspecies\t17\t\x01").unwrap();
        writeln!(file, "20\tPhenylobacterium\tgenus\t19\t\x01").unwrap();

        taxonomy_file
    }

    #[test]
    fn test_get_sequence() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_sequences").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(taxonomy_file.to_str().unwrap(), AggregationMethod::Lca).unwrap();
        let proteins = Proteins::try_from_database_file(database_file.to_str().unwrap(), &taxon_aggregator).unwrap();

        println!("{:?}", proteins);

        //assert_eq!(proteins.proteins.len(), 4);
        assert_eq!(proteins.get_sequence(&proteins[0]), "MLPGLALLLLAAWTARALEV");
        assert_eq!(proteins.get_sequence(&proteins[1]), "PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG");
        assert_eq!(proteins.get_sequence(&proteins[2]), "KWDSDPSGTKTCIDT");
        assert_eq!(proteins.get_sequence(&proteins[3]), "KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH");
    }

    #[test]
    fn test_get_taxon() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(taxonomy_file.to_str().unwrap(), AggregationMethod::Lca).unwrap();
        let proteins = Proteins::try_from_database_file(database_file.to_str().unwrap(), &taxon_aggregator).unwrap();

        let taxa = vec![1, 2, 6, 17];
        for (i, protein) in proteins.proteins.iter().enumerate() {
            assert_eq!(protein.taxon_id, taxa[i]);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(taxonomy_file.to_str().unwrap(), AggregationMethod::Lca).unwrap();
        let proteins = Proteins::try_from_database_file(database_file.to_str().unwrap(), &taxon_aggregator).unwrap();

        for protein in proteins.proteins.iter() {
            println!("{:?}", protein.functional_annotations);
            assert_eq!(decode(&protein.functional_annotations), "GO:0009279;IPR:IPR016364;IPR:IPR008816");
        }
    }
}
