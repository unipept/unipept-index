use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use sa_index::{KmerTable, ReadBinary, SuffixArray};
#[cfg(feature = "mmap")]
use sa_index::ReadBinaryMmap;
use sa_index::suffix_to_protein_index::SuffixToProteinMapping;
use sa_mappings::proteins::Proteins;

pub fn load_suffix_array_file(file: &str) -> Result<SuffixArray, Box<dyn Error>> {
    #[cfg(feature = "mmap")]
    return SuffixArray::read_binary_mmap(std::path::Path::new(file));

    #[cfg(not(feature = "mmap"))]
    {
        let f = File::open(file)?;
        let mut reader = BufReader::new(f);
        SuffixArray::read_binary(&mut reader)
    }
}

pub fn load_proteins_file(file: &str) -> Result<Proteins, Box<dyn Error>> {
    #[cfg(feature = "mmap")]
    return Proteins::read_binary_mmap(std::path::Path::new(file));

    #[cfg(not(feature = "mmap"))]
    {
        let f = File::open(file)?;
        let mut reader = BufReader::new(f);
        Proteins::read_binary(&mut reader)
    }
}

pub fn load_mapping_file(file: &str) -> Result<SuffixToProteinMapping, Box<dyn Error>> {
    #[cfg(feature = "mmap")]
    return SuffixToProteinMapping::read_binary_mmap(std::path::Path::new(file));

    #[cfg(not(feature = "mmap"))]
    {
        let f = File::open(file)?;
        let mut reader = BufReader::new(f);
        SuffixToProteinMapping::read_binary(&mut reader)
    }
}

pub fn load_kmer_table_file(file: &str) -> Result<KmerTable, Box<dyn Error>> {
    let f = File::open(file)?;
    let mut reader = BufReader::new(f);
    KmerTable::read_binary(&mut reader)
}
