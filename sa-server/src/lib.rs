use std::{error::Error, fs::File, io::BufReader};

use sa_index::{ReadBinary, ReadBinaryMmap, SuffixArray, suffix_to_protein_index::legacy::SuffixToProteinMapping};
use sa_mappings::proteins::Proteins;

pub fn load_suffix_array_file(file: &str, use_mmap: bool) -> Result<SuffixArray, Box<dyn Error>> {
    if use_mmap {
        return SuffixArray::read_binary_mmap(std::path::Path::new(file));
    }
    let f = File::open(file)?;
    let mut reader = BufReader::new(f);
    SuffixArray::read_binary(&mut reader)
}

pub fn load_proteins_file(file: &str, use_mmap: bool) -> Result<Proteins, Box<dyn Error>> {
    if use_mmap {
        return Proteins::read_binary_mmap(std::path::Path::new(file));
    }

    let proteins_file = File::open(file)?;
    let mut reader = BufReader::new(proteins_file);
    Proteins::read_binary(&mut reader)
}

pub fn load_mapping_file(file: &str, use_mmap: bool) -> Result<SuffixToProteinMapping, Box<dyn Error>> {
    if use_mmap {
        return SuffixToProteinMapping::read_binary_mmap(std::path::Path::new(file));
    }
    let f = File::open(file)?;
    let mut reader = BufReader::new(f);
    SuffixToProteinMapping::read_binary(&mut reader)
}
