use std::error::Error;
use std::fs::File;
use std::io::{BufReader, Read};
use sa_compression::load_compressed_suffix_array;
use sa_index::binary::{load_suffix_array, load_suffix_array_mmap};
use sa_index::suffix_to_protein_index::{SuffixToProteinIndex, load_mapping};
use sa_index::SuffixArray;
use sa_mappings::proteins::Proteins;

pub fn load_suffix_array_file(file: &str, use_mmap: bool) -> Result<SuffixArray, Box<dyn Error>> {
    if use_mmap {
        return load_suffix_array_mmap(std::path::Path::new(file));
    }

    let mut sa_file = File::open(file)?;
    let mut reader = BufReader::new(&mut sa_file);

    let mut bits_per_value_buffer = [0_u8; 1];
    reader
        .read_exact(&mut bits_per_value_buffer)
        .map_err(|_| "Could not read the flags from the binary file")?;
    let bits_per_value = bits_per_value_buffer[0];

    if bits_per_value == 64 {
        load_suffix_array(&mut reader)
    } else {
        load_compressed_suffix_array(&mut reader, bits_per_value as usize)
    }
}

pub fn load_proteins_file(file: &str) -> Result<Proteins, Box<dyn Error>> {
    Proteins::load_from_binary(file)
}

pub fn load_mapping_file(file: &str) -> Result<Box<dyn SuffixToProteinIndex>, Box<dyn Error>> {
    let f = File::open(file)?;
    let mut reader = BufReader::new(f);
    load_mapping(&mut reader)
}
