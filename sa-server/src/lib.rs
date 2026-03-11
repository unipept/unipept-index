use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read};
use sa_compression::load_compressed_suffix_array;
use sa_index::binary::{load_suffix_array, load_suffix_array_mmap};
use sa_index::SuffixArray;
use sa_mappings::proteins::Proteins;

pub fn load_suffix_array_file(file: &str, use_mmap: bool) -> Result<SuffixArray, Box<dyn Error>> {
    if use_mmap {
        load_suffix_array_mmap(std::path::Path::new(file))?;
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

pub fn load_proteins_file(file: &str, use_mmap: bool) -> Result<Proteins, Box<dyn Error>> {
    if use_mmap {
        Proteins::try_from_binary_file(&file, true)
    } else {
        let proteins_filename = "/mnt/data/uniprot-2025-04/suffix-array/proteins.bin";
        let proteins = Proteins::try_from_database_file(&file)?;
        let mut proteins_file = open_file_buffer(&proteins_filename, 100 * 1024 * 1024)?;
        proteins.write_binary(&mut proteins_file)?;
        Ok(proteins)
    }
}

fn open_file_buffer(file: &str, buffer_size: usize) -> std::io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true) // if the file already exists, empty the file
        .open(file)?;

    Ok(BufWriter::with_capacity(buffer_size, file))
}
