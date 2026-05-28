use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    time::{SystemTime, SystemTimeError, UNIX_EPOCH}
};

use clap::Parser;
use sa_builder::{Arguments, SuffixToProteinMappingStyle, build_ssa};
use sa_index::array::dump_compressed_suffix_array;
use sa_index::array::dump_suffix_array;
use sa_index::suffix_to_protein_index::{DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein};
use sa_index::{KmerTable, WriteBinary};
use sa_mappings::proteins::{InMemoryProteins, ProteinsBackend as _};
use text_compression::ProteinTextBackend as _;

fn main() {
    let Arguments {
        database_file,
        output_sa,
        output_proteins,
        sparseness_factor,
        construction_algorithm,
        compress_sa,
        output_mapping,
        mapping_style,
        output_kmer_table,
        kmer_size,
    } = Arguments::parse();

    let proteins = timed("loading the proteins", || {
        InMemoryProteins::load_from_tsv(&database_file)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()))
    });
    let bits_per_value = (proteins.text().len() as f64).log2().ceil() as usize;

    // Not the most efficient way to get the text, but still acceptable
    let text: Vec<u8> = proteins.text().iter().collect();
    let sa = timed("building the suffix array", || {
        build_ssa(text, &construction_algorithm, sparseness_factor)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()))
    });
    eprintln!("\tAmount of items: {}", sa.len());

    // Build the k-mer table while the SA Vec is still in memory (before it is consumed below).
    let kmer_table: Option<KmerTable> = output_kmer_table.as_ref().map(|_| {
        timed(&format!("building k-mer table (k={})", kmer_size), || {
            KmerTable::build_from_raw_sa(&sa, proteins.text().len(), |i| proteins.text().get(i), kmer_size)
        })
    });

    // open the output file
    let mut sa_file =
        open_file_buffer(&output_sa, 100 * 1024 * 1024).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("dumping the suffix array", || {
        eprintln!("\tAmount of bits per item: {}", if compress_sa { bits_per_value } else { 64 });
        if compress_sa {
            if let Err(err) = dump_compressed_suffix_array(sa, sparseness_factor, bits_per_value, &mut sa_file) {
                eprint_and_exit(err.to_string().as_str());
            };
        } else {
            if let Err(err) = dump_suffix_array(sa, sparseness_factor, &mut sa_file) {
                eprint_and_exit(err.to_string().as_str());
            }
        }
    });

    let mut mapping_file = open_file_buffer(&output_mapping, 100 * 1024 * 1024)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("writing suffix-to-protein mapping binary", || {
        let text_len = proteins.text().len();
        let result = match &mapping_style {
            SuffixToProteinMappingStyle::Dense =>
                DenseSuffixToProtein::from_text_parts(text_len, |i: usize| proteins.text().get(i)).write_binary(&mut mapping_file),
            SuffixToProteinMappingStyle::Sparse =>
                SparseSuffixToProtein::from_text_parts(text_len, |i: usize| proteins.text().get(i)).write_binary(&mut mapping_file),
            SuffixToProteinMappingStyle::BitVec =>
                BitVecSuffixToProtein::from_text_parts(text_len, |i: usize| proteins.text().get(i)).write_binary(&mut mapping_file),
        };
        result.unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
        mapping_file.flush().unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    });
    eprintln!("\tOutput: {}", output_mapping);

    let mut proteins_file = open_file_buffer(&output_proteins, 100 * 1024 * 1024)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("writing proteins binary", || {
        proteins.write_binary(&mut proteins_file).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    });
    eprintln!("\tOutput: {}", output_proteins);

    if let (Some(table), Some(kmer_table_path)) = (kmer_table, output_kmer_table) {
        let mut kmer_file = open_file_buffer(&kmer_table_path, 64 * 1024 * 1024)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

        timed("writing k-mer table", || {
            table.write_binary(&mut kmer_file)
                .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
        });
        eprintln!("\tOutput: {}", kmer_table_path);
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

fn eprint_and_exit(err: &str) -> ! {
    eprintln!("{}", err);
    std::process::exit(1);
}

fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}

fn timed<T, F: FnOnce() -> T>(msg: &str, f: F) -> T {
    eprintln!();
    eprintln!("📋 Started {}...", msg);
    let start = get_time_ms().unwrap();
    let result = f();
    eprintln!("✅ Successfully {} in {} seconds!", msg, (get_time_ms().unwrap() - start) / 1000.0);
    result
}
