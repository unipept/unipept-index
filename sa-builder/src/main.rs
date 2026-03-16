use std::{
    fs::{File, OpenOptions},
    io::BufWriter,
    time::{SystemTime, SystemTimeError, UNIX_EPOCH}
};

use clap::Parser;
use sa_builder::{Arguments, build_ssa};
use sa_index::array::dump_compressed_suffix_array;
use sa_index::array::dump_suffix_array;
use sa_index::suffix_to_protein_index::dump_mapping;
use sa_index::WriteBinary;
use sa_mappings::proteins::Proteins;

fn main() {
    let Arguments {
        database_file,
        output_sa,
        output_proteins,
        sparseness_factor,
        construction_algorithm,
        compress_sa,
        output_mapping,
        mapping_style
    } = Arguments::parse();
    eprintln!();
    eprintln!("📋 Started loading the proteins...");
    let start_proteins_time = get_time_ms().unwrap();
    let data = Proteins::raw_text_from_tsv(&database_file)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    let bits_per_value = (data.len() as f64).log2().ceil() as usize;
    eprintln!(
        "✅ Successfully loaded the proteins in {} seconds!",
        (get_time_ms().unwrap() - start_proteins_time) / 1000.0
    );

    eprintln!();
    eprintln!("📋 Started building the suffix array...");
    let start_ssa_time = get_time_ms().unwrap();
    let sa = build_ssa(data, &construction_algorithm, sparseness_factor)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!(
        "✅ Successfully built the suffix array in {} seconds!",
        (get_time_ms().unwrap() - start_ssa_time) / 1000.0
    );
    eprintln!("\tAmount of items: {}", sa.len());

    // open the output file
    let mut file =
        open_file_buffer(&output_sa, 100 * 1024 * 1024).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    eprintln!();
    eprintln!("📋 Started dumping the suffix array...");
    let start_dump_time = get_time_ms().unwrap();

    if compress_sa {
        if let Err(err) = dump_compressed_suffix_array(sa, sparseness_factor, bits_per_value, &mut file) {
            eprint_and_exit(err.to_string().as_str());
        };

        eprintln!(
            "✅ Successfully dumped the suffix array in {} seconds!",
            (get_time_ms().unwrap() - start_dump_time) / 1000.0
        );
        eprintln!("\tAmount of bits per item: {}", bits_per_value);
    } else {
        if let Err(err) = dump_suffix_array(&sa, sparseness_factor, &mut file) {
            eprint_and_exit(err.to_string().as_str());
        }

        eprintln!(
            "✅ Successfully dumped the suffix array in {} seconds!",
            (get_time_ms().unwrap() - start_dump_time) / 1000.0
        );
        eprintln!("\tAmount of bits per item: 64");
    }

    // Write the proteins binary file
    eprintln!();
    eprintln!("Started writing proteins binary...");
    let start_proteins_bin_time = get_time_ms().unwrap();
    let proteins = Proteins::load_from_tsv(&database_file)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    let mut proteins_file = open_file_buffer(&output_proteins, 100 * 1024 * 1024)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    proteins.write_binary(&mut proteins_file).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!(
        "Successfully wrote proteins binary in {} seconds!",
        (get_time_ms().unwrap() - start_proteins_bin_time) / 1000.0
    );
    eprintln!("\tOutput: {}", output_proteins);

    if let Some(mapping_path) = output_mapping {
        eprintln!();
        eprintln!("Started writing suffix-to-protein mapping binary...");
        let start_mapping_time = get_time_ms().unwrap();
        let mut mapping_file = open_file_buffer(&mapping_path, 100 * 1024 * 1024)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
        dump_mapping(&mapping_style, proteins.text(), &mut mapping_file)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
        eprintln!(
            "Successfully wrote mapping binary in {} seconds!",
            (get_time_ms().unwrap() - start_mapping_time) / 1000.0
        );
        eprintln!("\tOutput: {}", mapping_path);
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

pub fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}
