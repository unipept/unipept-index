use std::{fs::{File, OpenOptions}, io::Result};

use clap::Parser;
use sa_builder::{
    build_ssa,
    Arguments
};
use sa_index::binary::dump_suffix_array;
use sa_compression::dump_compressed_suffix_array;
use sa_mappings::{
    proteins::Proteins,
    taxonomy::{
        AggregationMethod,
        TaxonAggregator
    }
};

fn main() {
    let Arguments {
        database_file,
        taxonomy,
        output,
        sparseness_factor,
        construction_algorithm,
        compress_sa
    } = Arguments::parse();

    let taxon_id_calculator = TaxonAggregator::try_from_taxonomy_file(&taxonomy, AggregationMethod::LcaStar).unwrap_or_else(
        |err| eprint_and_exit(err.to_string().as_str())
    );

    // read input
    let mut data = Proteins::try_from_database_file_without_annotations(&database_file, &taxon_id_calculator).unwrap_or_else(
        |err| eprint_and_exit(err.to_string().as_str())
    );

    // calculate sparse suffix array
    let sa = build_ssa(&mut data, &construction_algorithm, sparseness_factor).unwrap_or_else(
        |err| eprint_and_exit(err.to_string().as_str())
    );

    // open the output file
    let mut file = open_file(&output).unwrap_or_else(
        |err| eprint_and_exit(err.to_string().as_str())
    );

    if compress_sa {
        if let Err(err) = dump_compressed_suffix_array::<37>(sa, sparseness_factor, &mut file) {
            eprint_and_exit(err.to_string().as_str());
        };
    } else {
        if let Err(err) = dump_suffix_array(&sa, sparseness_factor, &mut file) {
            eprint_and_exit(err.to_string().as_str());
        };
    }
}

fn open_file(file: &str) -> Result<File> {
    OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true) // if the file already exists, empty the file
        .open(file)
}

fn eprint_and_exit(err: &str) -> ! {
    eprintln!("{}", err);
    std::process::exit(1);
}
