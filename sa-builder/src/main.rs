use std::{
    fs::{
        File,
        OpenOptions
    },
    io::{BufWriter, Result}
};

use clap::Parser;
use sa_builder::{
    build_ssa,
    Arguments
};
use sa_compression::dump_compressed_suffix_array;
use sa_index::binary::dump_suffix_array;
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

    eprintln!();
    eprintln!("📋 Started loading the taxon file...");
    let taxon_id_calculator =
        TaxonAggregator::try_from_taxonomy_file(&taxonomy, AggregationMethod::LcaStar)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!("✅ Successfully loaded the taxon file!");
    eprintln!("\tAggregation method: LCA*");

    eprintln!();
    eprintln!("📋 Started loading the proteins...");
    let mut data =
        Proteins::try_from_database_file_without_annotations(&database_file, &taxon_id_calculator)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!("✅ Successfully loaded the proteins!");

    eprintln!();
    eprintln!("📋 Started building the suffix array...");
    let sa = build_ssa(&mut data, &construction_algorithm, sparseness_factor)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!("✅ Successfully built the suffix array!");
    eprintln!("\tAmount of items: {}", sa.len());
    eprintln!("\tSample rate: {}", sparseness_factor);

    // open the output file
    let mut file =
        open_file_buffer(&output, 100 * 1024 * 1024).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    eprintln!();
    eprintln!("📋 Started dumping the suffix array...");

    if compress_sa {
        let bits_per_value = (data.len() as f64).log2().ceil() as usize;

        if let Err(err) =
            dump_compressed_suffix_array(sa, sparseness_factor, bits_per_value, &mut file)
        {
            eprint_and_exit(err.to_string().as_str());
        };

        eprintln!("✅ Successfully dumped the suffix array!");
        eprintln!("\tAmount of bits per item: {}", bits_per_value);
    } else {
        if let Err(err) = dump_suffix_array(&sa, sparseness_factor, &mut file) {
            eprint_and_exit(err.to_string().as_str());
        }

        eprintln!("✅ Successfully dumped the suffix array!");
        eprintln!("\tAmount of bits per item: 64");
    }
}

fn open_file_buffer(file: &str, buffer_size: usize) -> Result<BufWriter<File>> {
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
