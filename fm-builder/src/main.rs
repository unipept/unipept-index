use std::path::Path;
use fm_builder::{Arguments, build_fm};
use clap::Parser;
use sa_mappings::proteins::Proteins;

use std::time::{SystemTime, SystemTimeError, UNIX_EPOCH};

fn main() {
    let Arguments {
        database_file,
        output,
        sparseness_factor,
    } = Arguments::parse();
    let base_output_path = Path::new(&output);

    eprintln!();
    eprintln!("📋 Started loading the proteins...");
    let start_proteins_time = get_time_ms().unwrap();
    let data = Proteins::try_from_database_file_uncompressed(&database_file)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!(
        "✅ Successfully loaded the proteins in {} seconds!",
        (get_time_ms().unwrap() - start_proteins_time) / 1000.0
    );

    eprintln!();
    eprintln!("📋 Started building the FM-index...");
    let start_fm_time = get_time_ms().unwrap();
    build_fm(data, sparseness_factor, base_output_path)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    eprintln!(
        "✅ Successfully built the FM-index in {} seconds!",
        (get_time_ms().unwrap() - start_fm_time) / 1000.0
    );

}

fn eprint_and_exit(err: &str) -> ! {
    eprintln!("{}", err);
    std::process::exit(1);
}

pub fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}