use std::path::Path;
use fm_builder::build_fm;
use std::error::Error;
use std::fs;

fn main() -> Result<(), Box<dyn Error>> {
    let input_path = Path::new("../test_files/uniprot_entries.100M.txt");
    let text = fs::read(input_path)?;
    // let text  = b"BANANA$".to_vec();

    let base_path = Path::new("../test_files/uniprot.100M");
    build_fm(text, 8, base_path)?;

    Ok(())
}