use std::path::Path;
use fm_index::{FMIndex, FMIndexRange};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let base_path = Path::new("../test_files/uniprot");
    let fm_index = FMIndex::from_files(base_path)?;

    let range = FMIndexRange { begin: 0, end: 1001, begin_rev: 0, end_rev: 1001 };
    let range = fm_index.left_extension(b'E', range);
    let range = fm_index.right_extension(b'P', range);
    let range = fm_index.left_extension(b'P', range);
    println!("'PEP' occurs in BWT range: {}..{}", range.begin, range.end);

    for i in range.begin..range.end {
        let pos = fm_index.locate(i);
        println!("'PEP' occurs in position: {}", pos);
    }

    Ok(())
}