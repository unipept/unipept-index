use std::path::Path;
use fm_index::FMIndex;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let base_path = Path::new("../test_files/uniprot");
    let fm_index = FMIndex::from_files(base_path)?;

    let (sp, ep) = (0, 1001);
    let (sp, ep) = fm_index.extend(b'A', sp, ep);
    println!("'A' occurs in BWT range: {}..{}", sp, ep);

    for i in sp..ep {
        let pos = fm_index.locate(i);
        println!("'A' occurs in position: {}", pos);
    }

    Ok(())
}