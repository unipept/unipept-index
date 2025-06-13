use std::path::Path;
use fm_index::search_scheme::SearchScheme;
use fm_index::fm_index::FMIndex;
use fm_index::search::approximate_search;
use std::error::Error;
use std::fs;

fn main() -> Result<(), Box<dyn Error>> {
    let base_path = Path::new("../test_files/uniprot");
    let fm_index = FMIndex::from_files(base_path)?;

    let searches_path = Path::new("../test_files/kuch_k+2_searches.txt");
    let search_scheme = SearchScheme::from_file(searches_path)?;
    search_scheme.validate()?;
    let matches = approximate_search(&fm_index, b"ICQQADTVLAKKRVDLHMTREEMLTER".to_vec(), search_scheme)?;

    let input_path = Path::new("../test_files/uniprot_entries.1000.txt");
    let text = String::from_utf8(fs::read(input_path)?)?;
    for m in matches {
        println!("{}", m);
        for t in m..(m+5) {
            print!("{}", text.chars().nth(t).unwrap());
        }
        println!();
    }

    Ok(())
}