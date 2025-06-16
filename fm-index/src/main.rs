use std::io::BufReader;
use std::path::Path;
use fm_index::search_scheme::SearchScheme;
use fm_index::fm_index::FMIndex;
use fm_index::search::search_multiple;
use std::error::Error;
use std::io::BufRead;
use std::time::Instant;
use std::fs;
use std::fs::File;

fn main() -> Result<(), Box<dyn Error>> {
    let base_path = Path::new("../test_files/uniprot");
    println!("Start loading index...");
    let fm_index = FMIndex::from_files(base_path)?;
    println!("Index loaded");

    let searches_path = Path::new("../search_schemes/kuch_k+1/1/searches.txt");
    let search_scheme = SearchScheme::from_file(searches_path)?;
    search_scheme.validate()?;

    let peptides = BufReader::new(File::open("../test_files/casanovo_peptides_10.txt")?)
        .lines()
        .map(|line| Ok(line?.into_bytes()))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;

    println!("Start searching...");
    let start = Instant::now();
    let matches = search_multiple(&fm_index, peptides, search_scheme).map_err(|e| e as Box<dyn Error>)?;
    let duration = start.elapsed();
    println!("Searching done in {:?}", duration);

    println!("{}", matches.len());

    let input_path = Path::new("../test_files/uniprot_entries.100k.txt");
    let text = String::from_utf8(fs::read(input_path)?)?;
    for m in matches {
        for t in m.start_position..(m.start_position+m.length) {
            print!("{}", text.chars().nth(t).unwrap());
        }
        println!();
    }

    Ok(())
}