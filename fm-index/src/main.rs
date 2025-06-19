use std::io::BufReader;
use std::path::Path;
use fm_index::search_scheme::SearchScheme;
use fm_index::fm_index::FMIndex;
use fm_index::search::search_multiple;
use std::error::Error;
use std::io::{BufRead, Write, self};
use std::time::Instant;
use std::fs;
use std::fs::File;

static FM_PATH: &str = "../test_files/uniprot";
static SEARCH_SCHEME_PATH: &str = "../search_schemes/kuch_k+1/1/searches.txt";
static PEPTIDES_PATH: &str = "../test_files/casanovo_peptides_1000.txt";
static TEXT_PATH: &str = "../test_files/uniprot_entries.100k.txt";

fn main() -> Result<(), Box<dyn Error>> {
    let base_path = Path::new(FM_PATH);
    eprintln!("Start loading index...");
    let fm_index = FMIndex::from_files(base_path)?;
    eprintln!("Index loaded");

    let searches_path = Path::new(SEARCH_SCHEME_PATH);
    let search_scheme = SearchScheme::from_file(searches_path)?;
    search_scheme.validate()?;

    let peptides = BufReader::new(File::open(PEPTIDES_PATH)?)
        .lines()
        .map(|line| Ok(line?.into_bytes()))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;

    eprintln!("Start searching...");
    let start = Instant::now();
    let matches = search_multiple(&fm_index, peptides, search_scheme).map_err(|e| e as Box<dyn Error>)?;
    let duration = start.elapsed();
    eprintln!("Searching done in {:?}", duration);

    eprintln!("There are {} matches", matches.len());

    eprintln!("Reporting matches...");
    let mut writer = io::stdout();
    let input_path = Path::new(TEXT_PATH);
    let text = String::from_utf8(fs::read(input_path)?)?;
    for m in matches {
        let end = m.start_position + m.length;
        writeln!(writer, "{}", &text[m.start_position..end])?;
    }

    Ok(())
}