use std::io::BufReader;
use std::path::Path;
use fm_index::search_scheme::SearchScheme;
use fm_index::fm_index::FMIndex;
use fm_index::search::{search_multiple, search_multiple_grouped, search_multiple_exact, search_multiple_exact_grouped};
use std::error::Error;
use std::io::{BufRead, Write, self};
use std::time::Instant;
use std::fs;
use std::fs::File;

static FM_PATH: &str = "../test_files/uniprot";
static SEARCH_SCHEME_PATH: &str = "../search_schemes/kuch_k+1/1/searches.txt";
static PEPTIDES_PATH: &str = "../test_files/casanovo_peptides_1000.txt";
static TEXT_PATH: &str = "../test_files/uniprot_entries.100k.txt";

fn approximate_search(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<(), Box<dyn Error>> {
    let searches_path = Path::new(SEARCH_SCHEME_PATH);
    let search_scheme = SearchScheme::from_file(searches_path)?;
    search_scheme.validate()?;

    eprintln!("Start searching...");
    let start = Instant::now();
    let matches = search_multiple(fm_index, patterns, &search_scheme).map_err(|e| e as Box<dyn Error>)?;
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

fn exact_search(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<(), Box<dyn Error>> {
    eprintln!("Start searching...");
    let start = Instant::now();
    let matches = search_multiple_exact(&fm_index, patterns).map_err(|e| e as Box<dyn Error>)?;
    let duration = start.elapsed();
    eprintln!("Searching done in {:?}", duration);

    eprintln!("There are {} matches", matches.len());

    eprintln!("Reporting matches...");
    let mut writer = io::stdout();
    for m in matches {
        writeln!(writer, "{}", m)?;
    }

    Ok(())
}


fn parameter_search(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<(), Box<dyn Error>> {

    let mut search_schemes = vec![None];
    for i in 0..3 {
        let path_str = format!("../search_schemes/kuch_k+1/{}/searches.txt", i);
        println!("{}", &path_str);
        let searches_path = Path::new(&path_str);
        search_schemes.push(Some(SearchScheme::from_file(searches_path)?));
    }

    let mut writer = io::stdout();
    writeln!(writer, "Pattern,Edit distance,matches")?;

    for (i, search_scheme) in search_schemes[0..3].iter().enumerate() {

        let patterns: Vec<Vec<u8>> = patterns
            .iter()
            .cloned()
            .filter(|v| v.len() >= 7+i)
            .collect();

        eprintln!("Start searching for ED {}...", i);
        let start = Instant::now();
        let matches: Vec<usize> = match search_scheme {
            Some(scheme) => search_multiple_grouped(fm_index, patterns.clone(), &scheme).map_err(|e| e as Box<dyn Error>)?
                .iter()
                .map(|match_set| match_set.len())
                .collect(),
            None => search_multiple_exact_grouped(fm_index, patterns.clone()).map_err(|e| e as Box<dyn Error>)?
                .iter()
                .map(|match_set| match_set.len())
                .collect()
        };
        let duration = start.elapsed();
        eprintln!("Searching done in {:?}", duration);

        eprintln!("Reporting matches...");
        for j in 0..matches.len() {
            writeln!(writer, "{},{},{}", String::from_utf8(patterns[j].clone()).expect("Invalid UTF-8"), i, matches[j])?;
        }

    }

    Ok(())
}


fn main() -> Result<(), Box<dyn Error>> {

    let base_path = Path::new(FM_PATH);
    eprintln!("Start loading index...");
    let fm_index = FMIndex::from_files(base_path)?;
    eprintln!("Index loaded");

    let patterns = BufReader::new(File::open(PEPTIDES_PATH)?)
        .lines()
        .map(|line| Ok(line?.into_bytes()))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;

    parameter_search(&fm_index, patterns)
}