use std::collections::HashSet;
use std::io::BufReader;
use std::path::Path;
use fm_index::search_scheme::SearchScheme;
use fm_index::fm_index::FMIndex;
use fm_index::search::{search_multiple, search_multiple_exact, search_multiple_exact_grouped, search_multiple_grouped, FMMatch};
use std::error::Error;
use std::io::{BufRead, Write, self};
use std::time::Instant;
use std::fs;
use std::fs::File;
use std::env;

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
        writeln!(writer, "{}", m.start_position)?;
    }

    Ok(())
}


fn parameter_search(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<(), Box<dyn Error>> {

    let max_ed = 3;

    let mut search_schemes = vec![None];
    for i in 1..max_ed {
        let path_str = format!("../search_schemes/kuch_k+1/{}/searches.txt", i);
        println!("{}", &path_str);
        let searches_path = Path::new(&path_str);
        search_schemes.push(Some(SearchScheme::from_file(searches_path)?));
    }

    let mut writer = io::stdout();
    writeln!(writer, "Pattern,Edit distance,matches")?;

    for (i, search_scheme) in search_schemes[0..max_ed].iter().enumerate() {

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

        eprintln!("Reporting match counts...");
        for j in 0..matches.len() {
            writeln!(writer, "{},{},{}", String::from_utf8(patterns[j].clone()).expect("Invalid UTF-8"), i, matches[j])?;
        }

    }

    Ok(())
}

fn parameter_search_in_text(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, ed: i32) -> Result<HashSet<FMMatch>, Box<dyn Error>> {

    let mut search_scheme = None;
    if ed > 0 {
        let path_str = format!("../search_schemes/kuch_k+1/{}/searches.txt", ed);
        let searches_path = Path::new(&path_str);
        search_scheme = Some(SearchScheme::from_file(searches_path)?);
    }

    let patterns: Vec<Vec<u8>> = patterns
        .iter()
        .cloned()
        .filter(|v| v.len() >= 7 + ed as usize)
        .collect();

    eprintln!("Start searching for ED {}...", ed);
    let start = Instant::now();
    let matches = match search_scheme {
        Some(scheme) => search_multiple(fm_index, patterns.clone(), &scheme).map_err(|e| e as Box<dyn Error>)?,
        None => search_multiple_exact(fm_index, patterns.clone()).map_err(|e| e as Box<dyn Error>)?
    };
    
    let duration = start.elapsed();
    eprintln!("Searching done in {:?}", duration);

    Ok(matches)
}


fn main() -> Result<(), Box<dyn Error>> {

    let args: Vec<String> = env::args().collect();
    let output_base = &args[1];

    let base_path = Path::new(FM_PATH);
    eprintln!("Start loading index...");
    let fm_index = FMIndex::from_files(base_path)?;
    eprintln!("Index loaded");

    let patterns = BufReader::new(File::open(PEPTIDES_PATH)?)
        .lines()
        .map(|line| Ok(line?.into_bytes()))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;

    let max_ed=2;
    let mut matches_per_ed: Vec<HashSet<FMMatch>> = Vec::new();
    for ed in 0..=max_ed {
        matches_per_ed.push(parameter_search_in_text(&fm_index, patterns.clone(), ed)?);
    }
    
    drop(fm_index);

    for ed in 0..=max_ed {
        eprintln!("Reporting matches for edit distance {}...", ed);
        let output_path = format!("{}_ed{}.txt", output_base, ed);
        let output_path = Path::new(&output_path);
        let mut writer = File::create(output_path)?;
        let input_path = Path::new(TEXT_PATH);
        let text = String::from_utf8(fs::read(input_path)?)?;

        let mut sorted_matches: Vec<_> = matches_per_ed[ed as usize].clone().into_iter().collect();
        sorted_matches.sort_by_key(|m| m.start_position);

        let mut matches_intext: HashSet<&str> = HashSet::new();
        for m in sorted_matches {
            let end = m.start_position + m.length;
            let _ = matches_intext.insert(&text[m.start_position..end]);
        }

        for m in matches_intext {
            writeln!(writer, "{}", m)?;
        }
    }
    

    Ok(())
}