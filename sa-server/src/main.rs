use std::{
    error::Error,
    sync::Arc
};
use axum::{
    Json, Router,
    extract::{DefaultBodyLimit, State},
    http::StatusCode,
    routing::post
};
use clap::Parser;
use sa_index::{
    peptide_search::{SearchResult, search_all_peptides},
    sa_searcher::Searcher,
    suffix_to_protein_index::SuffixToProteinMapping
};
use serde::Deserialize;
use sa_server::{load_mapping_file, load_proteins_file, load_suffix_array_file};

/// Enum that represents all possible commandline arguments
#[derive(Parser, Debug)]
pub struct Arguments {
    /// Path to the database file. If --mmap is set, this should point to the binary proteins file
    /// (.proteins.bin); otherwise it should point to the TSV database file.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    index_file: String,
    /// Path to the prebuilt suffix-to-protein mapping binary file.
    #[arg(long)]
    mapping_file: String,
    /// Use memory-mapped I/O to load the suffix array and ProteinText. When set, --database-file
    /// must point to a binary proteins file (.proteins.bin). Makes startup near-instant by letting
    /// the OS page in data on demand, at the cost of slower initial queries while pages are loaded.
    #[arg(short, long, default_value_t = false)]
    mmap: bool
}

/// Function used by serde to place a default value in the cutoff field of the input
fn default_cutoff() -> usize {
    10000
}

/// Function used by serde to use `true` as a default value
#[allow(dead_code)]
fn default_true() -> bool {
    true
}

/// Struct representing the input arguments accepted by the endpoints
///
/// # Arguments
/// * `peptides` - List of peptides we want to process
/// * `cutoff` - The maximum amount of matches to process, default value 10000
/// * `equate_il` - True if we want to equalize I and L during search
/// * `clean_taxa` - True if we only want to use proteins marked as "valid"
#[derive(Debug, Deserialize)]
struct InputData {
    peptides: Vec<String>,
    #[serde(default = "default_cutoff")] // default value is 10000
    cutoff: usize,
    #[serde(default = "bool::default")]
    // default value is false // TODO: maybe default should be true?
    equate_il: bool,
    #[serde(default = "bool::default")] // default false
    tryptic: bool
}

#[tokio::main]
async fn main() {
    let args = Arguments::parse();
    if let Err(err) = start_server(args).await {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}

/// Endpoint executed for peptide matching, without any analysis
///
/// # Arguments
/// * `state(searcher)` - The searcher object provided by the server
/// * `data` - InputData object provided by the user with the peptides to be searched and the config
///
/// # Returns
///
/// Returns the search results from the index as a JSON
async fn search(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>
) -> Result<Json<Vec<SearchResult>>, StatusCode> {
    let search_result = search_all_peptides(&searcher, &data.peptides, data.cutoff, data.equate_il, data.tryptic);

    Ok(Json(search_result))
}

/// Starts the server with the provided commandline arguments
///
/// # Arguments
/// * `args` - The provided commandline arguments
///
/// # Returns
///
/// Returns ()
///
/// # Errors
///
/// Returns any error occurring during the startup or uptime of the server
async fn start_server(args: Arguments) -> Result<(), Box<dyn Error>> {
    let Arguments { database_file, index_file, mmap, mapping_file } = args;

    eprintln!();
    eprintln!("Started loading the suffix array...");
    let suffix_array = load_suffix_array_file(&index_file, mmap)?;
    eprintln!("Successfully loaded the suffix array!");
    eprintln!("\tAmount of items: {}", suffix_array.len());
    eprintln!("\tAmount of bits per item: {}", suffix_array.bits_per_value());
    eprintln!("\tSample rate: {}", suffix_array.sample_rate());

    eprintln!();
    eprintln!("Started loading the proteins...");
    let proteins = load_proteins_file(&database_file, mmap)?;
    eprintln!("Successfully loaded the proteins!");

    eprintln!();
    eprintln!("Started loading the suffix-to-protein mapping...");
    let SuffixToProteinMapping(mapping) = load_mapping_file(&mapping_file, mmap)?;
    eprintln!("Successfully loaded the suffix-to-protein mapping!");

    let searcher = Arc::new(Searcher::new(suffix_array, proteins, mapping));

    // build our application with a route
    let app = Router::new()
        .route("/search", post(search))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;

    eprintln!();
    eprintln!("🚀 Server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}
