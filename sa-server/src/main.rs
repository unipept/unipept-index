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
use sa_index::{peptide_search::{SearchResult, search_all_peptides}, sa_searcher::Searcher, SuffixArray, SuffixArrayBackend};
use serde::Deserialize;
use sa_server::{load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};

/// Enum that represents all possible commandline arguments
#[derive(Parser, Debug)]
pub struct Arguments {
    /// Path to the database file. This should point to the binary proteins file (.proteins.bin)
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    index_file: String,
    /// Path to the prebuilt suffix-to-protein mapping binary file.
    #[arg(long)]
    mapping_file: String,
    /// Optional path to a pre-built k-mer bounds table file (produced by sa-builder
    /// --output-kmer-table). When provided, binary search is accelerated by ~60 %.
    #[arg(long)]
    kmer_table_file: Option<String>,
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
    State(searcher): State<Arc<Searcher<SuffixArray>>>,
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
    let Arguments { database_file, index_file, mapping_file, kmer_table_file } = args;

    eprintln!();
    eprintln!("Started loading the suffix array...");
    let suffix_array = load_suffix_array_file(&index_file)?;
    eprintln!("Successfully loaded the suffix array!");
    eprintln!("\tAmount of items: {}", suffix_array.len());
    eprintln!("\tAmount of bits per item: {}", suffix_array.bits_per_value());
    eprintln!("\tSample rate: {}", suffix_array.sample_rate());

    eprintln!();
    eprintln!("Started loading the proteins...");
    let proteins = load_proteins_file(&database_file)?;
    eprintln!("Successfully loaded the proteins!");

    eprintln!();
    eprintln!("Started loading the suffix-to-protein mapping...");
    let mapping = load_mapping_file(&mapping_file)?;
    eprintln!("Successfully loaded the suffix-to-protein mapping!");

    let mut searcher = Searcher::new(suffix_array, proteins, mapping);

    if let Some(ref path) = kmer_table_file {
        eprintln!();
        eprintln!("Started loading the k-mer table...");
        let table = load_kmer_table_file(path)?;
        eprintln!("Successfully loaded the k-mer table! (k={})", table.k);
        searcher = searcher.with_kmer_table(table);
    }

    let searcher = Arc::new(searcher);

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
