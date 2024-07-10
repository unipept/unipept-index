use std::{
    error::Error,
    fs::File,
    io::{
        BufReader,
        Read
    },
    sync::Arc
};

use axum::{
    extract::{
        DefaultBodyLimit,
        State
    },
    http::StatusCode,
    routing::{
        get,
        post
    },
    Json,
    Router
};
use clap::Parser;
use sa_compression::load_compressed_suffix_array;
use sa_index::{
    binary::load_suffix_array,
    peptide_search::{
        analyse_all_peptides,
        search_all_peptides,
        OutputData,
        SearchOnlyResult,
        SearchResultWithAnalysis
    },
    sa_searcher::Searcher,
    suffix_to_protein_index::SparseSuffixToProtein,
    SuffixArray
};
use sa_mappings::{
    functionality::FunctionAggregator,
    proteins::Proteins,
    taxonomy::{
        AggregationMethod,
        TaxonAggregator
    }
};
use serde::{
    Deserialize,
    Serialize
};

/// Enum that represents all possible commandline arguments
#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be
    /// concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    index_file:    String,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy:      String
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
#[derive(Debug, Deserialize, Serialize)]
#[allow(non_snake_case)]
struct InputData {
    peptides:         Vec<String>,
    #[serde(default = "default_cutoff")] // default value is 10000
    cutoff: usize,
    #[serde(default = "bool::default")]
    // default value is false // TODO: maybe default should be true?
    equate_il: bool
}

#[tokio::main]
async fn main() {
    let args = Arguments::parse();
    if let Err(err) = start_server(args).await {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}

/// Basic handler used to check the server status
async fn root() -> &'static str {
    "Server is online"
}

/// Endpoint executed for peptide matching and taxonomic and functional analysis
///
/// # Arguments
/// * `state(searcher)` - The searcher object provided by the server
/// * `data` - InputData object provided by the user with the peptides to be searched and the config
///
/// # Returns
///
/// Returns the search and analysis results from the index as a JSON
async fn analyse(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>
) -> Result<Json<OutputData<SearchResultWithAnalysis>>, StatusCode> {
    let search_result = analyse_all_peptides(
        &searcher,
        &data.peptides,
        data.cutoff,
        data.equate_il
    );

    Ok(Json(search_result))
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
) -> Result<Json<OutputData<SearchOnlyResult>>, StatusCode> {
    let search_result = search_all_peptides(
        &searcher,
        &data.peptides,
        data.cutoff,
        data.equate_il,
    );

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
    let Arguments {
        database_file,
        index_file,
        taxonomy
    } = args;

    eprintln!();
    eprintln!("📋 Started loading the suffix array...");
    let sa = load_suffix_array_file(&index_file)?;
    eprintln!("✅ Successfully loaded the suffix array!");
    eprintln!("\tAmount of items: {}", sa.len());
    eprintln!("\tAmount of bits per item: {}", sa.bits_per_value());
    eprintln!("\tSample rate: {}", sa.sample_rate());

    eprintln!();
    eprintln!("📋 Started loading the taxon file...");
    let taxon_id_calculator =
        TaxonAggregator::try_from_taxonomy_file(&taxonomy, AggregationMethod::LcaStar)?;
    eprintln!("✅ Successfully loaded the taxon file!");
    eprintln!("\tAggregation method: LCA*");

    eprintln!();
    eprintln!("📋 Started creating the function aggregator...");
    let function_aggregator = FunctionAggregator {};
    eprintln!("✅ Successfully created the function aggregator!");

    eprintln!();
    eprintln!("📋 Started loading the proteins...");
    let proteins = Proteins::try_from_database_file(&database_file, &taxon_id_calculator)?;
    let suffix_index_to_protein = Box::new(SparseSuffixToProtein::new(&proteins.input_string));
    eprintln!("✅ Successfully loaded the proteins!");

    let searcher = Arc::new(Searcher::new(
        sa,
        suffix_index_to_protein,
        proteins,
        taxon_id_calculator,
        function_aggregator
    ));

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `POST /analyse` goes to `analyse` and set max payload size to 5 MB
        .route("/analyse", post(analyse))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher.clone())
        // `POST /search` goes to `search` and set max payload size to 5 MB
        .route("/search", post(search))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;

    eprintln!();
    eprintln!("🚀 Server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}

fn load_suffix_array_file(file: &str) -> Result<SuffixArray, Box<dyn Error>> {
    // Open the suffix array file
    let mut sa_file = File::open(file)?;

    // Create a buffer reader for the file
    let mut reader = BufReader::new(&mut sa_file);

    // Read the bits per value from the binary file (1 byte)
    let mut bits_per_value_buffer = [0_u8; 1];
    reader
        .read_exact(&mut bits_per_value_buffer)
        .map_err(|_| "Could not read the flags from the binary file")?;
    let bits_per_value = bits_per_value_buffer[0];

    if bits_per_value == 64 {
        load_suffix_array(&mut reader)
    } else {
        load_compressed_suffix_array(&mut reader, bits_per_value as usize)
    }
}
