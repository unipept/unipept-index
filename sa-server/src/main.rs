//! The `sa-server` binary: loads a prebuilt index and serves peptide searches over HTTP.
//!
//! One route, `POST /search`, taking [`InputData`] and returning an array of search results. The
//! request and response fields are tabulated in the repository README, which is where a caller
//! looks; this file documents why the server is put together the way it is.
//!
//! **This is not how search reaches Unipept.** The API site calls `sa_index::peptide_search`
//! directly rather than proxying this binary, so this server is not the untrusted boundary and
//! hardening it protects nothing in production — the limits belong in the caller that exposes
//! those functions to the network. See the note on resource limits in `sa_index::peptide_search`.
//!
//! Which storage backend each index structure uses is a compile-time choice; see
//! [`sa_server::backends`]. Startup logs the combination because nothing else at runtime reveals
//! it.

use std::{error::Error, sync::Arc};

use axum::{
    Json, Router,
    body::Body,
    extract::{DefaultBodyLimit, State},
    http::{HeaderValue, Response, header},
    response::IntoResponse,
    routing::post
};
use clap::Parser;
use sa_index::{SuffixArrayBackend, peptide_search::search_all_peptides_json, sa_searcher::Searcher};
use sa_server::{ActiveSearcher, load_kmer_table_file, load_mapping_file, load_proteins_file, load_suffix_array_file};
use serde::Deserialize;

/// Serve peptide searches over a prebuilt suffix-array index.
///
/// All four files come from one `sa-builder` run and only mean anything as a set: the suffix array
/// indexes positions in the text stored in `proteins.bin`, and the mapping resolves those same
/// positions back to entries in it. Each file is loaded and validated on its own, so a mismatched
/// set gets past loading; `Searcher::try_new` is what refuses to start on one.
#[derive(Parser, Debug)]
pub struct Arguments {
    /// Path to the binary proteins file (proteins.bin) holding the protein table and text.
    #[arg(short, long)]
    database_file: String,
    /// Path to the binary suffix array (sa.bin).
    #[arg(short, long)]
    index_file: String,
    /// Path to the prebuilt suffix-to-protein mapping binary file.
    #[arg(long)]
    mapping_file: String,
    /// Optional path to a pre-built k-mer bounds table file (produced by sa-builder
    /// --output-kmer-table). When provided, the binary search starts from precomputed bounds
    /// instead of the whole array, which is a large saving on short peptides.
    #[arg(long)]
    kmer_table_file: Option<String>,
    /// Address to bind the HTTP server to.
    #[arg(long, default_value = "0.0.0.0:3000")]
    address: String
}

/// The `cutoff` a request that omits the field gets.
///
/// Serde has no way to express a literal default, so the value lives in a function. It is the
/// historical Unipept limit rather than a property of the index — a short peptide can match
/// millions of proteins, and a caller that wants all of them asks for them explicitly.
fn default_cutoff() -> usize {
    10000
}

/// The body of a `POST /search` request.
///
/// Every field but `peptides` is optional and defaults to the pre-existing behaviour, so an older
/// caller that sends only a peptide list keeps getting the same answers.
#[derive(Debug, Deserialize)]
struct InputData {
    /// The peptides to search for. Each one is answered independently, in the order given.
    peptides: Vec<String>,
    /// Stop collecting matches for a peptide after this many; the result then reports
    /// `cutoff_used`. Defaults to [`default_cutoff`].
    #[serde(default = "default_cutoff")]
    cutoff: usize,
    /// Treat I and L as interchangeable. Defaults to false: the caller decides, because the
    /// answer differs and mass spectrometry cannot always distinguish the two.
    #[serde(default = "bool::default")]
    equate_il: bool,
    /// Return only matches at tryptic boundaries. Defaults to false.
    #[serde(default = "bool::default")]
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

/// `POST /search`: matches every peptide in the request against the index.
///
/// Search only — no taxonomic or functional analysis is done here, and the results come back in
/// request order. The searcher is shared behind an `Arc` because it is read-only once loaded and
/// far too large to clone; `search_all_peptides_json` fans the peptides out over rayon internally,
/// so concurrency does not come from axum handing out one handler task per request.
///
/// # Why this is not `Json<Vec<SearchResult>>`
///
/// `search_all_peptides_json` hands back the body already serialised, one chunk per peptide, built
/// on the rayon workers that did the search. Serialising here instead would put the whole answer —
/// hundreds of megabytes on a large non-tryptic request — through a single-threaded
/// `serde_json::to_vec`, which is the largest serial stretch a request has. It also could not
/// borrow: a `SearchResult` holds references into the `Searcher`, so it never satisfies the
/// `'static` bound a response body needs.
///
/// The bytes and the headers are exactly what `Json` produced: `content-type: application/json`
/// and nothing else, with hyper deriving `content-length` from the body. Note that `Vec<u8>`'s own
/// `IntoResponse` would set `application/octet-stream`, which is why the response is built
/// explicitly rather than through a header tuple.
async fn search(State(searcher): State<Arc<ActiveSearcher>>, data: Json<InputData>) -> Response<Body> {
    let chunks = search_all_peptides_json(&searcher, &data.peptides, data.cutoff, data.equate_il, data.tryptic);

    json_response(chunks)
}

/// Joins pre-serialised JSON chunks into the response `Json` would have produced.
///
/// Split out from [`search`] so it can be tested without an index: the chunks are already proved
/// byte-identical to a single `serde_json` pass by `sa_index::peptide_search`, and this is the
/// other half — that the HTTP response around them did not change either.
fn json_response(chunks: Vec<Vec<u8>>) -> Response<Body> {
    // Sized exactly, so the concatenation neither grows by doubling nor over-allocates.
    let mut body = Vec::with_capacity(chunks.iter().map(Vec::len).sum());
    for chunk in &chunks {
        body.extend_from_slice(chunk);
    }
    drop(chunks);

    let mut response = body.into_response();
    response.headers_mut().insert(header::CONTENT_TYPE, HeaderValue::from_static("application/json"));
    response
}

/// Loads the index, then serves until the process is killed.
///
/// Loading is eager and strictly ordered: every file is read and the set is checked for
/// consistency *before* the listener is bound, so a bad invocation fails at startup with a message
/// instead of binding a port and failing per request. Each step is announced on stderr because on
/// the full database this takes minutes, and silence during it is indistinguishable from a hang.
///
/// # Errors
///
/// Any failure loading the index, any mismatch between the files, and any error from binding or
/// serving. All of them reach `main`, which prints and exits non-zero.
async fn start_server(args: Arguments) -> Result<(), Box<dyn Error>> {
    let Arguments {
        database_file,
        index_file,
        mapping_file,
        kmer_table_file,
        address
    } = args;

    // Storage is a compile-time choice per structure, with a large effect on memory use, and
    // nothing else at runtime reveals which combination this binary has.
    eprintln!();
    eprintln!("Storage backends: {}", sa_server::backend_summary());

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

    // Refuse to start on a mismatched set rather than serving wrong answers from one. The three
    // files are loaded independently and each is well-formed on its own, so nothing before this
    // point would have noticed.
    let mut searcher = Searcher::try_new(suffix_array, proteins, mapping)?;

    if let Some(ref path) = kmer_table_file {
        eprintln!();
        eprintln!("Started loading the k-mer table...");
        let table = load_kmer_table_file(path)?;
        eprintln!("Successfully loaded the k-mer table! (k={})", table.k);
        searcher = searcher.try_with_kmer_table(table)?;
    }

    let searcher = Arc::new(searcher);

    // 5 MB of peptides, well above axum's 2 MB default: a realistic batch is tens of thousands of
    // short strings, and the response is orders of magnitude larger than the request that asked for
    // it anyway. This bounds what one request can buffer, not what it can cost — see the note on
    // resource limits in the crate docs.
    let app = Router::new()
        .route("/search", post(search))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind(&address).await?;

    eprintln!();
    eprintln!("🚀 Server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use axum::body::to_bytes;
    use serde::Serialize;

    use super::*;

    /// The shape `Json<Vec<SearchResult>>` serialises to, kept here as the reference the streamed
    /// response is compared against. It only has to serialise the same way, not be the same type.
    #[derive(Serialize)]
    struct ReferenceResult {
        sequence: String,
        proteins: Vec<ReferenceProtein>,
        cutoff_used: bool
    }

    #[derive(Serialize)]
    struct ReferenceProtein {
        taxon: u32,
        uniprot_accession: String,
        functional_annotations: String
    }

    fn reference(count: usize) -> Vec<ReferenceResult> {
        (0..count)
            .map(|index| ReferenceResult {
                sequence: format!("PEPTIDE{index}"),
                proteins: (0..index)
                    .map(|protein| ReferenceProtein {
                        taxon: protein as u32,
                        uniprot_accession: format!("P{protein:05}"),
                        functional_annotations: "EC:1.1.1.-;GO:0009279".to_string()
                    })
                    .collect(),
                cutoff_used: index % 2 == 0
            })
            .collect()
    }

    /// Chunks in the shape `sa_index::peptide_search::json_chunk` produces them: each prefixed with
    /// the `,` that separates it, the first `,` overwritten with `[`, the last carrying the `]`.
    fn chunks_for(results: &[ReferenceResult]) -> Vec<Vec<u8>> {
        if results.is_empty() {
            return vec![b"[]".to_vec()];
        }
        let mut chunks: Vec<Vec<u8>> = results
            .iter()
            .map(|result| {
                let mut chunk = vec![b','];
                serde_json::to_writer(&mut chunk, result).unwrap();
                chunk
            })
            .collect();
        chunks[0][0] = b'[';
        chunks.last_mut().unwrap().push(b']');
        chunks
    }

    async fn body_bytes(response: Response<Body>) -> Vec<u8> {
        to_bytes(response.into_body(), usize::MAX).await.unwrap().to_vec()
    }

    /// The whole contract of moving serialisation out of `Json`: same bytes, same headers.
    #[tokio::test]
    async fn json_response_matches_what_json_would_have_returned() {
        for count in [0, 1, 2, 17] {
            let results = reference(count);

            let ours = json_response(chunks_for(&results));
            let theirs = Json(&results).into_response();

            assert_eq!(ours.status(), theirs.status(), "count={count}");
            assert_eq!(ours.headers(), theirs.headers(), "count={count}");
            assert_eq!(body_bytes(ours).await, body_bytes(theirs).await, "count={count}");
        }
    }

    /// `Vec<u8>`'s own `IntoResponse` sets `application/octet-stream`, so the content type has to be
    /// overwritten rather than merely added. A regression here would be invisible to a client that
    /// parses the body regardless.
    #[tokio::test]
    async fn json_response_is_labelled_as_json() {
        let response = json_response(chunks_for(&reference(3)));
        assert_eq!(response.headers().get(header::CONTENT_TYPE).unwrap(), "application/json");
        assert_eq!(response.headers().len(), 1, "Json sets exactly one header");
    }
}
