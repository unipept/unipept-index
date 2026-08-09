//! Loading the index files, in whichever way this build was compiled to.
//!
//! The three index structures each have two readers — one that copies into owned memory, one
//! that maps the file — and which is used is decided at compile time by the `mmap` feature. The
//! loaders below hide that choice so `main` never mentions it.
//!
//! See [`BACKEND`] for reporting which one is active, and the `sa-index` crate docs for the
//! trade-off.

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use sa_index::{KmerTable, ReadBinary, SuffixArray};
#[cfg(feature = "mmap")]
use sa_index::ReadBinaryMmap;
use sa_index::suffix_to_protein_index::SuffixToProteinMapping;
use sa_mappings::proteins::Proteins;

/// Which storage backend this binary was compiled with.
///
/// Baked in at compile time, so it cannot be inspected any other way at runtime — which is
/// exactly why the server logs it at startup. The two have very different memory profiles, and
/// telling them apart otherwise means checking how the binary was built.
pub const BACKEND: &str = if cfg!(feature = "mmap") { "mmap" } else { "preloaded" };

/// Reads a structure that has both an owned and a memory-mapped reader, choosing by feature.
///
/// The four loaders below were four copies of this three-line body.
macro_rules! load_by_backend {
    ($ty:ty, $file:expr) => {{
        #[cfg(feature = "mmap")]
        {
            <$ty>::read_binary_mmap(std::path::Path::new($file))
        }
        #[cfg(not(feature = "mmap"))]
        {
            let f = File::open($file)?;
            let mut reader = BufReader::new(f);
            <$ty>::read_binary(&mut reader)
        }
    }};
}

/// Loads the suffix array.
pub fn load_suffix_array_file(file: &str) -> Result<SuffixArray, Box<dyn Error>> {
    load_by_backend!(SuffixArray, file)
}

/// Loads the protein table and the concatenated protein text (both live in `proteins.bin`).
pub fn load_proteins_file(file: &str) -> Result<Proteins, Box<dyn Error>> {
    load_by_backend!(Proteins, file)
}

/// Loads the suffix-to-protein mapping.
pub fn load_mapping_file(file: &str) -> Result<SuffixToProteinMapping, Box<dyn Error>> {
    load_by_backend!(SuffixToProteinMapping, file)
}

/// Loads a pre-built k-mer bounds table.
///
/// Unlike the three above this has no memory-mapped variant: the table is small relative to the
/// index and is read into owned memory in both configurations.
pub fn load_kmer_table_file(file: &str) -> Result<KmerTable, Box<dyn Error>> {
    let f = File::open(file)?;
    let mut reader = BufReader::new(f);
    KmerTable::read_binary(&mut reader)
}
