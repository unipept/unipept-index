//! Loading the index files this build was compiled to use.
//!
//! Every index structure has two implementations — one holding owned memory, one borrowing a
//! memory mapping — and both are always compiled. [`backends`] names one of each for this build;
//! the loaders below are then just `LoadIndex::load` on those names, since each concrete type
//! already knows whether it is read or mapped.
//!
//! See [`backends`] for what the storage features select and why, and
//! [`backends::backend_summary`] for reporting it at startup.

use std::{error::Error, path::Path};

use sa_index::{KmerTable, LoadIndex};

pub mod backends;

pub use backends::{
    ActiveMapping, ActiveProteins, ActiveSa, ActiveSearcher, ActiveText, MAPPING_BACKEND, PROTEINS_BACKEND, SA_BACKEND,
    TEXT_BACKEND, backend_summary
};

/// Loads the suffix array.
pub fn load_suffix_array_file(file: &str) -> Result<ActiveSa, Box<dyn Error>> {
    ActiveSa::load(Path::new(file))
}

/// Loads the protein table and the concatenated protein text (both live in `proteins.bin`).
pub fn load_proteins_file(file: &str) -> Result<ActiveProteins, Box<dyn Error>> {
    ActiveProteins::load(Path::new(file))
}

/// Loads the suffix-to-protein mapping.
pub fn load_mapping_file(file: &str) -> Result<ActiveMapping, Box<dyn Error>> {
    ActiveMapping::load(Path::new(file))
}

/// Loads a pre-built k-mer bounds table.
///
/// Unlike the three above this has no memory-mapped variant: the table is small relative to the
/// index and is read into owned memory in every configuration.
pub fn load_kmer_table_file(file: &str) -> Result<KmerTable, Box<dyn Error>> {
    KmerTable::load(Path::new(file))
}
