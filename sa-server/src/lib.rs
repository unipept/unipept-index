//! Loading the index files, in whichever way this build was compiled to.
//!
//! Each index structure has two readers — one that copies into owned memory, one that maps the
//! file — and which is used is decided at compile time, per structure. The loaders below hide
//! that choice so `main` never mentions it.
//!
//! See [`backend_summary`] for reporting which ones are active, and the `sa-index` crate docs for
//! the trade-off.

use std::{error::Error, fs::File, io::BufReader};

#[cfg(feature = "mmap")]
use sa_index::ReadBinaryMmap;
use sa_index::{KmerTable, ReadBinary, SuffixArray, suffix_to_protein_index::SuffixToProteinMapping};
use sa_mappings::proteins::Proteins;

/// How the suffix array is stored in this build.
///
/// Baked in at compile time, so it cannot be inspected any other way at runtime — which is why
/// the server logs all four at startup. The configurations have very different memory profiles,
/// and telling them apart otherwise means checking how the binary was built.
pub const SA_BACKEND: &str = if cfg!(feature = "mmap") { "mmap" } else { "preloaded" };
/// How the concatenated protein text is stored in this build.
pub const TEXT_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-text"))) { "mmap" } else { "preloaded" };
/// How the protein metadata table is stored in this build.
pub const PROTEINS_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-proteins"))) { "mmap" } else { "preloaded" };
/// How the suffix-to-protein mapping is stored in this build.
pub const MAPPING_BACKEND: &str =
    if cfg!(all(feature = "mmap", not(feature = "preloaded-mapping"))) { "mmap" } else { "preloaded" };

/// One line naming the storage of every structure, for the startup log.
pub fn backend_summary() -> String {
    format!("sa={SA_BACKEND} text={TEXT_BACKEND} proteins={PROTEINS_BACKEND} mapping={MAPPING_BACKEND}")
}

/// Reads a structure that has both an owned and a memory-mapped reader, choosing by `$mapped_when`.
///
/// The predicate is passed in rather than hard-coded because the structures are no longer selected
/// together: `mmap` maps all of them, and a `preloaded-*` feature pulls one back without affecting
/// the others. The loaders below were four copies of this three-line body.
macro_rules! load_by_backend {
    ($ty:ty, $file:expr, $mapped_when:meta) => {{
        #[cfg($mapped_when)]
        {
            <$ty>::read_binary_mmap(std::path::Path::new($file))
        }
        #[cfg(not($mapped_when))]
        {
            let f = File::open($file)?;
            let mut reader = BufReader::new(f);
            <$ty>::read_binary(&mut reader)
        }
    }};
}

/// Loads the suffix array. It has no `preloaded-*` override — it follows `mmap`.
pub fn load_suffix_array_file(file: &str) -> Result<SuffixArray, Box<dyn Error>> {
    load_by_backend!(SuffixArray, file, feature = "mmap")
}

/// Loads the protein table and the concatenated protein text (both live in `proteins.bin`).
///
/// Note the predicate: the file has to be *mapped* whenever either of its two sections is mapped,
/// which is why this is not simply `not(preloaded-proteins)`. With only `preloaded-proteins` set
/// the metadata is copied out but the text still borrows the mapping, so the reader is still the
/// mmap one; only when both sections are preloaded does the plain buffered reader apply.
pub fn load_proteins_file(file: &str) -> Result<Proteins, Box<dyn Error>> {
    load_by_backend!(
        Proteins,
        file,
        all(feature = "mmap", not(all(feature = "preloaded-text", feature = "preloaded-proteins")))
    )
}

/// Loads the suffix-to-protein mapping.
pub fn load_mapping_file(file: &str) -> Result<SuffixToProteinMapping, Box<dyn Error>> {
    load_by_backend!(SuffixToProteinMapping, file, all(feature = "mmap", not(feature = "preloaded-mapping")))
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
