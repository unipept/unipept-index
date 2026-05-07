//! `Proteins` — per-build concrete backends with a type alias.
//!
//! Submodules:
//! - [`preloaded`] — [`InMemoryProteins`] (always available)
//! - [`mmap`]      — [`MmapBackedProteins`] (mmap builds only)
//!
//! The alias [`Proteins`] resolves to the active backend.

pub mod preloaded;
#[cfg(feature = "mmap")]
pub mod mmap;

pub use preloaded::InMemoryProteins;
#[cfg(feature = "mmap")]
pub use mmap::MmapBackedProteins;

/// Type alias — resolves to the active backend for this build.
#[cfg(feature = "mmap")]
pub type Proteins = MmapBackedProteins;
#[cfg(not(feature = "mmap"))]
pub type Proteins = InMemoryProteins;

// Re-export I/O traits used by callers.
pub use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};

// ── Shared types ──────────────────────────────────────────────────────────────

use fa_compression::algorithm1::decode;

pub static SEPARATION_CHARACTER: u8 = b'-';
pub static TERMINATION_CHARACTER: u8 = b'$';

pub struct Protein {
    pub uniprot_id: String,
    pub taxon_id: u32,
    pub functional_annotations: Vec<u8>,
}

impl Protein {
    pub fn get_functional_annotations(&self) -> String {
        decode(&self.functional_annotations)
    }
}

#[derive(Clone, Copy)]
pub struct ProteinRef<'a> {
    pub uniprot_id: &'a str,
    pub taxon_id: u32,
    pub functional_annotations: &'a [u8],
}

impl<'a> ProteinRef<'a> {
    pub fn get_functional_annotations(&self) -> String {
        decode(self.functional_annotations)
    }
}
