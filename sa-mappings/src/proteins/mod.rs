//! Protein metadata, in either of two backends.
//!
//! Submodules:
//! - [`preloaded`] — [`InMemoryProteins`], the metadata table in owned memory
//! - [`mmap`] — [`MmapBackedProteins`], the table borrowed from a mapping
//!
//! Both are always compiled; callers pick by naming a type. Everything that reads protein metadata
//! is written against [`ProteinsBackend`].
//!
//! # Two independent axes
//!
//! Both structs are generic over their text backend, so the metadata and the text are stored
//! independently: the metadata table is the biggest structure in the index and the text the
//! hottest, and the best place for one is not the best place for the other. All four pairings
//! load from the same `proteins.bin`, which holds both sections — and which reader each pairing
//! needs is its own `LoadIndex` impl, since the file must be mapped whenever *either* section is.
//!
//! Everything both backends share lives here rather than in a submodule: the [`ProteinsBackend`]
//! trait they implement, the [`Protein`] / [`ProteinRef`] pair a lookup returns, and the
//! [`SEPARATION_CHARACTER`] and [`TERMINATION_CHARACTER`] delimiters that give the concatenated
//! text its structure.

pub mod mmap;
pub mod preloaded;
#[cfg(test)]
pub(crate) mod test_fixtures;

// ── Shared types ──────────────────────────────────────────────────────────────
use fa_compression::algorithm1::decode;
pub use mmap::MmapBackedProteins;
pub use preloaded::InMemoryProteins;
// The I/O traits callers need, re-exported so they do not have to depend on `binary-traits`
// directly. `text-compression` re-exports them from there for the same reason.
pub use text_compression::{LoadIndex, ProteinTextBackend, ReadBinary, ReadBinaryMmap, WriteBinary};

/// Byte placed between consecutive protein sequences in the concatenated text.
///
/// A suffix that would run past the end of its protein hits this instead, which is what stops
/// matches spanning two proteins. Must not appear in any sequence.
pub static SEPARATION_CHARACTER: u8 = b'-';
/// Byte marking the end of the concatenated text.
///
/// It sorts below every residue, so the final suffix orders correctly against its own prefixes,
/// and it gives the tryptic C-terminus predicate a sentinel to read at the last position instead
/// of a special case. It is *not* what keeps the comparison loops in bounds — those check the
/// text length explicitly; see `Searcher::compare` in `sa-index`.
///
/// Must not appear in any sequence.
pub static TERMINATION_CHARACTER: u8 = b'$';

// ── ProteinsBackend trait ─────────────────────────────────────────────────────

/// Common interface for all proteins backends.
///
/// The associated type `Text` is the backend's own generic parameter, so the metadata storage and
/// the text storage are chosen independently — see the module docs. Each backend still has a
/// single, unconditional trait impl; no `#[cfg]` gates needed.
pub trait ProteinsBackend: Send + Sync {
    /// The concrete protein-text type this backend stores.
    type Text: ProteinTextBackend + Sync;
    /// Borrows the concatenated protein text the suffix array indexes.
    fn text(&self) -> &Self::Text;
    /// Number of proteins.
    fn len(&self) -> usize;
    /// Whether there are no proteins.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Returns the metadata for protein `index`.
    ///
    /// # Panics
    ///
    /// If `index` is out of range.
    fn get(&self, index: usize) -> ProteinRef<'_>;

    /// Reads every OS page in the backing store into the page cache.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn touch_all_pages(&self) {}

    /// Non-blocking hardware prefetch for the metadata record at `index`.
    /// What that record is differs per backend; each impl says which bytes it hints.
    /// Default is a no-op; backends override when a prefetch is meaningful.
    #[inline]
    fn prefetch(&self, _index: usize) {}

    /// Non-blocking hardware prefetch for the string data (UID + FA) at `index`.
    /// Default is a no-op; mmap-backed implementations override this.
    #[inline]
    fn prefetch_strings(&self, _index: usize) {}
}

/// One protein's metadata, owned.
///
/// Produced by the TSV loader and held by the preloaded backend. Reading code should prefer
/// [`ProteinRef`], which both backends can hand out without copying.
pub struct Protein {
    /// UniProt accession, e.g. `P12345`.
    pub uniprot_id: String,
    /// NCBI taxon id.
    pub taxon_id: u32,
    /// Functional annotations, encoded with `fa_compression`.
    pub functional_annotations: Vec<u8>
}

/// One protein's metadata, borrowed.
///
/// The common currency of the read path: the preloaded backend borrows from its `Vec<Protein>`
/// and the mmap backend borrows straight from the mapping, so neither copies per result.
#[derive(Clone, Copy)]
pub struct ProteinRef<'a> {
    /// UniProt accession, e.g. `P12345`.
    pub uniprot_id: &'a str,
    /// NCBI taxon id.
    pub taxon_id: u32,
    /// Functional annotations, still encoded; see [`Self::get_functional_annotations`].
    pub functional_annotations: &'a [u8]
}

impl ProteinRef<'_> {
    /// Decodes the functional annotations into their `GO:`/`EC:`/`IPR:` text form.
    ///
    /// Allocates, so it is done once per returned result rather than during search.
    pub fn get_functional_annotations(&self) -> String {
        decode(self.functional_annotations)
    }
}
