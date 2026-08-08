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

#[cfg(test)]
pub(crate) mod test_fixtures {
    //! The shared TSV fixture. Both backends' tests build an index from the same rows so that
    //! their assertions are comparable; each picks how many rows it wants.

    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;

    use tempdir::TempDir;

    /// `(uniprot_id, taxon_id, sequence, annotations)`, in the order they appear in the file.
    pub(crate) const TEST_PROTEINS: [(&str, u32, &str, &str); 4] = [
        ("P12345", 1, "MLPGLALLLLAAWTARALEV", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
        ("P54321", 2, "PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
        ("P67890", 6, "KWDSDPSGTKTCIDT", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
        ("P13579", 17, "KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ];

    /// Writes `proteins` as a UniProt TSV into `tmp_dir` and returns its path.
    pub(crate) fn write_database_file(
        tmp_dir: &TempDir,
        proteins: &[(&str, u32, &str, &str)]
    ) -> PathBuf {
        let path = tmp_dir.path().join("database.tsv");
        let mut f = File::create(&path).unwrap();
        for (uid, taxon, sequence, annotations) in proteins {
            writeln!(f, "{uid}\t{taxon}\t{sequence}\t{annotations}").unwrap();
        }
        path
    }
}

pub use preloaded::InMemoryProteins;
#[cfg(feature = "mmap")]
pub use mmap::MmapBackedProteins;

/// Type alias — resolves to the active backend for this build.
#[cfg(feature = "mmap")]
pub type Proteins = MmapBackedProteins;
/// Type alias — resolves to the active backend for this build.
#[cfg(not(feature = "mmap"))]
pub type Proteins = InMemoryProteins;

// Re-export I/O traits used by callers.
pub use text_compression::{WriteBinary, ReadBinary, ReadBinaryMmap};

// ── Shared types ──────────────────────────────────────────────────────────────

use fa_compression::algorithm1::decode;
pub use text_compression::ProteinTextBackend;

/// Byte placed between consecutive protein sequences in the concatenated text.
///
/// A suffix that would run past the end of its protein hits this instead, which is what stops
/// matches spanning two proteins. Must not appear in any sequence.
pub static SEPARATION_CHARACTER: u8 = b'-';
/// Byte marking the end of the concatenated text.
///
/// Every search relies on this: it guarantees a candidate comparison always terminates before
/// running off the end of the text, so the hot loops need no explicit end-of-text bounds check.
pub static TERMINATION_CHARACTER: u8 = b'$';

// ── ProteinsBackend trait ─────────────────────────────────────────────────────

/// Common interface for all proteins backends.
///
/// The associated type `Text` lets each backend declare its own concrete text
/// type (`InMemoryProteinText` or `MmapBackedProteinText`) while keeping a
/// single, unconditional trait impl per backend — no `#[cfg]` gates needed.
pub trait ProteinsBackend: Send + Sync {
    /// The concrete protein-text type this backend stores.
    type Text: ProteinTextBackend + Sync;
    /// Borrows the concatenated protein text the suffix array indexes.
    fn text(&self) -> &Self::Text;
    /// Number of proteins.
    fn len(&self) -> usize;
    /// Whether there are no proteins.
    fn is_empty(&self) -> bool { self.len() == 0 }
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

    /// Non-blocking hardware prefetch for the fixed-table entry at `index`.
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
    pub functional_annotations: Vec<u8>,
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
    pub functional_annotations: &'a [u8],
}

impl ProteinRef<'_> {
    /// Decodes the functional annotations into their `GO:`/`EC:`/`IPR:` text form.
    ///
    /// Allocates, so it is done once per returned result rather than during search.
    pub fn get_functional_annotations(&self) -> String {
        decode(self.functional_annotations)
    }
}
