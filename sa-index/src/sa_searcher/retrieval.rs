//! Protein retrieval: turning matched suffix positions into protein references.
//!
//! Self-contained phase (no shared search internals), kept in its own `impl` block to
//! separate it from the search machinery in the parent module.

use sa_mappings::proteins::{ProteinRef, ProteinsBackend};

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
use crate::Nullable;

use super::Searcher;

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Returns all the proteins that correspond with the provided suffixes.
    ///
    /// Two-pass prefetch pipeline (PREFETCH_DISTANCE = 32):
    /// Pass 1 — prefetch suffix_to_protein mapping entries D iterations ahead, collect protein_indices.
    /// Pass 2 — prefetch protein entries D iterations ahead, build ProteinRef result.
    ///
    /// Note: prefetch_strings is intentionally omitted — it reads the fixed-table entry to obtain
    /// string offsets, which causes a stall when the entry has not yet landed from the earlier
    /// prefetch hint (D/2 iterations × ~5 ns < ~80–100 ns DRAM latency).
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &[i64]) -> Vec<ProteinRef<'_>> {
        // D=32 → D/2 iterations × ~5 ns ≈ 80–100 ns gap before the protein read in
        // proteins.get(), giving the prefetch hint time to complete for most DRAM configs.
        const PREFETCH_DISTANCE: usize = 32;

        // Pass 1: prefetch suffix_to_protein mapping, collect protein_indices
        let mut protein_indices = Vec::with_capacity(suffixes.len());
        for (i, &suffix) in suffixes.iter().enumerate() {
            if let Some(&fs) = suffixes.get(i + PREFETCH_DISTANCE) {
                self.suffix_index_to_protein.prefetch_for_suffix(fs);
            }
            protein_indices.push(self.suffix_index_to_protein.suffix_to_protein(suffix));
        }

        // Pass 2: prefetch proteins (D ahead), build ProteinRefs
        let mut res = Vec::with_capacity(suffixes.len());
        for (i, &protein_index) in protein_indices.iter().enumerate() {
            if let Some(&fpi) = protein_indices.get(i + PREFETCH_DISTANCE) {
                if !fpi.is_null() { self.proteins.prefetch(fpi as usize); }
            }
            if !protein_index.is_null() {
                res.push(self.proteins.get(protein_index as usize));
            }
        }
        res
    }
}
