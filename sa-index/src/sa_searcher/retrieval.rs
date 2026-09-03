//! Protein retrieval: turning matched suffix positions into protein references.
//!
//! Self-contained phase (no shared search internals), kept in its own `impl` block to
//! separate it from the search machinery in the parent module.
//!
//! # Why this prefetches even in the preloaded build
//!
//! Prefetching reads like an mmap concern — hiding page faults — but the two-pass loop below is
//! just as valuable when everything is in owned RAM, and the reason is that "loaded" is not
//! "cached". Every structure this phase touches is orders of magnitude larger than any L3, on the
//! smaller reference database as much as on the full one:
//!
//! * suffix-to-protein mapping: 4 bytes per text position dense, ~1.25 bits as a bit vector
//! * protein metadata table: the largest structure after the suffix array
//! * protein text: 5 bits per residue
//!
//! Retrieval walks them at positions the suffix array happened to produce, which is effectively
//! random, so nearly every lookup is a last-level miss — ~80-100 ns, whether or not the page is
//! resident. The mmap build additionally risks a page fault, so it gains more, but the preloaded
//! build was leaving the same memory-level parallelism on the table by issuing one dependent load
//! at a time.
//!
//! Hence the shape: pass 1 issues hints for a whole batch, pass 2 consumes it once the loads are
//! in flight. The lookahead distance is `RETRIEVAL_PREFETCH_DISTANCE`.

use protein_metadata::{ProteinRef, ProteinsBackend};

use super::{RETRIEVAL_PREFETCH_DISTANCE, Searcher};
use crate::{Nullable, array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Returns all the proteins that correspond with the provided suffixes.
    ///
    /// Two-pass prefetch pipeline, look-ahead D = `RETRIEVAL_PREFETCH_DISTANCE`:
    /// Pass 1 — prefetch suffix_to_protein mapping entries D iterations ahead, collect protein_indices.
    /// Pass 2 — prefetch protein entries D iterations ahead, build ProteinRef result.
    ///
    /// D is 32 → D/2 iterations × ~5 ns ≈ 80–100 ns gap before the protein read in
    /// `proteins.get()`, giving the prefetch hint time to complete for most DRAM configs.
    ///
    /// Note: prefetch_strings is intentionally omitted — it reads the fixed-table entry to obtain
    /// string offsets, which causes a stall when the entry has not yet landed from the earlier
    /// prefetch hint (D/2 iterations × ~5 ns < ~80–100 ns DRAM latency).
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &[i64]) -> Vec<ProteinRef<'_>> {
        let distance = RETRIEVAL_PREFETCH_DISTANCE;

        // Pass 1: prefetch suffix_to_protein mapping, collect protein_indices
        let mut protein_indices = Vec::with_capacity(suffixes.len());
        for (i, &suffix) in suffixes.iter().enumerate() {
            if let Some(&fs) = suffixes.get(i + distance) {
                self.suffix_index_to_protein.prefetch_for_suffix(fs);
            }
            protein_indices.push(self.suffix_index_to_protein.suffix_to_protein(suffix));
        }

        // Pass 2: prefetch proteins (D ahead), build ProteinRefs
        let mut res = Vec::with_capacity(suffixes.len());
        for (i, &protein_index) in protein_indices.iter().enumerate() {
            if let Some(&fpi) = protein_indices.get(i + distance)
                && !fpi.is_null()
            {
                self.proteins.prefetch(fpi as usize);
            }
            if !protein_index.is_null() {
                res.push(self.proteins.get(protein_index as usize));
            }
        }
        res
    }

    /// Returns the taxon id of every protein that corresponds with the provided suffixes.
    ///
    /// Same two-pass pipeline as [`Self::retrieve_proteins`] and for the same reason — see the
    /// module docs; only the second pass differs, reading one field per protein instead of
    /// building a `ProteinRef`.
    ///
    /// The prefetch pays off more completely here than it does there. `proteins.prefetch` hints
    /// the metadata entry, which is the whole of what this reads; `retrieve_proteins` goes on to
    /// follow that entry's offsets into the UID and annotation blobs, which the same hint does not
    /// cover.
    ///
    /// Duplicates are kept. Suffixes commonly land in the same protein, and far more often in the
    /// same taxon, but deduplicating is the caller's business — `retrieve_proteins` does no
    /// post-processing either, and `peptide_search` is where a result is shaped.
    #[inline]
    pub fn retrieve_taxa(&self, suffixes: &[i64]) -> Vec<u32> {
        let distance = RETRIEVAL_PREFETCH_DISTANCE;

        // Pass 1: prefetch suffix_to_protein mapping, collect protein_indices
        let mut protein_indices = Vec::with_capacity(suffixes.len());
        for (i, &suffix) in suffixes.iter().enumerate() {
            if let Some(&fs) = suffixes.get(i + distance) {
                self.suffix_index_to_protein.prefetch_for_suffix(fs);
            }
            protein_indices.push(self.suffix_index_to_protein.suffix_to_protein(suffix));
        }

        // Pass 2: prefetch protein entries (D ahead), read taxon ids
        let mut res = Vec::with_capacity(suffixes.len());
        for (i, &protein_index) in protein_indices.iter().enumerate() {
            if let Some(&fpi) = protein_indices.get(i + distance)
                && !fpi.is_null()
            {
                self.proteins.prefetch(fpi as usize);
            }
            if !protein_index.is_null() {
                res.push(self.proteins.taxon_id(protein_index as usize));
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use crate::sa_searcher::{
        SearchAllSuffixesResult,
        test_utils::{example_searcher, repeated_residue_searcher}
    };

    #[test]
    fn test_retrieve_proteins() {
        // Proteins (distinct taxa): AI=10, CLACVAA=20, AC=30, KCRLY=40.
        let searcher = example_searcher();

        // "A" matches 5 suffixes: once in protein 0 (taxon 10), three times in protein 1
        // (taxon 20), once in protein 2 (taxon 30); never in protein 3 (KCRLY).
        let suffixes = match searcher.search_matching_suffixes_scalar(b"A", usize::MAX, false, false) {
            SearchAllSuffixesResult::SearchResult(s) | SearchAllSuffixesResult::MaxMatches(s) => s,
            SearchAllSuffixesResult::NoMatches => vec![]
        };

        let found = searcher.retrieve_proteins(&suffixes);
        assert_eq!(found.len(), suffixes.len(), "one protein per matched suffix");

        let mut taxa: Vec<u32> = found.iter().map(|p| p.taxon_id).collect();
        taxa.sort();
        assert_eq!(taxa, vec![10, 20, 20, 20, 30]);
    }

    #[test]
    fn test_retrieve_proteins_empty() {
        let searcher = example_searcher();

        assert!(searcher.retrieve_proteins(&[]).is_empty());
    }

    // > RETRIEVAL_PREFETCH_DISTANCE (32) suffixes, to exercise the two-pass prefetch-ahead loop.
    #[test]
    fn test_retrieve_proteins_many() {
        let n = 70usize;
        let searcher = repeated_residue_searcher('A', n);

        let suffixes = match searcher.search_matching_suffixes_scalar(b"A", usize::MAX, false, false) {
            SearchAllSuffixesResult::SearchResult(s) | SearchAllSuffixesResult::MaxMatches(s) => s,
            SearchAllSuffixesResult::NoMatches => vec![]
        };
        assert_eq!(suffixes.len(), n);
        let found = searcher.retrieve_proteins(&suffixes);
        assert_eq!(found.len(), n);
        assert!(found.iter().all(|p| p.taxon_id == 10));
    }

    // A separator position maps to u32::NULL and must be skipped (not returned).
    #[test]
    fn test_retrieve_proteins_skips_separators() {
        let searcher = example_searcher();

        // position 0 = protein 0 ('A', taxon 10); position 2 = separator ('-', null).
        let found = searcher.retrieve_proteins(&[0, 2]);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].taxon_id, 10);
    }
}
