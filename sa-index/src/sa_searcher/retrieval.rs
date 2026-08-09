//! Protein retrieval: turning matched suffix positions into protein references.
//!
//! Self-contained phase (no shared search internals), kept in its own `impl` block to
//! separate it from the search machinery in the parent module.
//!
//! # Why this prefetches even in the preloaded build
//!
//! Prefetching reads like an mmap concern — hiding page faults — but the two-pass loop below is
//! just as valuable when everything is in owned RAM, and the reason is that "loaded" is not
//! "cached". At UniProt scale the structures this phase touches are orders of magnitude larger
//! than any L3:
//!
//! * suffix-to-protein mapping: ~1.2 GB dense, a few hundred MB as a bit vector
//! * protein metadata table: gigabytes
//! * protein text: ~190 MB packed
//!
//! Retrieval walks them at positions the suffix array happened to produce, which is effectively
//! random, so nearly every lookup is a last-level miss — ~80-100 ns, whether or not the page is
//! resident. The mmap build additionally risks a page fault, so it gains more, but the preloaded
//! build was leaving the same memory-level parallelism on the table by issuing one dependent load
//! at a time.
//!
//! Hence the shape: pass 1 issues hints for a whole batch, pass 2 consumes it once the loads are
//! in flight. The lookahead distance is `SearchTuning::retrieval_prefetch_distance`.

use sa_mappings::proteins::{ProteinRef, ProteinsBackend};

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
use crate::Nullable;

use super::Searcher;

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Returns all the proteins that correspond with the provided suffixes.
    ///
    /// Two-pass prefetch pipeline, look-ahead D = `tuning.retrieval_prefetch_distance`:
    /// Pass 1 — prefetch suffix_to_protein mapping entries D iterations ahead, collect protein_indices.
    /// Pass 2 — prefetch protein entries D iterations ahead, build ProteinRef result.
    ///
    /// D defaults to 32 → D/2 iterations × ~5 ns ≈ 80–100 ns gap before the protein read in
    /// `proteins.get()`, giving the prefetch hint time to complete for most DRAM configs.
    ///
    /// Note: prefetch_strings is intentionally omitted — it reads the fixed-table entry to obtain
    /// string offsets, which causes a stall when the entry has not yet landed from the earlier
    /// prefetch hint (D/2 iterations × ~5 ns < ~80–100 ns DRAM latency).
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &[i64]) -> Vec<ProteinRef<'_>> {
        let distance = self.tuning.retrieval_prefetch_distance;

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
            if let Some(&fpi) = protein_indices.get(i + distance) {
                if !fpi.is_null() { self.proteins.prefetch(fpi as usize); }
            }
            if !protein_index.is_null() {
                res.push(self.proteins.get(protein_index as usize));
            }
        }
        res
    }

}

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
    use text_compression::ProteinText;

    use crate::{
        array::OriginalSA,
        sa_searcher::{test_helpers::get_example_proteins, SearchAllSuffixesResult, Searcher},
        suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping},
        SuffixArray,
    };

    #[test]
    fn test_retrieve_proteins() {
        // Proteins (distinct taxa): AI=10, CLACVAA=20, AC=30, KCRLY=40.
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18],
            1,
        ));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // "A" matches 5 suffixes: once in protein 0 (taxon 10), three times in protein 1
        // (taxon 20), once in protein 2 (taxon 30); never in protein 3 (KCRLY).
        let suffixes = match searcher.search_matching_suffixes(b"A", usize::MAX, false, false) {
            SearchAllSuffixesResult::SearchResult(s) | SearchAllSuffixesResult::MaxMatches(s) => s,
            SearchAllSuffixesResult::NoMatches => vec![],
        };

        let found = searcher.retrieve_proteins(&suffixes);
        assert_eq!(found.len(), suffixes.len(), "one protein per matched suffix");

        let mut taxa: Vec<u32> = found.iter().map(|p| p.taxon_id).collect();
        taxa.sort();
        assert_eq!(taxa, vec![10, 20, 20, 20, 30]);
    }

    #[test]
    fn test_retrieve_proteins_empty() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        assert!(searcher.retrieve_proteins(&[]).is_empty());
    }

    // > the default retrieval_prefetch_distance (32) suffixes to exercise the two-pass prefetch-ahead loop.
    #[test]
    fn test_retrieve_proteins_many() {
        let n = 70usize;
        let mut input = "A".repeat(n);
        input.push('$');
        let text = ProteinText::from_string(&input);
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 5,
            functional_annotations: vec![],
        }]);
        // SA of A^n$ is [n, n-1, …, 0].
        let sa = SuffixArray::Original(OriginalSA((0..=n as i64).rev().collect(), 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        let suffixes = match searcher.search_matching_suffixes(b"A", usize::MAX, false, false) {
            SearchAllSuffixesResult::SearchResult(s) | SearchAllSuffixesResult::MaxMatches(s) => s,
            SearchAllSuffixesResult::NoMatches => vec![],
        };
        assert_eq!(suffixes.len(), n);
        let found = searcher.retrieve_proteins(&suffixes);
        assert_eq!(found.len(), n);
        assert!(found.iter().all(|p| p.taxon_id == 5));
    }

    // A separator position maps to u32::NULL and must be skipped (not returned).
    #[test]
    fn test_retrieve_proteins_skips_separators() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // position 0 = protein 0 ('A', taxon 10); position 2 = separator ('-', null).
        let found = searcher.retrieve_proteins(&[0, 2]);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].taxon_id, 10);
    }
}
