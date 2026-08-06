//! Multi-peptide search orchestration.
//!
//! The searcher exposes two primitives: `search_matching_suffixes` (scalar, one peptide) and
//! `search_matching_suffixes_batched` (a chunk of peptides, with the binary-search phase
//! interleaved for memory-level parallelism). This module owns the one piece of glue every
//! consumer needs: take a whole list of peptides, split it into MLP batches, and run them in
//! parallel across rayon. Both the production path (`peptide_search::search_all_peptides`) and
//! the benchmark call this, so batching lives in one place and is measured as it ships.

use rayon::prelude::*;
use sa_mappings::proteins::ProteinsBackend;

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;

use super::{SearchAllSuffixesResult, Searcher};

/// Default cross-query MLP batch size. Chosen from the full-DB sweep: batching is a win on
/// long peptides (search-bound) and neutral on short/retrieval-bound ones; the gain plateaus
/// at the 8–16 knee, and 32/64 start regressing short peptides. 16 takes the whole knee
/// without hurting any length or either backend. A single peptide degrades to scalar work.
pub const DEFAULT_MLP_BATCH: usize = 16;

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Searches every peptide, returning one raw suffix result per input, in order.
    ///
    /// `batch_size` selects the strategy: `> 1` interleaves that many peptides per rayon task
    /// via `search_matching_suffixes_batched` (memory-level parallelism); `1` runs the scalar
    /// `search_matching_suffixes` one peptide per task. Both produce identical results
    /// (see `test_batched_matches_scalar`), so `batch_size` is a pure performance knob.
    ///
    /// Peptides are searched as given — callers that need length filtering or normalisation
    /// (e.g. skipping peptides shorter than the sample rate) must do it before calling.
    pub fn search_all_matching_suffixes(
        &self,
        peptides: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool,
        batch_size: usize,
    ) -> Vec<SearchAllSuffixesResult> {
        if batch_size > 1 {
            peptides
                .par_chunks(batch_size)
                .flat_map(|chunk| self.search_matching_suffixes_batched(chunk, max_matches, equate_il, tryptic))
                .collect()
        } else {
            peptides
                .par_iter()
                .map(|&p| self.search_matching_suffixes(p, max_matches, equate_il, tryptic))
                .collect()
        }
    }
}

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use sa_mappings::proteins::ProteinsBackend as _;

    use crate::{
        array::OriginalSA,
        sa_searcher::{test_helpers::get_example_proteins, SearchAllSuffixesResult, Searcher},
        suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping},
        SuffixArray,
    };

    // SA of the example text "AI-CLACVAA-AC-KCRLY$" (L→I normalised is irrelevant here: no
    // kmer table attached, so raw suffix search is used).
    fn example_searcher() -> Searcher<SuffixArray> {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp))
    }

    /// Orchestrated multi-peptide search must be identical whether scalar (batch 1) or batched,
    /// across the equate_il × tryptic combinations, and independent of batch size.
    #[test]
    fn test_search_all_matches_scalar() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"KCR", b"VAA", b"ZZZ", b"AI"];
        for equate_il in [true, false] {
            for tryptic in [true, false] {
                let scalar = s.search_all_matching_suffixes(&peptides, 1000, equate_il, tryptic, 1);
                for batch in [2usize, 3, 16] {
                    let batched = s.search_all_matching_suffixes(&peptides, 1000, equate_il, tryptic, batch);
                    assert_eq!(scalar.len(), peptides.len());
                    assert_eq!(batched.len(), peptides.len(), "one result per input, in order");
                    for (i, (a, b)) in scalar.iter().zip(batched.iter()).enumerate() {
                        assert_eq!(a, b, "peptide {} differs at batch {} (il={} tr={})", i, batch, equate_il, tryptic);
                    }
                }
            }
        }
    }

    /// max_matches cutoff must propagate through the orchestrator identically for both strategies.
    #[test]
    fn test_search_all_respects_cutoff() {
        let s = example_searcher();
        // "A" matches 5 suffixes; cap at 2 → MaxMatches with 2 entries on both paths.
        let peptides: Vec<&[u8]> = vec![b"A"];
        let scalar = s.search_all_matching_suffixes(&peptides, 2, true, false, 1);
        let batched = s.search_all_matching_suffixes(&peptides, 2, true, false, 16);
        assert!(matches!(scalar[0], SearchAllSuffixesResult::MaxMatches(_)));
        assert_eq!(scalar[0], batched[0]);
    }

    #[test]
    fn test_search_all_empty() {
        let s = example_searcher();
        let none: Vec<&[u8]> = vec![];
        assert!(s.search_all_matching_suffixes(&none, 1000, true, false, 1).is_empty());
        assert!(s.search_all_matching_suffixes(&none, 1000, true, false, 16).is_empty());
    }

    // A batch size that does not divide the input length must still return one result per
    // peptide in the original order (exercises the par_chunks remainder).
    #[test]
    fn test_search_all_ragged_batch() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"AI", b"VAA"]; // 5 peptides
        let scalar = s.search_all_matching_suffixes(&peptides, 1000, true, false, 1);
        let batched = s.search_all_matching_suffixes(&peptides, 1000, true, false, 4); // 4 + 1
        assert_eq!(scalar, batched);
    }
}
