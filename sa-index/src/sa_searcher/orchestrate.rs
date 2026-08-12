//! Multi-peptide search orchestration.
//!
//! The searcher exposes two search primitives: `search_matching_suffixes` (scalar, one
//! peptide) and `search_matching_suffixes_batched` (a chunk of peptides, with the
//! binary-search phase interleaved for memory-level parallelism). This module owns the one
//! piece of glue every consumer needs: take a whole list of peptides, split it into MLP
//! batches, and run them in parallel across rayon. Both the production path
//! (`peptide_search::search_all_peptides`) and the benchmark call this, so batching lives in
//! one place and is measured as it ships.
//!
//! Retrieval deliberately has no counterpart here. A cross-query batched retrieval was built
//! and measured (run3): the `retrieval_batch` sweep moved throughput by a median of +1.7%,
//! never clearing the 3.9% noise floor in any of 12 (file, backend, baseline) combinations —
//! including long peptides, the case it was designed for. It was reverted rather than carried.
//! `retrieve_proteins`' own two-pass prefetch already covers queries with many matches, and
//! queries with few matches are dominated by peptides where retrieval is a small share of the
//! work.

use rayon::prelude::*;
use sa_mappings::proteins::ProteinsBackend;

use super::{SearchAllSuffixesResult, Searcher};
use crate::{array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

/// Default cross-query MLP batch size. Chosen from the full-DB sweep: batching is a win on
/// long peptides (search-bound) and neutral on short/retrieval-bound ones; the gain plateaus
/// at the 8–16 knee, and 32/64 start regressing short peptides. 16 takes the whole knee
/// without hurting any length or either backend. A single peptide degrades to scalar work.
pub const DEFAULT_MLP_BATCH: usize = 16;

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Searches every peptide, returning one raw suffix result per input, in order.
    ///
    /// `tuning.mlp_batch` selects the strategy: `> 1` interleaves that many peptides per rayon
    /// task via `search_matching_suffixes_batched` (memory-level parallelism); `1` runs the scalar
    /// `search_matching_suffixes` one peptide per task. Both produce identical results
    /// (see `test_batched_matches_scalar`), so it is a pure performance knob — which is why it
    /// lives in [`SearchTuning`](super::SearchTuning) with the others rather than being threaded through as an
    /// argument every caller has to supply and no caller ever varied.
    ///
    /// Peptides are searched as given — callers that need length filtering or normalisation
    /// (e.g. skipping peptides shorter than the sample rate) must do it before calling.
    pub fn search_all_matching_suffixes(
        &self,
        peptides: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> Vec<SearchAllSuffixesResult> {
        let batch_size = self.tuning.mlp_batch;
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
    use crate::{
        SuffixArray,
        sa_searcher::{SearchAllSuffixesResult, Searcher, test_helpers::example_searcher}
    };

    /// A copy of `searcher` with `mlp_batch` set, so a test can compare batch sizes.
    ///
    /// The searcher is not cloneable and the tests only need to read it, so this takes it by
    /// reference and returns a view with one field of the tuning changed. `SearchTuning` is `Copy`,
    /// which is what makes that a two-line helper rather than a rebuild of the index.
    fn with_batch(searcher: &Searcher<SuffixArray>, mlp_batch: usize) -> Searcher<SuffixArray> {
        let mut copy = example_searcher();
        copy.tuning = searcher.tuning;
        copy.tuning.mlp_batch = mlp_batch;
        copy
    }

    /// Orchestrated multi-peptide search must be identical whether scalar (batch 1) or batched,
    /// across the equate_il × tryptic combinations, and independent of batch size.
    #[test]
    fn test_search_all_matches_scalar() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"KCR", b"VAA", b"ZZZ", b"AI"];
        for equate_il in [true, false] {
            for tryptic in [true, false] {
                let scalar = with_batch(&s, 1).search_all_matching_suffixes(&peptides, 1000, equate_il, tryptic);
                for batch in [2usize, 3, 16] {
                    let batched =
                        with_batch(&s, batch).search_all_matching_suffixes(&peptides, 1000, equate_il, tryptic);
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
        let scalar = with_batch(&s, 1).search_all_matching_suffixes(&peptides, 2, true, false);
        let batched = with_batch(&s, 16).search_all_matching_suffixes(&peptides, 2, true, false);
        assert!(matches!(scalar[0], SearchAllSuffixesResult::MaxMatches(_)));
        assert_eq!(scalar[0], batched[0]);
    }

    #[test]
    fn test_search_all_empty() {
        let s = example_searcher();
        let none: Vec<&[u8]> = vec![];
        assert!(with_batch(&s, 1).search_all_matching_suffixes(&none, 1000, true, false).is_empty());
        assert!(with_batch(&s, 16).search_all_matching_suffixes(&none, 1000, true, false).is_empty());
    }

    // A batch size that does not divide the input length must still return one result per
    // peptide in the original order (exercises the par_chunks remainder).
    #[test]
    fn test_search_all_ragged_batch() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"AI", b"VAA"]; // 5 peptides
        let scalar = with_batch(&s, 1).search_all_matching_suffixes(&peptides, 1000, true, false);
        let batched = with_batch(&s, 4).search_all_matching_suffixes(&peptides, 1000, true, false); // 4 + 1
        assert_eq!(scalar, batched);
    }
}
