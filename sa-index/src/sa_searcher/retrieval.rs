//! Protein retrieval: turning matched suffix positions into protein references.
//!
//! Self-contained phase (no shared search internals), kept in its own `impl` block to
//! separate it from the search machinery in the parent module.

use sa_mappings::proteins::{ProteinRef, ProteinsBackend};

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;
use crate::Nullable;

use super::Searcher;

/// Default number of queries whose retrieval is interleaved by `retrieve_proteins_batched`.
/// 16 queries × the observed average of a few matched suffixes each lands at or above the
/// default `retrieval_prefetch_distance` of 32, which is exactly the point at which the
/// look-ahead starts firing for the peptide lengths that never reached it before (26–50 aa
/// peptides match only a handful of suffixes, so single-query retrieval prefetched nothing at
/// all). Larger groups only add scratch-buffer pressure — the pipeline is already full at that
/// depth.
///
/// Mirrored by `SearchTuning::retrieval_batch`, which is what `retrieve_all_proteins` actually
/// reads; this constant is the default that field is initialised to.
pub const DEFAULT_RETRIEVAL_BATCH: usize = 16;

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

    /// Retrieves proteins for many queries at once, interleaving the random reads across
    /// queries so the prefetch pipeline stays full even when individual result lists are
    /// far shorter than the prefetch distance. Returns one Vec per input, in order —
    /// identical to calling `retrieve_proteins` on each input separately.
    ///
    /// Why this exists: `retrieve_proteins` guards its look-ahead on `suffixes.get(i + D)`,
    /// so a query matching fewer than D (32 by default) suffixes issues *no* prefetch at all.
    /// That is the common case for long peptides (26–50 aa match a handful of suffixes each)
    /// and exactly
    /// the case where the two random reads per suffix — the suffix→protein mapping entry and
    /// the protein entry — are pure ~80–100 ns DRAM latency with nothing to overlap them
    /// with. Flattening `batch` queries into one index makes the look-ahead cross query
    /// boundaries, so ~D misses are in flight regardless of how short each list is. Same
    /// trick as `search_matching_suffixes_batched`, applied to the retrieval phase.
    ///
    /// As in `retrieve_proteins`, `prefetch_strings` is deliberately not issued: it reads the
    /// fixed-table entry to obtain the string offsets, which stalls whenever that entry has
    /// not yet landed from the earlier `prefetch` hint.
    ///
    /// Null protein indices (separator positions) are skipped, not returned, so an output
    /// Vec can be shorter than its input slice — identical to the single-query path.
    pub fn retrieve_proteins_batched(&self, suffix_lists: &[&[i64]], batch: usize) -> Vec<Vec<ProteinRef<'_>>> {
        // batch == 0 would make `chunks` panic; a group of one degrades to the single-query
        // behaviour (look-ahead only fires for lists longer than D), which is still correct.
        let batch = batch.max(1);
        let distance = self.tuning.retrieval_prefetch_distance;

        let mut results: Vec<Vec<ProteinRef<'_>>> = Vec::with_capacity(suffix_lists.len());

        // Scratch buffers are hoisted out of the group loop: the flattened suffixes, the
        // exclusive-prefix offsets that map a flat index back to its query, and the protein
        // indices produced by pass 1. Reusing them means one allocation for the whole call
        // instead of three per group, and the buffers stay warm in cache between groups.
        let mut flat_suffixes: Vec<i64> = Vec::new();
        let mut offsets: Vec<usize> = Vec::new();
        let mut protein_indices: Vec<u32> = Vec::new();

        for group in suffix_lists.chunks(batch) {
            // Flatten, so pass 1 has one contiguous array to look ahead in. The copy is
            // sequential (~1 ns per i64, hardware-prefetched) against two random reads per
            // suffix at ~80–100 ns each, so it stays well below the noise floor of the reads
            // it enables.
            flat_suffixes.clear();
            offsets.clear();
            offsets.push(0);
            for suffixes in group {
                flat_suffixes.extend_from_slice(suffixes);
                offsets.push(flat_suffixes.len());
            }

            // Pass 1: prefetch suffix_to_protein D ahead across the whole group, collect
            // protein_indices. The look-ahead now spills over into the next query's suffixes.
            protein_indices.clear();
            protein_indices.reserve(flat_suffixes.len());
            for (i, &suffix) in flat_suffixes.iter().enumerate() {
                if let Some(&fs) = flat_suffixes.get(i + distance) {
                    self.suffix_index_to_protein.prefetch_for_suffix(fs);
                }
                protein_indices.push(self.suffix_index_to_protein.suffix_to_protein(suffix));
            }

            // Pass 2: prefetch proteins (D ahead, again across query boundaries) and split the
            // ProteinRefs back out per query using the offsets recorded above.
            for query in 0..group.len() {
                let (start, end) = (offsets[query], offsets[query + 1]);
                let mut res = Vec::with_capacity(end - start);
                for i in start..end {
                    if let Some(&fpi) = protein_indices.get(i + distance) {
                        if !fpi.is_null() { self.proteins.prefetch(fpi as usize); }
                    }
                    let protein_index = protein_indices[i];
                    if !protein_index.is_null() {
                        res.push(self.proteins.get(protein_index as usize));
                    }
                }
                results.push(res);
            }
        }

        results
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

    fn example_searcher() -> Searcher<SuffixArray> {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp))
    }

    // Compares the batched retrieval against the per-query reference. ProteinRef holds borrowed
    // fields, so identity is checked on the observable payload (taxon + accession).
    fn assert_batched_matches_single(searcher: &Searcher<SuffixArray>, suffix_lists: &[&[i64]], batch: usize) {
        let expected: Vec<Vec<(u32, String)>> = suffix_lists
            .iter()
            .map(|suffixes| {
                searcher.retrieve_proteins(suffixes).iter().map(|p| (p.taxon_id, p.uniprot_id.to_string())).collect()
            })
            .collect();

        let got: Vec<Vec<(u32, String)>> = searcher
            .retrieve_proteins_batched(suffix_lists, batch)
            .iter()
            .map(|proteins| proteins.iter().map(|p| (p.taxon_id, p.uniprot_id.to_string())).collect())
            .collect();

        assert_eq!(got.len(), suffix_lists.len(), "one Vec per input, in order (batch {})", batch);
        assert_eq!(got, expected, "batched retrieval differs from per-query retrieval at batch {}", batch);
    }

    // Result lists of wildly varying length — empty, length 1, well below the default retrieval_prefetch_distance (32)
    // and well above it — must give exactly the per-query result, independent of the batch size.
    #[test]
    fn test_retrieve_batched_matches_single() {
        let n = 70usize;
        let mut input = "A".repeat(n);
        input.push('$');
        let text = ProteinText::from_string(&input);
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: "P0".to_string(),
            taxon_id: 5,
            functional_annotations: vec![],
        }]);
        let sa = SuffixArray::Original(OriginalSA((0..=n as i64).rev().collect(), 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        let long: Vec<i64> = (0..n as i64).collect(); // 70 > default distance
        let medium: Vec<i64> = (0..20).collect(); // < default distance
        let short = vec![3i64];
        let empty: Vec<i64> = vec![];

        let suffix_lists: Vec<&[i64]> = vec![
            &empty, &short, &medium, &long, &empty, &empty, &short, &long, &medium, &short, &long, &empty,
        ];

        for batch in [1usize, 2, 16] {
            assert_batched_matches_single(&searcher, &suffix_lists, batch);
        }

        // Sanity check that the reference itself is non-trivial: the long list yields n proteins.
        let got = searcher.retrieve_proteins_batched(&suffix_lists, 16);
        assert!(got[0].is_empty());
        assert_eq!(got[1].len(), 1);
        assert_eq!(got[2].len(), 20);
        assert_eq!(got[3].len(), n);
    }

    // Separator positions map to u32::NULL and must be skipped by the batched path too, so an
    // output Vec is shorter than its input slice — exactly as in the single-query path.
    #[test]
    fn test_retrieve_batched_skips_separators() {
        let searcher = example_searcher();

        // Text "AI-CLACVAA-AC-KCRLY$": separators at 2, 10 and 13.
        let with_seps = vec![0i64, 2, 3, 10, 13, 14];
        let only_seps = vec![2i64, 10, 13];
        let no_seps = vec![0i64, 1, 4];
        let suffix_lists: Vec<&[i64]> = vec![&with_seps, &only_seps, &no_seps, &only_seps];

        for batch in [1usize, 2, 16] {
            assert_batched_matches_single(&searcher, &suffix_lists, batch);
        }

        let got = searcher.retrieve_proteins_batched(&suffix_lists, 16);
        // [0, 2, 3, 10, 13, 14]: positions 2, 10 and 13 are separators, so only 0, 3 and 14 return.
        assert_eq!(got[0].len(), 3);
        assert!(got[1].is_empty(), "all-separator list retrieves nothing");
        assert_eq!(got[2].len(), 3);
    }

    // A number of queries that is not a multiple of the batch size must still return one Vec per
    // input, in order (exercises the short trailing chunk).
    #[test]
    fn test_retrieve_batched_ragged() {
        let searcher = example_searcher();

        let a = vec![0i64, 1];
        let b = vec![3i64, 4, 5, 6];
        let c: Vec<i64> = vec![];
        let d = vec![14i64];
        let e = vec![11i64, 12];
        let suffix_lists: Vec<&[i64]> = vec![&a, &b, &c, &d, &e]; // 5 queries

        for batch in [2usize, 3, 4] {
            // 2+2+1, 3+2, 4+1
            assert_batched_matches_single(&searcher, &suffix_lists, batch);
        }

        assert!(searcher.retrieve_proteins_batched(&[], 16).is_empty());
    }
}
