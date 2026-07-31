//! Scalar (one-peptide-at-a-time) search pipeline: narrow the SA window with a binary
//! search, then collect matching suffixes. The memory-level-parallel counterpart lives in
//! `super::batched`; both call the shared primitives (`compare`, `iterate_sa_range`,
//! `prefetch_kmer_range`) that stay in the parent module.

use std::cmp::min;
use std::sync::atomic::Ordering;
use std::time::Instant;

use sa_mappings::proteins::ProteinsBackend;

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;

use super::BoundSearch::{Maximum, Minimum};
use super::{BoundSearch, BoundSearchResult, SearchAllSuffixesResult, Searcher};

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Binary search within the SA window `[left, right)` for the minimum or maximum bound
    /// of `search_string`.
    ///
    /// `lcp_skip` is the number of leading characters known to match for every entry in the
    /// window (e.g. from a k-mer table lookup). Both LCP accumulators start at `lcp_skip` so
    /// those characters are never re-compared.
    ///
    /// Returns `(found, bound_index)` where `found` is true if any suffix matched the full
    /// search string, and `bound_index` is the min (inclusive) or max (inclusive) SA index.
    fn binary_search_bound_in_range(
        &self,
        bound: BoundSearch,
        search_string: &[u8],
        left: usize,
        right: usize,
        lcp_skip: usize,
    ) -> (bool, usize) {
        let mut lo = left;
        let mut hi = right;
        let mut lcp_left: usize = lcp_skip;
        let mut lcp_right: usize = lcp_skip;
        let mut found = false;

        while hi - lo > 1 {
            let center = (lo + hi) / 2;
            self.prefetch_binary_search_pivots(lo, center, hi);
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(center), skip, bound);

            found |= lcp_center == search_string.len();

            if (retval && bound == Minimum) || (!retval && bound == Maximum) {
                hi = center;
                lcp_right = lcp_center;
            } else {
                lo = center;
                lcp_left = lcp_center;
            }
        }

        // Handle the edge case where the initial left boundary was never a center comparison:
        // when the window narrowed to [left, left+1), we must check whether the minimum bound
        // should be at `left` rather than `left+1`.
        if hi == left + 1 && lo == left {
            let (retval, lcp_center) = self.compare(
                search_string, self.sa.get(lo), min(lcp_left, lcp_right), bound,
            );

            found |= lcp_center == search_string.len();

            if bound == Minimum && retval {
                hi = lo;
            }
        }

        match bound {
            Minimum => (found, hi),
            Maximum => (found, lo),
        }
    }

    /// Searches for the minimum and maximum bound for a string in the suffix array.
    ///
    /// When a k-mer table is attached, the first `k` characters of `search_string` are used
    /// to narrow the binary search window before running the search, reducing random memory
    /// accesses by ~60 %.
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// Returns the minimum and maximum bound of all matches in the suffix array, or `NoMatches` if
    /// no matches were found
    pub fn search_bounds(&self, search_string: &[u8]) -> BoundSearchResult {
        let full_range = (0, self.sa.len(), 0);
        let (left, right, lcp_skip) = match &self.kmer_table {
            Some(table) if search_string.len() >= table.k => {
                match table.lookup(&search_string[..table.k]) {
                    Some((lo, hi)) => (lo, hi + 1, table.k),
                    None => return BoundSearchResult::NoMatches,
                }
            }
            _ => full_range,
        };

        let (found_min, min_bound) =
            self.binary_search_bound_in_range(Minimum, search_string, left, right, lcp_skip);

        if !found_min {
            return BoundSearchResult::NoMatches;
        }

        let (_, max_bound) =
            self.binary_search_bound_in_range(Maximum, search_string, left, right, lcp_skip);

        BoundSearchResult::SearchResult((min_bound, max_bound + 1))
    }

    /// Searches for the suffixes matching a search string
    /// During search I and L can be equated
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    /// * `max_matches` - The maximum amount of matches processed, if more matches are found we
    ///   don't process them
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    /// * `tryptic` - Boolean indicating if we only want tryptic matches.
    ///
    /// # Returns
    ///
    /// Returns all the matching suffixes
    #[inline]
    pub fn search_matching_suffixes(
        &self,
        search_string: &[u8],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> SearchAllSuffixesResult {
        self.prefetch_kmer_range(search_string);

        // Cap pre-allocation at 4096 entries (32 KB) so callers passing large max_matches
        // don't wastefully over-allocate for peptides that match rarely.
        let mut matching_suffixes: Vec<i64> = Vec::with_capacity(max_matches.min(4096));
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        let mut skip: usize = 0;
        while skip < self.sa.sample_rate() as usize {
            // il_locations is built in ascending index order, so partition_point gives us
            // the first position that is relevant for this skip value in O(log n).
            // These are still absolute positions within `search_string`, not within the suffix.
            let il_locations_from_skip = &il_locations[il_locations.partition_point(|&x| x < skip)..];
            let current_search_string_prefix = &search_string[..skip];
            let current_search_string_suffix = &search_string[skip..];
            let t_bounds = Instant::now();
            let search_bound_result = self.search_bounds(&search_string[skip..]);
            self.search_bounds_ns.fetch_add(t_bounds.elapsed().as_nanos() as u64, Ordering::Relaxed);

            // if the shorter part is matched, see if what goes before the matched suffix matches
            // the unmatched part of the prefix
            if let BoundSearchResult::SearchResult((min_bound, max_bound)) = search_bound_result {
                let t_iter = Instant::now();

                // Fast-path: when equate_il=true, !tryptic, and skip=0, every entry in the SA
                // range is a valid match — no per-entry filtering needed.
                if equate_il && !tryptic && skip == 0 {
                    let range_size = max_bound - min_bound;
                    if range_size >= max_matches {
                        // The range is larger than the cutoff: collect the first max_matches
                        // entries directly. The tight collect() loop lets the compiler emit
                        // efficient (potentially SIMD) code for the Vec fill.
                        let result: Vec<i64> = self.sa
                            .iter_range(min_bound, min_bound + max_matches)
                            .collect();
                        self.match_iter_ns.fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);
                        return SearchAllSuffixesResult::MaxMatches(result);
                    }
                    // range_size < max_matches: collect all entries and continue to the next
                    // skip value (rare for short peptides where the range is large).
                    for s in self.sa.iter_range(min_bound, max_bound) {
                        matching_suffixes.push(s);
                    }
                    self.match_iter_ns.fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);
                } else {
                    // Generic path: tryptic filtering, IL checking, or skip > 0.
                    let text = self.proteins.text();
                    let hit_max = self.iterate_sa_range(
                        self.sa.iter_range(min_bound, max_bound),
                        max_bound.saturating_sub(min_bound),
                        text,
                        skip,
                        search_string,
                        current_search_string_prefix,
                        current_search_string_suffix,
                        il_locations_from_skip,
                        equate_il,
                        tryptic,
                        &mut matching_suffixes,
                        max_matches,
                    );
                    self.match_iter_ns.fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    if hit_max {
                        return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                    }
                }
            }

            skip += 1;

            // Only prefetch skip+1's SA range after confirming we're not returning early.
            // Issuing madvise when skip+1 will never execute wastes a syscall
            if skip < self.sa.sample_rate() as usize {
                self.prefetch_kmer_range(&search_string[skip..]);
            }
        }

        if matching_suffixes.is_empty() {
            SearchAllSuffixesResult::NoMatches
        } else {
            SearchAllSuffixesResult::SearchResult(matching_suffixes)
        }
    }

    /// Prefetches both potential next binary search pivots before the blocking sa.get(center)
    /// call. One fetch will be wasted; both are free (single non-blocking CPU instruction).
    /// From iteration 2 onward the needed SA entry is already in L1/L2 cache.
    #[inline]
    fn prefetch_binary_search_pivots(&self, lo: usize, center: usize, hi: usize) {
        self.sa.prefetch_sa_index((lo + center) / 2);
        self.sa.prefetch_sa_index((center + hi) / 2);
    }
}

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
    use text_compression::ProteinText;

    use crate::{
        array::OriginalSA,
        sa_searcher::{test_helpers::get_example_proteins, BoundSearchResult, SearchAllSuffixesResult, Searcher},
        suffix_to_protein_index::{BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinMapping},
        SuffixArray,
    };

    #[test]
    fn test_search_simple() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));

        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        // search bounds 'A'
        let bounds_res = searcher.search_bounds(b"A");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((4, 9)));

        // search bounds '$'
        let bounds_res = searcher.search_bounds(b"$");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((0, 1)));

        // search bounds 'AC'
        let bounds_res = searcher.search_bounds(b"AC");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((6, 8)));
    }

    #[test]
    fn test_search_sparse() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(vec![9, 0, 3, 12, 15, 6, 18], 3));

        let suffix_index_to_protein = SparseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::Sparse(suffix_index_to_protein));

        // search suffix 'VAA'
        let found_suffixes = searcher.search_matching_suffixes(b"VAA", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![7]));

        // search suffix 'AC'
        let found_suffixes = searcher.search_matching_suffixes(b"AC", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![5, 11]));
    }

    #[test]
    fn test_search_dense() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(vec![9, 0, 3, 12, 15, 6, 18], 3));

        let suffix_index_to_protein = DenseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::Dense(suffix_index_to_protein));

        // search suffix 'VAA'
        let found_suffixes = searcher.search_matching_suffixes(b"VAA", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![7]));

        // search suffix 'AC'
        let found_suffixes = searcher.search_matching_suffixes(b"AC", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![5, 11]));
    }

    #[test]
    fn test_il_equality() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));

        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        let bounds_res = searcher.search_bounds(b"I");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((13, 16)));

        // search bounds 'RIZ' with equal I and L
        let bounds_res = searcher.search_bounds(b"RIY");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((17, 18)));
    }

    #[test]
    fn test_il_equality_sparse() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(vec![9, 0, 3, 12, 15, 6, 18], 3));

        let suffix_index_to_protein = SparseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::Sparse(suffix_index_to_protein));

        // search bounds 'RIZ' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(b"RIY", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![16]));

        // search bounds 'RIZ' without equal I and L
        let found_suffixes = searcher.search_matching_suffixes(b"RIY", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::NoMatches);
    }

    // test edge case where an I or L is the first index in the sparse SA.
    #[test]
    fn test_l_first_index_in_sa() {
        let input_string = "LMPYY$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![0, 2, 4], 2));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(b"IM", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0]));
    }

    #[test]
    fn test_il_missing_matches() {
        let input_string = "AAILLL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![6, 0, 1, 5, 4, 3, 2], 1));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        let found_suffixes = searcher.search_matching_suffixes(b"I", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![2, 3, 4, 5]));
    }

    #[test]
    fn test_il_duplication() {
        let input_string = "IIIILL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![6, 5, 4, 3, 2, 1, 0], 1));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        let found_suffixes = searcher.search_matching_suffixes(b"II", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4]));
    }

    #[test]
    fn test_il_suffix_check() {
        let input_string = "IIIILL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![6, 4, 2, 0], 2));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        // search all places where II is in the string IIIILL, but with a sparse SA
        // this way we check if filtering the suffixes works as expected
        let found_suffixes = searcher.search_matching_suffixes(b"II", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2]));
    }

    #[test]
    fn test_il_duplication2() {
        let input_string = "IILLLL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![6, 5, 4, 3, 2, 1, 0], 1));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(b"II", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4]));
    }

    #[test]
    fn test_tryptic_search() {
        let input_string = "PAA-AAKPKAPAA$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);

        let sparse_sa = SuffixArray::Original(OriginalSA(vec![13, 3, 12, 11, 1, 4, 2, 5, 9, 8, 6, 10, 0, 7], 1));
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, SuffixToProteinMapping::BitVec(suffix_index_to_protein));

        let found_suffixes_1 = searcher.search_matching_suffixes(b"PAA", usize::MAX, false, true);
        assert_eq!(found_suffixes_1, SearchAllSuffixesResult::SearchResult(vec![0]));

        let found_suffixes_2 = searcher.search_matching_suffixes(b"APAA", usize::MAX, false, true);
        assert_eq!(found_suffixes_2, SearchAllSuffixesResult::SearchResult(vec![9]));
    }

    // Attaching a k-mer bounds table must not change results, only narrow the search.
    //
    // The table's k-mer prefixes must be L/I-free here: these test fixtures use a raw
    // (non-L→I-normalized) SA, and the k-mer table only groups suffixes contiguously on a
    // normalized SA (which production builds). Chars beyond the prefix may contain L — those
    // are resolved by `compare` in the within-range search, so `KCRLY` is fine (prefix `KCR`).
    #[test]
    fn test_search_with_kmer_table() {
        let example_searcher = || {
            let proteins = get_example_proteins();
            let stp = BitVecSuffixToProtein::new(proteins.text());
            let sa = SuffixArray::Original(OriginalSA(
                vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
            Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp))
        };
        // mix: uses the table (len >= k), bypasses it (len < k), and misses entirely
        let peptides: Vec<&[u8]> = vec![b"VAA", b"CVAA", b"KCR", b"KCRLY", b"AC", b"A", b"ZZZ"];

        let plain = example_searcher();
        let mut kmered = example_searcher();
        kmered.build_kmer_table(3);

        for p in &peptides {
            assert_eq!(
                kmered.search_matching_suffixes(p, usize::MAX, false, false),
                plain.search_matching_suffixes(p, usize::MAX, false, false),
                "k-mer vs plain mismatch for {:?}", std::str::from_utf8(p).unwrap()
            );
        }
    }

    // `max_matches` cutoff returns MaxMatches with exactly that many entries.
    #[test]
    fn test_max_matches_cutoff() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // "A" has 5 matches; cap at 2.
        match searcher.search_matching_suffixes(b"A", 2, true, false) {
            SearchAllSuffixesResult::MaxMatches(v) => assert_eq!(v.len(), 2),
            other => panic!("expected MaxMatches, got {:?}", other),
        }
    }

    // Exercises iterate_sa_range's two-pass prefetch batching (range >= 32) and the batch
    // refill (> BATCH_SIZE=64), using equate_il=false to force the generic (non-fast) path.
    #[test]
    fn test_iterate_sa_range_two_pass() {
        let n = 70usize;
        let mut input = "A".repeat(n);
        input.push('$');
        let text = ProteinText::from_string(&input);
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![],
        }]);
        // SA of A^n$ is [n, n-1, …, 0] ("$" < "A$" < "AA$" < …).
        let sa = SuffixArray::Original(OriginalSA((0..=n as i64).rev().collect(), 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        let found = searcher.search_matching_suffixes(b"A", usize::MAX, false, false);
        let expected: Vec<i64> = (0..n as i64).collect();
        assert_eq!(found, SearchAllSuffixesResult::SearchResult(expected));
    }
}
