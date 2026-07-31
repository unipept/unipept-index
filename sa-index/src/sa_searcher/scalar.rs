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
