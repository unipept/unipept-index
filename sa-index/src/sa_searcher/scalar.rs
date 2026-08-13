//! Scalar (one-peptide-at-a-time) search pipeline: narrow the SA window with a binary
//! search, then collect matching suffixes. The memory-level-parallel counterpart lives in
//! `super::batched`; both call the shared primitives (`compare`, `iterate_sa_range`) that stay
//! in the parent module.

use std::cmp::min;

use sa_mappings::proteins::{ProteinsBackend, SEPARATION_CHARACTER};

use super::{
    BoundSearch,
    BoundSearch::{Maximum, Minimum},
    BoundSearchResult, MAX_RESULT_PREALLOC, SearchAllSuffixesResult, Searcher,
    metrics::Timer,
    tryptic::tryptic_extension_chars
};
use crate::{array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

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
        lcp_skip: usize
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
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(lo), min(lcp_left, lcp_right), bound);

            found |= lcp_center == search_string.len();

            if bound == Minimum && retval {
                hi = lo;
            }
        }

        match bound {
            Minimum => (found, hi),
            Maximum => (found, lo)
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
        let (left, right, lcp_skip) = match &self.kmer_table {
            Some(table) if search_string.len() >= table.k => match table.lookup(&search_string[..table.k]) {
                Some((lo, hi)) => (lo, hi + 1, table.k),
                None => return BoundSearchResult::NoMatches
            },
            _ => (0, self.sa.len(), 0)
        };
        self.search_bounds_within(search_string, left, right, lcp_skip)
    }

    /// `search_bounds` over the whole suffix array, never consulting the k-mer table.
    ///
    /// Required for search strings containing the protein separator, which the k-mer table
    /// cannot represent: `ALPHABET` in `kmer_table.rs` holds amino acids only, so
    /// `bytes_to_kmer_index` returns `None` for any k-mer containing `-`, and `search_bounds`
    /// would turn that into `NoMatches`. Routing the separator variant of the left-extended
    /// tryptic search through the table would therefore silently drop every protein-start
    /// match — roughly 3 % of all tryptic hits, and every protein's N-terminal peptide.
    ///
    /// Costs a full-height binary search (~35 random SA reads instead of ~13). Measured on the
    /// full DB that is affordable: even 26-50aa tryptic queries, where this fixed cost is the
    /// largest share of a ~4 µs budget, came out 1.40x *faster* overall (run5). Adding `-` to the
    /// k-mer table's ALPHABET would remove the special case, at +29 MB and a rebuild of every
    /// table file — not worth it while the bypass is free. Bounded by
    /// `test_extended_protein_start_with_kmer_table`.
    fn search_bounds_full_range(&self, search_string: &[u8]) -> BoundSearchResult {
        self.search_bounds_within(search_string, 0, self.sa.len(), 0)
    }

    /// Shared tail of `search_bounds` / `search_bounds_full_range`: the two bound searches over
    /// an already-chosen window.
    fn search_bounds_within(
        &self,
        search_string: &[u8],
        left: usize,
        right: usize,
        lcp_skip: usize
    ) -> BoundSearchResult {
        let (found_min, min_bound) = self.binary_search_bound_in_range(Minimum, search_string, left, right, lcp_skip);

        if !found_min {
            return BoundSearchResult::NoMatches;
        }

        let (_, max_bound) = self.binary_search_bound_in_range(Maximum, search_string, left, right, lcp_skip);

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
        // There is no k-mer prefetch pass here, and there is no longer one in the batched path
        // either. The call used to sit immediately before the `search_bounds` below that repeats
        // the identical k-mer table lookup — no intervening work, so no latency to hide — and
        // once the mmap `madvise` behind it was removed, the hint it issued became a no-op and
        // only the discarded probe remained. Removing it from this path measured +0.3% median
        // over 6 (bucket, backend) combos (run3, inside the 3.9% noise floor); removing it from
        // the batched path measured neutral on both backends. See `super::batched`.

        let mut matching_suffixes: Vec<i64> = Vec::with_capacity(max_matches.min(MAX_RESULT_PREALLOC));
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        let sample = self.sa.sample_rate() as usize;

        // Tryptic searches replace the most expensive skip pass with left-extended searches.
        //
        // A `skip = j` pass finds matches at `ms ≡ -j (mod s)` by searching a j-character-shorter
        // string, so its SA range grows ~20x per character dropped — `skip = s-1` alone is ~95 %
        // of all candidates examined, a ~40x overscan (measured: 195,293 candidates per 5-10aa
        // peptide to find 134 matches).
        //
        // But a tryptic match at `ms` requires `text[ms-1]` to be K, R or a separator — which is
        // exactly the character needed to extend the search string *leftward* instead of
        // truncating it. `X + peptide` is one character longer, so its range is ~20x *smaller*.
        // Its SA entry sits at `ms-1`, covering `ms ≡ 1 (mod s)` — the same positions
        // `skip = s-1` covers, since `-(s-1) ≡ 1 (mod s)`.
        //
        // Only for s >= 2: at s = 1 every position is sampled, `skip = 0` already covers
        // everything, and there is no truncated pass to replace (extending would drop ms = 0).
        //
        // Measured on the full DB (run5 vs run4, 20 reps, results bit-identical — the
        // `candidates_accepted` counters match exactly):
        //
        //   bucket   candidates examined      throughput        tryptic vs non-tryptic
        //   5-10aa   1.95e9 -> 1.87e8 (10.4x)  741 -> 7,632 qps   12.1x slower -> 1.2x
        //   11-25aa  3.96e7 -> 1.96e7 ( 2.0x)  29k  ->  60k  qps   6.6x slower -> 3.3x
        //   26-50aa  4.44e6 -> 2.29e6 ( 1.9x) 247k  -> 346k  qps   1.3x slower -> 0.9x
        //
        // The gain concentrates on short peptides because the ~20x-per-character range growth
        // only holds while the range is dominated by random-match statistics. A 26-50aa peptide
        // is already near its floor of roughly one occurrence per protein, so dropping a
        // character barely widens it — which is also why those buckets were never the problem.
        // Long peptides end up *faster* than their non-tryptic counterparts (0.9x) simply
        // because a non-tryptic query returns up to `max_matches` hits to retrieve, and a
        // tryptic one returns a handful.
        let use_extended = tryptic && sample >= 2;
        let skip_end = if use_extended { sample - 1 } else { sample };

        let mut skip: usize = 0;
        while skip < skip_end {
            // il_locations is built in ascending index order, so partition_point gives us
            // the first position that is relevant for this skip value in O(log n).
            // These are still absolute positions within `search_string`, not within the suffix.
            let il_locations_from_skip = &il_locations[il_locations.partition_point(|&x| x < skip)..];
            let current_search_string_prefix = &search_string[..skip];
            let current_search_string_suffix = &search_string[skip..];
            let t_bounds = Timer::start();
            let search_bound_result = self.search_bounds(&search_string[skip..]);
            self.metrics.search_bounds_ns.add(t_bounds.elapsed_ns());

            // if the shorter part is matched, see if what goes before the matched suffix matches
            // the unmatched part of the prefix
            if let BoundSearchResult::SearchResult((min_bound, max_bound)) = search_bound_result {
                let t_iter = Timer::start();

                // Fast-path: when equate_il=true, !tryptic, and skip=0, every entry in the SA
                // range is a valid match — no per-entry filtering needed.
                //
                // The same holds when equate_il=false but the peptide contains no I or L:
                // `compare` (used to narrow the SA range above) always normalizes L->I on both
                // sides, because the index itself is built with every L replaced by I. So an
                // in-range suffix can only fail to be an exact match if the text has an I where
                // the peptide has an L (or vice versa) — impossible when the peptide has
                // neither character. In that case equate_il is moot and the fast path is sound.
                if (equate_il || il_locations.is_empty()) && !tryptic && skip == 0 {
                    let range_size = max_bound - min_bound;
                    if range_size >= max_matches {
                        // The range is larger than the cutoff: collect the first max_matches
                        // entries directly. The tight collect() loop lets the compiler emit
                        // efficient (potentially SIMD) code for the Vec fill.
                        let result: Vec<i64> = self.sa.iter_range(min_bound, min_bound + max_matches).collect();
                        self.metrics.match_iter_ns.add(t_iter.elapsed_ns());
                        return SearchAllSuffixesResult::MaxMatches(result);
                    }
                    // range_size < max_matches: collect all entries and continue to the next
                    // skip value (rare for short peptides where the range is large).
                    for s in self.sa.iter_range(min_bound, max_bound) {
                        matching_suffixes.push(s);
                    }
                    self.metrics.match_iter_ns.add(t_iter.elapsed_ns());
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
                        max_matches
                    );
                    self.metrics.match_iter_ns.add(t_iter.elapsed_ns());
                    if hit_max {
                        return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                    }
                }
            }

            skip += 1;
        }

        // Left-extended phase, replacing the skip = sample-1 pass (see above).
        if use_extended {
            let text = self.proteins.text();
            let last_is_kr = matches!(search_string.last(), Some(b'K' | b'R'));
            // One reusable buffer for `X + search_string`: three allocations per peptide in this
            // path would cost more than the extra searches save.
            let mut extended = Vec::with_capacity(search_string.len() + 1);

            for &prefix_char in tryptic_extension_chars(search_string) {
                extended.clear();
                extended.push(prefix_char);
                extended.extend_from_slice(search_string);

                let t_bounds = Timer::start();
                // The separator cannot be represented in the k-mer table — see
                // `search_bounds_full_range` for why routing it there would silently drop every
                // protein-start match.
                let bounds = if prefix_char == SEPARATION_CHARACTER {
                    self.search_bounds_full_range(&extended)
                } else {
                    self.search_bounds(&extended)
                };
                self.metrics.search_bounds_ns.add(t_bounds.elapsed_ns());

                if let BoundSearchResult::SearchResult((min_bound, max_bound)) = bounds {
                    let t_iter = Timer::start();
                    let hit_max = self.iterate_extended_sa_range(
                        self.sa.iter_range(min_bound, max_bound),
                        max_bound.saturating_sub(min_bound),
                        text,
                        search_string,
                        &il_locations,
                        equate_il,
                        last_is_kr,
                        &mut matching_suffixes,
                        max_matches
                    );
                    self.metrics.match_iter_ns.add(t_iter.elapsed_ns());
                    if hit_max {
                        return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                    }
                }
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

#[cfg(test)]
mod tests {
    use crate::sa_searcher::{
        BoundSearchResult, SearchAllSuffixesResult,
        test_utils::{
            EXAMPLE_SA_FULL, EXAMPLE_SA_SPARSE3, Mapping, TRYPTIC_FIXTURE, example_searcher, example_searcher_with,
            repeated_residue_searcher, searcher_over_text, tryptic_fixture_peptides
        }
    };

    #[test]
    fn test_search_simple() {
        let searcher = example_searcher_with(&EXAMPLE_SA_FULL, 1, Mapping::BitVec);

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

    // A sparse SA needs the `skip` loop to recover matches at unsampled positions, and the
    // answer must not depend on which suffix-to-protein representation the index carries — the
    // three are interchangeable by construction, and search must not have grown a dependence on
    // one of them.
    #[test]
    fn test_search_sparse() {
        for mapping in [Mapping::Dense, Mapping::Sparse, Mapping::BitVec] {
            let searcher = example_searcher_with(&EXAMPLE_SA_SPARSE3, 3, mapping);

            // 'VAA' sits at 7, which sparseness 3 does not sample: only skip = 1 finds it.
            let found_suffixes = searcher.search_matching_suffixes(b"VAA", usize::MAX, false, false);
            assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![7]));

            // 'AC' matches twice, at a sampled (position 12 → skip 0) and an unsampled position.
            let found_suffixes = searcher.search_matching_suffixes(b"AC", usize::MAX, false, false);
            assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![5, 11]));
        }
    }

    #[test]
    fn test_il_equality() {
        let searcher = example_searcher_with(&EXAMPLE_SA_FULL, 1, Mapping::BitVec);

        let bounds_res = searcher.search_bounds(b"I");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((13, 16)));

        // search bounds 'RIZ' with equal I and L
        let bounds_res = searcher.search_bounds(b"RIY");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((17, 18)));
    }

    #[test]
    fn test_il_equality_sparse() {
        let searcher = example_searcher_with(&EXAMPLE_SA_SPARSE3, 3, Mapping::Sparse);

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
        let searcher = searcher_over_text("LMPYY$", 2);

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(b"IM", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0]));
    }

    #[test]
    fn test_il_missing_matches() {
        let searcher = searcher_over_text("AAILLL$", 1);

        let found_suffixes = searcher.search_matching_suffixes(b"I", usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![2, 3, 4, 5]));
    }

    // With `equate_il`, every I/L arrangement of the same length is one text: the index stores
    // both as I, so "IIIILL" and "IILLLL" must answer identically and each match must be
    // reported exactly once, not once per spelling.
    #[test]
    fn test_il_duplication() {
        for text in ["IIIILL$", "IILLLL$"] {
            let searcher = searcher_over_text(text, 1);

            let found_suffixes = searcher.search_matching_suffixes(b"II", usize::MAX, true, false);
            assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4]), "text {text}");
        }
    }

    #[test]
    fn test_il_suffix_check() {
        let searcher = searcher_over_text("IIIILL$", 2);

        // search all places where II is in the string IIIILL, but with a sparse SA
        // this way we check if filtering the suffixes works as expected
        let found_suffixes = searcher.search_matching_suffixes(b"II", usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2]));
    }

    #[test]
    fn test_tryptic_search() {
        let searcher = searcher_over_text("PAA-AAKPKAPAA$", 1);

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
        // mix: uses the table (len >= k), bypasses it (len < k), and misses entirely
        let peptides: Vec<&[u8]> = vec![b"VAA", b"CVAA", b"KCR", b"KCRLY", b"AC", b"A", b"ZZZ"];

        let plain = example_searcher();
        let mut kmered = example_searcher();
        kmered.build_kmer_table(3);

        for p in &peptides {
            assert_eq!(
                kmered.search_matching_suffixes(p, usize::MAX, false, false),
                plain.search_matching_suffixes(p, usize::MAX, false, false),
                "k-mer vs plain mismatch for {:?}",
                std::str::from_utf8(p).unwrap()
            );
        }
    }

    // `max_matches` cutoff returns MaxMatches with exactly that many entries.
    #[test]
    fn test_max_matches_cutoff() {
        let searcher = example_searcher();

        // "A" has 5 matches; cap at 2.
        match searcher.search_matching_suffixes(b"A", 2, true, false) {
            SearchAllSuffixesResult::MaxMatches(v) => assert_eq!(v.len(), 2),
            other => panic!("expected MaxMatches, got {:?}", other)
        }
    }

    // ── left-extended tryptic search ─────────────────────────────────────

    // The left-extended tryptic path (sparseness >= 2) must return exactly what the dense index
    // returns, where `use_extended` is false and the original skip loop runs. A dense searcher
    // over the same text is therefore an exact oracle.
    //
    // This is the load-bearing test for the whole transform: no pre-existing test reached it,
    // because `test_tryptic_search` uses sparseness 1.
    #[test]
    fn test_extended_tryptic_matches_dense() {
        let dense = searcher_over_text(TRYPTIC_FIXTURE, 1);
        let peptides = tryptic_fixture_peptides();
        assert!(peptides.len() > 30, "corpus too small to be meaningful");

        // Guard against a vacuous pass: the fixture must actually produce tryptic hits.
        let hits = peptides
            .iter()
            .filter(|p| {
                !matches!(dense.search_matching_suffixes(p, usize::MAX, true, true), SearchAllSuffixesResult::NoMatches)
            })
            .count();
        assert!(hits >= 5, "fixture yields only {hits} tryptic hits; comparison would be weak");

        for &sparseness in &[2u8, 3] {
            let sparse = searcher_over_text(TRYPTIC_FIXTURE, sparseness);
            for equate_il in [false, true] {
                for p in &peptides {
                    assert_eq!(
                        sparse.search_matching_suffixes(p, usize::MAX, equate_il, true),
                        dense.search_matching_suffixes(p, usize::MAX, equate_il, true),
                        "sparseness={sparseness} equate_il={equate_il} peptide={:?}",
                        std::str::from_utf8(p).unwrap()
                    );
                }
            }
        }
    }

    // A match at a protein start is found through the `'-' + peptide` variant, which must bypass
    // the k-mer table: `'-'` is not in the table's ALPHABET, so a table lookup returns None and
    // `search_bounds` would report NoMatches — silently losing every protein-start match.
    //
    // "PKTR" starts protein 2 at position 21, which is *odd*, so at sparseness 2 it is not in the
    // suffix array and can only be reached via the extended search (whose SA entry is the
    // separator at 20). It also starts with proline, so the K/R variants are correctly skipped
    // and `'-'` is the only search performed — making this fail if the bypass is wrong.
    #[test]
    fn test_extended_protein_start_with_kmer_table() {
        let plain = searcher_over_text(TRYPTIC_FIXTURE, 2);
        let mut kmered = searcher_over_text(TRYPTIC_FIXTURE, 2);
        kmered.build_kmer_table(3);

        // Non-vacuous: the protein-start match really is there, and really is at position 21.
        assert_eq!(
            plain.search_matching_suffixes(b"PKTR", usize::MAX, true, true),
            SearchAllSuffixesResult::SearchResult(vec![21]),
            "protein-start tryptic match missing from the un-tabled searcher"
        );

        for p in [&b"PKTR"[..], b"RIY", b"KTR", b"MKA", b"AKT", b"QST"] {
            assert_eq!(
                kmered.search_matching_suffixes(p, usize::MAX, true, true),
                plain.search_matching_suffixes(p, usize::MAX, true, true),
                "k-mer table changed the tryptic result for {:?}",
                std::str::from_utf8(p).unwrap()
            );
        }
    }

    // A peptide starting with proline can only match at a protein start: K|R followed by P is not
    // a trypsin cut site, so the K and R extension variants are skipped entirely.
    #[test]
    fn test_extended_proline_start_only_matches_protein_start() {
        let s = searcher_over_text(TRYPTIC_FIXTURE, 2);

        // "PKTR" at 21 is a protein start → found.
        assert_eq!(
            s.search_matching_suffixes(b"PKTR", usize::MAX, true, true),
            SearchAllSuffixesResult::SearchResult(vec![21])
        );
        // "PQST" at 16 is preceded by K (15), but proline blocks the cut → no tryptic match,
        // even though the same peptide matches fine without the tryptic filter.
        assert_eq!(s.search_matching_suffixes(b"PQST", usize::MAX, true, true), SearchAllSuffixesResult::NoMatches);
        assert_eq!(
            s.search_matching_suffixes(b"PQST", usize::MAX, true, false),
            SearchAllSuffixesResult::SearchResult(vec![16])
        );
    }

    // sparseness 1 must be untouched by the transform: every position is sampled, skip=0 already
    // covers everything, and extending would drop ms = 0. "MKAPTR" at position 0 is the canary —
    // a protein start with no preceding character at all, so it is reachable *only* through the
    // skip=0 pass and never through a left-extended search.
    //
    // (It is tryptic because it ends in R and text[6] = 'V' is not proline. "MKA" would not be:
    // it ends in A at position 3, which holds 'P'.)
    #[test]
    fn test_extended_guard_leaves_dense_index_alone() {
        let dense = searcher_over_text(TRYPTIC_FIXTURE, 1);
        assert_eq!(
            dense.search_matching_suffixes(b"MKAPTR", usize::MAX, true, true),
            SearchAllSuffixesResult::SearchResult(vec![0]),
            "match at ms=0 lost on the dense index"
        );
        // and it survives on a sparse index too, via the skip=0 pass (0 is always sampled)
        for &sparseness in &[2u8, 3] {
            assert_eq!(
                searcher_over_text(TRYPTIC_FIXTURE, sparseness).search_matching_suffixes(
                    b"MKAPTR",
                    usize::MAX,
                    true,
                    true
                ),
                SearchAllSuffixesResult::SearchResult(vec![0]),
                "match at ms=0 lost at sparseness {sparseness}"
            );
        }
    }

    // I/L-free peptides must produce identical results for equate_il=true and equate_il=false:
    // `compare` (used to narrow the SA range) always normalizes L->I on both sides, because the
    // index is built with every L replaced by I. So an in-range suffix can only fail to be an
    // exact match if the text has an I where the peptide has an L (or vice versa) — impossible
    // when the peptide has neither character. equate_il is therefore moot for such peptides, and
    // the fast path should trigger identically for both. Exercises both fast-path branches:
    // range_size (70) >= max_matches (10), and range_size < max_matches (usize::MAX).
    #[test]
    fn test_il_free_fast_path_matches_equate_il_true() {
        let n = 70usize;
        let searcher = repeated_residue_searcher('A', n);

        // range_size (70) < max_matches: both collect the full range via the fast path.
        assert_eq!(
            searcher.search_matching_suffixes(b"A", usize::MAX, false, false),
            searcher.search_matching_suffixes(b"A", usize::MAX, true, false)
        );
        // range_size (70) >= max_matches (10): both hit the MaxMatches cutoff identically.
        assert_eq!(
            searcher.search_matching_suffixes(b"A", 10, false, false),
            searcher.search_matching_suffixes(b"A", 10, true, false)
        );
    }

    // Negative counterpart to the above: a peptide that DOES contain I/L must keep taking the
    // slow (validating) path when equate_il=false, even over a large SA range, and must still
    // discriminate I from L correctly. Text is all 'L'; `compare` treats L==I during the bound
    // search (narrowing the range), so searching "I" matches every position there — but none of
    // them is an actual 'I' in the text, so equate_il=false must reject all of them while
    // equate_il=true (which does take the fast path) accepts them.
    #[test]
    fn test_il_present_large_range_still_slow_path() {
        let n = 70usize;
        let searcher = repeated_residue_searcher('L', n);

        // equate_il=true takes the fast path and happily returns 10 raw (unvalidated) matches...
        match searcher.search_matching_suffixes(b"I", 10, true, false) {
            SearchAllSuffixesResult::MaxMatches(v) => assert_eq!(v.len(), 10),
            other => panic!("expected MaxMatches, got {:?}", other)
        }
        // ...but equate_il=false must still validate and reject them all: no position actually
        // holds 'I'. Covers both the range_size >= max_matches and < max_matches branches.
        assert_eq!(searcher.search_matching_suffixes(b"I", 10, false, false), SearchAllSuffixesResult::NoMatches);
        assert_eq!(
            searcher.search_matching_suffixes(b"I", usize::MAX, false, false),
            SearchAllSuffixesResult::NoMatches
        );
    }
}
