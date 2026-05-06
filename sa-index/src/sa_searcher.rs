use std::{cmp::min, ops::Deref, sync::atomic::{AtomicU64, Ordering}, time::Instant};

use sa_mappings::proteins::{ProteinRef, Proteins, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinTextSlice;

use crate::{
    KmerTable, Nullable, SuffixArray,
    sa_searcher::BoundSearch::{Maximum, Minimum},
    suffix_to_protein_index::{DenseSuffixToProtein, SparseSuffixToProtein, BitVecSuffixToProtein, SuffixToProteinIndex}
};

/// Enum indicating if we are searching for the minimum, or maximum bound in the suffix array
#[derive(Clone, Copy, PartialEq)]
enum BoundSearch {
    Minimum,
    Maximum
}

/// Enum representing the minimum and maximum bound of the found matches in the suffix array
#[derive(PartialEq, Debug)]
pub enum BoundSearchResult {
    NoMatches,
    SearchResult((usize, usize))
}

/// Enum representing the matching suffixes after searching a peptide in the suffix array
/// Both the MaxMatches and SearchResult indicate found suffixes, but MaxMatches is used when the
/// cutoff is reached.
#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    NoMatches,
    MaxMatches(Vec<i64>),
    SearchResult(Vec<i64>)
}

/// Custom implementation of partialEq for SearchAllSuffixesResult
/// We consider 2 SearchAllSuffixesResult equal if they exist of the same key, and the Vec contains
/// the same values, but the order can be different
impl PartialEq for SearchAllSuffixesResult {
    fn eq(&self, other: &Self) -> bool {
        fn array_eq_unordered(arr1: &[i64], arr2: &[i64]) -> bool {
            let mut arr1_copy = arr1.to_owned();
            let mut arr2_copy = arr2.to_owned();

            arr1_copy.sort();
            arr2_copy.sort();

            arr1_copy == arr2_copy
        }

        match (self, other) {
            (SearchAllSuffixesResult::MaxMatches(arr1), SearchAllSuffixesResult::MaxMatches(arr2)) => {
                array_eq_unordered(arr1, arr2)
            }
            (SearchAllSuffixesResult::SearchResult(arr1), SearchAllSuffixesResult::SearchResult(arr2)) => {
                array_eq_unordered(arr1, arr2)
            }
            (SearchAllSuffixesResult::NoMatches, SearchAllSuffixesResult::NoMatches) => true,
            _ => false
        }
    }
}

pub struct SparseSearcher(Searcher);

impl SparseSearcher {
    pub fn new(sa: SuffixArray, proteins: Proteins) -> Self {
        let suffix_index_to_protein = SparseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));
        Self(searcher)
    }
}

impl Deref for SparseSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub struct BitVecSearcher(Searcher);

impl BitVecSearcher {
    pub fn new(sa: SuffixArray, proteins: Proteins) -> Self {
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));
        Self(searcher)
    }
}

impl Deref for BitVecSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub struct DenseSearcher(Searcher);

impl DenseSearcher {
    pub fn new(sa: SuffixArray, proteins: Proteins) -> Self {
        let suffix_index_to_protein = DenseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));
        Self(searcher)
    }
}

impl Deref for DenseSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Struct that contains all the elements needed to search a peptide in the suffix array
/// This struct also contains all the functions used for search
///
/// # Arguments
/// * `sa` - The sparse suffix array representing the protein database
/// * `sparseness_factor` - The sparseness factor used by the suffix array
/// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix is
///   part of
/// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the
///   taxonomic analysis provided by Unipept
/// * `function_aggregator` - Object used to retrieve the functional annotations and to calculate
///   the functional analysis provided by Unipept
pub struct Searcher {
    pub sa: SuffixArray,
    pub proteins: Proteins,
    pub suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
    pub kmer_table: Option<KmerTable>,
    /// Total nanoseconds spent inside `search_bounds()` across all queries (since last drain).
    pub search_bounds_ns: AtomicU64,
    /// Total nanoseconds spent iterating matches in `search_matching_suffixes()` (since last drain).
    pub match_iter_ns: AtomicU64,
}

impl Searcher {
    /// Creates a new Searcher object
    ///
    /// # Arguments
    /// * `sa` - The sparse suffix array representing the protein database
    /// * `sparseness_factor` - The sparseness factor used by the suffix array
    /// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix
    ///   is part of
    /// * `proteins` - List of all the proteins where the suffix array is build on
    /// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the
    ///   taxonomic analysis provided by Unipept
    /// * `function_aggregator` - Object used to retrieve the functional annotations and to
    ///   calculate the functional analysis provided by Unipept
    ///
    /// # Returns
    ///
    /// Returns a new Searcher object
    pub fn new(sa: SuffixArray, proteins: Proteins, suffix_index_to_protein: Box<dyn SuffixToProteinIndex>) -> Self {
        Self {
            sa,
            proteins,
            suffix_index_to_protein,
            kmer_table: None,
            search_bounds_ns: AtomicU64::new(0),
            match_iter_ns: AtomicU64::new(0),
        }
    }

    /// Returns `(search_bounds_ns, match_iter_ns)` accumulated since the last call and resets both
    /// counters to zero.  Safe to call concurrently with ongoing searches (relaxed ordering).
    pub fn drain_timing_ns(&self) -> (u64, u64) {
        let bounds = self.search_bounds_ns.swap(0, Ordering::Relaxed);
        let iter = self.match_iter_ns.swap(0, Ordering::Relaxed);
        (bounds, iter)
    }

    /// Attaches a pre-built k-mer bounds table to this searcher.
    pub fn with_kmer_table(mut self, table: KmerTable) -> Self {
        self.kmer_table = Some(table);
        self
    }

    /// Builds and attaches a k-mer table with the given `k` using the already-loaded index data.
    pub fn build_kmer_table(&mut self, k: usize) {
        self.kmer_table = Some(KmerTable::build_from_sa(&self.sa, self.proteins.text(), k));
    }

    /// Normalizes L to I so that both map to the same character during suffix array comparisons.
    /// The index is built with L replaced by I, so all order comparisons must do the same.
    #[inline]
    fn normalize_li(c: u8) -> u8 {
        if c == b'L' { b'I' } else { c }
    }

    /// Compares `search_string` against the suffix starting at `suffix` in the protein text,
    /// skipping the first `skip` characters (known to match already).
    ///
    /// Returns `(bound_satisfied, matched_len)` where:
    /// - `bound_satisfied` is `true` when the bound condition holds: for `Minimum`, the search
    ///   string is ≤ the suffix; for `Maximum`, it is ≥ the suffix.
    /// - `matched_len` is how many characters of `search_string` matched the suffix.
    ///
    /// L and I are treated as equal during matching because the index was built with L → I.
    fn compare(&self, search_string: &[u8], suffix: i64, skip: usize, bound: BoundSearch) -> (bool, usize) {
        let text = self.proteins.text();
        let mut i_text = (suffix as usize) + skip;
        let mut i_search = skip;
        let mut bound_satisfied = false;

        let condition_check = match bound {
            Minimum => |a: u8, b: u8| a < b,
            Maximum => |a: u8, b: u8| a > b,
        };

        // Advance while characters match (treating L == I).
        while i_search < search_string.len()
            && i_text < text.len()
            && Self::normalize_li(search_string[i_search]) == Self::normalize_li(text.get(i_text))
        {
            i_text += 1;
            i_search += 1;
        }

        if !search_string.is_empty() {
            if i_search == search_string.len() {
                bound_satisfied = true;
            } else if i_text < text.len() {
                // The index has L replaced by I, so normalize both sides before ordering.
                bound_satisfied = condition_check(
                    Self::normalize_li(search_string[i_search]),
                    Self::normalize_li(text.get(i_text)),
                );
            }
        }

        (bound_satisfied, i_search)
    }

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

    /// Check if a cut is the start of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the start of a protein.
    #[inline]
    fn check_start_of_protein(&self, cut_index: usize) -> bool {
        cut_index == 0 || self.proteins.text().get(cut_index - 1) == SEPARATION_CHARACTER
    }

    /// Check if a cut is the end of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the end of a protein.
    #[inline]
    fn check_end_of_protein(&self, cut_index: usize) -> bool {
        self.proteins.text().get(cut_index) == TERMINATION_CHARACTER
            || self.proteins.text().get(cut_index) == SEPARATION_CHARACTER
    }

    /// Check if a cut is a tryptic cut, so check if the amino acid preceding the cut is K or R and the amino acid at the cut is not P.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is a tryptic cut.
    #[inline]
    fn check_tryptic_cut(&self, cut_index: usize) -> bool {
        (self.proteins.text().get(cut_index - 1) == b'K' || self.proteins.text().get(cut_index - 1) == b'R')
            && self.proteins.text().get(cut_index) != b'P'
    }

    /// Returns true of the prefixes are the same
    /// if `equate_il` is set to true, L and I are considered the same
    ///
    /// # Arguments
    /// * `search_string_prefix` - The unchecked prefix of the string/peptide that is searched
    /// * `index_prefix` - The unchecked prefix from the protein from the suffix array
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string_prefix` and `index_prefix` are considered the same, otherwise
    /// false
    #[inline]
    fn check_prefix(search_string_prefix: &[u8], index_prefix: ProteinTextSlice, equate_il: bool) -> bool {
        index_prefix.equals_slice(search_string_prefix, equate_il)
    }

    /// Returns true of the search_string and index_string are equal
    /// This is automatically true if `equate_il` is set to true, since there matched during
    /// search where I = L If `equate_il` is set to false, we need to check if the I and
    /// L locations have the same character
    ///
    /// # Arguments
    /// * `skip` - The used skip factor during the search iteration
    /// * `il_locations` - The locations of the I's and L's in the **original** peptide
    /// * `search_string` - The peptide that is being searched, but already with the skipped prefix
    ///   removed from it
    /// * `index_string` - The suffix that search_string matches with when I and L were equalized
    ///   during search
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string` and `index_string` are considered the same, otherwise false
    fn check_suffix(
        skip: usize,
        il_locations: &[usize],
        search_string: &[u8],
        text_slice: ProteinTextSlice,
        equate_il: bool
    ) -> bool {
        if equate_il { true } else { text_slice.check_il_locations(skip, il_locations, search_string) }
    }

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

    /// Issues hardware prefetch hints for all text positions that will be accessed
    /// during validation of a candidate match at [ms, me). Called in Pass 1 of the
    /// two-pass batching loop to give DRAM latency hiding time before Pass 2 reads.
    #[inline]
    fn prefetch_match_positions(text: &text_compression::ProteinText, ms: usize, me: usize) {
        text.prefetch_at(ms.saturating_sub(1));
        text.prefetch_at(ms);
        text.prefetch_at(me.saturating_sub(1));
        text.prefetch_at(me);
    }

    /// Issues an early OS prefetch hint for the k-mer's SA range (skip=0 case), giving the OS
    /// more lead time to load those pages into the page cache before binary search starts.
    #[inline]
    #[cfg_attr(not(feature = "mmap"), allow(unused_variables))]
    fn prefetch_kmer_range(&self, search_string: &[u8]) {
        #[cfg(feature = "mmap")]
        if let Some(table) = &self.kmer_table {
            if search_string.len() >= table.k {
                if let Some((lo, hi)) = table.lookup(&search_string[..table.k]) {
                    self.sa.prefetch_sa_range(lo, hi + 1);
                }
            }
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

    /// Checks whether the candidate suffix `raw` is a valid match for the current search
    /// parameters. Returns `Some(match_start)` when valid, `None` otherwise.
    #[inline]
    fn validate_candidate(
        &self,
        text: &text_compression::ProteinText,
        raw: i64,
        skip: usize,
        search_string: &[u8],
        prefix: &[u8],
        suffix_str: &[u8],
        il_locations: &[usize],
        equate_il: bool,
        tryptic: bool,
    ) -> Option<i64> {
        let suffix = raw as usize;
        if suffix < skip { return None; }
        let match_start = suffix - skip;
        let match_end = suffix + search_string.len() - skip;
        let valid = (skip == 0
            || Self::check_prefix(prefix, ProteinTextSlice::new(text, match_start, suffix), equate_il))
            && Self::check_suffix(skip, il_locations, suffix_str, ProteinTextSlice::new(text, suffix, match_end), equate_il)
            && (!tryptic
                || ((self.check_start_of_protein(match_start) || self.check_tryptic_cut(match_start))
                    && (self.check_end_of_protein(match_end) || self.check_tryptic_cut(match_end))));
        if valid { Some(match_start as i64) } else { None }
    }

    /// Iterates the SA range and validates each candidate suffix.
    /// Returns `true` if `max_matches` was reached, `false` otherwise.
    ///
    /// Two-pass batching with hardware prefetch hints to hide DRAM latency:
    /// Pass 1 — fills a batch and issues prefetch hints for the text positions to validate.
    /// Pass 2 — validates candidates after prefetches have had time to complete.
    fn iterate_sa_range(
        &self,
        mut sa_iter: impl Iterator<Item = i64>,
        range_size: usize,
        text: &text_compression::ProteinText,
        skip: usize,
        search_string: &[u8],
        prefix: &[u8],
        suffix_str: &[u8],
        il_locations: &[usize],
        equate_il: bool,
        tryptic: bool,
        matching_suffixes: &mut Vec<i64>,
        max_matches: usize,
    ) -> bool {
        // Tuned on x86_64 Zen4/Intel Sapphire Rapids: DRAM latency ~80–100 ns, one SA entry
        // read per ~2–3 ns at that cache level → 64 entries ≈ 192 ns gap, comfortably above
        // the latency floor. Re-benchmark on ARM or NVMe-backed mmap.
        const BATCH_SIZE: usize = 64;
        // Minimum range for prefetch to resolve before use.
        // Below this threshold the two-pass overhead exceeds the latency-hiding benefit.
        const PREFETCH_THRESHOLD: usize = 32;

        if range_size < PREFETCH_THRESHOLD {
            for raw in sa_iter {
                if let Some(v) = self.validate_candidate(
                    text, raw, skip, search_string, prefix, suffix_str, il_locations, equate_il, tryptic,
                ) {
                    matching_suffixes.push(v);
                    if matching_suffixes.len() >= max_matches { return true; }
                }
            }
            return false;
        }

        let mut raw_batch = [0i64; BATCH_SIZE];

        loop {
            // --- Pass 1: fill batch and prefetch text positions ---
            let mut batch_len = 0usize;
            for s in &mut sa_iter {
                let su = s as usize;
                if su >= skip {
                    Self::prefetch_match_positions(text, su - skip, su + search_string.len() - skip);
                }
                raw_batch[batch_len] = s;
                batch_len += 1;
                if batch_len == BATCH_SIZE { break; }
            }
            if batch_len == 0 { break; }

            // --- Pass 2: validate (prefetches have had time to complete) ---
            for &raw in &raw_batch[..batch_len] {
                if let Some(v) = self.validate_candidate(
                    text, raw, skip, search_string, prefix, suffix_str, il_locations, equate_il, tryptic,
                ) {
                    matching_suffixes.push(v);
                    if matching_suffixes.len() >= max_matches { return true; }
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use sa_mappings::proteins::{Protein, Proteins};
    use text_compression::ProteinText;

    use crate::{
        sa_searcher::{BoundSearchResult, SearchAllSuffixesResult, Searcher}, suffix_to_protein_index::{BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein}, SuffixArray
    };

    #[test]
    fn test_partial_eq_search_all_suffixes_result() {
        let search_all_suffixes_result_1 = SearchAllSuffixesResult::SearchResult(vec![1, 2, 3]);
        let search_all_suffixes_result_2 = SearchAllSuffixesResult::SearchResult(vec![3, 2, 1]);
        let search_all_suffixes_result_3 = SearchAllSuffixesResult::SearchResult(vec![1, 2, 4]);
        let search_all_suffixes_result_4 = SearchAllSuffixesResult::MaxMatches(vec![1, 2, 3]);
        let search_all_suffixes_result_5 = SearchAllSuffixesResult::MaxMatches(vec![3, 2, 1]);
        let search_all_suffixes_result_6 = SearchAllSuffixesResult::MaxMatches(vec![1, 2, 4]);
        let search_all_suffixes_result_7 = SearchAllSuffixesResult::NoMatches;
        let search_all_suffixes_result_8 = SearchAllSuffixesResult::NoMatches;

        assert_eq!(search_all_suffixes_result_1, search_all_suffixes_result_2);
        assert_ne!(search_all_suffixes_result_1, search_all_suffixes_result_3);
        assert_eq!(search_all_suffixes_result_4, search_all_suffixes_result_5);
        assert_ne!(search_all_suffixes_result_4, search_all_suffixes_result_6);
        assert_eq!(search_all_suffixes_result_7, search_all_suffixes_result_8);
        assert_ne!(search_all_suffixes_result_1, search_all_suffixes_result_7);
        assert_ne!(search_all_suffixes_result_4, search_all_suffixes_result_7);
    }

    fn get_example_proteins() -> Proteins {
        let input_string = "AI-CLACVAA-AC-KCRLY$";
        let text = ProteinText::from_string(input_string);

        Proteins::new(text, vec![
            Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            },
        ])
    }

    #[test]
    fn test_search_simple() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1);

        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

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
        let sa = SuffixArray::Original(vec![9, 0, 3, 12, 15, 6, 18], 3);

        let suffix_index_to_protein = SparseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

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
        let sa = SuffixArray::Original(vec![9, 0, 3, 12, 15, 6, 18], 3);

        let suffix_index_to_protein = DenseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

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
        let sa = SuffixArray::Original(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1);

        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

        let bounds_res = searcher.search_bounds(b"I");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((13, 16)));

        // search bounds 'RIZ' with equal I and L
        let bounds_res = searcher.search_bounds(b"RIY");
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((17, 18)));
    }

    #[test]
    fn test_il_equality_sparse() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![9, 0, 3, 12, 15, 6, 18], 3);

        let suffix_index_to_protein = SparseSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![0, 2, 4], 2);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![6, 0, 1, 5, 4, 3, 2], 1);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![6, 5, 4, 3, 2, 1, 0], 1);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![6, 4, 2, 0], 2);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![6, 5, 4, 3, 2, 1, 0], 1);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

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

        let sparse_sa = SuffixArray::Original(vec![13, 3, 12, 11, 1, 4, 2, 5, 9, 8, 6, 10, 0, 7], 1);
        let suffix_index_to_protein = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        let found_suffixes_1 = searcher.search_matching_suffixes(b"PAA", usize::MAX, false, true);
        assert_eq!(found_suffixes_1, SearchAllSuffixesResult::SearchResult(vec![0]));

        let found_suffixes_2 = searcher.search_matching_suffixes(b"APAA", usize::MAX, false, true);
        assert_eq!(found_suffixes_2, SearchAllSuffixesResult::SearchResult(vec![9]));
    }
}
