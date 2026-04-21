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
        /// Returns true if `arr1` and `arr2` contains the same elements, the order of the elements
        /// is ignored
        ///
        /// # Arguments
        /// * `arr1` - The first array used in the comparison
        /// * `arr2` - The second array used in the comparison
        ///
        /// # Returns
        ///
        /// Returns true if arr1 and arr2 contains the same elements, the order of the elements is
        /// ignored
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

    /// Compares the `search_string` to the `suffix`
    /// During search this function performs extra logic since the suffix array is build with I ==
    /// L, while ` self.proteins.input_string` is the original text where I != L
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide being searched in the suffix array
    /// * `suffix` - The current suffix from the suffix array we are comparing with in the binary
    ///   search
    /// * `skip` - How many characters we can skip in the comparison because we already know these
    ///   match
    /// * `bound` - Indicates if we are searching for the min of max bound
    ///
    /// # Returns
    ///
    /// The first argument is true if `bound` == `Minimum` and `search_string` <= `suffix` or if
    /// `bound` == `Maximum` and `search_string` >= `suffix` The second argument indicates how
    /// far the `suffix` and `search_string` matched
    fn compare(&self, search_string: &[u8], suffix: i64, skip: usize, bound: BoundSearch) -> (bool, usize) {
        let text = self.proteins.text();
        let mut index_in_suffix = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;

        // Depending on if we are searching for the min of max bound our condition is different
        let condition_check = match bound {
            Minimum => |a: u8, b: u8| a < b,
            Maximum => |a: u8, b: u8| a > b
        };

        // match as long as possible
        while index_in_search_string < search_string.len()
            && index_in_suffix < text.len()
            && (search_string[index_in_search_string] == text.get(index_in_suffix)
                || (search_string[index_in_search_string] == b'L' && text.get(index_in_suffix) == b'I')
                || (search_string[index_in_search_string] == b'I' && text.get(index_in_suffix) == b'L'))
        {
            index_in_suffix += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty
        // search string should not be found)
        if !search_string.is_empty() {
            if index_in_search_string == search_string.len() {
                is_cond_or_equal = true
            } else if index_in_suffix < text.len() {
                // in our index every L was replaced by a I, so we need to replace them if we want
                // to search in the right direction
                let peptide_char = if search_string[index_in_search_string] == b'L' {
                    b'I'
                } else {
                    search_string[index_in_search_string]
                };

                let protein_char = if text.get(index_in_suffix) == b'L' {
                    b'I'
                } else {
                    text.get(index_in_suffix)
                };

                is_cond_or_equal = condition_check(peptide_char, protein_char);
            }
        }

        (is_cond_or_equal, index_in_search_string)
    }

    /// Searches for the minimum or maximum bound for a string in the suffix array
    ///
    /// # Arguments
    /// * `bound` - Indicates if we are searching the minimum or maximum bound
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// The first argument is true if a match was found
    /// The second argument indicates the index of the minimum or maximum bound for the match
    /// (depending on `bound`)
    /// Core binary search within the window `[left, right)` starting character comparisons
    /// at position `lcp_skip` (the first `lcp_skip` characters are known to match).
    ///
    /// All elements in `[left, right)` are guaranteed to share at least `lcp_skip` characters
    /// with `search_string`, so both LCP accumulators are initialised to `lcp_skip`.
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
            // Prefetch both potential next pivots before the blocking sa.get(center) call.
            // One fetch will be wasted; both are free (single non-blocking CPU instruction).
            // From iteration 2 onward the needed SA entry is already in L1/L2 cache.
            self.sa.prefetch_sa_index((lo + center) / 2);
            self.sa.prefetch_sa_index((center + hi) / 2);
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(center), skip, bound);

            found |= lcp_center == search_string.len();

            if retval && bound == Minimum || !retval && bound == Maximum {
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
        let (left, right, lcp_skip) = if let Some(table) = &self.kmer_table {
            if search_string.len() >= table.k {
                match table.lookup(&search_string[..table.k]) {
                    Some((lo, hi)) => (lo, hi + 1, table.k),
                    None => return BoundSearchResult::NoMatches,
                }
            } else {
                (0, self.sa.len(), 0)
            }
        } else {
            (0, self.sa.len(), 0)
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
        // Issue an early OS prefetch hint for the k-mer's SA range (skip=0 case).
        // This gives the OS more lead time to load those pages into the page cache
        // while we collect IL locations below, before binary search starts.
        if let Some(table) = &self.kmer_table {
            if search_string.len() >= table.k {
                if let Some((lo, hi)) = table.lookup(&search_string[..table.k]) {
                    self.sa.prefetch_sa_range(lo, hi + 1);
                }
            }
        }

        // Pre-allocate up to 1M entries (8 MB); fall back to empty Vec for very large or
        // unbounded max_matches (e.g. usize::MAX in tests) to avoid a capacity overflow.
        let mut matching_suffixes: Vec<i64> = if max_matches <= 1 << 20 {
            Vec::with_capacity(max_matches)
        } else {
            Vec::new()
        };
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        // Batch size for two-pass prefetch: large enough that all prefetches from pass 1
        // resolve before they are consumed in pass 2 (DRAM latency ~100 ns, each SA read ~3 ns,
        // so 64 entries × 3 ns = 192 ns gap is sufficient).
        // Allocated once outside the skip loop and reused (via .clear()) across iterations
        // to avoid repeated heap allocations for each skip value.
        const BATCH_SIZE: usize = 64;
        let mut raw_batch: Vec<i64> = Vec::with_capacity(BATCH_SIZE);

        let mut skip: usize = 0;
        while skip < self.sa.sample_rate() as usize {
            let mut il_locations_start = 0;
            while il_locations_start < il_locations.len() && il_locations[il_locations_start] < skip {
                il_locations_start += 1;
            }
            let il_locations_current_suffix = &il_locations[il_locations_start..];
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
                    // Generic path: iterate with two-pass prefetching to hide DRAM latency for
                    // protein-text accesses required by prefix checks, IL validation, or tryptic
                    // filtering.
                    //
                    // Pass 1 — collect a batch of raw SA entries from the iterator and issue
                    //          non-blocking hardware prefetches for the protein-text positions
                    //          each entry will access.
                    // Pass 2 — validate each buffered entry; the prefetches from pass 1 will
                    //          have had time to complete, so protein-text reads are cheap.
                    let text = self.proteins.text();
                    let mut sa_iter = self.sa.iter_range(min_bound, max_bound);

                    loop {
                        // --- Pass 1 ---
                        raw_batch.clear();
                        for s in &mut sa_iter {
                            let su = s as usize;
                            if su >= skip {
                                let ms = su - skip;
                                let me = su + search_string.len() - skip;
                                // Prefetch positions needed by prefix check, IL check, and
                                // tryptic checks respectively.
                                text.prefetch_at(ms.saturating_sub(1)); // tryptic start / K|R before
                                text.prefetch_at(ms);                   // prefix start
                                text.prefetch_at(me.saturating_sub(1)); // tryptic cut before end
                                text.prefetch_at(me);                   // tryptic end / separator
                            }
                            raw_batch.push(s);
                            if raw_batch.len() == BATCH_SIZE { break; }
                        }
                        if raw_batch.is_empty() { break; }

                        // --- Pass 2 ---
                        for &raw in &raw_batch {
                            let suffix = raw as usize;

                            if suffix >= skip {
                                let match_start = suffix - skip;
                                let match_end = suffix + search_string.len() - skip;

                                if (skip == 0
                                    || Self::check_prefix(
                                        current_search_string_prefix,
                                        ProteinTextSlice::new(text, match_start, suffix),
                                        equate_il
                                    ))
                                    && Self::check_suffix(
                                        skip,
                                        il_locations_current_suffix,
                                        current_search_string_suffix,
                                        ProteinTextSlice::new(text, suffix, match_end),
                                        equate_il
                                    )
                                    && (!tryptic
                                        || ((self.check_start_of_protein(match_start) || self.check_tryptic_cut(match_start))
                                            && (self.check_end_of_protein(match_end) || self.check_tryptic_cut(match_end))))
                                {
                                    matching_suffixes.push((suffix - skip) as i64);

                                    if matching_suffixes.len() >= max_matches {
                                        self.match_iter_ns.fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);
                                        return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                                    }
                                }
                            }
                        }
                    }
                    self.match_iter_ns.fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
            }
            skip += 1;
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

    /// Returns all the proteins that correspond with the provided suffixes
    ///
    /// # Arguments
    /// * `suffixes` - List of suffix indices
    ///
    /// # Returns
    ///
    /// Returns the proteins that every suffix is a part of
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &[i64]) -> Vec<ProteinRef<'_>> {
        // ── Pass 1: collect protein indices ──────────────────────────────────
        // Prefetch the mapping entry for suffix[i + N1] before reading suffix[i].
        // This hides the DRAM latency of the (random) mapping lookup behind N1
        // iterations of cheap computation.
        const N1: usize = 16;
        let mut protein_indices = Vec::with_capacity(suffixes.len());
        for (i, &suffix) in suffixes.iter().enumerate() {
            if let Some(&fs) = suffixes.get(i + N1) {
                self.suffix_index_to_protein.prefetch_for_suffix(fs);
            }
            protein_indices.push(self.suffix_index_to_protein.suffix_to_protein(suffix));
        }

        // ── Pass 2: collect ProteinRefs ──────────────────────────────────────
        // protein_indices is a small Vec (~135 × 4 bytes) that fits in L1 cache,
        // so protein_indices[i + N2] is free to read at iteration i. We use it
        // to prefetch the fixed-table entry for that future protein before the
        // blocking proteins.get() call for the current protein.
        const N2: usize = 16;
        let mut res = Vec::with_capacity(suffixes.len());
        for (i, &protein_index) in protein_indices.iter().enumerate() {
            if let Some(&fpi) = protein_indices.get(i + N2) {
                if !fpi.is_null() {
                    self.proteins.prefetch(fpi as usize);
                }
            }
            if !protein_index.is_null() {
                res.push(self.proteins.get(protein_index as usize));
            }
        }
        res
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
