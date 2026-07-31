mod batched;
mod retrieval;
mod scalar;
#[cfg(test)]
mod test_helpers;

use std::sync::atomic::{AtomicU64, Ordering};

use sa_mappings::proteins::{Proteins, ProteinsBackend, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::{ProteinTextBackend, ProteinTextSlice};

use crate::{
    KmerTable, array::SuffixArrayBackend,
    sa_searcher::BoundSearch::{Maximum, Minimum},
    suffix_to_protein_index::{
        SuffixToProteinMappingBackend, SuffixToProteinMapping,
    },
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
pub struct Searcher<SA: SuffixArrayBackend, P: ProteinsBackend = Proteins, STPM: SuffixToProteinMappingBackend = SuffixToProteinMapping> {
    pub sa: SA,
    pub proteins: P,
    pub suffix_index_to_protein: STPM,
    pub kmer_table: Option<KmerTable>,
    /// Total nanoseconds spent inside `search_bounds()` across all queries (since last drain).
    pub search_bounds_ns: AtomicU64,
    /// Total nanoseconds spent iterating matches in `search_matching_suffixes()` (since last drain).
    pub match_iter_ns: AtomicU64,
}

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
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
    pub fn new(sa: SA, proteins: P, suffix_index_to_protein: STPM) -> Self {
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
        let text_len = text.len();
        let mut i_text = (suffix as usize) + skip;
        let mut i_search = skip;
        let mut bound_satisfied = false;

        let condition_check = match bound {
            Minimum => |a: u8, b: u8| a < b,
            Maximum => |a: u8, b: u8| a > b,
        };

        // Advance while characters match (treating L == I).
        while i_search < search_string.len()
            && i_text < text_len
            && Self::normalize_li(search_string[i_search]) == Self::normalize_li(text.get(i_text))
        {
            i_text += 1;
            i_search += 1;
        }

        if !search_string.is_empty() {
            if i_search == search_string.len() {
                bound_satisfied = true;
            } else if i_text < text_len {
                // The index has L replaced by I, so normalize both sides before ordering.
                bound_satisfied = condition_check(
                    Self::normalize_li(search_string[i_search]),
                    Self::normalize_li(text.get(i_text)),
                );
            }
        }

        (bound_satisfied, i_search)
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
    fn check_prefix(search_string_prefix: &[u8], index_prefix: ProteinTextSlice<'_, P::Text>, equate_il: bool) -> bool {
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
    #[inline]
    fn check_suffix(
        skip: usize,
        il_locations: &[usize],
        search_string: &[u8],
        text_slice: ProteinTextSlice<'_, P::Text>,
        equate_il: bool
    ) -> bool {
        if equate_il { true } else { text_slice.check_il_locations(skip, il_locations, search_string) }
    }

    /// Issues hardware prefetch hints for all text positions that will be accessed
    /// during validation of a candidate match at [ms, me). Called in Pass 1 of the
    /// two-pass batching loop to give DRAM latency hiding time before Pass 2 reads.
    #[inline]
    fn prefetch_match_positions(text: &P::Text, ms: usize, me: usize) {
        text.prefetch_at(ms.saturating_sub(1));
        text.prefetch_at(ms);
        text.prefetch_at(me.saturating_sub(1));
        text.prefetch_at(me);
    }

    /// Issues an early OS prefetch hint for the k-mer's SA range (skip=0 case), giving the OS
    /// more lead time to load those pages into the page cache before binary search starts.
    #[inline]
    fn prefetch_kmer_range(&self, search_string: &[u8]) {
        if let Some(table) = &self.kmer_table {
            if search_string.len() >= table.k {
                if let Some((lo, hi)) = table.lookup(&search_string[..table.k]) {
                    self.sa.prefetch_sa_range(lo, hi + 1);
                }
            }
        }
    }

    /// Checks whether the candidate suffix `raw` is a valid match for the current search
    /// parameters. Returns `Some(match_start)` when valid, `None` otherwise.
    #[inline]
    fn validate_candidate(
        &self,
        text: &P::Text,
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
        text: &P::Text,
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

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use std::sync::atomic::Ordering;

    use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
    use text_compression::ProteinText;

    use super::SearchAllSuffixesResult;
    use crate::{
        array::OriginalSA,
        sa_searcher::{test_helpers::get_example_proteins, Searcher},
        suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping},
        SuffixArray,
    };

    // A full suffix array over the example proteins (positions never L/I-normalized here).
    fn example_searcher() -> Searcher<SuffixArray> {
        let proteins = get_example_proteins();
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1));
        Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp))
    }

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

    // Direct test of the protein-boundary checks used by tryptic filtering.
    // Text "AI-CLACVAA-AC-KCRLY$": separators at 2/10/13, termination at 19.
    #[test]
    fn test_check_protein_boundaries() {
        let searcher = example_searcher();

        // start of a protein: index 0, or immediately after a separator
        for i in [0usize, 3, 11, 14] {
            assert!(searcher.check_start_of_protein(i), "expected protein start at {i}");
        }
        for i in [1usize, 5, 12] {
            assert!(!searcher.check_start_of_protein(i), "unexpected protein start at {i}");
        }

        // end of a protein: the position itself holds a separator or termination char
        for i in [2usize, 10, 13, 19] {
            assert!(searcher.check_end_of_protein(i), "expected protein end at {i}");
        }
        for i in [0usize, 1, 5] {
            assert!(!searcher.check_end_of_protein(i), "unexpected protein end at {i}");
        }
    }

    // A tryptic cut is valid iff preceded by K or R and NOT followed by proline (P).
    #[test]
    fn test_check_tryptic_cut() {
        // K(0) A(1) R(2) C(3) K(4) P(5) D(6) $(7)
        let text = ProteinText::from_string("KARCKPD$");
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![],
        }]);
        let stp = BitVecSuffixToProtein::new(proteins.text());
        // check_tryptic_cut only reads the text, so the SA content is irrelevant here.
        let sa = SuffixArray::Original(OriginalSA(vec![0, 1, 2, 3, 4, 5, 6, 7], 1));
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        assert!(searcher.check_tryptic_cut(1)); // after K, C follows
        assert!(searcher.check_tryptic_cut(3)); // after R, C follows
        assert!(!searcher.check_tryptic_cut(5)); // after K but P follows (proline blocks)
        assert!(!searcher.check_tryptic_cut(2)); // preceded by A, not K/R
    }

    // drain_timing_ns returns the accumulated counters and resets them to zero.
    #[test]
    fn test_drain_timing_ns() {
        let searcher = example_searcher();
        assert_eq!(searcher.drain_timing_ns(), (0, 0));

        searcher.search_bounds_ns.store(123, Ordering::Relaxed);
        searcher.match_iter_ns.store(456, Ordering::Relaxed);
        assert_eq!(searcher.drain_timing_ns(), (123, 456));
        assert_eq!(searcher.drain_timing_ns(), (0, 0)); // reset after draining
    }
}
