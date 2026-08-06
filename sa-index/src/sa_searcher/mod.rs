mod batched;
mod orchestrate;
mod retrieval;
mod scalar;
#[cfg(test)]
mod test_helpers;

pub use orchestrate::DEFAULT_MLP_BATCH;

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

/// The query-invariant half of the tryptic-cut predicate, resolved once per peptide instead of
/// once per candidate match.
///
/// A candidate that reaches the tryptic check has already been confirmed to equal the search
/// string over `[match_start, match_end)` — either by the suffix array range itself or by the
/// `check_prefix` / `check_suffix` conjuncts that run before it. So the two text positions the
/// naive predicate reads *inside* the match are redundant: the peptide already knows them.
///
/// * `text[match_start]` is `search_string[0]`
/// * `text[match_end - 1]` is `search_string.last()`
///
/// The substitution holds under both `equate_il` settings: K, R and P are untouched by the
/// L → I normalization the index is built with, so a text position can only hold P/K/R when the
/// corresponding peptide character does.
#[derive(Clone, Copy)]
enum TrypticQuery {
    /// Not filtering for tryptic matches — every candidate passes the boundary check.
    Off,
    /// Filtering, with both query-invariant terms already resolved.
    On {
        /// `search_string[0] != b'P'`, standing in for `text[match_start] != b'P'`.
        /// Proline directly after a K/R blocks the trypsin cut at the N-terminus.
        first_not_proline: bool,
        /// `search_string.last()` is K or R, standing in for `text[match_end - 1]` being a
        /// trypsin cut site.
        last_is_kr: bool
    },
    /// Degenerate zero-length query: `match_start == match_end`, so neither substitution above
    /// exists and we fall back to the original four-read formulation. Unreachable in production
    /// (`search_proteins_for_peptide` drops peptides shorter than the sparseness factor), but
    /// `search_matching_suffixes` is public, so stay bit-exact for it anyway.
    ZeroLength
}

impl TrypticQuery {
    #[inline]
    fn new(tryptic: bool, search_string: &[u8]) -> Self {
        if !tryptic {
            return Self::Off;
        }
        match (search_string.first(), search_string.last()) {
            (Some(&first), Some(&last)) => Self::On {
                first_not_proline: first != b'P',
                last_is_kr: last == b'K' || last == b'R'
            },
            _ => Self::ZeroLength
        }
    }
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

    /// Tryptic-boundary check for a candidate match spanning `[match_start, match_end)`.
    ///
    /// Equivalent to
    /// `(check_start_of_protein(ms) || check_tryptic_cut(ms)) && (check_end_of_protein(me) || check_tryptic_cut(me))`
    /// but reads **two** text positions instead of four: only `ms - 1` and `me`. The other two
    /// (`ms` and `me - 1`) sit inside the match, where the text equals the peptide, so `query`
    /// supplies them without touching memory (see `TrypticQuery`).
    ///
    /// This is the hottest loop in the index: `tryptic` disables the fast path in `scalar.rs`,
    /// so the candidate scan has to run until `max_matches` *accepted* matches accumulate —
    /// work that scales as `max_matches / acceptance_rate`. Halving the reads also halves the
    /// prefetch hints `prefetch_match_positions` has to keep in flight, which matters more than
    /// the ALU saving (see the comment there).
    #[inline]
    fn check_tryptic_boundaries(
        &self,
        text: &P::Text,
        match_start: usize,
        match_end: usize,
        query: TrypticQuery
    ) -> bool {
        match query {
            TrypticQuery::Off => true,
            TrypticQuery::On { first_not_proline, last_is_kr } => {
                // N-terminus. The `match_start == 0` guard has to come first: without it the
                // `match_start - 1` read underflows. (The original formulation was only safe
                // because `check_start_of_protein` short-circuits the `||` on that same test.)
                let n_term_ok = match_start == 0 || {
                    let before = text.get(match_start - 1);
                    before == SEPARATION_CHARACTER
                        || ((before == b'K' || before == b'R') && first_not_proline)
                };

                // C-terminus. `text[match_end]` answers protein-end and proline-block at once.
                n_term_ok && {
                    let after = text.get(match_end);
                    after == TERMINATION_CHARACTER
                        || after == SEPARATION_CHARACTER
                        || (last_is_kr && after != b'P')
                }
            }
            // Zero-length query: `match_end == match_start`, so the peptide has no character to
            // stand in for either read and the original formulation is the only correct one.
            TrypticQuery::ZeroLength => {
                (self.check_start_of_protein(match_start) || self.check_tryptic_cut(match_start))
                    && (self.check_end_of_protein(match_end) || self.check_tryptic_cut(match_end))
            }
        }
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

    /// Issues hardware prefetch hints for the text positions that will be accessed during
    /// validation of a candidate match at [ms, me). Called in Pass 1 of the two-pass batching
    /// loop to give DRAM latency hiding time before Pass 2 reads.
    ///
    /// Two hints, not four. The tryptic boundary check now only reads `ms - 1` and `me`
    /// (`check_tryptic_boundaries`), and those two hints also cover what `check_prefix` and
    /// `check_suffix` read inside `[ms, me)`: the text is 5-bit packed, so one 64-byte cache
    /// line holds ~102 characters, and the hints for `ms - 1` and `me` pull in the lines at
    /// both ends of the span. Any peptide short enough to span at most two lines (~100 aa —
    /// every realistic query) is therefore fully covered.
    ///
    /// The count matters more than the arithmetic it saves: `iterate_sa_range` uses
    /// `BATCH_SIZE = 64`, so four hints per candidate meant 256 outstanding prefetches against
    /// the ~10-12 line-fill buffers a core has, ~20x oversubscribed — most hints were evicted
    /// before Pass 2 could use them. Two per candidate doubles the effective prefetch depth.
    #[inline]
    fn prefetch_match_positions(text: &P::Text, ms: usize, me: usize) {
        text.prefetch_at(ms.saturating_sub(1));
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
        tryptic: TrypticQuery,
    ) -> Option<i64> {
        let suffix = raw as usize;
        if suffix < skip { return None; }
        let match_start = suffix - skip;
        let match_end = suffix + search_string.len() - skip;
        // Order matters: the tryptic check stays last because it assumes the span
        // [match_start, match_end) has already been confirmed equal to `search_string`.
        let valid = (skip == 0
            || Self::check_prefix(prefix, ProteinTextSlice::new(text, match_start, suffix), equate_il))
            && Self::check_suffix(skip, il_locations, suffix_str, ProteinTextSlice::new(text, suffix, match_end), equate_il)
            && self.check_tryptic_boundaries(text, match_start, match_end, tryptic);
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

        // Resolve the query-invariant half of the tryptic predicate once per peptide rather
        // than once per candidate — this loop runs `max_matches / acceptance_rate` times.
        let tryptic = TrypticQuery::new(tryptic, search_string);

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
    use text_compression::{ProteinText, ProteinTextBackend as _};

    use super::{SearchAllSuffixesResult, TrypticQuery};
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

    // Builds a searcher over an arbitrary protein text. The tryptic boundary checks only read
    // the text, so the SA content is irrelevant — an identity permutation keeps it well-formed.
    fn searcher_over(text: &str) -> Searcher<SuffixArray> {
        let protein_text = ProteinText::from_string(text);
        let proteins = Proteins::new(protein_text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![],
        }]);
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let sa = SuffixArray::Original(OriginalSA((0..text.len() as i64).collect(), 1));
        Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp))
    }

    // The original four-text-read formulation, kept verbatim as the oracle for the equivalence
    // test below. `check_start_of_protein` MUST stay on the left of the `||`: it is the only
    // thing that stops `check_tryptic_cut(0)` from reading `text[-1]`.
    fn old_tryptic_predicate(searcher: &Searcher<SuffixArray>, match_start: usize, match_end: usize) -> bool {
        (searcher.check_start_of_protein(match_start) || searcher.check_tryptic_cut(match_start))
            && (searcher.check_end_of_protein(match_end) || searcher.check_tryptic_cut(match_end))
    }

    // A text exercising every boundary case the tryptic filter can hit:
    //
    //   index 0..: K  P  A  R  C  K  A  -  R  P  K  R  A  -  P  K  -  A  I  L  K  $
    //
    // - match_start == 0 (the underflow guard)
    // - protein starts right after the separators at 7, 13 and 16
    // - protein ends at the separators 7/13/16 and at the terminator 21
    // - K/R-preceded cuts (positions 1, 4, 6, 11, 12, 16, 21)
    // - P-blocked cuts ("KP" at 0-1, "RP" at 8-9, "KP"-like "K-" vs "AP" at 13-14)
    // - a single-residue protein ("K" between the separators at 13 and 16 is "PK")
    // - I and L adjacent at 18/19, so peptides can differ from the text at L/I positions
    const BOUNDARY_TEXT: &str = "KPARCKA-RPKRA-PK-AILK$";

    // Every candidate the filter can ever see: for each span of the text, the peptide IS the
    // text over that span, which is exactly the invariant the two dropped reads rely on.
    // The new two-read predicate must agree with the old four-read one on every single one.
    #[test]
    fn test_tryptic_boundaries_equivalent_to_original() {
        let searcher = searcher_over(BOUNDARY_TEXT);
        let text = searcher.proteins.text();
        let text_len = BOUNDARY_TEXT.len();

        let mut checked = 0usize;
        let mut accepted = 0usize;
        for match_start in 0..text_len {
            // `match_end == text_len` would read one past the end in BOTH formulations, and
            // cannot occur for a real match: the text always ends in the terminator, which no
            // peptide contains.
            for match_end in (match_start + 1)..text_len {
                let peptide: Vec<u8> = (match_start..match_end).map(|i| text.get(i)).collect();
                let query = TrypticQuery::new(true, &peptide);

                let expected = old_tryptic_predicate(&searcher, match_start, match_end);
                let actual = searcher.check_tryptic_boundaries(text, match_start, match_end, query);
                assert_eq!(
                    actual, expected,
                    "mismatch at [{match_start}, {match_end}) for peptide {:?}",
                    String::from_utf8_lossy(&peptide)
                );

                // L/I equating: the index stores L as I, so a peptide may carry an L where the
                // text holds an I. K, R and P survive that normalization untouched, so swapping
                // I <-> L in the peptide must not move the predicate either way.
                let swapped: Vec<u8> = peptide
                    .iter()
                    .map(|&c| match c {
                        b'I' => b'L',
                        b'L' => b'I',
                        other => other,
                    })
                    .collect();
                let swapped_actual = searcher.check_tryptic_boundaries(
                    text, match_start, match_end, TrypticQuery::new(true, &swapped),
                );
                assert_eq!(
                    swapped_actual, expected,
                    "I/L swap changed the verdict at [{match_start}, {match_end})"
                );

                checked += 1;
                accepted += expected as usize;
            }
        }

        // Guard against the test silently degenerating: the text must produce a meaningful mix
        // of accepted and rejected candidates, otherwise the equivalence proves nothing.
        assert_eq!(checked, text_len * (text_len - 1) / 2);
        assert!(accepted > 0 && accepted < checked, "degenerate fixture: {accepted}/{checked} accepted");
    }

    // The zero-length fallback keeps reading the text, so it must also match the oracle.
    //
    // `cut == 0` is skipped because it panics identically in both formulations: with
    // `match_start == match_end == 0`, `check_end_of_protein(0)` is false for a text that does
    // not start with a separator, and the `check_tryptic_cut(0)` that follows it underflows on
    // `text[-1]`. That is pre-existing behaviour of the degenerate zero-length path (only the
    // N-terminal `||` has a `match_start == 0` guard, the C-terminal one does not), and it is
    // unreachable in production, where empty peptides are dropped before the search.
    #[test]
    fn test_tryptic_boundaries_zero_length_matches_original() {
        let searcher = searcher_over(BOUNDARY_TEXT);
        let text = searcher.proteins.text();

        for cut in 1..BOUNDARY_TEXT.len() {
            assert_eq!(
                searcher.check_tryptic_boundaries(text, cut, cut, TrypticQuery::new(true, b"")),
                old_tryptic_predicate(&searcher, cut, cut),
                "zero-length mismatch at {cut}"
            );
        }
    }

    // `tryptic = false` short-circuits the whole check, without touching the text.
    #[test]
    fn test_tryptic_boundaries_off_accepts_everything() {
        let searcher = searcher_over(BOUNDARY_TEXT);
        let text = searcher.proteins.text();

        for match_start in 0..BOUNDARY_TEXT.len() {
            for match_end in match_start..BOUNDARY_TEXT.len() {
                let query = TrypticQuery::new(false, b"PAAP");
                assert!(searcher.check_tryptic_boundaries(text, match_start, match_end, query));
            }
        }
    }

    // The two query-invariant terms are read off the peptide itself, so they must be right for
    // every combination of first/last residue that changes a verdict.
    #[test]
    fn test_tryptic_query_hoisting() {
        assert!(matches!(TrypticQuery::new(false, b"PAAK"), TrypticQuery::Off));
        assert!(matches!(TrypticQuery::new(true, b""), TrypticQuery::ZeroLength));

        let flags = |peptide: &[u8]| match TrypticQuery::new(true, peptide) {
            TrypticQuery::On { first_not_proline, last_is_kr } => (first_not_proline, last_is_kr),
            _ => panic!("expected TrypticQuery::On for {:?}", String::from_utf8_lossy(peptide)),
        };

        // (first_not_proline, last_is_kr)
        assert_eq!(flags(b"PAAC"), (false, false)); // starts with P
        assert_eq!(flags(b"PAAK"), (false, true)); // starts with P, ends with K
        assert_eq!(flags(b"AAAP"), (true, false)); // ends with P
        assert_eq!(flags(b"AAAK"), (true, true)); // ends with K
        assert_eq!(flags(b"AAAR"), (true, true)); // ends with R
        assert_eq!(flags(b"AAAC"), (true, false)); // neither
        assert_eq!(flags(b"P"), (false, false)); // single residue: first == last
        assert_eq!(flags(b"K"), (true, true));
        assert_eq!(flags(b"R"), (true, true));
        // K/R are never L/I-normalized, so a lowercase-free peptide of I/L cannot be a cut site.
        assert_eq!(flags(b"ILI"), (true, false));

        // And the flags actually drive the verdict. Text "KPAAK-": the K at 0 precedes the
        // match, so the N-terminal cut hinges entirely on the peptide's own first residue.
        let searcher = searcher_over("KPAAK-AAAA$");
        let text = searcher.proteins.text();
        // [1, 5) is "PAAK": preceded by K, but the peptide starts with P -> N-term blocked.
        assert!(!searcher.check_tryptic_boundaries(text, 1, 5, TrypticQuery::new(true, b"PAAK")));
        // [2, 5) is "AAK": preceded by P, which is neither K/R nor a separator -> N-term fails.
        assert!(!searcher.check_tryptic_boundaries(text, 2, 5, TrypticQuery::new(true, b"AAK")));
        // [0, 5) is "KPAAK": match_start == 0 is a protein start, and the match ends on the
        // separator at 5 -> both ends valid.
        assert!(searcher.check_tryptic_boundaries(text, 0, 5, TrypticQuery::new(true, b"KPAAK")));
        // [0, 1) is "K": protein start, but text[1] is P, so the K cut is proline-blocked.
        assert!(!searcher.check_tryptic_boundaries(text, 0, 1, TrypticQuery::new(true, b"K")));
        // [6, 10) is "AAAA": protein start (separator at 5) but the peptide does not end in
        // K/R and text[10] is the terminator -> C-term valid via end-of-protein.
        assert!(searcher.check_tryptic_boundaries(text, 6, 10, TrypticQuery::new(true, b"AAAA")));
        // [6, 9) is "AAA": neither end-of-protein nor a K/R terminus -> rejected.
        assert!(!searcher.check_tryptic_boundaries(text, 6, 9, TrypticQuery::new(true, b"AAA")));
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
