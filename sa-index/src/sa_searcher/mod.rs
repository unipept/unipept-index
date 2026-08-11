mod batched;
pub(crate) mod metrics;
mod orchestrate;
mod retrieval;
mod scalar;
// The shared fixtures build a `SuffixArray::Original` / `SuffixToProteinMapping::BitVec`, which
// only exist in the preloaded configuration. Every consumer is gated the same way, so gating the
// module itself keeps `--features mmap` compiling.
#[cfg(all(test, not(feature = "mmap")))]
mod test_helpers;

pub use orchestrate::DEFAULT_MLP_BATCH;
use sa_mappings::proteins::{Proteins, ProteinsBackend, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::{ProteinTextBackend, ProteinTextSlice};

use crate::{
    KmerTable,
    array::SuffixArrayBackend,
    sa_searcher::{
        BoundSearch::{Maximum, Minimum},
        metrics::Counter
    },
    suffix_to_protein_index::{SuffixToProteinMapping, SuffixToProteinMappingBackend}
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

/// Characters that can legally precede a tryptic match: the two trypsin cut residues plus the
/// protein separator (a match at a protein start needs no cut).
///
/// Searching `X + peptide` for each of these — instead of truncating the peptide — is what makes
/// the tryptic path cheap; see `search_matching_suffixes` for the derivation. Ordered K, R,
/// separator so the two common cases run first.
const TRYPTIC_EXTENSION_CHARS: [u8; 3] = [b'K', b'R', b'-'];

/// The separator-only subset, used when the peptide starts with proline: `K|R` followed by P is
/// not a trypsin cut site, so those two searches are guaranteed empty and are skipped entirely.
const TRYPTIC_EXTENSION_CHARS_PROLINE: [u8; 1] = [b'-'];

/// Left-extension characters to search for this peptide.
#[inline]
fn tryptic_extension_chars(search_string: &[u8]) -> &'static [u8] {
    debug_assert_eq!(TRYPTIC_EXTENSION_CHARS[2], SEPARATION_CHARACTER);
    if search_string.first() == Some(&b'P') { &TRYPTIC_EXTENSION_CHARS_PROLINE } else { &TRYPTIC_EXTENSION_CHARS }
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

/// Upper bound on `validate_batch` — sizes the on-stack candidate buffer in
/// `iterate_sa_range`, so it must stay a compile-time constant.
///
/// 256 `i64`s is 2 KB of stack per call. The ceiling exists to bound that, not because the
/// sweep found anything at 256: the measured optimum is 64 and the curve is flat above it, so
/// this is headroom for re-tuning on other hardware rather than a value in use.
pub const MAX_VALIDATE_BATCH: usize = 256;

/// Cap on how much result space is pre-allocated per peptide, in suffixes.
///
/// Callers pass a `max_matches` cutoff that is an upper bound, not an estimate — the server's
/// default is 10 000 — while the overwhelming majority of peptides match a handful of times.
/// Reserving the full cutoff would allocate 80 KB per peptide and touch none of it; capping at
/// 4096 entries (32 KB) keeps the common case to one allocation while letting rare high-hit
/// peptides grow normally.
pub(crate) const MAX_RESULT_PREALLOC: usize = 4096;

/// Runtime-tunable batch and prefetch-lookahead sizes for the search and retrieval hot
/// paths. Every field is a pure performance knob: results are identical for any setting
/// (asserted by `test_tuning_does_not_change_results`).
///
/// All three defaults are confirmed by the run3 full-DB sweep (3 peptide-length buckets x
/// {preloaded, mmap} x {fast-path, validating} baselines, 20 reps, 3.9% noise floor). Two
/// knobs that the same sweep found dead were removed rather than left to be re-measured:
/// `retrieval_batch` (cross-query batched retrieval, median +1.7%) and `scalar_kmer_prefetch`
/// (+0.3%).
#[derive(Clone, Copy, Debug)]
pub struct SearchTuning {
    /// Candidates per two-pass validation batch in `iterate_sa_range`.
    /// Clamped to `1..=MAX_VALIDATE_BATCH`.
    ///
    /// The one knob measured to matter, and it is a cliff rather than a peak: 16 -> 32 gains
    /// ~10% in all 6 (bucket, backend) combos, then it plateaus. Against this default of 64,
    /// 32 wins in 1/6 and 128 wins in 5/6 but never above the noise floor, and 128 regresses
    /// long peptides on the preloaded backend by 2.7%. Do not lower it below 32.
    pub validate_batch: usize,
    /// Minimum SA range size before `iterate_sa_range` uses two-pass validation
    /// instead of a straight loop.
    ///
    /// Swept over {8, 16, 32, 64}: median full-range swing +0.9%, inside the noise floor
    /// everywhere. Left tunable for re-measurement on other hardware, not because 32 is
    /// known to be special.
    pub validate_prefetch_threshold: usize,
    /// Prefetch look-ahead distance (in suffixes) inside protein retrieval.
    ///
    /// Swept over {8, 16, 32, 64}: median full-range swing +1.2%, inside the noise floor
    /// everywhere. Same caveat as `validate_prefetch_threshold`.
    pub retrieval_prefetch_distance: usize,
    /// Issue `madvise(MADV_WILLNEED)` over an SA range before scanning it.
    ///
    /// **Off by default, and measured not worth enabling.** Under a memory ceiling it does remove
    /// 23-25% of major faults, but the throughput that buys decays to nothing as threads rise
    /// (+12.0% at the core count, ~0% at 96), because oversubscription already overlaps those
    /// faults; and it costs -3.7% with the index resident. Full numbers on
    /// `MmapBackedSA::advise_willneed_range`, which also records the -16.8% regression that got
    /// the first version of this removed.
    ///
    /// Kept as a knob for the two cases that could still favour it: slower storage, and running
    /// at the core count where `RAYON_NUM_THREADS` cannot be raised.
    ///
    /// A no-op on the preloaded backend, where the trait method is the default no-op.
    pub willneed: bool
}

impl Default for SearchTuning {
    fn default() -> Self {
        Self {
            validate_batch: 64,              // confirmed by run3; 16 costs ~10%, 128 gains nothing
            validate_prefetch_threshold: 32, // measured flat over 8..64
            retrieval_prefetch_distance: 32, // measured flat over 8..64
            willneed: false                  // regressed when resident; unmeasured under a cap
        }
    }
}

/// Everything needed to search a peptide against the index, plus the search itself.
///
/// The three generic parameters are the storage backends, which the `mmap` feature resolves at
/// compile time; see the crate docs. They default to the aliases for the active build, so most
/// callers write `Searcher<SuffixArray>`.
///
/// Construct with [`Searcher::new`], then optionally attach a k-mer table with
/// [`Searcher::with_kmer_table`]. The searcher is immutable during search and `Sync`, so one
/// instance serves every request.
pub struct Searcher<
    SA: SuffixArrayBackend,
    P: ProteinsBackend = Proteins,
    STPM: SuffixToProteinMappingBackend = SuffixToProteinMapping
> {
    pub sa: SA,
    pub proteins: P,
    pub suffix_index_to_protein: STPM,
    pub kmer_table: Option<KmerTable>,
    /// Batch and prefetch-lookahead sizes for the hot paths. Public so callers (the benchmark)
    /// can mutate it between configurations, exactly as they already do with `kmer_table`.
    pub tuning: SearchTuning,
    /// Total nanoseconds spent inside `search_bounds()` across all queries (since last drain).
    /// Only accumulated with the `metrics` feature enabled; a no-op ZST otherwise.
    pub(crate) search_bounds_ns: Counter,
    /// Total nanoseconds spent iterating matches in `search_matching_suffixes()` (since last drain).
    /// Only accumulated with the `metrics` feature enabled; a no-op ZST otherwise.
    pub(crate) match_iter_ns: Counter,
    /// Candidate suffixes inspected by `iterate_sa_range` (since last drain), i.e. every entry
    /// the SA-range scan looked at, accepted or not. `metrics` only.
    pub(crate) candidates_examined: Counter,
    /// Candidate suffixes `iterate_sa_range` accepted as real matches (since last drain).
    /// `metrics` only.
    ///
    /// Together with `candidates_examined` this settles why tryptic search is ~12.5x slower
    /// than non-tryptic on 5–10 aa peptides: a low accepted/examined ratio means the scan is
    /// simply sifting ~1/ratio times more candidates to reach `max_matches` (make each check
    /// cheaper), whereas a ratio near 1 with the cutoff rarely reached means whole SA ranges
    /// are being scanned to exhaustion (a `max_candidates` scan cap is the fix).
    pub(crate) candidates_accepted: Counter
}

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Builds a searcher over the three index structures.
    ///
    /// All three must come from the same index build: the suffix array indexes positions in
    /// `proteins`' text, and `suffix_index_to_protein` maps those same positions to indices into
    /// `proteins`. Mixing builds produces wrong answers rather than errors.
    ///
    /// The sparseness factor is not passed in — it is read from the suffix array, which records
    /// it in its file header.
    ///
    /// Starts with no k-mer table and default [`SearchTuning`].
    pub fn new(sa: SA, proteins: P, suffix_index_to_protein: STPM) -> Self {
        Self {
            sa,
            proteins,
            suffix_index_to_protein,
            kmer_table: None,
            tuning: SearchTuning::default(),
            search_bounds_ns: Counter::new(),
            match_iter_ns: Counter::new(),
            candidates_examined: Counter::new(),
            candidates_accepted: Counter::new()
        }
    }

    /// Returns `(search_bounds_ns, match_iter_ns)` accumulated since the last call and resets both
    /// counters to zero.  Safe to call concurrently with ongoing searches (relaxed ordering).
    ///
    /// Present in both feature configurations so callers need no `cfg`; without the `metrics`
    /// feature it always returns `(0, 0)`.
    pub fn drain_timing_ns(&self) -> (u64, u64) {
        (self.search_bounds_ns.drain(), self.match_iter_ns.drain())
    }

    /// Returns `(candidates_examined, candidates_accepted)` accumulated by `iterate_sa_range`
    /// since the last call and resets both counters to zero. Same contract as
    /// `drain_timing_ns`: always present, always `(0, 0)` without the `metrics` feature.
    ///
    /// The ratio is the SA-range scan's acceptance rate; the tryptic paths are the interesting
    /// ones, since a non-tryptic I/L-free query never enters `iterate_sa_range` at all.
    pub fn drain_candidate_counts(&self) -> (u64, u64) {
        (self.candidates_examined.drain(), self.candidates_accepted.drain())
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
            Maximum => |a: u8, b: u8| a > b
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
                bound_satisfied =
                    condition_check(Self::normalize_li(search_string[i_search]), Self::normalize_li(text.get(i_text)));
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
                    before == SEPARATION_CHARACTER || ((before == b'K' || before == b'R') && first_not_proline)
                };

                // C-terminus. `text[match_end]` answers protein-end and proline-block at once.
                n_term_ok && {
                    let after = text.get(match_end);
                    after == TERMINATION_CHARACTER || after == SEPARATION_CHARACTER || (last_is_kr && after != b'P')
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
    /// The count matters more than the arithmetic it saves: `iterate_sa_range` batches 64
    /// candidates by default, so four hints per candidate meant 256 outstanding prefetches against
    /// the ~10-12 line-fill buffers a core has, ~20x oversubscribed — most hints were evicted
    /// before Pass 2 could use them. Two per candidate doubles the effective prefetch depth.
    #[inline]
    fn prefetch_match_positions(text: &P::Text, ms: usize, me: usize) {
        text.prefetch_at(ms.saturating_sub(1));
        text.prefetch_at(me);
    }

    /// Checks whether the candidate suffix `raw` is a valid match for the current search
    /// parameters. Returns `Some(match_start)` when valid, `None` otherwise.
    ///
    /// The parameters are the query state hoisted out of the candidate loop; bundling them into a
    /// struct would put a field load on the innermost path, so they stay positional.
    #[allow(clippy::too_many_arguments)]
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
        tryptic: TrypticQuery
    ) -> Option<i64> {
        let suffix = raw as usize;
        if suffix < skip {
            return None;
        }
        let match_start = suffix - skip;
        let match_end = suffix + search_string.len() - skip;
        // Order matters: the tryptic check stays last because it assumes the span
        // [match_start, match_end) has already been confirmed equal to `search_string`.
        let valid = (skip == 0
            || Self::check_prefix(prefix, ProteinTextSlice::new(text, match_start, suffix), equate_il))
            && Self::check_suffix(
                skip,
                il_locations,
                suffix_str,
                ProteinTextSlice::new(text, suffix, match_end),
                equate_il
            )
            && self.check_tryptic_boundaries(text, match_start, match_end, tryptic);
        if valid { Some(match_start as i64) } else { None }
    }

    /// Checks only the C-terminal half of the tryptic predicate.
    ///
    /// Used by the left-extended search path, where the N-terminal half holds by construction:
    /// the search string was `X + peptide` for `X` in {K, R, separator}, so the character before
    /// the match *is* `X`. Nothing to read and nothing to check.
    ///
    /// `last_is_kr` is the query-invariant stand-in for `text[match_end - 1]` — see
    /// [`TrypticQuery`] for why the substitution is sound.
    ///
    /// # On reading `text[match_end]`
    ///
    /// `match_end` is one past the match, so this reads the character *after* it. For a real
    /// peptide that is always in bounds: the text ends with [`TERMINATION_CHARACTER`], and a
    /// match reaching `text.len()` would have to have consumed the terminator — that is, the
    /// query's own last character would be `$`, which no protein sequence contains.
    ///
    /// A caller *can* construct such a query, since `$` is in the alphabet. Both backends
    /// tolerate the resulting `text.get(text.len())`: the preloaded one reads zero padding inside
    /// the last allocated word, and the mmap one reads inside the page-rounded mapping. Neither
    /// faults, and the value read cannot make the predicate accept, because such a query ends in
    /// `$` and so has `last_is_kr == false`. It is nonetheless an unchecked read that happens to
    /// land somewhere harmless rather than a designed guarantee — worth knowing before changing
    /// how the text is stored or how far the query alphabet extends.
    #[inline]
    fn check_tryptic_c_term(&self, text: &P::Text, match_end: usize, last_is_kr: bool) -> bool {
        let after = text.get(match_end);
        after == TERMINATION_CHARACTER || after == SEPARATION_CHARACTER || (last_is_kr && after != b'P')
    }

    /// Validates a hit from a left-extended search (`X + search_string`).
    ///
    /// The SA entry sits at `match_start - 1` (the `X`), so `match_start = raw + 1` — the mirror
    /// image of `validate_candidate`'s `suffix - skip`, which is why this cannot reuse it.
    ///
    /// Two conjuncts of the normal predicate are gone:
    /// * the prefix check, because the SA range already matched `X + search_string` in full, so
    ///   `text[match_start..match_end]` equals `search_string` modulo L → I;
    /// * the N-terminal tryptic check, free by construction (see `check_tryptic_c_term`).
    ///
    /// What remains is the I/L check (only when `equate_il` is false) and one text read at
    /// `match_end`.
    #[inline]
    fn validate_extended_candidate(
        &self,
        text: &P::Text,
        raw: i64,
        search_string: &[u8],
        il_locations: &[usize],
        equate_il: bool,
        last_is_kr: bool
    ) -> Option<i64> {
        let match_start = raw as usize + 1;
        let match_end = match_start + search_string.len();
        // skip = 0: il_locations are absolute positions in `search_string`, and the slice starts
        // at match_start, so the two index spaces already line up.
        let valid = Self::check_suffix(
            0,
            il_locations,
            search_string,
            ProteinTextSlice::new(text, match_start, match_end),
            equate_il
        ) && self.check_tryptic_c_term(text, match_end, last_is_kr);
        if valid { Some(match_start as i64) } else { None }
    }

    /// Scans the SA range of a left-extended search, validating each hit.
    /// Returns `true` if `max_matches` was reached.
    ///
    /// Same two-pass prefetch shape as `iterate_sa_range`, but the hint set is smaller because
    /// the validation is: `match_end` always (the C-terminal read), plus `match_start` only when
    /// `equate_il` is false and the I/L check will actually read the span.
    #[allow(clippy::too_many_arguments)]
    fn iterate_extended_sa_range(
        &self,
        mut sa_iter: impl Iterator<Item = i64>,
        range_size: usize,
        text: &P::Text,
        search_string: &[u8],
        il_locations: &[usize],
        equate_il: bool,
        last_is_kr: bool,
        matching_suffixes: &mut Vec<i64>,
        max_matches: usize
    ) -> bool {
        let batch_size = self.tuning.validate_batch.clamp(1, MAX_VALIDATE_BATCH);
        let prefetch_threshold = self.tuning.validate_prefetch_threshold;
        let needs_span = !equate_il && !il_locations.is_empty();

        let mut examined = 0u64;
        let mut accepted = 0u64;

        let hit_max = 'scan: {
            if range_size < prefetch_threshold {
                for raw in sa_iter {
                    examined += 1;
                    if let Some(v) =
                        self.validate_extended_candidate(text, raw, search_string, il_locations, equate_il, last_is_kr)
                    {
                        accepted += 1;
                        matching_suffixes.push(v);
                        if matching_suffixes.len() >= max_matches {
                            break 'scan true;
                        }
                    }
                }
                break 'scan false;
            }

            let mut raw_batch = [0i64; MAX_VALIDATE_BATCH];

            loop {
                let mut batch_len = 0usize;
                for s in &mut sa_iter {
                    let ms = s as usize + 1;
                    text.prefetch_at(ms + search_string.len());
                    if needs_span {
                        text.prefetch_at(ms);
                    }
                    raw_batch[batch_len] = s;
                    batch_len += 1;
                    if batch_len == batch_size {
                        break;
                    }
                }
                if batch_len == 0 {
                    break;
                }

                for &raw in &raw_batch[..batch_len] {
                    examined += 1;
                    if let Some(v) =
                        self.validate_extended_candidate(text, raw, search_string, il_locations, equate_il, last_is_kr)
                    {
                        accepted += 1;
                        matching_suffixes.push(v);
                        if matching_suffixes.len() >= max_matches {
                            break 'scan true;
                        }
                    }
                }
            }
            false
        };

        self.candidates_examined.add(examined);
        self.candidates_accepted.add(accepted);
        hit_max
    }

    /// Iterates the SA range and validates each candidate suffix.
    /// Returns `true` if `max_matches` was reached, `false` otherwise.
    ///
    /// Two-pass batching with hardware prefetch hints to hide DRAM latency:
    /// Pass 1 — fills a batch and issues prefetch hints for the text positions to validate.
    /// Pass 2 — validates candidates after prefetches have had time to complete.
    ///
    /// Same rationale as `validate_candidate` for the positional parameter list.
    #[allow(clippy::too_many_arguments)]
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
        max_matches: usize
    ) -> bool {
        // Default 64, tuned on x86_64 Zen4/Intel Sapphire Rapids: DRAM latency ~80–100 ns, one
        // SA entry read per ~2–3 ns at that cache level → 64 entries ≈ 192 ns gap, comfortably
        // above the latency floor. Runtime-settable because that reasoning does not transfer to
        // ARM or NVMe-backed mmap; the buffer stays on the stack (a heap Vec here would cost
        // more than the batching saves), so the fill length is clamped to the array's size.
        let batch_size = self.tuning.validate_batch.clamp(1, MAX_VALIDATE_BATCH);
        // Default 32: the minimum range for a prefetch to resolve before use. Below this the
        // two-pass overhead exceeds the latency-hiding benefit.
        let prefetch_threshold = self.tuning.validate_prefetch_threshold;

        // Resolve the query-invariant half of the tryptic predicate once per peptide rather
        // than once per candidate — this loop runs `max_matches / acceptance_rate` times.
        let tryptic = TrypticQuery::new(tryptic, search_string);

        // Local counters, folded into the shared atomics exactly once below. A per-candidate
        // RMW on a cache line every rayon worker shares would dominate the loop it measures.
        let mut examined = 0u64;
        let mut accepted = 0u64;

        let hit_max = 'scan: {
            if range_size < prefetch_threshold {
                for raw in sa_iter {
                    examined += 1;
                    if let Some(v) = self.validate_candidate(
                        text,
                        raw,
                        skip,
                        search_string,
                        prefix,
                        suffix_str,
                        il_locations,
                        equate_il,
                        tryptic
                    ) {
                        accepted += 1;
                        matching_suffixes.push(v);
                        if matching_suffixes.len() >= max_matches {
                            break 'scan true;
                        }
                    }
                }
                break 'scan false;
            }

            let mut raw_batch = [0i64; MAX_VALIDATE_BATCH];

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
                    if batch_len == batch_size {
                        break;
                    }
                }
                if batch_len == 0 {
                    break;
                }

                // --- Pass 2: validate (prefetches have had time to complete) ---
                for &raw in &raw_batch[..batch_len] {
                    examined += 1;
                    if let Some(v) = self.validate_candidate(
                        text,
                        raw,
                        skip,
                        search_string,
                        prefix,
                        suffix_str,
                        il_locations,
                        equate_il,
                        tryptic
                    ) {
                        accepted += 1;
                        matching_suffixes.push(v);
                        if matching_suffixes.len() >= max_matches {
                            break 'scan true;
                        }
                    }
                }
            }
            false
        };

        self.candidates_examined.add(examined);
        self.candidates_accepted.add(accepted);
        hit_max
    }
}

#[cfg(all(test, not(feature = "mmap")))]
mod tests {
    use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
    use text_compression::{ProteinText, ProteinTextBackend as _};

    use super::{MAX_VALIDATE_BATCH, SearchAllSuffixesResult, SearchTuning, TrypticQuery};
    use crate::{
        SuffixArray,
        array::OriginalSA,
        sa_searcher::{
            BoundSearchResult, Searcher,
            test_helpers::{example_searcher, searcher_over_text}
        },
        suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping}
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
            functional_annotations: vec![]
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
            functional_annotations: vec![]
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
                    actual,
                    expected,
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
                        other => other
                    })
                    .collect();
                let swapped_actual =
                    searcher.check_tryptic_boundaries(text, match_start, match_end, TrypticQuery::new(true, &swapped));
                assert_eq!(swapped_actual, expected, "I/L swap changed the verdict at [{match_start}, {match_end})");

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
            _ => panic!("expected TrypticQuery::On for {:?}", String::from_utf8_lossy(peptide))
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
    // Without the `metrics` feature the counters are no-op ZSTs, so the drain is always (0, 0)
    // — the API stays present either way so callers need no `cfg`.
    #[test]
    fn test_drain_timing_ns() {
        let searcher = example_searcher();
        assert_eq!(searcher.drain_timing_ns(), (0, 0));

        searcher.search_bounds_ns.store(123);
        searcher.match_iter_ns.store(456);
        let expected = if cfg!(feature = "metrics") { (123, 456) } else { (0, 0) };
        assert_eq!(searcher.drain_timing_ns(), expected);
        assert_eq!(searcher.drain_timing_ns(), (0, 0)); // reset after draining
    }

    // Same contract for the candidate counters.
    #[test]
    fn test_drain_candidate_counts() {
        let searcher = example_searcher();
        assert_eq!(searcher.drain_candidate_counts(), (0, 0));

        searcher.candidates_examined.store(70);
        searcher.candidates_accepted.store(7);
        let expected = if cfg!(feature = "metrics") { (70, 7) } else { (0, 0) };
        assert_eq!(searcher.drain_candidate_counts(), expected);
        assert_eq!(searcher.drain_candidate_counts(), (0, 0)); // reset after draining
    }

    // With `metrics` on, `iterate_sa_range` must actually count what it scans — and only what
    // *it* scans: the counters exist to measure the acceptance rate of the validating path, so
    // the fast path (which accepts a whole SA range without inspecting entries) must not
    // contribute. A text of all 'L' searched for "I" with equate_il=false enters the validating
    // path and rejects every candidate: 70 examined, 0 accepted.
    #[cfg(feature = "metrics")]
    #[test]
    fn test_candidate_counts_are_accumulated() {
        let n = 70usize;
        let mut input = "L".repeat(n);
        input.push('$');
        let text = ProteinText::from_string(&input);
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);
        let sa = SuffixArray::Original(OriginalSA((0..=n as i64).rev().collect(), 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // "I" matches all 70 'L' positions during the bound search (L is normalized to I), but
        // equate_il=false rejects every one of them: acceptance rate 0.
        searcher.drain_candidate_counts();
        assert_eq!(
            searcher.search_matching_suffixes(b"I", usize::MAX, false, false),
            SearchAllSuffixesResult::NoMatches
        );
        assert_eq!(searcher.drain_candidate_counts(), (n as u64, 0));

        // Same range with equate_il=true accepts everything it examines — but takes the fast
        // path, which bypasses iterate_sa_range entirely, so nothing is counted at all.
        searcher.search_matching_suffixes(b"I", usize::MAX, true, false);
        assert_eq!(searcher.drain_candidate_counts(), (0, 0));
    }

    // The two-pass path's stack buffer is sized by the compile-time MAX_VALIDATE_BATCH; a
    // larger runtime `validate_batch` must clamp to it, not index past the array.
    #[test]
    fn test_validate_batch_clamps_to_max() {
        let n = 300usize; // > MAX_VALIDATE_BATCH, so the batch loop refills at least once
        let mut input = "L".repeat(n);
        input.push('$');
        let text = ProteinText::from_string(&input);
        let proteins = Proteins::new(text, vec![Protein {
            uniprot_id: String::new(),
            taxon_id: 0,
            functional_annotations: vec![]
        }]);
        let sa = SuffixArray::Original(OriginalSA((0..=n as i64).rev().collect(), 1));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let mut searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // "L" enters the validating path (equate_il=false, peptide holds an L) and matches
        // every position, so the scan runs over the whole 300-entry range.
        let expected = searcher.search_matching_suffixes(b"L", usize::MAX, false, false);

        for validate_batch in [MAX_VALIDATE_BATCH, MAX_VALIDATE_BATCH + 1, usize::MAX, 0] {
            searcher.tuning = SearchTuning { validate_batch, ..SearchTuning::default() };
            assert_eq!(
                searcher.search_matching_suffixes(b"L", usize::MAX, false, false),
                expected,
                "validate_batch = {validate_batch} changed the result"
            );
        }
    }

    // ── SearchTuning equivalence ────────────────────────────────────────────────────────
    //
    // Every SearchTuning field is a pure performance knob, so the entire grid must produce
    // byte-identical search *and* retrieval output. This is the guard on that claim: without
    // it, a sweep that changed results would look like a speed-up.

    /// 500 short proteins built from a rotating set of motifs. The repetition is what makes the
    /// fixture useful: a 2-mer query lands an SA range well above the default two-pass
    /// threshold (32) while the longer, rarer peptides stay below it, so one run covers both
    /// branches of `iterate_sa_range`. K/R/P give the tryptic filter a real mix of accepts and
    /// rejects, and the I/L pairs exercise both `equate_il` settings.
    fn tuning_fixture_text() -> String {
        let motifs = ["AAKR", "AILK", "PAAK", "CVAAR", "KCRLY", "AAIP"];
        let mut s = String::new();
        for i in 0..500 {
            s.push_str(motifs[i % motifs.len()]);
            s.push_str(motifs[(i * 3 + 1) % motifs.len()]);
            s.push('-');
        }
        // Every motif recurs 20-40 times, so nothing built from them alone yields a *small* SA
        // range. One unique protein, from residues that appear nowhere else, supplies the
        // below-threshold ranges that exercise the straight (non-batched) scan.
        s.push_str(RARE_PROTEIN);
        s.push('$');
        s
    }

    /// The one non-repeated protein in `tuning_fixture_text`. Ends in K so it also has a valid
    /// tryptic C-terminus, and carries an L so `equate_il` is not moot for it.
    const RARE_PROTEIN: &str = "MWYFQDNTLHK";

    /// Peptides spanning the interesting shapes: very common 1–3-mers (large SA ranges, so the
    /// two-pass path and the `max_matches` cutoff), unique ones from `RARE_PROTEIN` (ranges of
    /// 1–2 entries, so the straight scan), I/L-carrying and I/L-free ones, K/R- and
    /// P-terminated ones, and a total miss.
    const TUNING_PEPTIDES: &[&[u8]] = &[
        b"A",
        b"AA",
        b"AAK",
        b"K",
        b"KR",
        b"IL",
        b"LI",
        b"AILK",
        b"CVAAR",
        b"KCRLY",
        b"AAKRAAIP",
        b"PAAK",
        b"AAIP",
        b"MWYFQ",
        b"MWYFQDNTLHK",
        b"NTLHK",
        b"NTIHK",
        b"ZZ"
    ];

    /// Runs the full pipeline — search then retrieval — over `peptides` on both the scalar
    /// (MLP batch 1) and batched (16) orchestration paths, and reduces the outcome to a plainly
    /// comparable value: the result discriminant, the matched suffixes in the order they were
    /// produced, and each retrieved protein's taxon + accession (`ProteinRef` borrows, so it is
    /// not `PartialEq` itself).
    ///
    /// `max_matches` is deliberately finite and larger than the smallest swept `validate_batch`
    /// so the cutoff fires mid-batch for the common peptides.
    /// One row per (mlp_batch, peptide): the batch size, the result-variant tag, the sorted
    /// matching suffixes, and the retrieved `(taxon_id, uniprot_id)` pairs. Comparing these
    /// rows across tuning settings is what proves the knobs are behaviour-neutral.
    type TuningRow = (usize, &'static str, Vec<i64>, Vec<(u32, String)>);

    fn tuning_run(
        searcher: &Searcher<SuffixArray>,
        peptides: &[&[u8]],
        equate_il: bool,
        tryptic: bool
    ) -> Vec<TuningRow> {
        let mut out = Vec::new();
        for mlp_batch in [1usize, 16] {
            let results = searcher.search_all_matching_suffixes(peptides, 64, equate_il, tryptic, mlp_batch);
            let tags: Vec<&'static str> = results
                .iter()
                .map(|r| match r {
                    SearchAllSuffixesResult::NoMatches => "none",
                    SearchAllSuffixesResult::MaxMatches(_) => "max",
                    SearchAllSuffixesResult::SearchResult(_) => "hit"
                })
                .collect();
            let suffixes: Vec<Vec<i64>> = results
                .iter()
                .map(|r| match r {
                    SearchAllSuffixesResult::NoMatches => vec![],
                    SearchAllSuffixesResult::MaxMatches(v) | SearchAllSuffixesResult::SearchResult(v) => v.clone()
                })
                .collect();

            let proteins: Vec<Vec<(u32, String)>> = suffixes
                .iter()
                .map(|v| searcher.retrieve_proteins(v).iter().map(|p| (p.taxon_id, p.uniprot_id.to_string())).collect())
                .collect();

            for i in 0..peptides.len() {
                out.push((mlp_batch, tags[i], suffixes[i].clone(), proteins[i].clone()));
            }
        }
        out
    }

    #[test]
    fn test_tuning_does_not_change_results() {
        let text = tuning_fixture_text();

        // Dense SA (skip is always 0), sparse SA (exercises skip = 1, 2 and therefore the
        // prefix/suffix validation), and a k-mer-narrowed one — which is the only
        // configuration in which `scalar_kmer_prefetch` has anything to switch off.
        let mut kmered = searcher_over_text(&text, 1);
        kmered.build_kmer_table(3);
        // The third element is the sparseness factor: `search_matching_suffixes` indexes
        // `search_string[skip..]` for every skip below it, so peptides shorter than that are
        // out of contract (production drops them before the search) and must be filtered out.
        let mut fixtures = [
            ("dense", searcher_over_text(&text, 1), 1usize),
            ("sparse", searcher_over_text(&text, 3), 3usize),
            ("kmer", kmered, 1usize)
        ];

        // The fixture must actually reach both scan paths, otherwise sweeping validate_batch
        // and validate_prefetch_threshold proves nothing. "AA" must exceed the largest swept
        // threshold (and the largest swept batch, so the batch refills); "MWYFQ" must fall
        // below the smallest non-zero one.
        match fixtures[0].1.search_bounds(b"AA") {
            BoundSearchResult::SearchResult((lo, hi)) => {
                assert!(hi - lo > MAX_VALIDATE_BATCH, "fixture too small: 'AA' range is only {}", hi - lo)
            }
            other => panic!("expected 'AA' to match, got {other:?}")
        }
        match fixtures[0].1.search_bounds(b"MWYFQ") {
            BoundSearchResult::SearchResult((lo, hi)) => {
                assert!(hi - lo < 8, "'MWYFQ' range {} is not below the threshold", hi - lo)
            }
            other => panic!("expected 'MWYFQ' to match, got {other:?}")
        }

        for (name, searcher, min_len) in fixtures.iter_mut() {
            let peptides: Vec<&[u8]> = TUNING_PEPTIDES.iter().copied().filter(|p| p.len() >= *min_len).collect();

            for equate_il in [false, true] {
                for tryptic in [false, true] {
                    searcher.tuning = SearchTuning::default();
                    let expected = tuning_run(searcher, &peptides, equate_il, tryptic);

                    // Every peptide matching nothing would make the comparison vacuous.
                    assert!(
                        expected.iter().any(|(_, tag, _, _)| *tag != "none"),
                        "{name}: fixture matches nothing (il={equate_il} tr={tryptic})"
                    );

                    for validate_batch in [1usize, 16, 64, MAX_VALIDATE_BATCH] {
                        for validate_prefetch_threshold in [0usize, 8, 32] {
                            for retrieval_prefetch_distance in [1usize, 8, 32] {
                                // `willneed` only decides whether readahead is *requested* for a
                                // range about to be read anyway, so it must never change which
                                // suffixes come back. It is swept here for exactly that reason.
                                for willneed in [false, true] {
                                    let tuning = SearchTuning {
                                        validate_batch,
                                        validate_prefetch_threshold,
                                        retrieval_prefetch_distance,
                                        willneed
                                    };
                                    searcher.tuning = tuning;
                                    assert_eq!(
                                        tuning_run(searcher, &peptides, equate_il, tryptic),
                                        expected,
                                        "{name}: results changed (il={equate_il} tr={tryptic}) \
                                         for {tuning:?}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
