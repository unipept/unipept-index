//! Searching peptides against the index.
//!
//! This file holds the [`Searcher`] itself and the primitives every search path shares: `compare`,
//! which orders a peptide against a suffix; `opening_window`, which picks the SA window a bound
//! search starts from and is the one place the k-mer table is consulted; and the candidate scan —
//! `iterate_sa_range` and its left-extended twin — which walks an SA range and validates what it
//! finds. Everything a search needs beyond that lives in a sibling:
//!
//! * `scalar` / `batched` — one peptide at a time versus B interleaved for memory-level
//!   parallelism. `batched` also owns `search_all_matching_suffixes_batched`, the entry point every
//!   multi-peptide caller goes through.
//! * `retrieval` — turning matched suffix positions into proteins.
//! * `tryptic` — the trypsin-cut filter the candidate scan applies.
//! * `tuning` — the constants the hot paths are fixed at, and what measured them;
//!   `results` — what a search returns.
//! * `measure` — the instrumentation, which costs nothing unless the `measure` feature is on.

mod batched;
pub(crate) mod measure;
mod results;
mod retrieval;
mod scalar;
#[cfg(test)]
mod test_utils;
mod tryptic;
mod tuning;

use protein_metadata::ProteinsBackend;
use protein_text::{ProteinTextBackend, ProteinTextSlice};
pub use results::{BoundSearchResult, SearchAllSuffixesResult};
use tryptic::TrypticQuery;
pub(crate) use tuning::{
    MAX_RESULT_PREALLOC, MLP_BATCH, PREFETCH_THRESHOLD, RETRIEVAL_PREFETCH_DISTANCE, VALIDATE_BATCH
};

use crate::{
    KmerTable,
    array::SuffixArrayBackend,
    kmer_table::KmerLookup,
    sa_searcher::{
        BoundSearch::{Maximum, Minimum},
        measure::SearchMeasurements
    },
    suffix_to_protein_index::SuffixToProteinMappingBackend
};

/// Enum indicating if we are searching for the minimum, or maximum bound in the suffix array
#[derive(Clone, Copy, PartialEq)]
enum BoundSearch {
    Minimum,
    Maximum
}

/// Everything needed to search a peptide against the index, plus the search itself.
///
/// The three generic parameters are the storage backends. Both implementations of each are always
/// compiled and nothing here names one, so a caller is free to combine them however it likes;
/// `sa-server` picks one combination per build in its `backends` module.
///
/// Construct with [`Searcher::new`], then optionally attach a k-mer table with
/// [`Searcher::with_kmer_table`]. The searcher is immutable during search and `Sync`, so one
/// instance serves every request.
pub struct Searcher<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> {
    /// The suffix array, which also carries the sparseness factor the index was built with.
    pub sa: SA,
    /// The protein metadata table and, through it, the concatenated protein text.
    pub proteins: P,
    /// Resolves a matched text position to an index into `proteins`.
    pub suffix_index_to_protein: STPM,
    /// Optional accelerator narrowing the opening binary search; see [`KmerTable`]. Attached
    /// through [`Searcher::with_kmer_table`] or [`Searcher::build_kmer_table`], never at
    /// construction.
    pub kmer_table: Option<KmerTable>,
    /// Instrumentation counters, drained through [`Searcher::drain_timing_ns`] and
    /// [`Searcher::drain_candidate_counts`]. Zero-sized without the `measure` feature.
    pub(crate) measurements: SearchMeasurements
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
    /// Starts with no k-mer table. Everything else about how a search runs is a constant; see
    /// the `tuning` module.
    pub fn new(sa: SA, proteins: P, suffix_index_to_protein: STPM) -> Self {
        Self::try_new(sa, proteins, suffix_index_to_protein).unwrap_or_else(|err| panic!("{err}"))
    }

    /// Builds a searcher, rejecting three structures that did not come from the same index build.
    ///
    /// The three files are loaded independently, each well-formed on its own, and nothing used to
    /// relate them: point `--mapping-file` at an older build and every load succeeds, the server
    /// reports itself ready, and the answers are silently wrong — protein indices resolved against
    /// the wrong text.
    ///
    /// There is no build identifier in the format to compare, and adding one would invalidate every
    /// existing index. Two relationships already hold implicitly, though, and both are exact:
    ///
    /// * the suffix array samples the text, so it holds `ceil(text_len / sample_rate)` entries;
    /// * a dense or bitvec mapping records one entry (or bit) per text position, so its length *is*
    ///   the text length. Sparse records protein starts and cannot be compared, so it is skipped
    ///   rather than guessed at.
    ///
    /// A different database changes the text length, which moves both. That does not prove three
    /// files share a build — two different databases could coincide in length — but it catches the
    /// realistic mistake, at no cost and without touching the format.
    pub fn try_new(sa: SA, proteins: P, suffix_index_to_protein: STPM) -> Result<Self, String> {
        let text_len = proteins.text().len();
        let sample_rate = sa.sample_rate() as usize;

        if sample_rate == 0 {
            return Err("the suffix array declares a sample rate of 0".to_string());
        }

        let expected = text_len.div_ceil(sample_rate);
        if sa.len() != expected {
            return Err(format!(
                "the suffix array holds {} entries, but the protein text is {text_len} characters at a \
                 sample rate of {sample_rate}, which needs {expected}; either these files are from \
                 different builds, or sa.bin was built at a sparseness factor an older sa-builder \
                 sampled incorrectly (8, 9, 16, 25 and 27 produced a smaller stride than their header \
                 declared) and needs rebuilding",
                sa.len()
            ));
        }

        if matches!(suffix_index_to_protein.implied_text_len(), Some(n) if n != text_len) {
            let mapping_len = suffix_index_to_protein.implied_text_len().unwrap_or(text_len);
            return Err(format!(
                "the suffix-to-protein mapping was built for a text of {mapping_len} characters, but the \
                 protein text is {text_len}; these files are from different builds"
            ));
        }

        Ok(Self {
            sa,
            proteins,
            suffix_index_to_protein,
            kmer_table: None,
            measurements: SearchMeasurements::new()
        })
    }

    /// Attaches a pre-built k-mer bounds table to this searcher.
    pub fn with_kmer_table(self, table: KmerTable) -> Self {
        self.try_with_kmer_table(table).unwrap_or_else(|err| panic!("{err}"))
    }

    /// Attaches a k-mer table, rejecting one that was not built from this suffix array.
    ///
    /// A table records SA index ranges, and the searcher feeds them straight into the binary search
    /// as bounds. A table from a *different* build therefore points past the end of a smaller array:
    /// the preloaded backend panics on the first such lookup and the mmap one reads past its data
    /// region and returns a fabricated suffix position, which then indexes the suffix-to-protein
    /// map. Two different wrong behaviours, both arriving mid-query rather than at startup.
    ///
    /// Mixing builds cannot be detected in general — `sa-builder` writes no build identifier — but
    /// this one mismatch is both the likeliest (rebuild the index, forget the table) and cheap to
    /// catch, since the bound is accumulated while the table is read.
    pub fn try_with_kmer_table(mut self, table: KmerTable) -> Result<Self, String> {
        if table.highest_bound() >= self.sa.len() {
            return Err(format!(
                "the k-mer table points at suffix-array index {} but the array holds {} entries; \
                 it was built from a different index",
                table.highest_bound(),
                self.sa.len()
            ));
        }
        self.kmer_table = Some(table);
        Ok(self)
    }

    /// Builds a k-mer table from this searcher's own suffix array and text, and attaches it.
    ///
    /// The build-mismatch check [`Self::try_with_kmer_table`] performs is unnecessary here and is
    /// not run: a table built from `self.sa` cannot point past `self.sa`. Use that method instead
    /// when the table comes from a file.
    ///
    /// Costs a full parallel pass over the suffix array plus a dense `AMINO_ACID_COUNT^k` table —
    /// ~127 MB at k=5. `sa-builder` normally writes the table once and the server loads it; this
    /// is for callers holding an index with no table file beside it.
    ///
    /// # Panics
    ///
    /// If `k` exceeds [`MAX_KMER_K`](crate::kmer_table::MAX_KMER_K).
    pub fn build_kmer_table(&mut self, k: usize) {
        self.kmer_table = Some(KmerTable::build_from_sa(&self.sa, self.proteins.text(), k));
    }

    /// The SA window a bound search for `search_string` should open with.
    ///
    /// `Some((left, right, lcp_skip))` is a half-open window plus the number of leading characters
    /// every entry in it is already known to share, so the bound search never re-compares them.
    /// `None` means the k-mer table settled the question: no suffix carries this prefix.
    ///
    /// This is the *only* place the k-mer table is consulted, and it is shared by both search
    /// paths — `scalar::search_bounds_scalar` and `batched::search_bounds_batched` — so the
    /// two cannot disagree about which queries the table may answer for.
    ///
    /// # The three answers, and why the table needs all three
    ///
    /// [`KmerLookup`] separates "this k-mer occurs nowhere" from "this k-mer is not in my
    /// alphabet", and only the first is a result:
    ///
    /// * [`Absent`](KmerLookup::Absent) → `None`. The table covers the k-mer and the index does not
    ///   hold it.
    /// * [`NotRepresentable`](KmerLookup::NotRepresentable) → the full array. `ALPHABET` in
    ///   `kmer_table.rs` holds amino acids only, so a k-mer containing the protein separator `-` or
    ///   the terminator `$` has no bucket and the table abstains rather than answering.
    /// * shorter than `k`, or no table attached → the full array, for the same reason: nothing to
    ///   look up.
    ///
    /// The abstention is not a corner case. The left-extended tryptic search looks up
    /// `'-' + peptide` for every tryptic query (see [`TrypticQuery`] and `tryptic_extension_chars`),
    /// so it lands here once per peptide. Reading that as "no matches" would silently drop every
    /// protein-start match — roughly 3 % of all tryptic hits, and every protein's N-terminal
    /// peptide. Both search paths used to carry a second, table-free bounds function called from a
    /// `prefix_char == SEPARATION_CHARACTER` branch at the call site; those four pieces are what
    /// this replaces. Bounded by `scalar::tests::test_extended_protein_start_with_kmer_table` and
    /// its batched twin.
    ///
    /// # Cost of the fallback
    ///
    /// A full-height binary search: ~35 random SA reads instead of ~12 at the default k=5 on the
    /// full UniProt index — see [`KmerTable`]. Measured on the full DB that is affordable: even
    /// 26-50aa tryptic queries, where this fixed cost is the largest share of a ~4 µs budget, came
    /// out 1.40x *faster* overall (run5), because the left-extension replaced a far more expensive
    /// truncated pass. Adding `-` to `ALPHABET` would remove the fallback at +29 MB and a rebuild
    /// of every table file — not worth it while the fallback is this cheap.
    #[inline]
    fn opening_window(&self, search_string: &[u8]) -> Option<(usize, usize, usize)> {
        let full_range = (0, self.sa.len(), 0);
        let Some(table) = &self.kmer_table else {
            return Some(full_range);
        };
        if search_string.len() < table.k {
            return Some(full_range);
        }
        match table.lookup(&search_string[..table.k]) {
            KmerLookup::Range(lo, hi) => Some((lo, hi + 1, table.k)),
            KmerLookup::Absent => None,
            KmerLookup::NotRepresentable => Some(full_range)
        }
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
            } else {
                // The *suffix* ran out first, so it is a proper prefix of the query and therefore
                // sorts below it. Both bounds have to agree on that: `Minimum` wants
                // `query <= suffix`, which is false, and `Maximum` wants `query >= suffix`, which
                // is true. Leaving this case to fall through gave `false` for both, so the two
                // bounds drew opposite conclusions about the same pair.
                //
                // It is not reachable today — a query can only outlive the text by containing `$`
                // at a non-final position, and the maximum pass is gated on the query being a full
                // prefix of some suffix, which forces `$` to be last. Those are mutually exclusive
                // while the text holds exactly one `$`. Written out anyway, so the correctness of
                // the bound search stops depending on that coincidence.
                bound_satisfied = matches!(bound, Maximum);
            }
        }

        (bound_satisfied, i_search)
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
        let needs_span = !equate_il && !il_locations.is_empty();

        let mut examined = 0u64;
        let mut accepted = 0u64;

        let hit_max = 'scan: {
            if range_size < PREFETCH_THRESHOLD {
                for raw in sa_iter {
                    examined += 1;
                    if let Some(v) =
                        self.validate_extended_candidate(text, raw, search_string, il_locations, equate_il, last_is_kr)
                    {
                        accepted += 1;
                        matching_suffixes.push(v);
                        if matching_suffixes.len() > max_matches {
                            break 'scan true;
                        }
                    }
                }
                break 'scan false;
            }

            let mut raw_batch = [0i64; VALIDATE_BATCH];

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
                    if batch_len == VALIDATE_BATCH {
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
                        if matching_suffixes.len() > max_matches {
                            break 'scan true;
                        }
                    }
                }
            }
            false
        };

        self.measurements.candidates_examined.add(examined);
        self.measurements.candidates_accepted.add(accepted);
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
        // `VALIDATE_BATCH` sizes the on-stack buffer below (a heap Vec here would cost more than
        // the batching saves) and `PREFETCH_THRESHOLD` decides whether this range is worth
        // batching at all. Both are fixed; see `tuning` for the sweep that could not move them.

        // Resolve the query-invariant half of the tryptic predicate once per peptide rather
        // than once per candidate — this loop runs `max_matches / acceptance_rate` times.
        let tryptic = TrypticQuery::new(tryptic, search_string);

        // Local counters, folded into the shared atomics exactly once below. A per-candidate
        // RMW on a cache line every rayon worker shares would dominate the loop it measures.
        let mut examined = 0u64;
        let mut accepted = 0u64;

        let hit_max = 'scan: {
            if range_size < PREFETCH_THRESHOLD {
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
                        if matching_suffixes.len() > max_matches {
                            break 'scan true;
                        }
                    }
                }
                break 'scan false;
            }

            let mut raw_batch = [0i64; VALIDATE_BATCH];

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
                    if batch_len == VALIDATE_BATCH {
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
                        if matching_suffixes.len() > max_matches {
                            break 'scan true;
                        }
                    }
                }
            }
            false
        };

        self.measurements.candidates_examined.add(examined);
        self.measurements.candidates_accepted.add(accepted);
        hit_max
    }
}

#[cfg(test)]
mod tests {
    use protein_metadata::ProteinsBackend as _;

    use super::{
        BoundSearch::{Maximum, Minimum},
        SearchAllSuffixesResult
    };
    use crate::{
        KmerTable,
        array::{InMemorySA, MmapBackedSA},
        sa_searcher::test_utils::{
            Fingerprint, MappedMetaOwnedText, MappedProteins, OwnedMetaMappedText, OwnedProteins, example_searcher,
            fingerprint, repeated_residue_searcher, searcher_over_text
        },
        suffix_to_protein_index::{InMemorySuffixToProteinMapping, MmapBackedSuffixToProteinMapping}
    };

    // `compare` is the ordering primitive the whole binary search rests on. Every other test
    // reaches it only through a bound search, which reports a range and hides both halves of
    // what it returns: whether the bound holds, and how far the two strings agreed.
    #[test]
    fn test_compare_reports_bound_and_match_length() {
        // A(0) C(1) G(2) K(3) L(4) $(5)
        let searcher = searcher_over_text("ACGKL$", 1);

        // When the peptide is a prefix of the suffix it satisfies *both* bounds, even though
        // neither strict inequality holds: it belongs inside the range from either side, which
        // is what lets one comparison drive both ends of the binary search.
        assert_eq!(searcher.compare(b"ACG", 0, 0, Minimum), (true, 3));
        assert_eq!(searcher.compare(b"ACG", 0, 0, Maximum), (true, 3));

        // "ACG" < "CGKL$": below the suffix, so the Minimum bound holds and Maximum does not.
        // They agree on nothing, hence a match length of 0.
        assert_eq!(searcher.compare(b"ACG", 1, 0, Minimum), (true, 0));
        assert_eq!(searcher.compare(b"ACG", 1, 0, Maximum), (false, 0));

        // "ACK" > "ACG..": above the suffix, and they agree on the first two characters.
        assert_eq!(searcher.compare(b"ACK", 0, 0, Maximum), (true, 2));
        assert_eq!(searcher.compare(b"ACK", 0, 0, Minimum), (false, 2));

        // `skip` asserts the first n characters already matched, so the scan resumes there and
        // the reported length still counts from the start of the peptide.
        assert_eq!(searcher.compare(b"ACK", 0, 2, Maximum), (true, 2));

        // L is normalised to I on both sides, so a peptide spelling it either way compares equal
        // — both are full matches against the 'L' at position 4.
        assert_eq!(searcher.compare(b"L", 4, 0, Minimum), (true, 1));
        assert_eq!(searcher.compare(b"I", 4, 0, Minimum), (true, 1));
    }

    // The two-pass candidate scan: a range above `PREFETCH_THRESHOLD` (32) with more candidates
    // than `VALIDATE_BATCH` (64), so the batch loop refills at least once. `equate_il=false` on
    // an I/L-free peptide is what keeps the fast path from swallowing the whole range.
    #[test]
    fn test_iterate_sa_range_two_pass() {
        let n = 70usize;
        let searcher = repeated_residue_searcher('A', n);

        let found = searcher.search_matching_suffixes_scalar(b"A", usize::MAX, false, false);
        let expected: Vec<i64> = (0..n as i64).collect();
        assert_eq!(found, SearchAllSuffixesResult::SearchResult(expected));
    }

    // `with_kmer_table` is the builder half of the k-mer API — `build_kmer_table` is what every
    // other test uses — and attaching a table must narrow the search without changing what it
    // finds.
    #[test]
    fn test_with_kmer_table_attaches_without_changing_results() {
        let plain = example_searcher();
        let table = KmerTable::build_from_sa(&plain.sa, plain.proteins.text(), 3);
        let tabled = example_searcher().with_kmer_table(table);

        assert!(tabled.kmer_table.is_some(), "table not attached");
        for peptide in [&b"VAA"[..], b"KCR", b"AC", b"A", b"ZZZ"] {
            assert_eq!(
                tabled.search_matching_suffixes_scalar(peptide, usize::MAX, false, false),
                plain.search_matching_suffixes_scalar(peptide, usize::MAX, false, false),
                "with_kmer_table changed the result for {:?}",
                std::str::from_utf8(peptide).unwrap()
            );
        }
    }

    #[test]
    fn every_backend_combination_returns_identical_results() {
        let expected = fingerprint::<InMemorySA, OwnedProteins, InMemorySuffixToProteinMapping>();

        // An agreement test over a fingerprint that stopped distinguishing anything would pass
        // forever, so check the reference has both kinds of row before comparing against it.
        assert!(
            expected.iter().any(|row| !row.6.is_empty()),
            "the fixture retrieved no proteins — it can no longer tell the backends apart"
        );
        assert!(expected.iter().any(|row| row.4 == "none"), "the fixture has no miss to check");

        // Spelled out rather than generated: sixteen lines that a failure can name, against a macro
        // whose expansion nothing could read.
        let combinations: Vec<(&str, Fingerprint)> = vec![
            (
                "sa=owned  proteins=owned      mapping=mapped",
                fingerprint::<InMemorySA, OwnedProteins, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=owned/map  mapping=owned ",
                fingerprint::<InMemorySA, OwnedMetaMappedText, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=owned/map  mapping=mapped",
                fingerprint::<InMemorySA, OwnedMetaMappedText, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=map/owned  mapping=owned ",
                fingerprint::<InMemorySA, MappedMetaOwnedText, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=map/owned  mapping=mapped",
                fingerprint::<InMemorySA, MappedMetaOwnedText, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=mapped     mapping=owned ",
                fingerprint::<InMemorySA, MappedProteins, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=owned  proteins=mapped     mapping=mapped",
                fingerprint::<InMemorySA, MappedProteins, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=owned      mapping=owned ",
                fingerprint::<MmapBackedSA, OwnedProteins, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=owned      mapping=mapped",
                fingerprint::<MmapBackedSA, OwnedProteins, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=owned/map  mapping=owned ",
                fingerprint::<MmapBackedSA, OwnedMetaMappedText, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=owned/map  mapping=mapped",
                fingerprint::<MmapBackedSA, OwnedMetaMappedText, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=map/owned  mapping=owned ",
                fingerprint::<MmapBackedSA, MappedMetaOwnedText, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=map/owned  mapping=mapped",
                fingerprint::<MmapBackedSA, MappedMetaOwnedText, MmapBackedSuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=mapped     mapping=owned ",
                fingerprint::<MmapBackedSA, MappedProteins, InMemorySuffixToProteinMapping>()
            ),
            (
                "sa=mapped proteins=mapped     mapping=mapped",
                fingerprint::<MmapBackedSA, MappedProteins, MmapBackedSuffixToProteinMapping>()
            ),
        ];

        for (name, rows) in combinations {
            assert_eq!(rows, expected, "{name} disagrees with the fully-owned combination");
        }
    }
}
