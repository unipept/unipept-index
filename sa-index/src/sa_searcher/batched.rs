//! Cross-query batching for memory-level parallelism (MLP), and the entry point that applies it.
//!
//! Interleaves B independent peptide searches so ~B random memory misses are in flight
//! per core, hiding the DRAM latency the scalar dependent-chain binary search cannot.
//! Each stream's per-step logic is identical to the scalar path in the parent module, so
//! results are identical (see `test_batched_matches_scalar`).
//!
//! Two layers, both here. `search_matching_suffixes_batched` is the kernel: single-threaded,
//! rayon-free, one chunk of peptides in lockstep. `search_all_matching_suffixes` is the entry
//! point every caller uses: it splits a whole peptide list into [`MLP_BATCH`]-sized chunks and
//! runs them across rayon. Both the production path (`peptide_search::search_all_peptides`) and
//! the benchmark go through it, so the decomposition is chosen in one place and measured as it
//! ships — the kernel is `pub(crate)` precisely so nothing outside can choose a different one.
//!
//! There is no scalar/batched choice to make any more. That used to be a runtime knob, and this
//! module was two: `orchestrate` held the branch and `batched` the kernel. `MLP_BATCH` is a
//! constant now, so batching is unconditional and the branch is gone. The scalar primitive is not
//! dead — `peptide_search::search_proteins_for_peptide` is the single-query path and calls
//! `search_matching_suffixes_scalar` directly — it simply is not something this path selects between.
//!
//! Retrieval deliberately has no counterpart here. A cross-query batched retrieval was built
//! and measured (run3): the `retrieval_batch` sweep moved throughput by a median of +1.7%,
//! never clearing the 3.9% noise floor in any of 12 (file, backend, baseline) combinations —
//! including long peptides, the case it was designed for. It was reverted rather than carried.
//! `retrieve_proteins`' own two-pass prefetch already covers queries with many matches, and
//! queries with few matches are dominated by peptides where retrieval is a small share of the
//! work.

use std::cmp::min;

use protein_metadata::{ProteinsBackend, SEPARATION_CHARACTER};
use rayon::prelude::*;
use text_compression::ProteinTextBackend;

use super::{
    BoundSearch,
    BoundSearch::{Maximum, Minimum},
    BoundSearchResult, MAX_RESULT_PREALLOC, MLP_BATCH, SearchAllSuffixesResult, Searcher,
    measure::Timer,
    tryptic::{TRYPTIC_EXTENSION_CHARS, tryptic_extension_chars}
};
use crate::{array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// The `_batched` twin of `binary_search_bound_scalar`, advancing many independent streams in
    /// lockstep. Three stages per level create memory-level parallelism:
    ///   A: compute `center` + prefetch `SA[center]` for all active streams
    ///   B: read `SA[center]` + prefetch `text[suffix]` for all active streams
    ///   C: `compare` + update `lo/hi/lcp` for all active streams
    fn binary_search_bound_batched(&self, bound: BoundSearch, streams: &mut [BsStream]) {
        loop {
            // Stage A: center + SA prefetch
            let mut any = false;
            for s in streams.iter_mut() {
                if s.hi - s.lo > 1 {
                    s.active = true;
                    any = true;
                    s.center = (s.lo + s.hi) / 2;
                    self.sa.prefetch_sa_index(s.center);
                } else {
                    s.active = false;
                }
            }
            if !any {
                break;
            }

            // Stage B: read SA + prefetch the (random) text position
            let text = self.proteins.text();
            let text_len = text.len();
            for s in streams.iter_mut() {
                if s.active {
                    s.suffix = self.sa.get(s.center);
                    let skip = min(s.lcp_left, s.lcp_right);
                    let pos = s.suffix as usize + skip;
                    if pos < text_len {
                        text.prefetch_at(pos);
                    }
                }
            }

            // Stage C: compare + update (prefetches from A/B have now landed)
            for s in streams.iter_mut() {
                if s.active {
                    let skip = min(s.lcp_left, s.lcp_right);
                    let (retval, lcp_center) = self.compare(s.ss, s.suffix, skip, bound);
                    s.found |= lcp_center == s.ss.len();
                    if (retval && bound == Minimum) || (!retval && bound == Maximum) {
                        s.hi = s.center;
                        s.lcp_right = lcp_center;
                    } else {
                        s.lo = s.center;
                        s.lcp_left = lcp_center;
                    }
                }
            }
        }

        // Edge-case tail (mirrors the scalar version): when the window narrowed to
        // [left, left+1) with `lo` still at the original left, `left` itself was never
        // a center and must be checked.
        for s in streams.iter_mut() {
            if s.hi == s.left0 + 1 && s.lo == s.left0 {
                let skip = min(s.lcp_left, s.lcp_right);
                let (retval, lcp_center) = self.compare(s.ss, self.sa.get(s.lo), skip, bound);
                s.found |= lcp_center == s.ss.len();
                if bound == Minimum && retval {
                    s.hi = s.lo;
                }
            }
            s.res_found = s.found;
            s.res_bound = match bound {
                Minimum => s.hi,
                Maximum => s.lo
            };
        }
    }

    /// The `_batched` twin of `search_bounds_scalar`: many search strings at once.
    fn search_bounds_batched(&self, strings: &[&[u8]]) -> Vec<BoundSearchResult> {
        self.search_bounds_inner_batched(strings, true)
    }

    /// `search_bounds_batched` with the k-mer table forced off.
    ///
    /// Needed for the separator variant of the left-extended tryptic search: `'-'` is absent from
    /// the k-mer table's ALPHABET, so a lookup returns `None` and the stream would be dropped as
    /// `NoMatches`, silently losing every protein-start match. Mirrors
    /// `search_bounds_full_range_scalar`.
    fn search_bounds_full_range_batched(&self, strings: &[&[u8]]) -> Vec<BoundSearchResult> {
        self.search_bounds_inner_batched(strings, false)
    }

    fn search_bounds_inner_batched(&self, strings: &[&[u8]], use_kmer_table: bool) -> Vec<BoundSearchResult> {
        let n = strings.len();
        let mut out: Vec<BoundSearchResult> = (0..n).map(|_| BoundSearchResult::NoMatches).collect();

        // Stage 1: initial range per stream (k-mer narrowing or full range).
        let mut ranges: Vec<Option<(usize, usize, usize)>> = vec![None; n];
        let mut min_streams: Vec<BsStream> = Vec::with_capacity(n);
        for (i, &ss) in strings.iter().enumerate() {
            let range = match &self.kmer_table {
                Some(table) if use_kmer_table && ss.len() >= table.k => {
                    match table.lookup(&ss[..table.k]) {
                        Some((lo, hi)) => (lo, hi + 1, table.k),
                        None => continue // out[i] stays NoMatches
                    }
                }
                _ => (0, self.sa.len(), 0)
            };
            ranges[i] = Some(range);
            min_streams.push(BsStream::new(i, ss, range.0, range.1, range.2));
        }

        // Stage 2: batched minimum bound.
        self.binary_search_bound_batched(Minimum, &mut min_streams);

        let mut min_bounds: Vec<usize> = vec![0; n];
        let mut max_streams: Vec<BsStream> = Vec::with_capacity(min_streams.len());
        for s in &min_streams {
            if s.res_found {
                min_bounds[s.idx] = s.res_bound;
                let (left, right, lcp_skip) = ranges[s.idx].unwrap();
                max_streams.push(BsStream::new(s.idx, s.ss, left, right, lcp_skip));
            }
            // else: out[s.idx] stays NoMatches (mirrors `if !found_min return NoMatches`)
        }

        // Stage 3: batched maximum bound.
        self.binary_search_bound_batched(Maximum, &mut max_streams);
        for s in &max_streams {
            out[s.idx] = BoundSearchResult::SearchResult((min_bounds[s.idx], s.res_bound + 1));
        }

        out
    }

    /// The `_batched` twin of `search_matching_suffixes_scalar`, over one chunk of peptides.
    /// Returns one result per input, in order — identical to calling the `_scalar` form on each
    /// peptide, but with the
    /// binary-search phase interleaved for memory-level parallelism. Single-threaded and
    /// rayon-free: the MLP kernel. Cross-thread parallelism and chunking live in
    /// [`Searcher::search_all_matching_suffixes`] below, which is the entry point callers use.
    pub(crate) fn search_matching_suffixes_batched(
        &self,
        strings: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> Vec<SearchAllSuffixesResult> {
        let n = strings.len();
        let sample = self.sa.sample_rate() as usize;

        // Per peptide, same cap as the scalar path — see `MAX_RESULT_PREALLOC`. It matters more
        // here: this allocates one result vector per peptide in the batch up front.
        let mut matching: Vec<Vec<i64>> =
            strings.iter().map(|_| Vec::with_capacity(max_matches.min(MAX_RESULT_PREALLOC))).collect();
        let il_locs: Vec<Vec<usize>> = strings
            .iter()
            .map(|ss| {
                ss.iter()
                    .enumerate()
                    .filter_map(|(i, &c)| if c == b'I' || c == b'L' { Some(i) } else { None })
                    .collect()
            })
            .collect();
        let mut done: Vec<Option<SearchAllSuffixesResult>> = (0..n).map(|_| None).collect();

        // Mirrors the scalar path: for tryptic searches the skip = sample-1 pass is replaced by
        // left-extended searches, which cover the same positions with a ~20x smaller SA range.
        // See `search_matching_suffixes_scalar` in scalar.rs for the derivation.
        let use_extended = tryptic && sample >= 2;
        let skip_end = if use_extended { sample - 1 } else { sample };

        for skip in 0..skip_end {
            let active: Vec<usize> = (0..n).filter(|&i| done[i].is_none()).collect();
            if active.is_empty() {
                break;
            }

            let sub: Vec<&[u8]> = active.iter().map(|&i| &strings[i][skip..]).collect();
            let t_bounds = Timer::start();
            let bounds = self.search_bounds_batched(&sub);
            self.measurements.search_bounds_ns.add(t_bounds.elapsed_ns());

            let t_iter = Timer::start();
            for (ai, &i) in active.iter().enumerate() {
                let (min_bound, max_bound) = match &bounds[ai] {
                    BoundSearchResult::SearchResult((lo, hi)) => (*lo, *hi),
                    BoundSearchResult::NoMatches => continue
                };
                let search_string = strings[i];
                // See scalar.rs's search_matching_suffixes_scalar for the soundness argument: `compare`
                // always normalizes L->I on both sides (the index is built with L replaced by
                // I), so an I/L-free peptide cannot mismatch an in-range suffix regardless of
                // equate_il — the fast path applies whenever equate_il is true, or the peptide
                // has no I/L to begin with.
                if (equate_il || il_locs[i].is_empty()) && !tryptic && skip == 0 {
                    let range_size = max_bound - min_bound;
                    // Strictly larger: a range of exactly max_matches is a complete set. Kept in
                    // step with `scalar.rs`, which makes the same distinction the same way.
                    if range_size > max_matches {
                        let result: Vec<i64> = self.sa.iter_range(min_bound, min_bound + max_matches).collect();
                        done[i] = Some(SearchAllSuffixesResult::MaxMatches(result));
                    } else {
                        for s in self.sa.iter_range(min_bound, max_bound) {
                            matching[i].push(s);
                        }
                    }
                } else {
                    let il_from_skip = &il_locs[i][il_locs[i].partition_point(|&x| x < skip)..];
                    let prefix = &search_string[..skip];
                    let suffix_str = &search_string[skip..];
                    let text = self.proteins.text();
                    let hit_max = self.iterate_sa_range(
                        self.sa.iter_range(min_bound, max_bound),
                        max_bound.saturating_sub(min_bound),
                        text,
                        skip,
                        search_string,
                        prefix,
                        suffix_str,
                        il_from_skip,
                        equate_il,
                        tryptic,
                        &mut matching[i],
                        max_matches
                    );
                    if hit_max {
                        done[i] =
                            Some(SearchAllSuffixesResult::truncated(std::mem::take(&mut matching[i]), max_matches));
                    }
                }
            }
            self.measurements.match_iter_ns.add(t_iter.elapsed_ns());
        }

        // Left-extended phase. One batched bound search per extension character, so the streams
        // stay interleaved exactly as in the skip loop above.
        if use_extended {
            let text = self.proteins.text();
            // Flat buffer holding every `X + peptide` for this variant, with per-stream offsets —
            // one allocation per variant instead of one per peptide.
            let mut ext_buf: Vec<u8> = Vec::new();
            let mut ext_spans: Vec<(usize, usize)> = Vec::new();

            for &prefix_char in TRYPTIC_EXTENSION_CHARS.iter() {
                let active: Vec<usize> = (0..n)
                    .filter(|&i| done[i].is_none() && tryptic_extension_chars(strings[i]).contains(&prefix_char))
                    .collect();
                if active.is_empty() {
                    continue;
                }

                ext_buf.clear();
                ext_spans.clear();
                for &i in &active {
                    let start = ext_buf.len();
                    ext_buf.push(prefix_char);
                    ext_buf.extend_from_slice(strings[i]);
                    ext_spans.push((start, ext_buf.len()));
                }
                let extended: Vec<&[u8]> = ext_spans.iter().map(|&(s, e)| &ext_buf[s..e]).collect();

                let t_bounds = Timer::start();
                // The separator is not representable in the k-mer table — routing it there would
                // silently drop every protein-start match. See `search_bounds_full_range_batched`.
                let bounds = if prefix_char == SEPARATION_CHARACTER {
                    self.search_bounds_full_range_batched(&extended)
                } else {
                    self.search_bounds_batched(&extended)
                };
                self.measurements.search_bounds_ns.add(t_bounds.elapsed_ns());

                let t_iter = Timer::start();
                for (ai, &i) in active.iter().enumerate() {
                    let (min_bound, max_bound) = match &bounds[ai] {
                        BoundSearchResult::SearchResult((lo, hi)) => (*lo, *hi),
                        BoundSearchResult::NoMatches => continue
                    };
                    let search_string = strings[i];
                    let last_is_kr = matches!(search_string.last(), Some(b'K' | b'R'));
                    let hit_max = self.iterate_extended_sa_range(
                        self.sa.iter_range(min_bound, max_bound),
                        max_bound.saturating_sub(min_bound),
                        text,
                        search_string,
                        &il_locs[i],
                        equate_il,
                        last_is_kr,
                        &mut matching[i],
                        max_matches
                    );
                    if hit_max {
                        done[i] =
                            Some(SearchAllSuffixesResult::truncated(std::mem::take(&mut matching[i]), max_matches));
                    }
                }
                self.measurements.match_iter_ns.add(t_iter.elapsed_ns());
            }
        }

        (0..n)
            .map(|i| {
                if let Some(r) = done[i].take() {
                    r
                } else if matching[i].is_empty() {
                    SearchAllSuffixesResult::NoMatches
                } else {
                    SearchAllSuffixesResult::SearchResult(std::mem::take(&mut matching[i]))
                }
            })
            .collect()
    }

    /// Searches every peptide, returning one raw suffix result per input, in order.
    ///
    /// Peptides are split into `MLP_BATCH`-sized chunks, one rayon task each, and each chunk goes
    /// through `search_matching_suffixes_batched` so its binary searches interleave. The output is
    /// identical to running `search_matching_suffixes_scalar` over the list one peptide at a time — that
    /// is what `test_orchestrated_matches_scalar` asserts, over three suffix-array shapes and
    /// both search options, retrieval included.
    ///
    /// Peptides are searched as given — callers that need length filtering or normalisation
    /// (e.g. skipping peptides shorter than the sample rate) must do it before calling.
    pub fn search_all_matching_suffixes(
        &self,
        peptides: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> Vec<SearchAllSuffixesResult> {
        peptides
            .par_chunks(MLP_BATCH)
            .flat_map(|chunk| self.search_matching_suffixes_batched(chunk, max_matches, equate_il, tryptic))
            .collect()
    }
}

/// Per-stream state for the batched binary search (`binary_search_bound_batched`).
struct BsStream<'a> {
    idx: usize,
    ss: &'a [u8],
    left0: usize,
    lo: usize,
    hi: usize,
    lcp_left: usize,
    lcp_right: usize,
    found: bool,
    center: usize,
    suffix: i64,
    active: bool,
    res_found: bool,
    res_bound: usize
}

impl<'a> BsStream<'a> {
    #[inline]
    fn new(idx: usize, ss: &'a [u8], left: usize, right: usize, lcp_skip: usize) -> Self {
        BsStream {
            idx,
            ss,
            left0: left,
            lo: left,
            hi: right,
            lcp_left: lcp_skip,
            lcp_right: lcp_skip,
            found: false,
            center: 0,
            suffix: 0,
            active: false,
            res_found: false,
            res_bound: 0
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sa_searcher::{
        BoundSearchResult, PREFETCH_THRESHOLD, SearchAllSuffixesResult, VALIDATE_BATCH,
        test_utils::{
            EXAMPLE_SA_SPARSE3, Mapping, PreloadedSearcher, TRYPTIC_FIXTURE, example_searcher, example_searcher_with,
            repeated_residue_searcher, searcher_over_text, tryptic_fixture_peptides
        }
    };

    #[test]
    fn test_batched_matches_scalar() {
        // Assert the batched search returns per-peptide results identical to the scalar
        // search across equate_il and max_matches settings. (macro avoids naming the
        // `ProteinsBackend as _` trait bound in a generic helper.)
        macro_rules! check_batched {
            ($searcher:expr, $peptides:expr) => {{
                for &eq in &[false, true] {
                    for &tr in &[false, true] {
                        for &mm in &[usize::MAX, 1usize, 2usize] {
                            let scalar: Vec<_> = $peptides
                                .iter()
                                .map(|p| $searcher.search_matching_suffixes_scalar(p, mm, eq, tr))
                                .collect();
                            let batched = $searcher.search_matching_suffixes_batched($peptides, mm, eq, tr);
                            for i in 0..$peptides.len() {
                                assert_eq!(
                                    batched[i],
                                    scalar[i],
                                    "mismatch: peptide={:?} equate_il={} tryptic={} max_matches={}",
                                    std::str::from_utf8($peptides[i]).unwrap(),
                                    eq,
                                    tr,
                                    mm
                                );
                            }
                        }
                    }
                }
            }};
        }

        // Dense/original SA (sample_rate 1)
        let searcher = example_searcher();
        let peptides: Vec<&[u8]> =
            vec![b"A", b"AC", b"AI", b"CLA", b"KCRLY", b"VAA", b"CVAA", b"LACVAA", b"C", b"ZZ", b"$"];
        check_batched!(&searcher, &peptides);

        // Sparse SA (sample_rate 3) — exercises skip = 0, 1, 2
        let searcher = example_searcher_with(&EXAMPLE_SA_SPARSE3, 3, Mapping::Sparse);
        let peptides: Vec<&[u8]> = vec![b"CLA", b"ACVAA", b"KCRLY", b"VAA", b"LACVAA", b"CVAA", b"CLACVAA", b"ZZZ"];
        check_batched!(&searcher, &peptides);
    }

    // The batched left-extended tryptic path must agree with the dense scalar index, which does
    // not use the transform at all (sparseness 1 keeps the original skip loop) and is therefore
    // an exact oracle.
    //
    // `test_batched_matches_scalar` above passes tryptic=false throughout, so without this the
    // batched extended path would have no coverage at all.
    #[test]
    fn test_batched_extended_tryptic_matches_dense_scalar() {
        let dense = searcher_over_text(TRYPTIC_FIXTURE, 1);
        let owned = tryptic_fixture_peptides();
        let peptides: Vec<&[u8]> = owned.iter().map(|p| p.as_slice()).collect();

        // Guard against a vacuous pass.
        let hits = peptides
            .iter()
            .filter(|p| {
                !matches!(
                    dense.search_matching_suffixes_scalar(p, usize::MAX, true, true),
                    SearchAllSuffixesResult::NoMatches
                )
            })
            .count();
        assert!(hits >= 5, "fixture yields only {hits} tryptic hits");

        for &sparseness in &[2u8, 3] {
            let sparse = searcher_over_text(TRYPTIC_FIXTURE, sparseness);
            for equate_il in [false, true] {
                let batched = sparse.search_matching_suffixes_batched(&peptides, usize::MAX, equate_il, true);
                for (i, p) in peptides.iter().enumerate() {
                    assert_eq!(
                        batched[i],
                        dense.search_matching_suffixes_scalar(p, usize::MAX, equate_il, true),
                        "sparseness={sparseness} equate_il={equate_il} peptide={:?}",
                        std::str::from_utf8(p).unwrap()
                    );
                }
            }
        }
    }

    // The separator variant must bypass the k-mer table in the batched path too, exactly as in
    // the scalar one — otherwise every protein-start match disappears once a table is attached.
    // "PKTR" starts protein 2 at position 21, which is unsampled at sparseness 2, so it is
    // reachable only through the '-' extended search.
    #[test]
    fn test_batched_extended_protein_start_with_kmer_table() {
        let plain = searcher_over_text(TRYPTIC_FIXTURE, 2);
        let mut kmered = searcher_over_text(TRYPTIC_FIXTURE, 2);
        kmered.build_kmer_table(3);

        let peptides: Vec<&[u8]> = vec![b"PKTR", b"RIY", b"KTR", b"MKAPTR", b"QST"];

        let plain_res = plain.search_matching_suffixes_batched(&peptides, usize::MAX, true, true);
        assert_eq!(
            plain_res[0],
            SearchAllSuffixesResult::SearchResult(vec![21]),
            "protein-start tryptic match missing from the batched un-tabled searcher"
        );

        let kmer_res = kmered.search_matching_suffixes_batched(&peptides, usize::MAX, true, true);
        for (i, p) in peptides.iter().enumerate() {
            assert_eq!(
                kmer_res[i],
                plain_res[i],
                "k-mer table changed the batched tryptic result for {:?}",
                std::str::from_utf8(p).unwrap()
            );
        }
    }

    #[test]
    fn test_batched_empty() {
        let searcher = example_searcher();
        assert!(searcher.search_matching_suffixes_batched(&[], usize::MAX, false, false).is_empty());
    }

    // The batched search must give the same result with a k-mer table as the plain scalar
    // search (covers search_bounds_batched's k-mer branch). L/I-free prefixes only, since the
    // raw test SA is not L->I normalized (see the scalar k-mer test for the reason).
    #[test]
    fn test_batched_with_kmer_table() {
        let reference = example_searcher();
        let mut kmered = example_searcher();
        kmered.build_kmer_table(3);

        let peptides: Vec<&[u8]> = vec![b"VAA", b"CVAA", b"KCR", b"KCRLY", b"AC", b"ZZZ"];
        let batched = kmered.search_matching_suffixes_batched(&peptides, usize::MAX, false, false);
        for (i, p) in peptides.iter().enumerate() {
            assert_eq!(
                batched[i],
                reference.search_matching_suffixes_scalar(p, usize::MAX, false, false),
                "batched+kmer vs plain scalar mismatch for {:?}",
                std::str::from_utf8(p).unwrap()
            );
        }
    }

    // Regression guard for the I/L-free fast-path extension on the batched path (mirrors
    // scalar.rs's `test_il_free_fast_path_matches_equate_il_true`). With peptides that contain
    // neither I nor L, equate_il=false must behave identically to equate_il=true, and the
    // batched kernel must agree with the scalar one — at a scale that exercises both the
    // range_size >= max_matches and range_size < max_matches fast-path branches.
    #[test]
    fn test_batched_il_free_fast_path_at_scale() {
        let searcher = repeated_residue_searcher('A', 70);
        let peptides: Vec<&[u8]> = vec![b"A", b"A", b"A"];
        for &mm in &[usize::MAX, 10usize] {
            let scalar: Vec<_> =
                peptides.iter().map(|p| searcher.search_matching_suffixes_scalar(p, mm, false, false)).collect();
            let batched_false = searcher.search_matching_suffixes_batched(&peptides, mm, false, false);
            let batched_true = searcher.search_matching_suffixes_batched(&peptides, mm, true, false);
            assert_eq!(batched_false, scalar, "batched(equate_il=false) vs scalar mismatch at max_matches={mm}");
            assert_eq!(batched_false, batched_true, "batched equate_il=false vs true mismatch at max_matches={mm}");
        }
    }

    /// What the entry point has to reproduce: the scalar primitive, applied one peptide at a time.
    ///
    /// This is the reference every test below compares against. It used to be
    /// `search_all_matching_suffixes` at `mlp_batch = 1` — the same function under a different
    /// setting, which could only ever catch a batching bug that the batch size happened to expose.
    /// Comparing against the primitive the batched path is built from is the stronger check, and it
    /// is the only one still possible now that the batch size is a constant.
    fn scalar_reference(
        searcher: &PreloadedSearcher,
        peptides: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> Vec<SearchAllSuffixesResult> {
        peptides
            .iter()
            .map(|&p| searcher.search_matching_suffixes_scalar(p, max_matches, equate_il, tryptic))
            .collect()
    }

    /// The whole-list entry point must equal the scalar primitive, across the
    /// equate_il × tryptic combinations.
    #[test]
    fn test_search_all_matches_scalar() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"KCR", b"VAA", b"ZZZ", b"AI"];
        for equate_il in [true, false] {
            for tryptic in [true, false] {
                let scalar = scalar_reference(&s, &peptides, 1000, equate_il, tryptic);
                let batched = s.search_all_matching_suffixes(&peptides, 1000, equate_il, tryptic);
                assert_eq!(scalar.len(), peptides.len());
                assert_eq!(batched.len(), peptides.len(), "one result per input, in order");
                for (i, (a, b)) in scalar.iter().zip(batched.iter()).enumerate() {
                    assert_eq!(a, b, "peptide {} differs (il={} tr={})", i, equate_il, tryptic);
                }
            }
        }
    }

    /// max_matches cutoff must propagate through the entry point exactly as through the primitive.
    #[test]
    fn test_search_all_respects_cutoff() {
        let s = example_searcher();
        // "A" matches 5 suffixes; cap at 2 → MaxMatches with 2 entries on both paths.
        let peptides: Vec<&[u8]> = vec![b"A"];
        let scalar = scalar_reference(&s, &peptides, 2, true, false);
        let batched = s.search_all_matching_suffixes(&peptides, 2, true, false);
        assert!(matches!(scalar[0], SearchAllSuffixesResult::MaxMatches(_)));
        assert_eq!(scalar[0], batched[0]);
    }

    #[test]
    fn test_search_all_empty() {
        let s = example_searcher();
        let none: Vec<&[u8]> = vec![];
        assert!(s.search_all_matching_suffixes(&none, 1000, true, false).is_empty());
    }

    // A peptide count that does not divide MLP_BATCH must still return one result per peptide in
    // the original order (exercises the par_chunks remainder). Five peptides against a batch of 16
    // is one short chunk; the fixture below is what covers a full chunk plus a remainder.
    #[test]
    fn test_search_all_ragged_batch() {
        let s = example_searcher();
        let peptides: Vec<&[u8]> = vec![b"A", b"AC", b"CLA", b"AI", b"VAA"];
        assert_eq!(
            scalar_reference(&s, &peptides, 1000, true, false),
            s.search_all_matching_suffixes(&peptides, 1000, true, false)
        );
    }

    /// 500 short proteins built from a rotating set of motifs. The repetition is what makes the
    /// fixture useful: a 2-mer query lands an SA range well above the two-pass threshold
    /// (`PREFETCH_THRESHOLD`) while the longer, rarer peptides stay below it, so one run covers
    /// both branches of `iterate_sa_range`. K/R/P give the tryptic filter a real mix of accepts and
    /// rejects, and the I/L pairs exercise both `equate_il` settings.
    fn pipeline_fixture_text() -> String {
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

    /// The one non-repeated protein in `pipeline_fixture_text`. Ends in K so it also has a valid
    /// tryptic C-terminus, and carries an L so `equate_il` is not moot for it.
    const RARE_PROTEIN: &str = "MWYFQDNTLHK";

    /// Peptides spanning the interesting shapes: very common 1–3-mers (large SA ranges, so the
    /// two-pass path and the `max_matches` cutoff), unique ones from `RARE_PROTEIN` (ranges of
    /// 1–2 entries, so the straight scan), I/L-carrying and I/L-free ones, K/R- and
    /// P-terminated ones, and a total miss.
    ///
    /// Eighteen of them, against an `MLP_BATCH` of 16: one full chunk and a remainder, so the
    /// ragged tail is covered here too and at the size that actually ships.
    const PIPELINE_PEPTIDES: &[&[u8]] = &[
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

    /// One row per peptide: the result-variant tag, the matched suffixes in the order they were
    /// produced, and the retrieved `(taxon_id, uniprot_id)` pairs (`ProteinRef` borrows, so it is
    /// not `PartialEq` itself).
    type PipelineRow = (&'static str, Vec<i64>, Vec<(u32, String)>);

    fn rows_from(searcher: &PreloadedSearcher, results: &[SearchAllSuffixesResult]) -> Vec<PipelineRow> {
        results
            .iter()
            .map(|result| {
                let tag = match result {
                    SearchAllSuffixesResult::NoMatches => "none",
                    SearchAllSuffixesResult::MaxMatches(_) => "max",
                    SearchAllSuffixesResult::SearchResult(_) => "hit"
                };
                let suffixes = match result {
                    SearchAllSuffixesResult::NoMatches => vec![],
                    SearchAllSuffixesResult::MaxMatches(v) | SearchAllSuffixesResult::SearchResult(v) => v.clone()
                };
                let proteins = searcher
                    .retrieve_proteins(&suffixes)
                    .iter()
                    .map(|p| (p.taxon_id, p.uniprot_id.to_string()))
                    .collect();
                (tag, suffixes, proteins)
            })
            .collect()
    }

    #[test]
    fn test_orchestrated_matches_scalar() {
        let text = pipeline_fixture_text();

        // Dense SA (skip is always 0), sparse SA (exercises skip = 1, 2 and therefore the
        // prefix/suffix validation), and a k-mer-narrowed one.
        let mut kmered = searcher_over_text(&text, 1);
        kmered.build_kmer_table(3);
        // The third element is the sparseness factor: `search_matching_suffixes_scalar` indexes
        // `search_string[skip..]` for every skip below it, so peptides shorter than that are
        // out of contract (production drops them before the search) and must be filtered out.
        let fixtures = [
            ("dense", searcher_over_text(&text, 1), 1usize),
            ("sparse", searcher_over_text(&text, 3), 3usize),
            ("kmer", kmered, 1usize)
        ];

        // The fixture must actually reach both scan paths, otherwise the comparison proves
        // nothing about `iterate_sa_range`. "AA" must exceed the validation batch (so the batch
        // refills at least once, which also puts it past the threshold); "MWYFQ" must fall below
        // the threshold, so the straight scan runs.
        match fixtures[0].1.search_bounds_scalar(b"AA") {
            BoundSearchResult::SearchResult((lo, hi)) => {
                assert!(hi - lo > VALIDATE_BATCH, "fixture too small: 'AA' range is only {}", hi - lo)
            }
            other => panic!("expected 'AA' to match, got {other:?}")
        }
        match fixtures[0].1.search_bounds_scalar(b"MWYFQ") {
            BoundSearchResult::SearchResult((lo, hi)) => {
                assert!(hi - lo < PREFETCH_THRESHOLD, "'MWYFQ' range {} is not below the threshold", hi - lo)
            }
            other => panic!("expected 'MWYFQ' to match, got {other:?}")
        }

        for (name, searcher, min_len) in fixtures.iter() {
            let peptides: Vec<&[u8]> = PIPELINE_PEPTIDES.iter().copied().filter(|p| p.len() >= *min_len).collect();

            for equate_il in [false, true] {
                for tryptic in [false, true] {
                    // Finite and well below the number of matches the common peptides have, so
                    // the cutoff fires mid-scan rather than never.
                    let expected = rows_from(searcher, &scalar_reference(searcher, &peptides, 64, equate_il, tryptic));

                    // Every peptide matching nothing would make the comparison vacuous.
                    assert!(
                        expected.iter().any(|(tag, _, _)| *tag != "none"),
                        "{name}: fixture matches nothing (il={equate_il} tr={tryptic})"
                    );

                    let actual =
                        rows_from(searcher, &searcher.search_all_matching_suffixes(&peptides, 64, equate_il, tryptic));
                    assert_eq!(actual, expected, "{name}: orchestrated results differ (il={equate_il} tr={tryptic})");
                }
            }
        }
    }
}
