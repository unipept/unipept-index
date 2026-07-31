//! Cross-query batching for memory-level parallelism (MLP).
//!
//! Interleaves B independent peptide searches so ~B random memory misses are in flight
//! per core, hiding the DRAM latency the scalar dependent-chain binary search cannot.
//! Each stream's per-step logic is identical to the scalar path in the parent module, so
//! results are identical (see `test_batched_matches_scalar`). Lives in its own `impl`
//! block to keep the batched pipeline out of the scalar `Searcher` code.

use std::cmp::min;
use std::sync::atomic::Ordering;
use std::time::Instant;

use sa_mappings::proteins::ProteinsBackend;
use text_compression::ProteinTextBackend;

use crate::array::SuffixArrayBackend;
use crate::suffix_to_protein_index::SuffixToProteinMappingBackend;

use super::BoundSearch::{Maximum, Minimum};
use super::{BoundSearch, BoundSearchResult, SearchAllSuffixesResult, Searcher};

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Batched `binary_search_bound_in_range`, advancing many independent streams in
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
                Maximum => s.lo,
            };
        }
    }

    /// Batched `search_bounds` for many search strings at once.
    fn search_bounds_batched(&self, strings: &[&[u8]]) -> Vec<BoundSearchResult> {
        let n = strings.len();
        let mut out: Vec<BoundSearchResult> =
            (0..n).map(|_| BoundSearchResult::NoMatches).collect();

        // Stage 1: initial range per stream (k-mer narrowing or full range).
        let mut ranges: Vec<Option<(usize, usize, usize)>> = vec![None; n];
        let mut min_streams: Vec<BsStream> = Vec::with_capacity(n);
        for (i, &ss) in strings.iter().enumerate() {
            let range = match &self.kmer_table {
                Some(table) if ss.len() >= table.k => match table.lookup(&ss[..table.k]) {
                    Some((lo, hi)) => (lo, hi + 1, table.k),
                    None => continue, // out[i] stays NoMatches
                },
                _ => (0, self.sa.len(), 0),
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

    /// Batched `search_matching_suffixes` over many peptides. Returns one result per
    /// input, in order — identical to calling the scalar version on each peptide, but
    /// with the binary-search phase interleaved for memory-level parallelism.
    pub fn search_matching_suffixes_batched(
        &self,
        strings: &[&[u8]],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool,
    ) -> Vec<SearchAllSuffixesResult> {
        let n = strings.len();
        let sample = self.sa.sample_rate() as usize;

        let mut matching: Vec<Vec<i64>> =
            strings.iter().map(|_| Vec::with_capacity(max_matches.min(4096))).collect();
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

        for &ss in strings {
            self.prefetch_kmer_range(ss);
        }

        for skip in 0..sample {
            let active: Vec<usize> = (0..n).filter(|&i| done[i].is_none()).collect();
            if active.is_empty() {
                break;
            }

            let sub: Vec<&[u8]> = active.iter().map(|&i| &strings[i][skip..]).collect();
            let t_bounds = Instant::now();
            let bounds = self.search_bounds_batched(&sub);
            self.search_bounds_ns
                .fetch_add(t_bounds.elapsed().as_nanos() as u64, Ordering::Relaxed);

            let t_iter = Instant::now();
            for (ai, &i) in active.iter().enumerate() {
                let (min_bound, max_bound) = match &bounds[ai] {
                    BoundSearchResult::SearchResult((lo, hi)) => (*lo, *hi),
                    BoundSearchResult::NoMatches => continue,
                };
                let search_string = strings[i];
                if equate_il && !tryptic && skip == 0 {
                    let range_size = max_bound - min_bound;
                    if range_size >= max_matches {
                        let result: Vec<i64> =
                            self.sa.iter_range(min_bound, min_bound + max_matches).collect();
                        done[i] = Some(SearchAllSuffixesResult::MaxMatches(result));
                    } else {
                        for s in self.sa.iter_range(min_bound, max_bound) {
                            matching[i].push(s);
                        }
                    }
                } else {
                    let il_from_skip =
                        &il_locs[i][il_locs[i].partition_point(|&x| x < skip)..];
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
                        max_matches,
                    );
                    if hit_max {
                        done[i] = Some(SearchAllSuffixesResult::MaxMatches(std::mem::take(
                            &mut matching[i],
                        )));
                    }
                }
            }
            self.match_iter_ns
                .fetch_add(t_iter.elapsed().as_nanos() as u64, Ordering::Relaxed);

            if skip + 1 < sample {
                for &i in &active {
                    if done[i].is_none() {
                        self.prefetch_kmer_range(&strings[i][skip + 1..]);
                    }
                }
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
    res_bound: usize,
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
            res_bound: 0,
        }
    }
}
