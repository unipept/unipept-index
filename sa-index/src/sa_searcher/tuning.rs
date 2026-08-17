//! The runtime performance knobs, and the two compile-time ceilings they work against.
//!
//! Nothing here changes an answer, only how long it takes to produce one — a claim the
//! equivalence test at the bottom of this file exists to keep true.

use serde::{Deserialize, Serialize};

use super::DEFAULT_MLP_BATCH;

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

/// Every runtime performance knob the searcher has, in one place.
///
/// Each field is a *pure* performance knob: results are identical for any setting (asserted by
/// `test_tuning_does_not_change_results`). Nothing here changes an answer, only how long it takes
/// to produce one.
///
/// Defaults are confirmed by the run3 full-DB sweep (3 peptide-length buckets x {preloaded, mmap} x
/// {fast-path, validating} baselines, 20 reps, 3.9% noise floor). Two knobs that the same sweep
/// found dead were removed rather than left to be re-measured: `retrieval_batch` (cross-query
/// batched retrieval, median +1.7%) and `scalar_kmer_prefetch` (+0.3%).
///
/// # Adding a knob
///
/// Add a field here, give it a default, and document what measurement justified that default. That
/// is the whole change: the struct is `Serialize`/`Deserialize`, so the benchmark harness records
/// every field into its output and accepts `--tune <field>=<value>` for any of them without being
/// modified, and the benchmark report picks the new knob up on its own — as a swept axis if a suite
/// declares `tune.<field>`, and otherwise in the line stating what was held and at what value.
///
/// `deny_unknown_fields` is what makes a typo in `--tune` an error instead of a silently ignored
/// setting, which would otherwise read as "this knob does nothing".
///
/// # What does not belong here
///
/// Only things the *searcher* reads. `RAYON_NUM_THREADS` is the notable exclusion: rayon's global
/// pool is built once per process, before any searcher exists, so a field here could not affect it.
/// It stays an environment variable, set per benchmark cell by the driver and recorded alongside
/// these values in the report.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchTuning {
    /// Peptides interleaved per rayon task, for cross-query memory-level parallelism.
    ///
    /// `> 1` runs `search_matching_suffixes_batched` over chunks of this size; `1` runs the scalar
    /// path, one peptide per task. See `DEFAULT_MLP_BATCH` for why the default is 16.
    ///
    /// Unlike the other fields this one is read at the top of the search rather than in a hot loop:
    /// it selects how the work is decomposed, not how memory is walked. It lives here anyway so
    /// that every knob is recorded, swept and reported through one mechanism.
    pub mlp_batch: usize,
    /// Candidates per two-pass validation batch in `iterate_sa_range`.
    /// Clamped to `1..=MAX_VALIDATE_BATCH`.
    ///
    /// The one knob measured to matter, and it is a cliff rather than a peak: 16 -> 32 gains
    /// ~10% in all 6 (bucket, backend) combos, then it plateaus. Against this default of 64,
    /// 32 wins in 1/6 and 128 wins in 5/6 but never above the noise floor, and 128 regresses
    /// long peptides on the preloaded backend by 2.7%. Do not lower it below 32.
    pub validate_batch: usize,
    /// Minimum SA range size, in entries, before `iterate_sa_range` and
    /// `iterate_extended_sa_range` use two-pass validation instead of a straight loop.
    ///
    /// This is the gate on the *candidate-validation* path, and the partner of
    /// `validate_batch`: that one sets the batch size, this one decides whether a range is big
    /// enough to batch at all. It does not reach retrieval, which prefetches unconditionally at
    /// `retrieval_prefetch_distance`. Below the threshold the two-pass overhead exceeds the
    /// latency it hides, so the scalar loop runs and `validate_batch` has no effect.
    ///
    /// Swept over {8, 16, 32, 64}: median full-range swing +0.9%, inside the noise floor
    /// everywhere. Note what that sweep could and could not see — every value in it leaves the
    /// two-pass path on for ranges above 64 and off for ranges below 8, so it priced the
    /// crossover, never the mechanism. Left tunable for re-measurement on other hardware, not
    /// because 32 is known to be special.
    ///
    /// The full-database 4x4 cross with `retrieval_prefetch_distance` (2dfa6517b7, three storage
    /// backends, 20 reps x 10,000 peptides) found **0 of 48 pairs** clearing their own noise floor
    /// against this default. Nothing in that run argues for moving it, and the argmaxes it produced
    /// disagree between backends — the sweep cannot tell 8 from 64 here. Both fields are kept for
    /// re-measurement rather than deleted; `docs/design/README.md` records the evidence separately
    /// from that decision.
    pub prefetch_threshold: usize,
    /// Prefetch look-ahead distance (in suffixes) inside protein retrieval.
    ///
    /// Swept over {8, 16, 32, 64}: median full-range swing +1.2%, inside the noise floor
    /// everywhere. Same caveat as `prefetch_threshold`, including the 0-of-48 full-database cross,
    /// and one of its own: a query matching fewer suffixes than the distance issues no prefetches
    /// at all, at any value in that sweep. Raising it widens that exclusion — tryptic queries match
    /// almost nothing and are mostly already in it.
    pub retrieval_prefetch_distance: usize
}

impl Default for SearchTuning {
    fn default() -> Self {
        Self {
            mlp_batch: DEFAULT_MLP_BATCH,    // the 8-16 knee; 32/64 regress short peptides
            validate_batch: 64,              // confirmed by run3; 16 costs ~10%, 128 gains nothing
            prefetch_threshold: 32,          // measured flat over 8..64
            retrieval_prefetch_distance: 32  // measured flat over 8..64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{MAX_VALIDATE_BATCH, SearchTuning};
    use crate::sa_searcher::{
        BoundSearchResult, SearchAllSuffixesResult,
        test_utils::{PreloadedSearcher, searcher_over_text}
    };

    // The two-pass path's stack buffer is sized by the compile-time MAX_VALIDATE_BATCH; a
    // larger runtime `validate_batch` must clamp to it, not index past the array.
    #[test]
    fn test_validate_batch_clamps_to_max() {
        let n = 300usize; // > MAX_VALIDATE_BATCH, so the batch loop refills at least once
        let mut searcher = searcher_over_text(&format!("{}$", "L".repeat(n)), 1);

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
        searcher: &mut PreloadedSearcher,
        peptides: &[&[u8]],
        equate_il: bool,
        tryptic: bool
    ) -> Vec<TuningRow> {
        let mut out = Vec::new();
        for mlp_batch in [1usize, 16] {
            // `mlp_batch` is part of the tuning now, so it is set the same way as every other knob.
            // The caller's outer loop reassigns the whole struct, so this override does not leak.
            searcher.tuning.mlp_batch = mlp_batch;
            let results = searcher.search_all_matching_suffixes(peptides, 64, equate_il, tryptic);
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
        // and prefetch_threshold proves nothing. "AA" must exceed the largest swept
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
                        for prefetch_threshold in [0usize, 8, 32] {
                            for retrieval_prefetch_distance in [1usize, 8, 32] {
                                // Functional update, so adding a knob to `SearchTuning` does not
                                // break this test — the new field simply starts at its default
                                // until someone sweeps it here too.
                                let tuning = SearchTuning {
                                    validate_batch,
                                    prefetch_threshold,
                                    retrieval_prefetch_distance,
                                    ..SearchTuning::default()
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
