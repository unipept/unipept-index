//! Opt-in instrumentation of the search path, behind the `measure` feature.
//!
//! # What it measures, and why the benchmark cannot
//!
//! `sa-benchmarks` brackets the search and retrieval phases with its own clock, so the wall time of
//! each is available in every build and owes nothing to this module. What it cannot reach is
//! anything *inside* one `search_all_matching_suffixes` call:
//!
//! * the split between the binary search and the range scan. Those alternate per peptide — three
//!   times per peptide at sparseness 3 — across every rayon worker, so there are tens of thousands
//!   of phase boundaries inside the one call an outside clock can bracket.
//! * how many SA entries the scan examined and how many it accepted. These are not timings at all;
//!   no profiler or clock produces them, and their ratio is what separates "the acceptance rate is
//!   low" from "whole ranges are being scanned to exhaustion" — see [`SearchMeasurements`].
//!
//! # Why it is a feature and not always on
//!
//! The counters live on `Searcher`, which every rayon worker shares, so each `fetch_add` is a
//! contended RMW on one cache line and each `Instant::now()` is a serializing clock read. Left
//! always on, that is instrumentation perturbing the very measurement it exists to produce: about
//! 2% of runtime for the scalar path and well under 1% at the batch size that ships. So the
//! benchmark carries exactly one instrumented arm, and its throughput is reported for scale only.
//!
//! Rather than sprinkling `#[cfg]` through the hot loops, both primitives are ZSTs with
//! no-op inherent methods when the feature is off: `Timer::start()` becomes nothing,
//! `elapsed_ns()` folds to the constant 0, and `Counter::add(0)` folds away with it. The call
//! sites read identically in both configurations, so a `Timer::start()` in `scalar.rs` is not
//! evidence that an uninstrumented build reads the clock — it does not.

#[cfg(feature = "measure")]
mod imp {
    use std::{
        sync::atomic::{AtomicU64, Ordering},
        time::Instant
    };

    /// A shared, monotonically increasing counter. Relaxed ordering throughout: the values are
    /// diagnostics, never used to synchronize anything.
    #[derive(Debug, Default)]
    pub struct Counter(AtomicU64);

    impl Counter {
        pub const fn new() -> Self {
            Self(AtomicU64::new(0))
        }

        /// Adds `v`. Callers in hot loops must accumulate locally and call this once per
        /// outer-loop invocation — a per-candidate RMW on a line shared by every worker would
        /// cost more than the code being measured.
        #[inline]
        pub fn add(&self, v: u64) {
            self.0.fetch_add(v, Ordering::Relaxed);
        }

        /// Overwrites the counter. Test/setup helper, not for the hot path.
        ///
        /// Only the `drain_*` tests below call it, so it is dead in a non-test build.
        #[allow(dead_code)]
        #[inline]
        pub fn store(&self, v: u64) {
            self.0.store(v, Ordering::Relaxed);
        }

        /// Returns the accumulated value and resets to zero.
        #[inline]
        pub fn drain(&self) -> u64 {
            self.0.swap(0, Ordering::Relaxed)
        }
    }

    /// A wall-clock stopwatch.
    #[derive(Clone, Copy)]
    pub struct Timer(Instant);

    impl Timer {
        #[inline]
        pub fn start() -> Self {
            Self(Instant::now())
        }

        #[inline]
        pub fn elapsed_ns(&self) -> u64 {
            self.0.elapsed().as_nanos() as u64
        }
    }
}

#[cfg(not(feature = "measure"))]
mod imp {
    /// Zero-sized stand-in: every method is a no-op, so the counter costs no memory, no
    /// atomics and no cache line.
    #[derive(Debug, Default)]
    pub struct Counter;

    impl Counter {
        pub const fn new() -> Self {
            Self
        }

        #[inline]
        pub fn add(&self, _v: u64) {}

        /// Mirrors the real `Counter::store`; see there for why this reads as dead code.
        #[allow(dead_code)]
        #[inline]
        pub fn store(&self, _v: u64) {}

        /// Always zero — the counter was never written.
        #[inline]
        pub fn drain(&self) -> u64 {
            0
        }
    }

    /// Zero-sized stand-in: `start()` emits no clock read and `elapsed_ns()` is the constant
    /// 0, which lets the optimizer delete the `Counter::add` that consumes it as well.
    #[derive(Clone, Copy)]
    pub struct Timer;

    impl Timer {
        #[inline]
        pub fn start() -> Self {
            Self
        }

        #[inline]
        pub fn elapsed_ns(&self) -> u64 {
            0
        }
    }
}

pub use imp::{Counter, Timer};
use sa_mappings::proteins::ProteinsBackend;

use super::Searcher;
use crate::{array::SuffixArrayBackend, suffix_to_protein_index::SuffixToProteinMappingBackend};

/// Every counter the search accumulates, as one field on [`Searcher`] rather than four.
///
/// Zero-sized without the `measure` feature, so the whole struct costs nothing and every `add`
/// below folds away with it.
#[derive(Debug, Default)]
pub struct SearchMeasurements {
    /// Total nanoseconds spent inside `search_bounds_scalar()` across all queries (since last drain).
    pub(super) search_bounds_ns: Counter,
    /// Total nanoseconds spent iterating matches in `search_matching_suffixes_scalar()` (since last drain).
    pub(super) match_iter_ns: Counter,
    /// Candidate suffixes inspected by `iterate_sa_range` (since last drain), i.e. every entry
    /// the SA-range scan looked at, accepted or not.
    pub(super) candidates_examined: Counter,
    /// Candidate suffixes `iterate_sa_range` accepted as real matches (since last drain).
    ///
    /// Together with `candidates_examined` this settles why tryptic search is ~12.5x slower
    /// than non-tryptic on 5–10 aa peptides: a low accepted/examined ratio means the scan is
    /// simply sifting ~1/ratio times more candidates to reach `max_matches` (make each check
    /// cheaper), whereas a ratio near 1 with the cutoff rarely reached means whole SA ranges
    /// are being scanned to exhaustion (a `max_candidates` scan cap is the fix).
    pub(super) candidates_accepted: Counter
}

impl SearchMeasurements {
    pub(super) const fn new() -> Self {
        Self {
            search_bounds_ns: Counter::new(),
            match_iter_ns: Counter::new(),
            candidates_examined: Counter::new(),
            candidates_accepted: Counter::new()
        }
    }
}

impl<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend> Searcher<SA, P, STPM> {
    /// Returns `(search_bounds_ns, match_iter_ns)` accumulated since the last call and resets both
    /// counters to zero.  Safe to call concurrently with ongoing searches (relaxed ordering).
    ///
    /// Present in both feature configurations so callers need no `cfg`; without the `measure`
    /// feature it always returns `(0, 0)`.
    pub fn drain_timing_ns(&self) -> (u64, u64) {
        (self.measurements.search_bounds_ns.drain(), self.measurements.match_iter_ns.drain())
    }

    /// Returns `(candidates_examined, candidates_accepted)` accumulated by `iterate_sa_range`
    /// since the last call and resets both counters to zero. Same contract as
    /// `drain_timing_ns`: always present, always `(0, 0)` without the `measure` feature.
    ///
    /// The ratio is the SA-range scan's acceptance rate; the tryptic paths are the interesting
    /// ones, since a non-tryptic I/L-free query never enters `iterate_sa_range` at all.
    pub fn drain_candidate_counts(&self) -> (u64, u64) {
        (self.measurements.candidates_examined.drain(), self.measurements.candidates_accepted.drain())
    }
}

#[cfg(test)]
mod tests {
    use crate::sa_searcher::test_utils::example_searcher;

    // drain_timing_ns returns the accumulated counters and resets them to zero.
    // Without the `measure` feature the counters are no-op ZSTs, so the drain is always (0, 0)
    // — the API stays present either way so callers need no `cfg`.
    #[test]
    fn test_drain_timing_ns() {
        let searcher = example_searcher();
        assert_eq!(searcher.drain_timing_ns(), (0, 0));

        searcher.measurements.search_bounds_ns.store(123);
        searcher.measurements.match_iter_ns.store(456);
        let expected = if cfg!(feature = "measure") { (123, 456) } else { (0, 0) };
        assert_eq!(searcher.drain_timing_ns(), expected);
        assert_eq!(searcher.drain_timing_ns(), (0, 0)); // reset after draining
    }

    // Same contract for the candidate counters.
    #[test]
    fn test_drain_candidate_counts() {
        let searcher = example_searcher();
        assert_eq!(searcher.drain_candidate_counts(), (0, 0));

        searcher.measurements.candidates_examined.store(70);
        searcher.measurements.candidates_accepted.store(7);
        let expected = if cfg!(feature = "measure") { (70, 7) } else { (0, 0) };
        assert_eq!(searcher.drain_candidate_counts(), expected);
        assert_eq!(searcher.drain_candidate_counts(), (0, 0)); // reset after draining
    }

    // With `measure` on, `iterate_sa_range` must actually count what it scans — and only what
    // *it* scans: the counters exist to measure the acceptance rate of the validating path, so
    // the fast path (which accepts a whole SA range without inspecting entries) must not
    // contribute. A text of all 'L' searched for "I" with equate_il=false enters the validating
    // path and rejects every candidate: 70 examined, 0 accepted.
    #[cfg(feature = "measure")]
    #[test]
    fn test_candidate_counts_are_accumulated() {
        use crate::sa_searcher::{SearchAllSuffixesResult, test_utils::searcher_over_text};

        let n = 70usize;
        let searcher = searcher_over_text(&format!("{}$", "L".repeat(n)), 1);

        // "I" matches all 70 'L' positions during the bound search (L is normalized to I), but
        // equate_il=false rejects every one of them: acceptance rate 0.
        searcher.drain_candidate_counts();
        assert_eq!(
            searcher.search_matching_suffixes_scalar(b"I", usize::MAX, false, false),
            SearchAllSuffixesResult::NoMatches
        );
        assert_eq!(searcher.drain_candidate_counts(), (n as u64, 0));

        // Same range with equate_il=true accepts everything it examines — but takes the fast
        // path, which bypasses iterate_sa_range entirely, so nothing is counted at all.
        searcher.search_matching_suffixes_scalar(b"I", usize::MAX, true, false);
        assert_eq!(searcher.drain_candidate_counts(), (0, 0));
    }
}
