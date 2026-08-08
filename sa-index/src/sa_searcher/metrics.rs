//! Instrumentation shims that compile away when the `metrics` feature is off.
//!
//! The counters live on `Searcher`, which every rayon worker shares, so each `fetch_add` is a
//! contended RMW on one cache line and each `Instant::now()` is a serializing clock read. In
//! `scalar.rs` that is 4 clock reads + 4 atomic RMWs per peptide at sample_rate 2 (~2% of
//! runtime at `mlp_batch=1`, well under 1% at the default batch of 16) — small, but it is
//! always-on production code perturbing the very measurement it exists to produce.
//!
//! Rather than sprinkling `#[cfg]` through the hot loops, both primitives are ZSTs with
//! no-op inherent methods when the feature is off: `Timer::start()` becomes nothing,
//! `elapsed_ns()` folds to the constant 0, and `Counter::add(0)` folds away with it. The call
//! sites read identically in both configurations.

#[cfg(feature = "metrics")]
mod imp {
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::Instant;

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
        /// Its only callers are the metrics tests, which live behind
        /// `#[cfg(all(test, not(feature = "mmap")))]`, so it reads as dead under `--all-features`.
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

#[cfg(not(feature = "metrics"))]
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
