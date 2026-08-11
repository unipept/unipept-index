//! The suffix array itself: sorted positions into the concatenated protein text.
//!
//! Four storage types, in two axes.
//!
//! **Packing.** [`OriginalSA`] stores one `i64` per entry; [`CompressedSA`] packs entries at the
//! minimum width the text length requires (29 bits for a 300 M-residue text), roughly halving the
//! file. Which one a file holds is recorded in its header, so the reader picks at runtime.
//!
//! **Sparseness.** The array may index only every n-th text position, trading search work for
//! size. The factor is likewise a header field, and search compensates by trying each of the n
//! possible alignments. A peptide shorter than the factor cannot be searched at all.
//!
//! **Storage.** [`InMemorySA`] holds owned memory and dispatches over the packing at runtime;
//! `MmapBackedSA` decodes straight out of a memory mapping and knows its packing from the same
//! header. The `mmap` feature picks which one [`SuffixArray`] means; see the crate docs.
//!
//! All four implement [`SuffixArrayBackend`], which is what the searcher is written against.

pub mod compressed;
#[cfg(feature = "mmap")]
pub mod mmap;
pub mod original;
pub mod preloaded;

pub use compressed::{CompressedSA, dump_compressed_suffix_array, load_compressed_suffix_array};
#[cfg(feature = "mmap")]
pub use mmap::MmapBackedSA;
pub use original::{OriginalRangeIter, OriginalSA, dump_suffix_array};
pub use preloaded::{InMemoryRangeIter, InMemorySA};

/// Type alias so existing call-sites can keep using `SuffixArray` unchanged.
#[cfg(feature = "mmap")]
pub type SuffixArray = MmapBackedSA;
/// Type alias so existing call-sites can keep using `SuffixArray` unchanged.
#[cfg(not(feature = "mmap"))]
pub type SuffixArray = InMemorySA;

/// Common interface implemented by every SA storage backend.
///
/// - [`OriginalSA`] and [`CompressedSA`] — owned memory, one packing each (always available).
/// - `MmapBackedSA` — memory-mapped, either packing (mmap feature only).
/// - [`InMemorySA`] — dispatches over `OriginalSA`/`CompressedSA` at runtime.
pub trait SuffixArrayBackend: Send + Sync {
    /// The concrete iterator type returned by [`Self::iter_range`].
    type RangeIter<'a>: Iterator<Item = i64> + ExactSizeIterator
    where
        Self: 'a;

    /// Number of entries in the array.
    ///
    /// With a sparseness factor above 1 this is smaller than the text length.
    fn len(&self) -> usize;

    /// Bits used to store one entry: 64 when unpacked, otherwise the compressed width.
    fn bits_per_value(&self) -> usize;

    /// The sparseness factor the array was built with — the array indexes every n-th text
    /// position. 1 means every position.
    ///
    /// Search reads this to decide how many alignments to try, and to reject peptides too short
    /// to be searchable.
    fn sample_rate(&self) -> u8;

    /// Returns the text position stored at `index`.
    ///
    /// # Panics
    ///
    /// If `index >= len()`.
    fn get(&self, index: usize) -> i64;

    /// Iterates entries in `start..end` (half-open).
    ///
    /// Preferred over a `get` loop when the range is walked in order: implementations amortise
    /// the unpacking across consecutive entries instead of re-deriving it per call.
    fn iter_range(&self, start: usize, end: usize) -> Self::RangeIter<'_>;

    /// Issues a hardware prefetch hint for the storage holding `index`, without reading it.
    ///
    /// Used to overlap the DRAM latency of the next binary-search probe with the current
    /// comparison. Out-of-range indices are ignored rather than panicking.
    fn prefetch_sa_index(&self, index: usize);

    /// Asks the kernel to start reading the pages holding entries `start..end`, without waiting.
    ///
    /// The counterpart to [`Self::prefetch_sa_index`] for pages that are not resident at all. A
    /// CPU prefetch hint cannot fault — see the `prefetch` crate docs — so it does nothing for an
    /// absent page, and the load that follows blocks the thread for a whole disk read.
    /// `madvise(MADV_WILLNEED)` queues real asynchronous readahead and returns, which is the only
    /// way to overlap misses from inside one thread.
    ///
    /// Measured worth attacking: at 75% residency the search phase splits roughly half into the
    /// dependent binary-search chain and half into the contiguous range scan this covers, and
    /// every mapping is opened `MADV_RANDOM`, so that scan takes one fault per page with no
    /// readahead at all.
    ///
    /// Default: no-op — only a mapping has anything to advise. Callers gate it on
    /// `SearchTuning::willneed` because each call is a syscall: worth it against a ~100 µs fault,
    /// pure overhead when the page was already resident.
    fn advise_willneed_range(&self, _start: usize, _end: usize) {}

    /// Reads every mapped page into the page cache.
    ///
    /// Default: no-op. Only the mmap backend implements it, as a warmup so the first requests do
    /// not pay the page faults.
    fn touch_all_pages(&self) {}

    /// Whether the array is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
