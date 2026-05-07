pub mod original;
pub mod compressed;
#[cfg(feature = "mmap")]
pub mod mmap;
pub mod preloaded;

pub use original::{OriginalSA, OriginalRangeIter, dump_suffix_array};
pub use compressed::{CompressedSA, dump_compressed_suffix_array, load_compressed_suffix_array};
#[cfg(feature = "mmap")]
pub use mmap::MmapBackedSA;
pub use preloaded::{InMemorySA, InMemoryRangeIter};

/// Type alias so existing call-sites can keep using `SuffixArray` unchanged.
#[cfg(feature = "mmap")]
pub type SuffixArray = MmapBackedSA;
#[cfg(not(feature = "mmap"))]
pub type SuffixArray = InMemorySA;

/// Common interface implemented by every SA storage backend.
///
/// - [`OriginalSA`] and [`CompressedSA`] — in-memory backends (always available).
/// - [`MmapBackedSA`] — mmap backend (mmap feature only).
/// - [`InMemorySA`] — runtime-selected wrapper over Original/Compressed (non-mmap only).
pub trait SuffixArrayBackend: Send + Sync {
    /// The concrete iterator type returned by [`iter_range`].
    type RangeIter<'a>: Iterator<Item = i64> + ExactSizeIterator where Self: 'a;

    fn len(&self) -> usize;
    fn bits_per_value(&self) -> usize;
    fn sample_rate(&self) -> u8;
    fn get(&self, index: usize) -> i64;
    fn iter_range(&self, start: usize, end: usize) -> Self::RangeIter<'_>;
    fn prefetch_sa_index(&self, index: usize);

    /// Touch every mapped page. Default: no-op (non-mmap backends).
    fn touch_all_pages(&self) {}
    /// Issue `MADV_WILLNEED` for the SA pages covering `lo..hi_exclusive`. Default: no-op.
    fn prefetch_sa_range(&self, _lo: usize, _hi_exclusive: usize) {}

    fn is_empty(&self) -> bool { self.len() == 0 }
}
