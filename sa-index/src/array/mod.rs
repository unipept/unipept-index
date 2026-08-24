//! The suffix array itself: sorted positions into the concatenated protein text.
//!
//! Two independent choices, both recorded in the file header so a reader learns them at load time
//! rather than from the build it happens to be.
//!
//! **Packing.** [`OriginalSA`] stores one `i64` per entry; [`CompressedSA`] packs entries at the
//! minimum width the text length requires (29 bits for a 300 M-residue text), roughly halving the
//! file.
//!
//! **Storage.** [`InMemorySA`] holds owned memory and dispatches over the packing at runtime;
//! [`MmapBackedSA`] decodes straight out of a memory mapping and handles either packing itself.
//! Both are always compiled and the searcher is generic over them; see the crate docs.
//!
//! A third property, **sparseness**, is a build parameter rather than a storage choice: the array
//! may index only every n-th text position, trading search work for size. Every type above carries
//! the factor and reports it through [`SuffixArrayBackend::sample_rate`]; search compensates by
//! trying each of the n possible alignments, and a peptide shorter than the factor cannot be
//! searched at all.
//!
//! All four types implement [`SuffixArrayBackend`], which is what the searcher is written against.
//!
//! # On-disk format for `sa.bin`
//!
//! ```text
//! [ bits_per_value: u8 ]   64 for the uncompressed packing; the compressed width otherwise
//! [ sparseness_factor: u8 ]
//! [ item_count: u64 little-endian ]
//! [ data ]                 item_count entries at bits_per_value each
//! ```
//!
//! At 64 bits the data is plain little-endian `i64`s, written by [`dump_suffix_array`]. Below 64
//! it is `bitarray`'s packing — most-significant-bit first within each little-endian `u64` word,
//! entries may straddle words — written by [`dump_compressed_suffix_array`]. Both emit the header
//! through the same private helper, so the two packings cannot drift apart.
//!
//! Readers: [`InMemorySA::read_binary`](text_compression::ReadBinary::read_binary), which
//! dispatches on the first byte, and `MmapBackedSA::read_binary_mmap`.

use std::{error::Error, io::Write};

pub mod mmap;
pub mod preloaded;
#[cfg(test)]
mod test_utils;

pub use mmap::MmapBackedSA;
pub use preloaded::{
    CompressedSA, InMemoryRangeIter, InMemorySA, OriginalRangeIter, OriginalSA, dump_compressed_suffix_array,
    dump_suffix_array, load_compressed_suffix_array
};

/// Writes the 10-byte header both packings share; see the module docs for the layout.
pub(super) fn write_sa_header(
    bits_per_value: usize,
    sparseness_factor: u8,
    item_count: usize,
    writer: &mut impl Write
) -> Result<(), Box<dyn Error>> {
    writer
        .write_all(&[bits_per_value as u8])
        .map_err(|_| "Could not write the required bits to the writer")?;
    writer
        .write_all(&[sparseness_factor])
        .map_err(|_| "Could not write the sparseness factor to the writer")?;
    // As `u64`, not `usize`: the readers always take 8 bytes here.
    writer
        .write_all(&(item_count as u64).to_le_bytes())
        .map_err(|_| "Could not write the size of the suffix array to the writer")?;
    Ok(())
}

/// Rejects a width no packing can use, before either reader trusts it.
///
/// `bits_per_value` is the first byte of an untrusted header, and both readers go on to shift by
/// it and to derive byte offsets from it. Outside `1..=64` that is a shift overflow — a panic in
/// debug, a wrapped value in release — and zero is separately useless: it makes the declared body
/// zero-length, so the file passes every size check and every lookup then reads off the end.
/// Neither reader used to check, and both duly panicked at the first `get` on a crafted header
/// instead of erroring at load, which is what [`ReadBinaryMmap`](text_compression::ReadBinaryMmap)
/// forbids.
///
/// Both call this, so the two agree on what a loadable file is; a width one accepted and the other
/// did not would be worse than either.
pub(super) fn check_bits_per_value(bits_per_value: usize) -> Result<(), Box<dyn Error>> {
    if !(1..=64).contains(&bits_per_value) {
        return Err(format!(
            "The SA header declares {bits_per_value} bits per value, outside the supported range of 1..=64"
        )
        .into());
    }
    Ok(())
}

/// Common interface implemented by every SA storage backend.
///
/// - [`OriginalSA`] and [`CompressedSA`] — owned memory, one packing each.
/// - [`MmapBackedSA`] — memory-mapped, either packing.
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
    /// `index` must be below [`len`](Self::len). Beyond it the behaviour is unspecified and
    /// differs per backend: the owned ones panic, while the mmap one may read a compressed entry
    /// out of the file's trailing slack and return a value that was never written.
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

    /// Reads every mapped page into the page cache, returning the number of bytes swept so a
    /// caller can report a bandwidth.
    ///
    /// Default: 0, nothing to fault. Only the mmap backend implements it, as a warmup so the first
    /// requests do not pay the page faults.
    fn touch_all_pages(&self) -> u64 {
        0
    }

    /// Whether the array is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
