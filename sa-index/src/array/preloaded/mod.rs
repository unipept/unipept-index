//! The two owned-memory suffix-array packings, and the runtime dispatch over them.
//!
//! These types own the `WriteBinary` implementations that produce the files *either* backend
//! reads, which is why `sa-builder` names only them. The header-sniffing `read_binary` below is
//! the single place that decides whether a file holds a 64-bit or a compressed array.
use std::{error::Error, io::BufRead};

use bitarray::DynBitArrayRangeIter;
use text_compression::{LoadIndex, ReadBinary, WriteBinary};

pub mod compressed;
pub mod original;
#[cfg(test)]
pub(super) mod test_utils;

pub use compressed::{CompressedSA, dump_compressed_suffix_array, load_compressed_suffix_array};
pub use original::{OriginalRangeIter, OriginalSA, dump_suffix_array};

use self::{compressed::load_compressed, original::load_original};
use super::SuffixArrayBackend;

/// Dispatch a method call to whichever backend variant is active.
macro_rules! dispatch {
    ($self:expr, $method:ident ( $($arg:expr),* )) => {
        match $self {
            Self::Original(b)   => b.$method($($arg),*),
            Self::Compressed(b) => b.$method($($arg),*),
        }
    };
}

// ── InMemorySA ───────────────────────────────────────────────────────────────

/// Wraps whichever packing a file holds, picked at load time from its `bits_per_value` header
/// field.
pub enum InMemorySA {
    Original(OriginalSA),
    Compressed(CompressedSA)
}

// ── InMemoryRangeIter ────────────────────────────────────────────────────────

/// [`InMemorySA::iter_range`](SuffixArrayBackend::iter_range)'s iterator: whichever of the two
/// packings' native iterators the loaded variant uses.
///
/// The enum is matched once per `next`, not once per range — the variant is fixed for the life of
/// the index, so the branch predicts perfectly, but it is a branch, and the mmap iterator has no
/// equivalent. Anything that widens this into real work per entry is a regression.
pub enum InMemoryRangeIter<'a> {
    Original(OriginalRangeIter<'a>),
    Compressed(DynBitArrayRangeIter<'a>)
}

impl Iterator for InMemoryRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> {
        dispatch!(self, next())
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = dispatch!(self, len());
        (n, Some(n))
    }
}

impl ExactSizeIterator for InMemoryRangeIter<'_> {}

// ── SuffixArrayBackend for InMemorySA ────────────────────────────────────────

impl SuffixArrayBackend for InMemorySA {
    type RangeIter<'a> = InMemoryRangeIter<'a>;

    fn len(&self) -> usize {
        dispatch!(self, len())
    }
    fn bits_per_value(&self) -> usize {
        dispatch!(self, bits_per_value())
    }
    fn sample_rate(&self) -> u8 {
        dispatch!(self, sample_rate())
    }
    #[inline]
    fn get(&self, index: usize) -> i64 {
        dispatch!(self, get(index))
    }
    // iter_range needs a manual match: each arm wraps its backend's native iterator
    // type into the appropriate InMemoryRangeIter variant.
    fn iter_range(&self, start: usize, end: usize) -> InMemoryRangeIter<'_> {
        match self {
            Self::Original(b) => InMemoryRangeIter::Original(b.iter_range(start, end)),
            Self::Compressed(b) => InMemoryRangeIter::Compressed(b.iter_range(start, end))
        }
    }

    #[inline]
    fn prefetch_sa_index(&self, index: usize) {
        dispatch!(self, prefetch_sa_index(index))
    }
}

// ── ReadBinary / WriteBinary ──────────────────────────────────────────────────

impl ReadBinary for InMemorySA {
    /// Reads the shared header and picks the variant its width names: 64 bits is the uncompressed
    /// packing, anything else the compressed one. This is the single place that decision is made
    /// for owned memory, and it is what `sa-server` calls in a preloaded build.
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf1 = [0u8; 1];
        reader.read_exact(&mut buf1).map_err(|_| "Could not read the required bits from the binary file")?;
        let bits_per_value = buf1[0] as usize;
        super::check_bits_per_value(bits_per_value)?;

        reader.read_exact(&mut buf1).map_err(|_| "Could not read the sample rate from the binary file")?;
        let sample_rate = buf1[0];
        super::check_sample_rate(sample_rate)?;

        let mut buf8 = [0u8; 8];
        reader
            .read_exact(&mut buf8)
            .map_err(|_| "Could not read the size of the suffix array from the binary file")?;
        let size = u64::from_le_bytes(buf8) as usize;

        if bits_per_value == 64 {
            let sa = load_original(reader, size)?;
            Ok(InMemorySA::Original(OriginalSA(sa, sample_rate)))
        } else {
            let sa = load_compressed(reader, bits_per_value, size)?;
            Ok(InMemorySA::Compressed(CompressedSA(sa, sample_rate)))
        }
    }
}

impl LoadIndex for InMemorySA {
    fn load(path: &std::path::Path) -> Result<Self, Box<dyn Error>> {
        text_compression::load_owned(path)
    }
}

impl WriteBinary for InMemorySA {
    fn write_binary<W: std::io::Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        dispatch!(self, write_binary(writer))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use text_compression::ReadBinary;

    use super::{
        super::test_utils::{
            assert_backend_holds, fit_to_width, owned_compressed, sample_sa, to_binary, to_file_bytes
        },
        InMemorySA, OriginalSA, SuffixArrayBackend
    };

    /// Loads a serialised array and asserts `read_binary` chose the variant the width calls for.
    fn load(bytes: &[u8], expect_compressed: bool) -> InMemorySA {
        let loaded = InMemorySA::read_binary(&mut Cursor::new(bytes)).unwrap();
        assert_eq!(matches!(loaded, InMemorySA::Compressed(_)), expect_compressed, "wrong variant chosen");
        loaded
    }

    /// A 64-bit header must come back as `Original`, entries intact.
    #[test]
    fn roundtrip_original() {
        let sa = sample_sa(300);
        let loaded = load(&to_binary(OriginalSA(sa.clone(), 3)), false);
        assert_backend_holds(&loaded, &sa, 3, 64);
    }

    /// Anything narrower must come back as `Compressed`. The widths straddle a `u64` boundary
    /// differently, which is where the header-sniffing branch would show a packing mismatch.
    #[test]
    fn roundtrip_compressed() {
        for bits in [8usize, 29, 40] {
            let sa = fit_to_width(&sample_sa(300), bits);
            let loaded = load(&to_binary(owned_compressed(&sa, 2, bits)), true);
            assert_backend_holds(&loaded, &sa, 2, bits);
        }
    }

    /// `CompressedSA::write_binary` and `dump_compressed_suffix_array` reach the body by different
    /// routes; a file from either must load the same way, since `sa-builder` uses the latter.
    #[test]
    fn both_compressed_writers_produce_the_same_file() {
        for bits in [8usize, 29, 40] {
            let sa = fit_to_width(&sample_sa(120), bits);
            let through_backend = to_binary(owned_compressed(&sa, 1, bits));
            let through_dump = to_file_bytes(&sa, 1, Some(bits));
            assert_eq!(through_backend, through_dump, "the two compressed writers disagree at {bits} bits");
        }
    }

    /// `iter_range` caches decoding state that `get` re-derives per call, so the two must agree at
    /// every start offset — including ranges opening mid-word — for both variants.
    #[test]
    fn iter_range_agrees_with_get() {
        let original = load(&to_binary(OriginalSA(sample_sa(200), 1)), false);
        let compressed_sa = fit_to_width(&sample_sa(200), 29);
        let compressed = load(&to_binary(owned_compressed(&compressed_sa, 1, 29)), true);

        for loaded in [&original, &compressed] {
            for start in [0usize, 1, 7, 63, 64, 65, 130] {
                for end in [start, start + 1, start + 13, start + 64, 200] {
                    if end > 200 {
                        continue;
                    }
                    let by_iter: Vec<i64> = loaded.iter_range(start, end).collect();
                    let by_get: Vec<i64> = (start..end).map(|i| loaded.get(i)).collect();
                    assert_eq!(by_iter, by_get, "iter_range({start}, {end})");
                    assert_eq!(loaded.iter_range(start, end).len(), end - start);
                }
            }
        }
    }

    /// This is the loader `sa-server` calls in the preloaded configuration. A file that is not an
    /// index — or one cut short of what its header promises — must come back as an error rather
    /// than an array that answers nonsense.
    #[test]
    fn read_binary_rejects_malformed_input() {
        let short_body = {
            let mut bytes = to_file_bytes(&sample_sa(5), 1, None);
            bytes.truncate(bytes.len() - 8);
            bytes
        };
        let cases: [(&str, Vec<u8>); 5] = [
            ("empty input", vec![]),
            ("width only", vec![64u8]),
            ("no count field", vec![64u8, 1]),
            ("count field cut short", vec![64u8, 1, 5, 0, 0, 0, 0, 0, 0]),
            ("body one entry short of the declared count", short_body)
        ];

        for (case, bytes) in cases {
            let result = InMemorySA::read_binary(&mut Cursor::new(bytes));
            assert!(result.is_err(), "{case} was accepted");
        }
    }

    /// `prefetch_sa_index` was covered for the mmap backend and for neither owned one, though both
    /// compute the address they hint by hand — `CompressedSA` scales the index by its bit width and
    /// takes a one-word slice, which is the arithmetic most likely to run off the end. Exercised
    /// through `InMemorySA` as well, since that forwards to whichever variant is live.
    #[test]
    fn prefetch_hints_are_harmless() {
        use super::super::test_utils::assert_prefetch_is_harmless;

        let sa = sample_sa(300);
        assert_prefetch_is_harmless(&OriginalSA(sa.clone(), 1), &sa);
        assert_prefetch_is_harmless(&load(&to_binary(OriginalSA(sa.clone(), 1)), false), &sa);

        for bits in [8usize, 29, 40] {
            let sa = fit_to_width(&sample_sa(300), bits);
            assert_prefetch_is_harmless(&owned_compressed(&sa, 2, bits), &sa);
            assert_prefetch_is_harmless(&load(&to_binary(owned_compressed(&sa, 2, bits)), true), &sa);
        }
    }

    /// A width outside `1..=64` must be refused by *both* readers, identically.
    ///
    /// Neither used to check. A `0` header made the declared body zero-length, so the file passed
    /// every size check and the first `get` read off the end of it; a `200` header overflowed the
    /// shift. Both cases panicked at lookup rather than erroring at load, on both backends.
    #[test]
    fn both_readers_reject_an_impossible_width() {
        use text_compression::ReadBinaryMmap;

        use crate::array::{MmapBackedSA, mmap::test_utils::write_to_tempfile};

        for bad in [0u8, 65, 200, 255] {
            let mut bytes = to_file_bytes(&sample_sa(5), 1, None);
            bytes[0] = bad;

            assert!(
                InMemorySA::read_binary(&mut Cursor::new(bytes.clone())).is_err(),
                "preloaded accepted a width of {bad}"
            );

            let tmp = write_to_tempfile(&bytes);
            assert!(MmapBackedSA::read_binary_mmap(tmp.path()).is_err(), "mmap accepted a width of {bad}");
        }
    }

    /// A truncated body must be refused by *both* readers, identically.
    ///
    /// The compressed preloaded path was the one sibling that never related the body it decoded to
    /// the count its header declared: `read_binary` refills the backing store with however many
    /// words the reader yielded, and nothing compared that to `len`. A short `sa.bin` therefore
    /// loaded cleanly, reported the declared entry count, and panicked on the first probe past the
    /// real data — inside a request handler, since this is a server startup path.
    #[test]
    fn both_readers_reject_a_truncated_body() {
        use text_compression::ReadBinaryMmap;

        use crate::array::{MmapBackedSA, mmap::test_utils::write_to_tempfile};

        for bits in [None, Some(29), Some(40)] {
            let sa = match bits {
                None => sample_sa(100),
                Some(b) => fit_to_width(&sample_sa(100), b)
            };
            let full = to_file_bytes(&sa, 1, bits);

            // A quarter of the body removed: unambiguously short, whatever the rounding.
            let cut = 10 + (full.len() - 10) * 3 / 4;

            let preloaded = InMemorySA::read_binary(&mut Cursor::new(full[..cut].to_vec()));
            assert!(preloaded.is_err(), "preloaded accepted a truncated body at bits={bits:?}");

            let tmp = write_to_tempfile(&full[..cut]);
            assert!(
                MmapBackedSA::read_binary_mmap(tmp.path()).is_err(),
                "mmap accepted a truncated body at bits={bits:?}"
            );

            // The intact file still loads and reads back, so the check above rejects for the right
            // reason rather than rejecting everything.
            let ok = InMemorySA::read_binary(&mut Cursor::new(full.clone())).expect("intact file must load");
            assert_backend_holds(&ok, &sa, 1, bits.unwrap_or(64));
        }
    }

    /// A sample rate of zero must be refused by *both* readers.
    ///
    /// It is the producer half of a pair: `sa-builder --sparseness-factor 0` used to be accepted,
    /// wrote `0` here, and produced an index that loaded cleanly and then matched no peptide at
    /// all, because the sample rate is the minimum searchable peptide length. The flag is
    /// range-checked now, but an index built before that still loads, so the readers reject it too.
    #[test]
    fn both_readers_reject_a_sample_rate_of_zero() {
        use text_compression::ReadBinaryMmap;

        use crate::array::{MmapBackedSA, mmap::test_utils::write_to_tempfile};

        for bits in [None, Some(29)] {
            let sa = match bits {
                None => sample_sa(20),
                Some(b) => fit_to_width(&sample_sa(20), b)
            };
            let mut bytes = to_file_bytes(&sa, 1, bits);
            bytes[1] = 0; // the sample-rate byte

            assert!(
                InMemorySA::read_binary(&mut Cursor::new(bytes.clone())).is_err(),
                "preloaded accepted a sample rate of 0 at bits={bits:?}"
            );
            let tmp = write_to_tempfile(&bytes);
            assert!(
                MmapBackedSA::read_binary_mmap(tmp.path()).is_err(),
                "mmap accepted a sample rate of 0 at bits={bits:?}"
            );

            // A rate of 1 is the ordinary case and must still load, so the check above is specific.
            bytes[1] = 1;
            assert!(InMemorySA::read_binary(&mut Cursor::new(bytes.clone())).is_ok());
            let tmp = write_to_tempfile(&bytes);
            assert!(MmapBackedSA::read_binary_mmap(tmp.path()).is_ok());
        }
    }

    /// The trait's `is_empty` default, which nothing else reaches.
    #[test]
    fn empty_arrays_report_empty() {
        assert!(InMemorySA::Original(OriginalSA(vec![], 1)).is_empty());
        assert!(!InMemorySA::Original(OriginalSA(vec![7], 1)).is_empty());
    }

    /// The two storage backends read the same bytes and must answer identically; the searcher is
    /// written against whichever one it is handed. `sa_searcher::backend_agreement` makes the same
    /// check end to end, but this one localises a failure to the array.
    #[test]
    fn agrees_with_the_mmap_backend() {
        use crate::array::mmap::test_utils::write_and_map;

        for bits in [None, Some(29), Some(40)] {
            let sa = match bits {
                None => sample_sa(250),
                Some(b) => fit_to_width(&sample_sa(250), b)
            };
            let preloaded = InMemorySA::read_binary(&mut Cursor::new(to_file_bytes(&sa, 2, bits))).unwrap();
            let (mapped, _tmp) = write_and_map(&sa, 2, bits);

            assert_backend_holds(&preloaded, &sa, 2, bits.unwrap_or(64));
            assert_backend_holds(&mapped, &sa, 2, bits.unwrap_or(64));
        }
    }
}
