// Non-mmap builds only — handles the runtime Original-vs-Compressed decision.
use std::{error::Error, io::BufRead};

use bitarray::BitArrayRangeIter;
use text_compression::{ReadBinary, WriteBinary};

use super::{SuffixArrayBackend, OriginalSA, OriginalRangeIter, CompressedSA};
use super::original::load_original;
use super::compressed::load_compressed;

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

/// Wraps either an Original or Compressed SA loaded from disk.
/// The variant is determined at runtime by the `bits_per_value` field in the binary header.
pub enum InMemorySA {
    Original(OriginalSA),
    Compressed(CompressedSA),
}

// ── InMemoryRangeIter ────────────────────────────────────────────────────────

pub enum InMemoryRangeIter<'a> {
    Original(OriginalRangeIter<'a>),
    Compressed(BitArrayRangeIter<'a>),
}

impl Iterator for InMemoryRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<i64> { dispatch!(self, next()) }

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

    fn len(&self) -> usize            { dispatch!(self, len()) }
    fn bits_per_value(&self) -> usize { dispatch!(self, bits_per_value()) }
    fn sample_rate(&self) -> u8       { dispatch!(self, sample_rate()) }
    #[inline]
    fn get(&self, index: usize) -> i64 { dispatch!(self, get(index)) }
    // iter_range needs a manual match: each arm wraps its backend's native iterator
    // type into the appropriate InMemoryRangeIter variant.
    fn iter_range(&self, start: usize, end: usize) -> InMemoryRangeIter<'_> {
        match self {
            Self::Original(b)   => InMemoryRangeIter::Original(b.iter_range(start, end)),
            Self::Compressed(b) => InMemoryRangeIter::Compressed(b.iter_range(start, end)),
        }
    }

    #[inline]
    fn prefetch_sa_index(&self, index: usize) { dispatch!(self, prefetch_sa_index(index)) }
}

// ── ReadBinary / WriteBinary ──────────────────────────────────────────────────

impl ReadBinary for InMemorySA {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf1 = [0u8; 1];
        reader.read_exact(&mut buf1).map_err(|_| "Could not read the required bits from the binary file")?;
        let bits_per_value = buf1[0] as usize;

        reader.read_exact(&mut buf1).map_err(|_| "Could not read the sample rate from the binary file")?;
        let sample_rate = buf1[0];

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8).map_err(|_| "Could not read the size of the suffix array from the binary file")?;
        let size = u64::from_le_bytes(buf8) as usize;

        if bits_per_value == 64 {
            let sa = load_original(reader, sample_rate, size)?;
            Ok(InMemorySA::Original(OriginalSA(sa, sample_rate)))
        } else {
            let sa = load_compressed(reader, bits_per_value, size)?;
            Ok(InMemorySA::Compressed(CompressedSA(sa, sample_rate)))
        }
    }
}

impl WriteBinary for InMemorySA {
    fn write_binary<W: std::io::Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        dispatch!(self, write_binary(writer))
    }
}
