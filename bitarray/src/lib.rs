//! Dense arrays of fixed-width values packed into `u64` words.
//!
//! The index stores hundreds of millions of values that need far fewer than 64 bits each — a
//! suffix array over a 300 M-residue text needs 29, and the protein text itself needs 5 for its
//! 24-letter alphabet. Packing them at their natural width rather than rounding up to a byte is
//! what keeps the compressed index small enough to be worth loading.
//!
//! Values are packed **most-significant-bit first within each little-endian `u64` word**, and a
//! value may straddle a word boundary. Both implementations below and the mmap readers in
//! `text-compression` and `sa-index` depend on that layout matching exactly.
//!
//! # Two implementations
//!
//! [`BitArray<BITS>`] fixes the width at compile time; [`DynBitArray`] takes it at runtime. They
//! are otherwise interchangeable, and `test_suite.rs` asserts they pack identically at every
//! width from 1 to 64. Which to use:
//!
//! * **[`BitArray<BITS>`] when the width is a property of the data.** The protein text is always
//!   5 bits per residue, so `text-compression` uses `BitArray<5>`. See [`constant`] for the
//!   const-folding this buys — a mechanism, not a measurement: nothing in the workspace benchmarks
//!   the two implementations against each other.
//! * **[`DynBitArray`] when the width comes from a file header.** The compressed suffix array
//!   chooses its width from the text length at build time and records it in the file, so the
//!   reader cannot know it until runtime.
//!
//! # The rest of the crate
//!
//! * [`Binary`] serialises a bit array's backing words — headerless, and read to EOF; the contract
//!   is worth reading before embedding one in a larger file.
//! * [`data_to_writer`] packs and writes a `Vec<i64>` in chunks, so building a file does not need
//!   a second full copy of the index in memory.
//! * [`hugepages`] explains why both constructors advise transparent huge pages *before* they
//!   touch the allocation, and why doing it afterwards would be pointless.
#![warn(missing_docs)]

#[cfg(test)]
#[macro_use]
mod test_suite;

mod binary;
pub mod constant;
pub mod dynamic;
pub mod hugepages;

use std::{
    cmp::max,
    io::{Result, Write}
};

pub use binary::Binary;
pub use constant::{BitArray, BitArrayRangeIter};
pub use dynamic::{DynBitArray, DynBitArrayRangeIter};

/// Writes packed bit data to a writer in chunks, minimising peak memory.
///
/// Builds and serialises `max_capacity` values at a time instead of packing the whole of `data`
/// into one array first. At index-build scale `data` is already several gigabytes, so the
/// difference is between one transient buffer and a second full copy of the index.
///
/// # Panics in the caller's future, if `max_capacity` is chosen badly
///
/// Each chunk is written as a whole number of `u64` words, so for the chunks to concatenate into
/// one continuous bit stream every chunk must occupy a whole number of words — that is,
/// `max_capacity * bits_per_value` must be a multiple of 64. The one caller in the workspace —
/// `sa_index::array::preloaded::compressed::dump_compressed_suffix_array` — passes `8 * 1024`,
/// which satisfies this for every width. See the note on the chunk size below.
pub fn data_to_writer(
    data: Vec<i64>,
    bits_per_value: usize,
    max_capacity: usize,
    writer: &mut impl Write
) -> Result<()> {
    // Round the requested chunk size down to a multiple of gcd(bits_per_value, 64).
    //
    // CAUTION: this is *not* the invariant the chunking actually needs. A chunk only lands on a
    // word boundary when `capacity * bits_per_value % 64 == 0`, i.e. when `capacity` is a
    // multiple of `64 / gcd`, not of `gcd`. For `bits_per_value = 5` the gcd is 1 and this line
    // rounds to nothing at all; the code is correct today only because the single caller passes
    // `8 * 1024`, which is a multiple of 64 and therefore satisfies the real invariant for every
    // width. A `max_capacity` of, say, 100 would silently emit a partly-filled trailing word per
    // chunk and corrupt the stream — and the `debug_assert` below only catches it in debug builds.
    // Tracked as a known issue; not changed here because doing so would alter the bytes this
    // function emits.
    let gcd = gcd(bits_per_value, 64);
    let capacity = max(gcd, max_capacity / gcd * gcd);
    debug_assert!(
        (capacity * bits_per_value).is_multiple_of(64),
        "chunk of {capacity} values at {bits_per_value} bits does not fill whole u64 words"
    );

    if data.len() <= capacity {
        let mut ba = DynBitArray::with_capacity(data.len(), bits_per_value);
        for (i, &v) in data.iter().enumerate() {
            ba.set(i, v as u64);
        }
        ba.write_binary(writer)?;
        return Ok(());
    }

    let mut ba = DynBitArray::with_capacity(capacity, bits_per_value);
    let chunks = data.chunks_exact(capacity);
    let remainder = chunks.remainder();

    for chunk in chunks {
        for (i, &v) in chunk.iter().enumerate() {
            ba.set(i, v as u64);
        }
        ba.write_binary(writer)?;
        ba.clear();
    }

    ba = DynBitArray::with_capacity(remainder.len(), bits_per_value);
    for (i, &v) in remainder.iter().enumerate() {
        ba.set(i, v as u64);
    }
    ba.write_binary(writer)?;

    Ok(())
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        if b < a {
            std::mem::swap(&mut b, &mut a);
        }
        b %= a;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(40, 64), 8);
        assert_eq!(gcd(64, 40), 8);
        assert_eq!(gcd(64, 64), 64);
        assert_eq!(gcd(32, 64), 32);
    }

    #[test]
    fn test_data_to_writer_no_chunks_needed() {
        let data = vec![0x1234567890, 0xabcdef0123, 0x4567890abc, 0xdef0123456];
        let mut writer = Vec::new();
        data_to_writer(data, 40, 2, &mut writer).unwrap();
        assert_eq!(writer, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01, 0x00, 0x00,
            0x00, 0x00, 0x56, 0x34, 0x12, 0xf0,
        ]);
    }

    #[test]
    fn test_data_to_writer_chunks_needed_no_remainder() {
        let data = vec![
            0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999,
            0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222,
            0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb,
            0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333, 0x44444444,
            0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd,
            0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666,
            0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff,
            0x00000000,
        ];
        let mut writer = Vec::new();
        data_to_writer(data, 32, 8, &mut writer).unwrap();
        assert_eq!(writer, vec![
            0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66,
            0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa,
            0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd,
            0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11,
            0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88,
            0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc,
            0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff,
            0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33,
            0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa,
            0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee,
            0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11,
            0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55,
            0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc,
            0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00,
            0xff, 0xff, 0xff, 0xff,
        ]);
    }

    #[test]
    fn test_data_to_writer_chunks_needed_plus_remainder() {
        let data = vec![
            0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999,
            0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222,
            0x33333333, 0x44444444, 0x55555555, 0x66666666, 0x77777777, 0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb,
            0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff, 0x00000000, 0x11111111, 0x22222222, 0x33333333,
        ];
        let mut writer = Vec::new();
        data_to_writer(data, 32, 8, &mut writer).unwrap();
        assert_eq!(writer, vec![
            0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66,
            0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88, 0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa,
            0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc, 0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd,
            0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11,
            0x44, 0x44, 0x44, 0x44, 0x33, 0x33, 0x33, 0x33, 0x66, 0x66, 0x66, 0x66, 0x55, 0x55, 0x55, 0x55, 0x88, 0x88,
            0x88, 0x88, 0x77, 0x77, 0x77, 0x77, 0xaa, 0xaa, 0xaa, 0xaa, 0x99, 0x99, 0x99, 0x99, 0xcc, 0xcc, 0xcc, 0xcc,
            0xbb, 0xbb, 0xbb, 0xbb, 0xee, 0xee, 0xee, 0xee, 0xdd, 0xdd, 0xdd, 0xdd, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff,
            0xff, 0xff, 0x22, 0x22, 0x22, 0x22, 0x11, 0x11, 0x11, 0x11, 0x00, 0x00, 0x00, 0x00, 0x33, 0x33, 0x33, 0x33,
        ]);
    }
}
