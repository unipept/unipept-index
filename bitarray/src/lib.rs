//! Dense arrays of fixed-width values packed into `u64` words.
//!
//! The index stores hundreds of millions of values that need far fewer than 64 bits each — a
//! suffix array over a 300 M-residue text needs 29, and the protein text itself needs 5 for its
//! 24-letter alphabet. Packing them at their natural width rather than rounding up to a byte is
//! what keeps the compressed index small enough to be worth loading.
//!
//! Values are packed **most-significant-bit first within each little-endian `u64` word**, and a
//! value may straddle a word boundary. Both implementations below and the mmap readers in
//! `protein-text` and `sa-index` depend on that layout matching exactly.
//!
//! # Two implementations
//!
//! [`BitArray<BITS>`] fixes the width at compile time; [`DynBitArray`] takes it at runtime. They
//! are otherwise interchangeable, and `test_suite.rs` asserts they pack identically at every
//! width from 1 to 64. Which to use:
//!
//! * **[`BitArray<BITS>`] when the width is a property of the data.** The protein text is always
//!   5 bits per residue, so `protein-text` uses `BitArray<5>`. See [`constant`] for the
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
//!
//! Both constructors advise transparent huge pages over the allocation before they zero it, via
//! [`memory_hints::hugepages`]. That ordering is the whole point and is argued there; callers do
//! not need to repeat it.
#![warn(missing_docs)]

#[cfg(test)]
#[macro_use]
mod test_suite;

pub mod binary;
pub mod constant;
pub mod dynamic;

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
/// # The chunk size is rounded, not trusted
///
/// Each chunk is written as a whole number of `u64` words, so for the chunks to concatenate into
/// one continuous bit stream every chunk must occupy a whole number of words — that is,
/// `capacity * bits_per_value` must be a multiple of 64. `max_capacity` is therefore a *maximum*
/// rather than the size used: it is rounded down to the nearest multiple of `64 / gcd`, which is
/// the smallest number of values that fills whole words at this width. Any `max_capacity` is
/// consequently safe, and the caller does not have to know the rule.
pub fn data_to_writer(
    data: Vec<i64>,
    bits_per_value: usize,
    max_capacity: usize,
    writer: &mut impl Write
) -> Result<()> {
    // Round the requested chunk size down to a multiple of `64 / gcd(bits_per_value, 64)`.
    //
    // That step is the smallest number of values whose packed width is a whole number of `u64`
    // words: `capacity * bits_per_value % 64 == 0` holds exactly when `capacity` is a multiple of
    // it. Rounding to `gcd` itself — which is what this did — is a different quantity and does not
    // imply the invariant: at `bits_per_value = 5` the gcd is 1, so it rounded to nothing at all,
    // and a `max_capacity` of 100 emitted a partly-filled trailing word per chunk, corrupting the
    // stream in release builds where the `debug_assert` below is compiled out.
    //
    // Both step sizes are powers of two no greater than 64, so for the `8 * 1024` both callers
    // pass this rounds to the same chunk size at every width from 1 to 64 and the bytes this
    // function emits are unchanged.
    let step = 64 / gcd(bits_per_value, 64);
    let capacity = max(step, max_capacity / step * step);
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

    /// A chunked write must produce the same bytes as packing everything in one array, at every
    /// width and for a `max_capacity` that is not a multiple of 64.
    ///
    /// `max_capacity = 100` is the case the old `gcd` rounding got wrong: it rounded to a multiple
    /// of `gcd(bits, 64)` rather than of `64 / gcd`, so at most widths the chunk did not fill whole
    /// words and every chunk boundary shifted the rest of the stream by a few bits. The
    /// `debug_assert` in `data_to_writer` catches it in a debug build; this catches it in any
    /// build, and pins the output to the unchunked packing rather than merely to itself.
    #[test]
    fn chunked_and_unchunked_writes_agree_at_every_width() {
        for bits in 1..=64 {
            let mask = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };
            let values: Vec<i64> = (0..250u64).map(|i| (i.wrapping_mul(0x9E37_79B9_7F4A_7C15) & mask) as i64).collect();

            let mut chunked = Vec::new();
            data_to_writer(values.clone(), bits, 100, &mut chunked).unwrap();

            let mut whole = DynBitArray::with_capacity(values.len(), bits);
            for (i, &v) in values.iter().enumerate() {
                whole.set(i, v as u64);
            }
            let mut unchunked = Vec::new();
            whole.write_binary(&mut unchunked).unwrap();

            assert_eq!(chunked, unchunked, "chunked write differs from unchunked at {bits} bits");
        }
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
