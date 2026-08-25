//! The shared test suite for the two bit-array implementations.
//!
//! `BitArray<BITS>` and `DynBitArray` are deliberately separate types (see the module docs on
//! each for why), but they must behave identically. Rather than maintain two copies of the same
//! ~150-line test module — which is what this crate did before — the suite is written once here
//! and instantiated by each implementation.
//!
//! It is instantiated *inside* the implementation's own module, so the tests keep direct access
//! to the private `data` and `len` fields and nothing has to be made `pub(crate)` for testing.
//!
//! Each implementation supplies a `$new!(capacity, BITS)` macro that hides the one real
//! difference between them: `BitArray` takes its width as a const generic parameter,
//! `DynBitArray` as a constructor argument.

/// One width's worth of the "overwriting replaces" property, for [`bitarray_test_suite`].
///
/// A separate macro rather than a loop because `BitArray` takes its width as a *const generic*, so
/// the constructor needs a literal. Widths are therefore enumerated at the call site.
macro_rules! assert_overwrite_replaces {
    ($new:ident, $bits:literal) => {{
        let bits: usize = $bits;
        let max = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };
        // 40 values reaches a word boundary at every width in the list.
        let mut ba = $new!(40, $bits);

        // Fill with all-ones, then overwrite each index with a distinct pattern. Any bit left
        // over from the first pass shows up as a wrong read.
        for i in 0..40 {
            ba.set(i, max);
        }
        for i in 0..40 {
            ba.set(i, (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) & max);
        }
        for i in 0..40 {
            let want = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) & max;
            assert_eq!(ba.get(i), want, "bits={bits} index={i}: stale bits after overwrite");
        }

        // And back to zero, where every leftover bit is visible.
        for i in 0..40 {
            ba.set(i, 0);
        }
        for i in 0..40 {
            assert_eq!(ba.get(i), 0, "bits={bits} index={i}: not cleared");
        }
    }};
}

/// One width's worth of the "an over-range value stays in its field" property.
macro_rules! assert_over_range_is_contained {
    ($new:ident, $bits:literal) => {{
        let bits: usize = $bits;
        let max = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };
        let mut ba = $new!(8, $bits);

        // Neighbours start at zero, so any bit that spills into them is visible. Filling with
        // `max` instead would hide the bug entirely: the spilled bits are ones, and ORing ones
        // into a field that is already all ones changes nothing.
        ba.set(4, u64::MAX); // deliberately far too wide for the field

        assert_eq!(ba.get(4), max, "bits={bits}: an over-range value should be truncated to the field");
        for i in 0..8 {
            if i == 4 {
                continue;
            }
            assert_eq!(ba.get(i), 0, "bits={bits}: index {i} disturbed by an over-range write at 4");
        }
    }};
}

/// One width's worth of the degenerate-length behaviour.
macro_rules! assert_degenerate_lengths {
    ($new:ident, $bits:literal) => {{
        let bits: usize = $bits;
        let max = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };

        let empty = $new!(0, $bits);
        assert_eq!(empty.len(), 0, "bits={bits}: empty length");
        assert!(empty.is_empty(), "bits={bits}: empty is_empty");

        let mut one = $new!(1, $bits);
        assert_eq!(one.len(), 1, "bits={bits}: single length");
        assert_eq!(one.get(0), 0, "bits={bits}: single starts zeroed");
        one.set(0, max);
        assert_eq!(one.get(0), max, "bits={bits}: single round trip");
        one.set(0, 0);
        assert_eq!(one.get(0), 0, "bits={bits}: single clear");
    }};
}

/// Generates the full behavioural test suite against a `$new!(capacity, BITS)` constructor macro.
macro_rules! bitarray_test_suite {
    ($new:ident) => {
        #[cfg(test)]
        mod tests {
            use super::*;

            #[test]
            fn test_with_capacity() {
                let ba = $new!(4, 40);
                assert_eq!(ba.data, vec![0, 0, 0]);
                assert_eq!(ba.len, 4);
            }

            #[test]
            fn test_get() {
                let mut ba = $new!(4, 40);
                ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];

                assert_eq!(ba.get(0), 0b0001110011111010110001000111111100110010);
                assert_eq!(ba.get(1), 0b1100001001010010011000010100110111001001);
                assert_eq!(ba.get(2), 0b1111001101001101101101101011101001010001);
                assert_eq!(ba.get(3), 0b0000100010010001010001001110101110011100);
            }

            /// Overwriting an index must replace the value, not OR the new one over the old.
            ///
            /// The straddling branch cleared `MASK >> start_block_offset` bits of the first word,
            /// which for a field starting near the top of a word is zero — so it cleared nothing
            /// and the previous value stayed underneath the new one. The parity test suite never
            /// wrote an index twice, which is why both implementations carried this identically.
            ///
            /// The widths below are chosen to cover fields that divide 64 evenly, fields that do
            /// not, and both ends of the range.
            #[test]
            fn overwriting_an_index_replaces_the_value() {
                assert_overwrite_replaces!($new, 1);
                assert_overwrite_replaces!($new, 3);
                assert_overwrite_replaces!($new, 5);
                assert_overwrite_replaces!($new, 8);
                assert_overwrite_replaces!($new, 17);
                assert_overwrite_replaces!($new, 28);
                assert_overwrite_replaces!($new, 32);
                assert_overwrite_replaces!($new, 40);
                assert_overwrite_replaces!($new, 63);
                assert_overwrite_replaces!($new, 64);
            }

            /// A value wider than the field must not disturb its neighbours.
            ///
            /// `set` ORs the value into place, so an unmasked over-range value spilled its high
            /// bits into the adjacent entry — corrupting a value the caller never wrote to.
            #[test]
            fn an_over_range_value_cannot_touch_its_neighbours() {
                assert_over_range_is_contained!($new, 1);
                assert_over_range_is_contained!($new, 3);
                assert_over_range_is_contained!($new, 5);
                assert_over_range_is_contained!($new, 17);
                assert_over_range_is_contained!($new, 28);
                assert_over_range_is_contained!($new, 40);
                assert_over_range_is_contained!($new, 63);
            }

            /// The degenerate lengths, which the rest of the suite never builds.
            #[test]
            fn zero_and_one_element_arrays_behave() {
                assert_degenerate_lengths!($new, 1);
                assert_degenerate_lengths!($new, 5);
                assert_degenerate_lengths!($new, 17);
                assert_degenerate_lengths!($new, 32);
                assert_degenerate_lengths!($new, 63);
                assert_degenerate_lengths!($new, 64);
            }

            #[test]
            fn test_set() {
                let mut ba = $new!(4, 40);

                ba.set(0, 0b0001110011111010110001000111111100110010_u64);
                ba.set(1, 0b1100001001010010011000010100110111001001_u64);
                ba.set(2, 0b1111001101001101101101101011101001010001_u64);
                ba.set(3, 0b0000100010010001010001001110101110011100_u64);

                assert_eq!(ba.data, vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000]);
            }

            #[test]
            fn test_bits_per_value() {
                assert_eq!($new!(4, 40).bits_per_value(), 40);
            }

            #[test]
            fn test_len_and_empty() {
                assert_eq!($new!(4, 40).len(), 4);
                assert!($new!(0, 40).is_empty());
                assert!(!$new!(4, 40).is_empty());
            }

            #[test]
            fn test_clear() {
                let mut ba = $new!(4, 40);
                ba.data = vec![0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4];
                ba.clear();
                assert_eq!(ba.data, vec![0, 0, 0]);
            }

            #[test]
            fn test_write_binary() {
                let mut ba = $new!(4, 40);
                ba.set(0, 0x1234567890_u64);
                ba.set(1, 0xabcdef0123_u64);
                ba.set(2, 0x4567890abc_u64);
                ba.set(3, 0xdef0123456_u64);

                let mut buf = Vec::new();
                ba.write_binary(&mut buf).unwrap();

                assert_eq!(buf, vec![
                    0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01,
                    0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0,
                ]);
            }

            #[test]
            fn test_read_binary() {
                let buf = [
                    0xef_u8, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12, 0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01,
                    0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
                ];
                let mut ba = $new!(4, 40);
                ba.read_binary(&buf[..]).unwrap();

                assert_eq!(ba.get(0), 0x1234567890);
                assert_eq!(ba.get(1), 0xabcdef0123);
                assert_eq!(ba.get(2), 0x4567890abc);
                assert_eq!(ba.get(3), 0xdef0123456);
            }

            #[test]
            fn test_iter_range_empty() {
                let ba = $new!(8, 32);
                assert!(ba.iter_range(3, 3).next().is_none());
                assert!(ba.iter_range(5, 3).next().is_none(), "reversed range yields nothing");
            }

            #[test]
            fn test_iter_range_single_entry() {
                let mut ba = $new!(4, 40);
                ba.set(2, 0xABCDEF1234_u64);
                assert_eq!(ba.iter_range(2, 3).collect::<Vec<i64>>(), vec![0xABCDEF1234_i64]);
            }

            #[test]
            fn test_iter_range_mid_block_start() {
                let mut ba = $new!(8, 32);
                for i in 0..8 {
                    ba.set(i, (i as u64) * 111 + 7);
                }
                assert_eq!(
                    ba.iter_range(1, 6).collect::<Vec<i64>>(),
                    (1..6).map(|i| ba.get(i) as i64).collect::<Vec<i64>>()
                );
            }

            #[test]
            fn test_iter_range_crosses_block_boundary() {
                let mut ba = $new!(16, 40);
                for i in 0..16 {
                    ba.set(i, i as u64 * 0x100000001 + 3);
                }
                for (start, end) in [(0, 16), (3, 13)] {
                    assert_eq!(
                        ba.iter_range(start, end).collect::<Vec<i64>>(),
                        (start..end).map(|i| ba.get(i) as i64).collect::<Vec<i64>>()
                    );
                }
            }

            #[test]
            fn test_iter_range_bits_per_value_64() {
                let mut ba = $new!(8, 64);
                for i in 0..8 {
                    ba.set(i, i as u64 * 0xDEAD_BEEF + 1);
                }
                for (start, end) in [(0, 8), (2, 6)] {
                    assert_eq!(
                        ba.iter_range(start, end).collect::<Vec<i64>>(),
                        (start..end).map(|i| ba.get(i) as i64).collect::<Vec<i64>>()
                    );
                }
            }

            #[test]
            fn test_iter_range_bits_per_value_1() {
                let mut ba = $new!(128, 1);
                for i in (0..128).step_by(3) {
                    ba.set(i, 1);
                }
                for (start, end) in [(0, 128), (60, 70)] {
                    assert_eq!(
                        ba.iter_range(start, end).collect::<Vec<i64>>(),
                        (start..end).map(|i| ba.get(i) as i64).collect::<Vec<i64>>()
                    );
                }
            }

            #[test]
            fn test_iter_range_exact_size() {
                let mut ba = $new!(10, 40);
                for i in 0..10 {
                    ba.set(i, i as u64 * 99);
                }
                assert_eq!(ba.iter_range(2, 8).len(), 6);
            }

            /// Every width, not just the handful the cases above pick. Catches off-by-ones in
            /// the `start_bit + BITS <= 64` split that only appear at particular alignments.
            #[test]
            fn test_iter_range_matches_get_at_every_offset() {
                let mut ba = $new!(64, 7);
                for i in 0..64 {
                    ba.set(i, (i as u64 * 37) & 0x7f);
                }
                for start in 0..64 {
                    for end in start..64 {
                        assert_eq!(
                            ba.iter_range(start, end).collect::<Vec<i64>>(),
                            (start..end).map(|i| ba.get(i) as i64).collect::<Vec<i64>>(),
                            "iter_range({start}, {end}) disagreed with get"
                        );
                    }
                }
            }
        }
    };
}

/// Asserts the two implementations agree bit-for-bit at a given width.
///
/// This is the test the duplicated suites could not express, and it is the one that matters:
/// nothing else pins `BitArray<N>` and `DynBitArray` to the same packing. They are separate types
/// purely so that one can const-fold its shifts, and a divergence would mean an index written by
/// one and read by the other decodes to garbage — which is exactly what the compressed suffix
/// array does, writing through `DynBitArray` and reading through either.
#[cfg(test)]
macro_rules! assert_implementations_agree {
    ($bits:literal) => {{
        const BITS: usize = $bits;
        const N: usize = 200;

        let mask: u64 = if BITS == 64 { u64::MAX } else { (1u64 << BITS) - 1 };
        // Deterministic pseudo-random values; a plain counter would never exercise the high bits.
        let values: Vec<u64> = (0..N).map(|i| (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) & mask).collect();

        let mut constant = $crate::BitArray::<BITS>::with_capacity(N);
        let mut dynamic = $crate::DynBitArray::with_capacity(N, BITS);
        for (i, &v) in values.iter().enumerate() {
            constant.set(i, v);
            dynamic.set(i, v);
        }

        for i in 0..N {
            assert_eq!(constant.get(i), values[i], "BitArray<{BITS}> lost value at {i}");
            assert_eq!(dynamic.get(i), values[i], "DynBitArray({BITS}) lost value at {i}");
        }

        // Identical packing, not merely identical values: the serialised bytes must match, since
        // one implementation writes the files the other reads.
        let (mut a, mut b) = (Vec::new(), Vec::new());
        $crate::Binary::write_binary(&constant, &mut a).unwrap();
        $crate::Binary::write_binary(&dynamic, &mut b).unwrap();
        assert_eq!(a, b, "BitArray<{BITS}> and DynBitArray({BITS}) packed differently");

        // And the sequential fast paths agree with each other over a block-crossing range.
        assert_eq!(
            constant.iter_range(7, N - 7).collect::<Vec<i64>>(),
            dynamic.iter_range(7, N - 7).collect::<Vec<i64>>(),
            "range iterators disagreed at {BITS} bits"
        );
    }};
}

/// Expands `assert_implementations_agree!` once per listed width.
#[cfg(test)]
macro_rules! seq_bits {
    ($($bits:literal),+ $(,)?) => { $( assert_implementations_agree!($bits); )+ };
}

#[cfg(test)]
mod parity {
    /// Sweeps every width the two implementations can be built at. `BitArray`'s width is a const
    /// generic, so the list has to be spelled out rather than looped over.
    #[test]
    fn implementations_agree_at_every_width() {
        seq_bits!(
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
            57, 58, 59, 60, 61, 62, 63, 64
        );
    }
}
