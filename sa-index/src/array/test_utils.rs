//! Fixtures shared by the suffix-array tests, of every packing and both storage backends.
//!
//! Fixtures that need a mapped file live in [`mmap::test_utils`](super::mmap::test_utils), and the
//! failing reader and writer that only the streaming backend has error paths for live in
//! [`preloaded::test_utils`](super::preloaded::test_utils).

use binary_traits::WriteBinary;
use bitarray::DynBitArray;

use super::{CompressedSA, SuffixArrayBackend};

/// A deterministic suffix array of `len` entries, spread across the value range so the high bits
/// of an entry are exercised rather than a run of small numbers that fits in a byte.
pub fn sample_sa(len: usize) -> Vec<i64> {
    (0..len).map(|i| ((i as i64).wrapping_mul(2_654_435_761)) & 0x0fff_ffff).collect()
}

/// Masks `sa` down to what `bits` can hold, for the widths the builder may choose.
pub fn fit_to_width(sa: &[i64], bits: usize) -> Vec<i64> {
    let mask = (u64::MAX >> (64 - bits)) as i64;
    sa.iter().map(|v| v & mask).collect()
}

/// Builds the owned compressed backend directly, bypassing the file format.
pub fn owned_compressed(sa: &[i64], sparseness: u8, bits: usize) -> CompressedSA {
    let mut bit_array = DynBitArray::with_capacity(sa.len(), bits);
    for (i, &v) in sa.iter().enumerate() {
        bit_array.set(i, v as u64);
    }
    CompressedSA(bit_array, sparseness)
}

/// Asserts a backend reports the metadata it was built with and reproduces `sa` entry for entry,
/// through both access paths — `get`, and `iter_range`, which decodes by its own route.
pub fn assert_backend_holds(backend: &impl SuffixArrayBackend, sa: &[i64], sparseness: u8, bits_per_value: usize) {
    assert_eq!(backend.len(), sa.len(), "length differs");
    assert_eq!(backend.bits_per_value(), bits_per_value, "width differs");
    assert_eq!(backend.sample_rate(), sparseness, "sparseness factor differs");
    assert_eq!(backend.is_empty(), sa.is_empty(), "emptiness differs");

    for (index, &expected) in sa.iter().enumerate() {
        assert_eq!(backend.get(index), expected, "entry {index} differs");
    }
    assert_eq!(backend.iter_range(0, sa.len()).collect::<Vec<i64>>(), sa, "iter_range disagrees with get");
}

/// Walks `prefetch_sa_index` over the whole array and well past its end, then re-reads every entry.
///
/// Every backend computes the address it hints by hand — the compressed one scales the index by its
/// bit width, the mmap one adds a file offset — so a wrong bound is an out-of-range index rather
/// than a wrong answer, and a hint may never disturb the entry it precedes. The trait promises
/// out-of-range indices are ignored rather than panicking, which is what the second half of the
/// walk exercises.
pub fn assert_prefetch_is_harmless(backend: &impl SuffixArrayBackend, sa: &[i64]) {
    for index in 0..sa.len() * 2 {
        backend.prefetch_sa_index(index);
    }
    for (index, &expected) in sa.iter().enumerate() {
        assert_eq!(backend.get(index), expected, "entry {index} differs after the hints");
    }
}

/// Serialises a backend through its [`WriteBinary`] impl. The first byte is the width the loaders
/// dispatch on.
pub fn to_binary(sa: impl WriteBinary) -> Vec<u8> {
    let mut buf = Vec::new();
    sa.write_binary(&mut buf).unwrap();
    buf
}

/// Serialises raw entries through the production writers `sa-builder` calls, which reach the body
/// by a different route than [`to_binary`] does. `bits` of `None` selects the uncompressed
/// packing.
pub fn to_file_bytes(sa: &[i64], sparseness: u8, bits: Option<usize>) -> Vec<u8> {
    let mut buf = Vec::new();
    match bits {
        None => super::dump_suffix_array(sa.to_vec(), sparseness, &mut buf).unwrap(),
        Some(bits) => super::dump_compressed_suffix_array(sa.to_vec(), sparseness, bits, &mut buf).unwrap()
    }
    buf
}
