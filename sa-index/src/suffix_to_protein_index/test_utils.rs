//! Fixtures shared by the mapping tests, of both backends.
//!
//! All three representations answer the same questions about the same text, so the text and the
//! expected answers live here once instead of in every variant's test module. Fixtures that need
//! a mapped file are in [`mmap::test_utils`](super::mmap::test_utils) instead.

use binary_traits::WriteBinary;
use protein_metadata::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use protein_text::{InMemoryProteinText, ProteinTextBackend};

use super::SuffixToProteinMappingBackend;
use crate::Nullable;

/// Three proteins — `ACG`, `CG` and `AAA` — separated and terminated, 11 bytes in total. Positions
/// 3, 6 and 10 hold a separator or the terminator and so belong to no protein.
pub fn sample_text() -> InMemoryProteinText {
    let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
    text.push(TERMINATION_CHARACTER as char);
    InMemoryProteinText::from_string(&text)
}

/// `protein_count` proteins of `protein_len` residues each, separated and terminated. Long enough
/// to span several of the bitvec rank structure's 512-bit superblocks, which the 11-byte
/// [`sample_text`] does not reach.
pub fn many_proteins_text(protein_count: usize, protein_len: usize) -> InMemoryProteinText {
    let protein = "ACGKLMNPQR"[..protein_len].to_string();
    let mut text = vec![protein; protein_count].join(&format!("{}", SEPARATION_CHARACTER as char));
    text.push(TERMINATION_CHARACTER as char);
    InMemoryProteinText::from_string(&text)
}

/// Serialises a mapping to a buffer, whose first byte is the type tag the loaders dispatch on.
pub fn to_binary(mapping: impl WriteBinary) -> Vec<u8> {
    let mut buf = Vec::new();
    mapping.write_binary(&mut buf).unwrap();
    buf
}

/// Asserts every answer [`sample_text`] fixes: one position inside each of the three proteins, and
/// all three positions that belong to no protein — both separators and the terminator.
pub fn assert_sample_lookups(mapping: &impl SuffixToProteinMappingBackend) {
    assert_eq!(mapping.suffix_to_protein(0), 0);
    assert_eq!(mapping.suffix_to_protein(5), 1);
    assert_eq!(mapping.suffix_to_protein(7), 2);
    assert_eq!(mapping.suffix_to_protein(3), u32::NULL);
    assert_eq!(mapping.suffix_to_protein(6), u32::NULL);
    assert_eq!(mapping.suffix_to_protein(10), u32::NULL);
}

/// Walks `prefetch_for_suffix` over the whole mapping and well past its end, then asserts the
/// lookups still hold. A backend that implements the hint computes the index it touches by hand,
/// so a wrong bound is an out-of-range panic rather than a wrong answer, and the hint may never
/// disturb the lookup it precedes. Expects a mapping built from [`sample_text`].
pub fn assert_prefetch_is_harmless(mapping: &impl SuffixToProteinMappingBackend) {
    for suffix in 0..sample_text().len() as i64 * 2 {
        mapping.prefetch_for_suffix(suffix);
    }
    assert_sample_lookups(mapping);
}

/// Asserts two mappings answer identically for every position of a `text_len` byte text.
pub fn assert_agree(
    expected: &impl SuffixToProteinMappingBackend,
    actual: &impl SuffixToProteinMappingBackend,
    text_len: usize
) {
    for suffix in 0..text_len as i64 {
        let (expected, actual) = (expected.suffix_to_protein(suffix), actual.suffix_to_protein(suffix));
        assert_eq!(expected, actual, "mismatch at suffix {}", suffix);
    }
}
