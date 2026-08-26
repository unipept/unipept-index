//! Fixtures shared by the two protein-text backends' tests.
//!
//! Both backends decode the same file and must answer identically, so what is asserted about one
//! belongs here rather than in either module — a helper that lives on one side tends to grow a
//! second, slightly different copy on the other.

use crate::ProteinTextBackend;

/// Walks `prefetch_at` over the whole text and well past its end, then re-reads every residue.
///
/// Both backends compute the address they hint by hand — the preloaded one scales the residue
/// index to a `u64` word, the mmap one to a byte offset into the mapping — so a wrong bound shows
/// up as an out-of-range index rather than a wrong answer, and a hint may never disturb the read
/// it precedes.
///
/// The walk deliberately runs to twice the length: [`ProteinTextBackend::prefetch_at`] promises
/// that an out-of-range index is ignored rather than panicking, because the search path hints a
/// lookahead position that may sit past the end and checking it at every call site would cost
/// more than the hint saves. That promise had no test on either backend.
pub fn assert_prefetch_is_harmless(text: &impl ProteinTextBackend, expected: &str) {
    for index in 0..expected.len() * 2 {
        text.prefetch_at(index);
    }

    assert_eq!(text.len(), expected.len(), "length differs after the hints");
    for (index, c) in expected.chars().enumerate() {
        assert_eq!(text.get(index), c as u8, "residue {index} differs after the hints");
    }
}
