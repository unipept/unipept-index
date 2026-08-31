//! Fixtures for the preloaded mapping tests, on top of the backend-agnostic ones in
//! [`suffix_to_protein_index::test_utils`](crate::suffix_to_protein_index::test_utils).

use binary_traits::{ReadBinary, WriteBinary};

use super::InMemorySuffixToProteinMapping;
use crate::suffix_to_protein_index::test_utils::{assert_sample_lookups, to_binary};

/// Asserts the tag byte `mapping` writes, and that `read_binary` picks the variant it names.
/// Expects a mapping built from [`sample_text`](crate::suffix_to_protein_index::test_utils::sample_text).
pub fn assert_dump_and_load(mapping: impl WriteBinary, expected_tag: u8) {
    let buf = to_binary(mapping);
    assert_eq!(buf[0], expected_tag);
    let loaded = InMemorySuffixToProteinMapping::read_binary(&mut std::io::Cursor::new(buf)).unwrap();
    assert_sample_lookups(&loaded);
}
