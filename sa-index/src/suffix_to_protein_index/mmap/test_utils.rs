//! Fixtures for the mmap mapping tests, on top of the backend-agnostic ones in
//! [`suffix_to_protein_index::test_utils`](crate::suffix_to_protein_index::test_utils).
//!
//! There is nothing to map until something has written a file, so every test here starts from a
//! preloaded mapping serialised to a temporary file — the same path `sa-builder` writes and
//! `sa-server` reads.

use tempfile::NamedTempFile;
use text_compression::ProteinTextBackend;

use super::MmapBackedSuffixToProteinMapping;
use crate::{
    ReadBinaryMmap, WriteBinary,
    suffix_to_protein_index::{
        SuffixToProteinMappingBackend,
        test_utils::{assert_sample_lookups, sample_text, to_binary}
    }
};

/// Writes a serialised mapping to a temporary file. The caller keeps the returned guard alive for
/// as long as the file stays mapped.
pub fn write_to_tempfile(buf: &[u8]) -> NamedTempFile {
    use std::io::Write;

    let mut tmp = NamedTempFile::new().unwrap();
    tmp.write_all(buf).unwrap();
    tmp.flush().unwrap();
    tmp
}

/// Serialises a preloaded mapping and maps the file back in. The returned file guard must outlive
/// the mapping.
pub fn write_and_map(mapping: impl WriteBinary) -> (MmapBackedSuffixToProteinMapping, NamedTempFile) {
    let tmp = write_to_tempfile(&to_binary(mapping));
    let mapped = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();
    (mapped, tmp)
}

/// Asserts the tag byte `mapping` writes, and that `read_binary_mmap` picks the variant it names.
/// Expects a mapping built from [`sample_text`].
pub fn assert_load_mmap(mapping: impl WriteBinary, expected_tag: u8) {
    let buf = to_binary(mapping);
    assert_eq!(buf[0], expected_tag);
    let tmp = write_to_tempfile(&buf);
    let loaded = MmapBackedSuffixToProteinMapping::read_binary_mmap(tmp.path()).unwrap();
    assert_sample_lookups(&loaded);
}

/// Walks the two hint methods over the whole mapping and past its end. Both compute file offsets
/// by hand, so a wrong bound is an out-of-range index on the mapping rather than a wrong answer,
/// and neither may disturb the lookups it precedes. Expects a mapping built from [`sample_text`].
pub fn assert_hints_are_harmless(mapping: impl WriteBinary) {
    let (loaded, _tmp) = write_and_map(mapping);
    loaded.touch_all_pages();
    for suffix in 0..sample_text().len() as i64 * 2 {
        loaded.prefetch_for_suffix(suffix);
    }
    assert_sample_lookups(&loaded);
}
