//! Fixtures for the mmap suffix-array tests, on top of the backend-agnostic ones in
//! [`array::test_utils`](crate::array::test_utils).
//!
//! There is nothing to map until something has written a file, so these start from
//! [`to_file_bytes`](crate::array::test_utils::to_file_bytes) — the writers `sa-builder` calls.

use std::io::Write;

use tempfile::NamedTempFile;
use text_compression::ReadBinaryMmap;

use super::MmapBackedSA;
use crate::array::{SuffixArrayBackend, test_utils::to_file_bytes};

/// Writes `buf` to a temporary file. The caller keeps the returned guard alive for as long as the
/// file stays mapped.
pub fn write_to_tempfile(buf: &[u8]) -> NamedTempFile {
    let mut tmp = NamedTempFile::new().unwrap();
    tmp.write_all(buf).unwrap();
    tmp.flush().unwrap();
    tmp
}

/// Writes a suffix array through the production writer and maps the file back in. `bits` of `None`
/// selects the uncompressed packing. The returned file guard must outlive the mapping.
pub fn write_and_map(sa: &[i64], sparseness: u8, bits: Option<usize>) -> (MmapBackedSA, NamedTempFile) {
    let tmp = write_to_tempfile(&to_file_bytes(sa, sparseness, bits));
    let mapped = MmapBackedSA::read_binary_mmap(tmp.path()).unwrap();
    (mapped, tmp)
}

/// Walks both hint methods over the whole array and past its end, then re-reads every entry. Both
/// compute byte offsets by hand, so a wrong bound is an out-of-range index on the mapping rather
/// than a wrong answer, and neither may disturb the entries it precedes.
pub fn assert_hints_are_harmless(sa: &[i64], sparseness: u8, bits: Option<usize>) {
    let (mapped, _tmp) = write_and_map(sa, sparseness, bits);
    mapped.touch_all_pages();
    for index in 0..sa.len() * 2 {
        mapped.prefetch_sa_index(index);
    }
    for (index, &expected) in sa.iter().enumerate() {
        assert_eq!(mapped.get(index), expected, "entry {index} differs after the hints");
    }
}
