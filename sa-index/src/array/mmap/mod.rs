use std::{error::Error, fs::File, path::Path};

use memmap2::Mmap;
use text_compression::{LoadIndex, ReadBinaryMmap};

#[cfg(test)]
pub(super) mod test_utils;

/// Suffix array read straight out of a memory mapping, in either packing.
///
/// Holds no entries of its own: every lookup decodes from `mmap`, so it may fault a page in. The
/// four fields after it are what [`read_binary_mmap`](ReadBinaryMmap::read_binary_mmap) took from
/// the file header.
pub struct MmapBackedSA {
    pub mmap: Mmap,
    /// Where the entries start — 10 bytes in, past the header.
    pub(crate) data_offset: usize,
    pub(crate) len: usize,
    pub(crate) bits_per_value: usize,
    pub(crate) sample_rate: u8
}

impl super::SuffixArrayBackend for MmapBackedSA {
    type RangeIter<'a> = MmapSaRangeIter<'a>;

    fn len(&self) -> usize {
        self.len
    }
    fn bits_per_value(&self) -> usize {
        self.bits_per_value
    }
    fn sample_rate(&self) -> u8 {
        self.sample_rate
    }

    #[inline]
    fn get(&self, index: usize) -> i64 {
        if self.bits_per_value == 64 {
            let offset = self.data_offset + index * 8;
            let bytes: [u8; 8] = self.mmap[offset..offset + 8].try_into().unwrap();
            i64::from_le_bytes(bytes)
        } else {
            let mask: u64 = (1u64 << self.bits_per_value) - 1;
            let bit_offset = index * self.bits_per_value;
            let start_block = bit_offset / 64;
            let start_block_offset = bit_offset % 64;
            let block_byte_offset = self.data_offset + start_block * 8;
            let start_val = read_u64_le(&self.mmap, block_byte_offset);
            if start_block_offset + self.bits_per_value <= 64 {
                ((start_val >> (64 - start_block_offset - self.bits_per_value)) & mask) as i64
            } else {
                let end_block_offset = (index + 1) * self.bits_per_value % 64;
                let end_val = read_u64_le(&self.mmap, block_byte_offset + 8);
                (((start_val << end_block_offset) | (end_val >> (64 - end_block_offset))) & mask) as i64
            }
        }
    }

    fn iter_range(&self, start: usize, end: usize) -> MmapSaRangeIter<'_> {
        MmapSaRangeIter::new(&self.mmap, self.data_offset, self.bits_per_value, start, end)
    }

    #[inline]
    fn prefetch_sa_index(&self, index: usize) {
        let byte_offset = self.data_offset + (index * self.bits_per_value) / 8;
        if byte_offset < self.mmap.len() {
            let ptr: *const u8 = &self.mmap[byte_offset];
            prefetch::prefetch_read(ptr);
        }
    }

    /// Warms the page cache over the SA data. See `text_compression::mmap::touch_all_pages` for
    /// why the sweep is shaped the way it is.
    fn touch_all_pages(&self) -> u64 {
        // Only the SA data, not the whole mapping: the 10-byte header would otherwise skew the
        // range for no benefit.
        let byte_len = (self.len * self.bits_per_value).div_ceil(8);
        text_compression::mmap::touch_all_pages(&self.mmap, self.data_offset..self.data_offset + byte_len)
    }

    // Do not add a `MADV_WILLNEED` over the SA range before scanning it. Tried twice, removed
    // twice. Resident it is pure cost — -16.8% qps with a 5-mer table, -3.7% with a 6-mer —
    // because every rayon thread contends on the same VMA's `mmap_lock` and the penalty scales
    // with range size (~54 KB per 5-mer range vs ~2.7 KB per 6-mer). Under a memory ceiling the
    // advice does land: major faults fall 23-25% at both a 167 GB and a 112 GB cap. But the
    // throughput that buys decays as threads rise (+12.0% at the core count, ~0% at 48-96),
    // because oversubscription and readahead are substitutes — with ~55 faults already in flight
    // across 96 threads, removing a quarter of them changes little — and it does not let the
    // thread count come down. If it is ever retried, fix the syscall count first with
    // `process_madvise` (Linux 5.10+).
}

impl ReadBinaryMmap for MmapBackedSA {
    /// Maps an SA file, checking that it is long enough for both the header and the entries the
    /// header declares. That check is what lets every lookup below index the mapping without
    /// bounds-checking against the file length.
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        // SAFETY: see the note in `text_compression::mmap` — an index file is written once by
        // sa-builder and is read-only for the lifetime of the process, so the mapping cannot be
        // truncated or written underneath us.
        let mmap = unsafe { Mmap::map(&file)? };

        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random)?;

        if mmap.len() < 10 {
            return Err("The binary file is too small to contain the SA header".into());
        }

        let bits_per_value = mmap[0] as usize;
        super::check_bits_per_value(bits_per_value)?;
        let sample_rate = mmap[1];
        let amount_of_items = u64::from_le_bytes(mmap[2..10].try_into()?) as usize;

        let header_bytes = 10usize;
        let total_bits = amount_of_items
            .checked_mul(bits_per_value)
            .ok_or("The SA header declares too many items or bits per value")?;
        let data_bytes = total_bits.div_ceil(8);

        if mmap.len() < header_bytes + data_bytes {
            return Err("The binary file is too small to contain the SA data".into());
        }

        Ok(MmapBackedSA {
            mmap,
            data_offset: 10,
            len: amount_of_items,
            bits_per_value,
            sample_rate
        })
    }
}

impl LoadIndex for MmapBackedSA {
    fn load(path: &Path) -> Result<Self, Box<dyn Error>> {
        Self::read_binary_mmap(path)
    }
}

// ── range iterator ────────────────────────────────────────────────────────────

/// Reads a u64 value in little-endian byte order from the given mmap at the given byte offset.
#[inline]
pub(super) fn read_u64_le(mmap: &Mmap, byte_offset: usize) -> u64 {
    let bytes: [u8; 8] = mmap[byte_offset..byte_offset + 8].try_into().unwrap();
    u64::from_le_bytes(bytes)
}

/// Sequential reader over a range of SA entries, decoding straight from the mapping. Handles
/// either packing.
///
/// # Bit layout
///
/// Entries are packed **most-significant-bit first within each little-endian `u64` word**, and an
/// entry may straddle a word boundary. That combination is unusual enough to be worth stating
/// plainly, because the two conventions pull in opposite directions: the *bytes* of a word are
/// read little-endian, but the *values* inside it are laid out from the top bit down. Hence
/// `word >> (64 - bit_offset - bits)` rather than the `word >> bit_offset` a purely
/// little-endian packing would use.
///
/// This must match `bitarray`'s packing exactly — the file is written through `DynBitArray` —
/// and `MmapBackedSA::get`, which decodes the same layout for random access.
///
/// # Why this exists alongside `get`
///
/// It caches the current and next words, so consecutive entries sharing a word cost no reload and
/// a straddling entry already has its second word to hand. Scanning a candidate range through
/// `get` would re-derive and re-load both per entry.
pub struct MmapSaRangeIter<'a> {
    mmap: &'a Mmap,
    data_offset: usize,
    bits_per_value: usize,
    mask: u64,
    current_word: u64,
    next_word: u64,
    block_idx: usize,
    bit_off: usize,
    remaining: usize
}

impl<'a> MmapSaRangeIter<'a> {
    pub fn new(mmap: &'a Mmap, data_offset: usize, bits_per_value: usize, start: usize, end: usize) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                mmap,
                data_offset,
                bits_per_value,
                mask: 0,
                current_word: 0,
                next_word: 0,
                block_idx: 0,
                bit_off: 0,
                remaining: 0
            };
        }

        let mask = if bits_per_value == 64 { u64::MAX } else { (1u64 << bits_per_value) - 1 };

        let bit_pos = start * bits_per_value;
        let block_idx = bit_pos / 64;
        let bit_off = bit_pos % 64;

        let current_word = read_u64_le(mmap, data_offset + block_idx * 8);
        let next_off = data_offset + (block_idx + 1) * 8;
        let next_word = if next_off + 8 <= mmap.len() { read_u64_le(mmap, next_off) } else { 0 };

        Self {
            mmap,
            data_offset,
            bits_per_value,
            mask,
            current_word,
            next_word,
            block_idx,
            bit_off,
            remaining
        }
    }
}

impl Iterator for MmapSaRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;

        let val = if self.bit_off + self.bits_per_value <= 64 {
            (self.current_word >> (64 - self.bit_off - self.bits_per_value)) & self.mask
        } else {
            let end_off = (self.bit_off + self.bits_per_value) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & self.mask
        };

        self.bit_off += self.bits_per_value;
        if self.bit_off >= 64 {
            self.bit_off -= 64;
            self.block_idx += 1;
            self.current_word = self.next_word;
            let next_off = self.data_offset + (self.block_idx + 1) * 8;
            self.next_word = if next_off + 8 <= self.mmap.len() { read_u64_le(self.mmap, next_off) } else { 0 };
        }

        Some(val as i64)
    }
}

impl ExactSizeIterator for MmapSaRangeIter<'_> {}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::{
        test_utils::{assert_hints_are_harmless, write_and_map, write_to_tempfile},
        *
    };
    use crate::array::{
        OriginalSA, SuffixArrayBackend, dump_compressed_suffix_array,
        test_utils::{fit_to_width, owned_compressed, sample_sa}
    };

    /// The mmap backend is the production storage layer. It must agree with the owned-memory
    /// backends that wrote the file, entry for entry.
    #[test]
    fn matches_the_uncompressed_backend() {
        let sa = sample_sa(500);
        let (mapped, _tmp) = write_and_map(&sa, 1, None);
        let owned = OriginalSA(sa.clone(), 1);

        assert_eq!(mapped.len(), owned.len());
        assert_eq!(mapped.bits_per_value(), 64);
        assert_eq!(mapped.sample_rate(), 1);
        for i in 0..sa.len() {
            assert_eq!(mapped.get(i), owned.get(i), "entry {i} differs");
        }
    }

    /// Same, for every compressed width the builder can choose. This is the case that exercises
    /// entries straddling a `u64` boundary, which is where a packing bug would hide.
    #[test]
    fn matches_the_compressed_backend_at_every_width() {
        for bits in [8usize, 13, 17, 28, 29, 31, 32, 33, 40, 63] {
            let sa = fit_to_width(&sample_sa(400), bits);

            let (mapped, _tmp) = write_and_map(&sa, 3, Some(bits));
            let owned = owned_compressed(&sa, 3, bits);

            assert_eq!(mapped.len(), sa.len(), "length differs at {bits} bits");
            assert_eq!(mapped.bits_per_value(), bits);
            assert_eq!(mapped.sample_rate(), 3);
            for (i, &expected) in sa.iter().enumerate() {
                assert_eq!(mapped.get(i), owned.get(i), "entry {i} differs at {bits} bits");
                assert_eq!(mapped.get(i), expected, "entry {i} lost at {bits} bits");
            }
        }
    }

    /// `iter_range` duplicates `get`'s unpacking with a cached word pair, so the two must agree
    /// at every start offset and length — including ranges that begin mid-word.
    #[test]
    fn iter_range_agrees_with_get() {
        for bits in [29usize, 32, 40] {
            let sa = fit_to_width(&sample_sa(200), bits);
            let (mapped, _tmp) = write_and_map(&sa, 1, Some(bits));

            for start in [0usize, 1, 7, 12, 63, 64, 65, 130] {
                for end in [start, start + 1, start + 13, start + 64, 200] {
                    if end > 200 || end < start {
                        continue;
                    }
                    let by_iter: Vec<i64> = mapped.iter_range(start, end).collect();
                    let by_get: Vec<i64> = (start..end).map(|i| mapped.get(i)).collect();
                    assert_eq!(by_iter, by_get, "iter_range({start}, {end}) at {bits} bits");
                    assert_eq!(mapped.iter_range(start, end).len(), end - start);
                }
            }
        }
    }

    #[test]
    fn truncated_files_error_rather_than_panicking() {
        let sa = sample_sa(100);
        let mut buf: Vec<u8> = Vec::new();
        dump_compressed_suffix_array(sa.clone(), 1, 29, &mut buf).unwrap();

        // The writer pads the body out to whole `u64` words, but the reader only requires
        // `ceil(items * bits / 8)` bytes, so the last few bytes of the file are slack that a
        // truncation can legitimately remove. Sweep up to the reader's actual requirement.
        let required = 10 + (sa.len() * 29).div_ceil(8);
        assert!(required <= buf.len());

        for cut in 0..required {
            let tmp = write_to_tempfile(&buf[..cut]);
            let err = MmapBackedSA::read_binary_mmap(tmp.path())
                .err()
                .unwrap_or_else(|| panic!("{cut} of {required} required bytes should not load"));
            assert!(err.to_string().contains("too small"), "unexpected error at {cut}: {err}");
        }
    }

    /// A header claiming more entries than the file holds must be rejected, not trusted.
    #[test]
    fn overlong_declared_length_is_rejected() {
        let sa = sample_sa(50);
        let mut buf: Vec<u8> = Vec::new();
        dump_compressed_suffix_array(sa, 1, 29, &mut buf).unwrap();
        buf[2..10].copy_from_slice(&1_000_000_u64.to_le_bytes());

        let tmp = write_to_tempfile(&buf);
        assert!(MmapBackedSA::read_binary_mmap(tmp.path()).is_err());
    }

    /// `prefetch_sa_index` and `touch_all_pages` are hints nothing else in the suite calls, yet
    /// both index the mapping by offsets they compute themselves.
    #[test]
    fn hints_stay_within_the_mapping() {
        assert_hints_are_harmless(&sample_sa(300), 1, None);
        for bits in [8usize, 29, 40] {
            assert_hints_are_harmless(&fit_to_width(&sample_sa(300), bits), 2, Some(bits));
        }
    }
}
