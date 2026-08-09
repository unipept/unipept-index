// This entire file is mmap-only.
use std::{fs::File, path::Path};
use std::error::Error;
use memmap2::Mmap;
use text_compression::ReadBinaryMmap;

/// Owned, mmap-backed suffix array backend.
pub struct MmapBackedSA {
    pub mmap: Mmap,
    pub(crate) data_offset: usize,
    pub(crate) len: usize,
    pub(crate) bits_per_value: usize,
    pub(crate) sample_rate: u8,
}

impl super::SuffixArrayBackend for MmapBackedSA {
    type RangeIter<'a> = MmapSaRangeIter<'a>;

    fn len(&self) -> usize { self.len }
    fn bits_per_value(&self) -> usize { self.bits_per_value }
    fn sample_rate(&self) -> u8 { self.sample_rate }

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

    /// Warms the page cache by reading one byte from every page of the SA data.
    ///
    /// The three steps are all load-bearing:
    ///
    /// 1. `Advice::Sequential` tells the kernel to read far ahead, so the sweep below faults in
    ///    long runs instead of one page at a time.
    /// 2. Touching one byte per 4 KiB page is what actually forces the fault; the read must be
    ///    laundered through `black_box`, or the optimizer deletes a loop whose result is unused
    ///    and the warmup silently does nothing.
    /// 3. `Advice::Random` restores the steady-state pattern. Search probes the array in an order
    ///    the kernel cannot predict, so leaving readahead on would make every later miss drag in
    ///    neighbouring pages that will not be used.
    ///
    /// Only the SA data is swept, not the whole mapping, so the 10-byte header does not skew the
    /// chunking.
    fn touch_all_pages(&self) {
        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Sequential);

        let byte_len = (self.len * self.bits_per_value).div_ceil(8);
        let data = &self.mmap[self.data_offset..self.data_offset + byte_len];
        for chunk in data.chunks(4096) {
            std::hint::black_box(chunk[0]);
        }

        #[cfg(unix)]
        let _ = self.mmap.advise(memmap2::Advice::Random);
    }

    // `prefetch_sa_range` is intentionally not overridden here: the default no-op in
    // `SuffixArrayBackend` (array/mod.rs) applies. A per-query `MADV_WILLNEED` over the
    // k-mer SA range was measured to *regress* the mmap backend by -16.8% qps (5-mer,
    // 26-50aa, batch=16) even though pages are already resident. Average range size for a
    // 5-mer is ~54 KB vs ~2.7 KB for a 6-mer (37.2e9 SA entries / 20^5 vs 20^6), and every
    // rayon thread contends on the same VMA's mmap_lock to issue the advice — the penalty
    // scales with range size, matching the observed 5-mer-hurts / 6-mer-roughly-even split.
    // On an index whose pages are already touched (see `touch_all_pages`), the syscall buys
    // nothing and the lock contention is pure cost.
}

impl ReadBinaryMmap for MmapBackedSA {
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

        Ok(MmapBackedSA { mmap, data_offset: 10, len: amount_of_items, bits_per_value, sample_rate })
    }
}

// ── range iterator ────────────────────────────────────────────────────────────

/// Reads a u64 value in little-endian byte order from the given mmap at the given byte offset.
#[inline]
pub(super) fn read_u64_le(mmap: &Mmap, byte_offset: usize) -> u64 {
    let bytes: [u8; 8] = mmap[byte_offset..byte_offset + 8].try_into().unwrap();
    u64::from_le_bytes(bytes)
}

/// Streaming sequential iterator over a contiguous range of a compressed or uncompressed
/// mmap-backed suffix array.
/// Sequential reader over a range of SA entries, decoding straight from the mapping.
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
    remaining: usize,
}

impl<'a> MmapSaRangeIter<'a> {
    pub fn new(
        mmap: &'a Mmap,
        data_offset: usize,
        bits_per_value: usize,
        start: usize,
        end: usize,
    ) -> Self {
        let remaining = end.saturating_sub(start);
        if remaining == 0 {
            return Self {
                mmap, data_offset, bits_per_value,
                mask: 0, current_word: 0, next_word: 0,
                block_idx: 0, bit_off: 0, remaining: 0,
            };
        }

        let mask = if bits_per_value == 64 { u64::MAX } else { (1u64 << bits_per_value) - 1 };

        let bit_pos   = start * bits_per_value;
        let block_idx = bit_pos / 64;
        let bit_off   = bit_pos % 64;

        let current_word = read_u64_le(mmap, data_offset + block_idx * 8);
        let next_off     = data_offset + (block_idx + 1) * 8;
        let next_word    = if next_off + 8 <= mmap.len() { read_u64_le(mmap, next_off) } else { 0 };

        Self { mmap, data_offset, bits_per_value, mask, current_word, next_word, block_idx, bit_off, remaining }
    }
}

impl Iterator for MmapSaRangeIter<'_> {
    type Item = i64;

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) { (self.remaining, Some(self.remaining)) }

    #[inline]
    fn next(&mut self) -> Option<i64> {
        if self.remaining == 0 { return None; }
        self.remaining -= 1;

        let val = if self.bit_off + self.bits_per_value <= 64 {
            (self.current_word >> (64 - self.bit_off - self.bits_per_value)) & self.mask
        } else {
            let end_off = (self.bit_off + self.bits_per_value) % 64;
            ((self.current_word << end_off) | (self.next_word >> (64 - end_off))) & self.mask
        };

        self.bit_off += self.bits_per_value;
        if self.bit_off >= 64 {
            self.bit_off   -= 64;
            self.block_idx += 1;
            self.current_word = self.next_word;
            let next_off = self.data_offset + (self.block_idx + 1) * 8;
            self.next_word = if next_off + 8 <= self.mmap.len() {
                read_u64_le(self.mmap, next_off)
            } else { 0 };
        }

        Some(val as i64)
    }
}

impl ExactSizeIterator for MmapSaRangeIter<'_> {}

