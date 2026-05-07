// This entire file is mmap-only.
use std::{fs::File, path::Path};
use std::error::Error;
use memmap2::Mmap;
use text_compression::ReadBinaryMmap;

/// Owned, mmap-backed suffix array backend.
pub struct MmapBackedSA {
    pub mmap: Mmap,
    pub data_offset: usize,
    pub len: usize,
    pub bits_per_value: usize,
    pub sample_rate: u8,
}

impl super::SuffixArrayBackend for MmapBackedSA {
    type RangeIter<'a> = MmapSaRangeIter<'a>;

    fn len(&self) -> usize { self.len }
    fn bits_per_value(&self) -> usize { self.bits_per_value }
    fn sample_rate(&self) -> u8 { self.sample_rate }

    fn get(&self, index: usize) -> i64 {
        get_mmap(&self.mmap, self.data_offset, self.bits_per_value, index)
    }

    fn iter_range(&self, start: usize, end: usize) -> MmapSaRangeIter<'_> {
        MmapSaRangeIter::new(&self.mmap, self.data_offset, self.bits_per_value, start, end)
    }

    fn prefetch_sa_index(&self, index: usize) {
        let byte_offset = self.data_offset + (index * self.bits_per_value) / 8;
        if byte_offset < self.mmap.len() {
            let ptr: *const u8 = &self.mmap[byte_offset];
            prefetch::prefetch_read(ptr);
        }
    }

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

    fn prefetch_sa_range(&self, lo: usize, hi_exclusive: usize) {
        #[cfg(unix)]
        {
            let byte_lo = self.data_offset + (lo * self.bits_per_value) / 8;
            let byte_hi = self.data_offset + (hi_exclusive * self.bits_per_value).div_ceil(8);
            let len = byte_hi.saturating_sub(byte_lo);
            if len > 0 && byte_hi <= self.mmap.len() {
                let _ = self.mmap.advise_range(memmap2::Advice::WillNeed, byte_lo, len);
            }
        }
    }
}

impl ReadBinaryMmap for MmapBackedSA {
    fn read_binary_mmap(path: &Path) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
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
pub(super) fn read_u64_le(mmap: &Mmap, byte_offset: usize) -> u64 {
    let bytes: [u8; 8] = mmap[byte_offset..byte_offset + 8].try_into().unwrap();
    u64::from_le_bytes(bytes)
}

/// Streaming sequential iterator over a contiguous range of a compressed or uncompressed
/// mmap-backed suffix array.
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

/// Returns the suffix array value at the given index from a memory-mapped file.
pub(super) fn get_mmap(mmap: &Mmap, data_offset: usize, bits_per_value: usize, index: usize) -> i64 {
    if bits_per_value == 64 {
        let offset = data_offset + index * 8;
        let bytes: [u8; 8] = mmap[offset..offset + 8].try_into().unwrap();
        i64::from_le_bytes(bytes)
    } else {
        let mask: u64 = (1u64 << bits_per_value) - 1;
        let bit_offset = index * bits_per_value;
        let start_block = bit_offset / 64;
        let start_block_offset = bit_offset % 64;
        let block_byte_offset = data_offset + start_block * 8;
        let start_val = read_u64_le(mmap, block_byte_offset);
        if start_block_offset + bits_per_value <= 64 {
            ((start_val >> (64 - start_block_offset - bits_per_value)) & mask) as i64
        } else {
            let end_block_offset = (index + 1) * bits_per_value % 64;
            let end_val = read_u64_le(mmap, block_byte_offset + 8);
            let a = start_val << end_block_offset;
            let b = end_val >> (64 - end_block_offset);
            ((a | b) & mask) as i64
        }
    }
}
