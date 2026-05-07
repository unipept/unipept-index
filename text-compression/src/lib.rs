use std::io::Write;

pub mod traits;
pub use traits::{WriteBinary, ReadBinary, ReadBinaryMmap};

pub mod preloaded;
#[cfg(feature = "mmap")]
pub mod mmap;

pub use preloaded::{InMemoryProteinText, dump_compressed_text, load_compressed_text};
#[cfg(feature = "mmap")]
pub use mmap::MmapProteinText;

/// Type alias — resolves to the single active backend for this build.
#[cfg(feature = "mmap")]
pub type ProteinText = MmapProteinText;
#[cfg(not(feature = "mmap"))]
pub type ProteinText = InMemoryProteinText;

/// Returns the number of bytes the BitArray data occupies for a given text length at 5 bits/value.
pub fn bit_array_byte_size(text_length: usize) -> usize {
    let extra = if (text_length * 5).is_multiple_of(64) { 0 } else { 1 };
    (text_length * 5 / 64 + extra) * 8
}

// ── ProteinTextSlice ──────────────────────────────────────────────────────────

pub struct ProteinTextSlice<'a> {
    text: &'a ProteinText,
    start: usize,
    end: usize,
}

impl<'a> ProteinTextSlice<'a> {
    pub fn new(text: &'a ProteinText, start: usize, end: usize) -> Self {
        Self { text, start, end }
    }

    pub fn get(&self, index: usize) -> u8 { self.text.get(self.start + index) }
    pub fn len(&self) -> usize { self.end - self.start }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    #[inline]
    pub fn equals_slice(&self, other: &[u8], equate_il: bool) -> bool {
        if equate_il {
            other.iter().zip(self.iter()).all(|(&s, t)| {
                s == t || (s == b'I' && t == b'L') || (s == b'L' && t == b'I')
            })
        } else {
            other.iter().zip(self.iter()).all(|(&s, t)| s == t)
        }
    }

    pub fn check_il_locations(&self, skip: usize, il_locations: &[usize], search_string: &[u8]) -> bool {
        for &il_location in il_locations {
            let index = il_location - skip;
            if search_string[index] != self.get(index) { return false; }
        }
        true
    }

    pub fn iter(&self) -> ProteinTextSliceIterator<'_> {
        ProteinTextSliceIterator { text_slice: self, index: 0 }
    }
}

// ── ProteinTextIterator ───────────────────────────────────────────────────────

pub struct ProteinTextIterator<'a> {
    pub protein_text: &'a ProteinText,
    pub index: usize,
}

pub struct ProteinTextSliceIterator<'a> {
    text_slice: &'a ProteinTextSlice<'a>,
    index: usize,
}

impl Iterator for ProteinTextSliceIterator<'_> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.text_slice.len() { return None; }
        self.index += 1;
        Some(self.text_slice.get(self.index - 1))
    }
}

impl Iterator for ProteinTextIterator<'_> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.protein_text.len() { return None; }
        self.index += 1;
        Some(self.protein_text.get(self.index - 1))
    }
}
