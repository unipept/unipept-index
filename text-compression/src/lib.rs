pub mod traits;
pub use traits::{WriteBinary, ReadBinary, ReadBinaryMmap};

pub mod preloaded;
#[cfg(feature = "mmap")]
pub mod mmap;

pub use preloaded::{InMemoryProteinText, dump_compressed_text, load_compressed_text};
#[cfg(feature = "mmap")]
pub use mmap::MmapBackedProteinText;

/// Decode table shared by both backends: 5-bit index → ASCII amino acid byte.
pub const BIT5_TO_CHAR: &[u8; 27] = b"ABCDEFGHIKLMNOPQRSTUVWXYZ-$";

/// Type alias — resolves to the single active backend for this build.
#[cfg(feature = "mmap")]
pub type ProteinText = MmapBackedProteinText;
#[cfg(not(feature = "mmap"))]
pub type ProteinText = InMemoryProteinText;

/// Returns the number of bytes the BitArray data occupies for a given text length at 5 bits/value.
pub fn bit_array_byte_size(text_length: usize) -> usize {
    let extra = if (text_length * 5).is_multiple_of(64) { 0 } else { 1 };
    (text_length * 5 / 64 + extra) * 8
}

// ── ProteinTextBackend ─────────────────────────────────────────────────────────

/// Full access interface for protein text backends.
pub trait ProteinTextBackend {
    fn get(&self, index: usize) -> u8;
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool { self.len() == 0 }

    fn prefetch_at(&self, index: usize);

    fn iter(&self) -> ProteinTextIterator<'_, Self> where Self: Sized {
        ProteinTextIterator { protein_text: self, index: 0 }
    }

    fn slice(&self, start: usize, end: usize) -> ProteinTextSlice<'_, Self> where Self: Sized {
        ProteinTextSlice::new(self, start, end)
    }
}

// ── ProteinTextSlice ──────────────────────────────────────────────────────────

pub struct ProteinTextSlice<'a, T: ProteinTextBackend> {
    text: &'a T,
    start: usize,
    end: usize,
}

impl<'a, T: ProteinTextBackend> ProteinTextSlice<'a, T> {
    pub fn new(text: &'a T, start: usize, end: usize) -> Self {
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

    pub fn iter(&self) -> ProteinTextSliceIterator<'_, T> {
        ProteinTextSliceIterator { text: self.text, pos: self.start, end: self.end }
    }
}

// ── ProteinTextIterator ───────────────────────────────────────────────────────

pub struct ProteinTextIterator<'a, T: ProteinTextBackend> {
    pub(crate) protein_text: &'a T,
    pub(crate) index: usize,
}

pub struct ProteinTextSliceIterator<'a, T: ProteinTextBackend> {
    text: &'a T,
    pos: usize,
    end: usize,
}

impl<T: ProteinTextBackend> Iterator for ProteinTextSliceIterator<'_, T> {
    type Item = u8;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.end { return None; }
        let c = self.text.get(self.pos);
        self.pos += 1;
        Some(c)
    }
}

impl<T: ProteinTextBackend> Iterator for ProteinTextIterator<'_, T> {
    type Item = u8;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.protein_text.len() { return None; }
        self.index += 1;
        Some(self.protein_text.get(self.index - 1))
    }
}
