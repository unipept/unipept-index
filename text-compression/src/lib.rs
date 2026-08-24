#![warn(missing_docs)]
//! The concatenated protein text, packed at 5 bits per residue.
//!
//! Every protein sequence in the database is concatenated into one text, separated by `-` and
//! terminated by `$`. The suffix array indexes positions in *this* text, so essentially every
//! search operation ends up reading it: comparing a candidate suffix against the query touches
//! one residue per character compared. It is the hottest data structure in the index.
//!
//! The alphabet is 25 amino-acid letters plus the two delimiters, so a residue needs 5 bits
//! rather than 8. At UniProt scale that is the difference between roughly 290 MB and 190 MB, and
//! the smaller footprint is worth more than the unpacking costs.
//!
//! # Two backends
//!
//! Both are always compiled; callers pick by naming a type, and everything here that touches the
//! text is generic over [`ProteinTextBackend`].
//!
//! * [`preloaded::InMemoryProteinText`] — the text decompressed into owned RAM. Faster per
//!   access, but the process pays the full resident size.
//! * [`mmap::MmapBackedProteinText`] — decoded straight out of a memory mapping, so the kernel
//!   decides what stays resident.
//!
//! Storing this structure differently from the rest of the index is worth a knob of its own
//! (`sa-server`'s `preloaded-text`): the text is the hottest structure in the index while the
//! protein metadata sharing its file is the biggest.
//!
//! `preloaded` owns the `WriteBinary` implementation that produces the file both backends read,
//! which is why `sa-builder` needs only that half.

pub use binary_traits::{LoadIndex, ReadBinary, ReadBinaryMmap, WriteBinary, load_owned};

pub mod mmap;
pub mod preloaded;
#[cfg(test)]
mod test_utils;

pub use mmap::MmapBackedProteinText;
pub use preloaded::InMemoryProteinText;

/// Decode table shared by both backends: 5-bit index → ASCII amino acid byte.
///
/// The inverse (ASCII → 5-bit) is built from this table at load time, so this array is the single
/// definition of the alphabet and its encoding. Changing the order changes the on-disk format.
///
/// Note it has 27 entries but is indexed by a 5-bit value, so codes 27..=31 are out of bounds. No
/// encoder can emit them, but a corrupt or truncated index file can, which would panic here.
/// Tracked as a known issue; padding the table to 32 entries would remove both the panic and the
/// bounds check from the hot path.
pub const BIT5_TO_CHAR: &[u8; 27] = b"ABCDEFGHIKLMNOPQRSTUVWXYZ-$";

/// Returns the number of bytes the BitArray data occupies for a given text length at 5 bits/value.
///
/// Rounded up to whole `u64` words, matching how `bitarray` allocates. Both mmap readers use this
/// to bounds-check a declared text length against the actual file size before mapping it, so it
/// must stay in step with `BitArray::<5>::with_capacity`.
pub fn bit_array_byte_size(text_length: usize) -> usize {
    let extra = if (text_length * 5).is_multiple_of(64) { 0 } else { 1 };
    (text_length * 5 / 64 + extra) * 8
}

// ── ProteinTextBackend ─────────────────────────────────────────────────────────

/// Full access interface for protein text backends.
pub trait ProteinTextBackend {
    /// Returns the ASCII residue at `index`, decoded from its 5-bit code.
    ///
    /// # Panics
    ///
    /// If `index` is past the end of the text.
    fn get(&self, index: usize) -> u8;

    /// Number of residues in the text, including separators and the terminator.
    fn len(&self) -> usize;

    /// Whether the text is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Reads at least one byte from every OS page backing this text, faulting it into the page
    /// cache.
    ///
    /// Returns the number of bytes swept, so a caller can report a bandwidth; see
    /// [`crate::mmap::touch_all_pages`].
    ///
    /// Default returns 0: an owned text is already resident, so there is nothing to fault.
    /// [`crate::MmapBackedProteinText`] overrides it. This lives on the *text* rather than only on
    /// the protein store because the two are independent axes — a build may own the metadata while
    /// the text stays mapped, and in that build nothing else can reach these pages.
    #[inline]
    fn touch_all_pages(&self) -> u64 {
        0
    }

    /// Issues a prefetch hint for the storage holding `index`, without reading it.
    ///
    /// Out-of-range indices are ignored rather than panicking: callers prefetch a lookahead
    /// position that may sit past the end, and checking it at every call site would cost more
    /// than the hint saves.
    fn prefetch_at(&self, index: usize);

    /// Iterates the whole text, one residue at a time.
    fn iter(&self) -> ProteinTextIterator<'_, Self>
    where
        Self: Sized
    {
        ProteinTextIterator { protein_text: self, index: 0 }
    }

    /// Borrows `start..end` (half-open) as a comparable slice.
    fn slice(&self, start: usize, end: usize) -> ProteinTextSlice<'_, Self>
    where
        Self: Sized
    {
        ProteinTextSlice::new(self, start, end)
    }
}

// ── ProteinTextSlice ──────────────────────────────────────────────────────────

/// A borrowed window onto the text, used to compare a candidate match against a query.
pub struct ProteinTextSlice<'a, T: ProteinTextBackend> {
    text: &'a T,
    start: usize,
    end: usize
}

impl<'a, T: ProteinTextBackend> ProteinTextSlice<'a, T> {
    /// Borrows `start..end` (half-open) of `text`. Bounds are not checked here; an out-of-range
    /// window panics when read.
    pub fn new(text: &'a T, start: usize, end: usize) -> Self {
        Self { text, start, end }
    }

    /// Returns the residue at `index`, counted from the start of the slice.
    pub fn get(&self, index: usize) -> u8 {
        self.text.get(self.start + index)
    }
    /// Number of residues in the window.
    pub fn len(&self) -> usize {
        self.end - self.start
    }
    /// Whether the window is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Compares this window against `other`, optionally treating I and L as equal.
    ///
    /// # Contract
    ///
    /// Comparison stops at the shorter of the two, so a longer `other` compares equal on its
    /// prefix. Callers are expected to have established the lengths already — the search path
    /// derives the window from the query length — so this does not re-check them.
    ///
    /// `equate_il` exists because the index may be built with L normalised to I: when the caller
    /// asks for I/L equivalence, a query I or L matches either residue in the text.
    #[inline]
    pub fn equals_slice(&self, other: &[u8], equate_il: bool) -> bool {
        if equate_il {
            other
                .iter()
                .zip(self.iter())
                .all(|(&s, t)| s == t || (s == b'I' && t == b'L') || (s == b'L' && t == b'I'))
        } else {
            other.iter().zip(self.iter()).all(|(&s, t)| s == t)
        }
    }

    /// Re-checks only the positions where the query holds an I or an L.
    ///
    /// The fast path for a query with no I/L is a plain equality test against an index built with
    /// L normalised to I. When the query *does* contain I or L, that test is too permissive at
    /// exactly those positions, so they are re-compared against the original characters here
    /// instead of falling back to a slower general comparison over the whole window.
    ///
    /// # Contract
    ///
    /// `il_locations` are positions in the unskipped query, and every one must be at least
    /// `skip`; the caller filters them when it computes them. A location below `skip` underflows
    /// `il_location - skip`.
    pub fn check_il_locations(&self, skip: usize, il_locations: &[usize], search_string: &[u8]) -> bool {
        for &il_location in il_locations {
            debug_assert!(il_location >= skip, "il_location {il_location} is below skip {skip}");
            let index = il_location - skip;
            if search_string[index] != self.get(index) {
                return false;
            }
        }
        true
    }

    /// Iterates the residues in this window.
    pub fn iter(&self) -> ProteinTextSliceIterator<'_, T> {
        ProteinTextSliceIterator { text: self.text, pos: self.start, end: self.end }
    }
}

// ── ProteinTextIterator ───────────────────────────────────────────────────────

/// Iterator over an entire [`ProteinTextBackend`].
pub struct ProteinTextIterator<'a, T: ProteinTextBackend> {
    pub(crate) protein_text: &'a T,
    pub(crate) index: usize
}

/// Iterator over a [`ProteinTextSlice`].
pub struct ProteinTextSliceIterator<'a, T: ProteinTextBackend> {
    text: &'a T,
    pos: usize,
    end: usize
}

impl<T: ProteinTextBackend> Iterator for ProteinTextSliceIterator<'_, T> {
    type Item = u8;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.end {
            return None;
        }
        let c = self.text.get(self.pos);
        self.pos += 1;
        Some(c)
    }
}

impl<T: ProteinTextBackend> Iterator for ProteinTextIterator<'_, T> {
    type Item = u8;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.protein_text.len() {
            return None;
        }
        self.index += 1;
        Some(self.protein_text.get(self.index - 1))
    }
}
