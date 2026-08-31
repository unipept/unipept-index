//! Turning a text position into a protein index.
//!
//! Search yields positions in the concatenated text; results are about proteins. This module
//! answers "which protein contains position i?", and returns `u32::NULL` for the separator and
//! terminator bytes, which belong to no protein.
//!
//! Three representations, trading space against lookup cost, chosen at build time
//! (`sa-builder --mapping-style`) and recorded in the file:
//!
//! * **Dense** — one `u32` per text position. One load per lookup, but 4 bytes per residue: ~1.2
//!   GB over the ~300 M-residue reference text, and ~300 GB at the scale of the full UniProt
//!   index, where it would be larger than the 160 GB suffix array and larger than everything else
//!   in that index put together. A choice for small databases only.
//! * **Sparse** — the start position of each protein, binary-searched. Smallest, but with m
//!   proteins it costs O(log m) dependent loads per lookup, each likely a cache miss.
//! * **BitVec** — a bit per text position marking the separators and the terminator, with a rank
//!   structure over it. Near-dense speed at ~1.25 bits per position; the default.
//!
//! Each has a preloaded and an mmap implementation, both always compiled and both wrapped by a
//! dispatch enum — [`InMemorySuffixToProteinMapping`] and [`MmapBackedSuffixToProteinMapping`] —
//! which is what a caller names to pick one. See the crate docs.

/// Generates a [`SuffixToProteinMappingBackend`] impl that forwards each listed method to the
/// active enum variant (`Dense`, `Sparse` or `BitVec`).
///
/// Both backend enums are the same three-way match repeated per method; they differ only in which
/// methods they list, since the preloaded one keeps the default no-op `touch_all_pages`. Declared
/// before the module declarations below so that its textual scope reaches both.
///
/// This is the one thing the two halves share, and it is only the dispatch shell: the lookups it
/// forwards to stay separate per the crate docs, so tuning one backend cannot perturb the other.
macro_rules! delegate_suffix_to_protein_mapping {
    ($mapping:ty { $(fn $method:ident(&self $(, $arg:ident: $arg_ty:ty)*) $(-> $ret:ty)?;)* }) => {
        impl $crate::suffix_to_protein_index::SuffixToProteinMappingBackend for $mapping {
            $(
                #[inline]
                fn $method(&self $(, $arg: $arg_ty)*) $(-> $ret)? {
                    match self {
                        Self::Dense(m) => m.$method($($arg),*),
                        Self::Sparse(m) => m.$method($($arg),*),
                        Self::BitVec(m) => m.$method($($arg),*)
                    }
                }
            )*
        }
    };
}

// TEMPORARY, removed in the searcher PR: the pre-split implementation, moved here unchanged.
// Not re-exported at this module root -- the new backends use the same type names, so callers
// name `legacy` explicitly until they switch.
pub mod legacy;

pub mod mmap;
pub mod preloaded;
#[cfg(test)]
mod test_utils;

pub use mmap::{
    MmapBackedSuffixToProteinMapping, MmapBitVecSuffixToProtein, MmapDenseSuffixToProtein, MmapSparseSuffixToProtein
};
pub use preloaded::{
    BitVecSuffixToProtein, DenseSuffixToProtein, InMemorySuffixToProteinMapping, SparseSuffixToProtein
};

/// Trait implemented by the SuffixToProtein mappings
pub trait SuffixToProteinMappingBackend: Send + Sync {
    /// Returns the index of the protein in the protein list for the given suffix, or `u32::NULL`
    /// if the position holds a separator or the terminator and so belongs to no protein.
    ///
    /// `suffix` must be a position within the text. Positions past the end are not a checked
    /// error: an implementation may return anything or panic.
    fn suffix_to_protein(&self, suffix: i64) -> u32;

    /// Non-blocking hardware prefetch hint for the data that `suffix_to_protein(suffix)` will access.
    ///
    /// Default is a no-op, kept by the two implementations that cannot name the address the
    /// lookup will touch: both sparse mappings, which binary-search. It accepts any `suffix`, in
    /// range or not, since a hint that misses is only a wasted load.
    #[inline]
    fn prefetch_for_suffix(&self, _suffix: i64) {}

    /// Reads at least one byte from every OS page in the mmap backing this mapping, ensuring all
    /// pages are resident in the page cache. Returns the number of bytes swept, so a caller can
    /// report a bandwidth.
    /// Default is 0; mmap-backed implementations override this.
    #[inline]
    fn touch_all_pages(&self) -> u64 {
        0
    }

    /// The text length this mapping was built for, when the representation records it.
    ///
    /// Used to check that the three index files came from the same `sa-builder` run — see
    /// `Searcher::try_new`. Dense stores one entry per
    /// text position and BitVec one bit, so both know the length exactly. Sparse stores protein
    /// start positions, which says nothing about the text length, so it keeps the `None` default:
    /// a mapping that cannot be compared is not the same as one that disagrees.
    #[inline]
    fn implied_text_len(&self) -> Option<usize> {
        None
    }
}
