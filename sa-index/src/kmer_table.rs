//! A precomputed suffix-array bounds table for every k-mer.
//!
//! Binary searching the suffix array costs O(log n) dependent, cache-missing probes before the
//! search has even narrowed to the right neighbourhood. This table skips that opening phase: for
//! every possible k-mer it stores the SA range whose suffixes start with it, so a search for a
//! peptide of at least k residues starts from those bounds instead of from the whole array.
//!
//! The table is dense — `alphabet_size^k` entries — so it is built once by `sa-builder` and
//! loaded by the server. That density is why `k` is small; see [`MAX_KMER_K`].
//!
//! It is an accelerator only: results are identical with and without it, which the golden
//! configuration matrix asserts.

use std::{
    error::Error,
    io::{BufRead, Write},
    sync::atomic::{AtomicUsize, Ordering}
};

use rayon::prelude::*;
use text_compression::ProteinTextBackend;

use crate::{LoadIndex, ReadBinary, WriteBinary, array::SuffixArrayBackend};

/// Amino acid alphabet used for k-mer indexing (no J; L is treated as I).
/// Index in this slice + 1 gives the 1-based `ascii_array` value for each character.
const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWYXBUZO";

/// Number of distinct amino acid values after normalizing L → I. No J; L shares I's slot → 24.
pub const AMINO_ACID_COUNT: usize = 24;

/// Maximum allowed k for a k-mer table.
///
/// Memory cost is `AMINO_ACID_COUNT^k × 16 bytes`:
///   k=4 →   ~5 MB (L3-cache friendly, recommended)
///   k=5 → ~128 MB
///   k=6 →   ~3 GB
///
/// Values above 7 are almost certainly a misconfiguration.
pub const MAX_KMER_K: usize = 7;

/// Builds the `ascii_array` lookup table at compile time: maps *any* byte → 1-based amino acid
/// index (0 = not in alphabet). L is mapped to the same slot as I so L→I normalization is free.
///
/// Sized for the whole `u8` range, not just ASCII, so that indexing it with an arbitrary query
/// byte is total. A peptide arrives from the network and is only uppercased before it reaches
/// [`KmerTable::lookup`]; a 128-entry table made every byte >= 128 an out-of-bounds index, one
/// line before the `char_idx == 0` test that exists to reject exactly those characters.
fn build_ascii_array() -> [u8; 256] {
    let mut array = [0u8; 256];
    let mut next_index: u8 = 1;
    for &c in ALPHABET {
        if c == b'L' {
            // L is not assigned its own slot; it maps to I's index instead (set below).
            continue;
        }
        array[c as usize] = next_index;
        next_index += 1;
    }
    // Map L to the same 1-based index as I so queries with L hit the same table entry.
    array[b'L' as usize] = array[b'I' as usize];
    array
}

/// Pre-computed k-mer SA bounds lookup table for accelerating suffix array binary search.
///
/// Stores the inclusive `(min_bound, max_bound)` SA range for every k-character amino acid
/// prefix. At query time the first `k` characters narrow the binary search window from the
/// full SA length (~32 iterations) to the k-mer's range (~13 iterations for k=4), reducing
/// random memory accesses and TLB pressure by ~60 %.
///
/// Memory: `AMINO_ACID_COUNT^k × 16 bytes` — ~5.3 MB for k=4 (fits in L3 cache).
pub struct KmerTable {
    /// Length of the k-mer prefix.
    pub k: usize,
    /// Maps any byte → 1-based amino acid index (0 = not in alphabet).
    /// L maps to the same index as I for transparent L→I normalization.
    ascii_array: [u8; 256],
    /// Flat `(min_bound, max_bound)` pairs indexed by `kmer_to_index(kmer)`.
    /// Absent k-mers are represented by `min_bound > max_bound`
    /// (sentinel: `(usize::MAX, 0)`).
    bounds: Vec<(usize, usize)>,
    /// The largest `max_bound` any present k-mer carries, or 0 for a table with no entries.
    ///
    /// Not part of the on-disk format — it is accumulated while reading, which the loader was
    /// already iterating anyway, so an existing `kmer_table.bin` still loads unchanged. It exists
    /// so [`Searcher::with_kmer_table`](crate::sa_searcher::Searcher::with_kmer_table) can reject a
    /// table built against a *different* suffix array: the searcher feeds these bounds straight
    /// into the binary search, so a table from a larger index sends `get` past the end of a smaller
    /// one — a panic on the preloaded backend, and a fabricated suffix position on the mmap one.
    highest_bound: usize
}

impl KmerTable {
    /// Builds the k-mer table via a single O(n) linear scan of the suffix array.
    ///
    /// Because the SA is sorted, each k-mer's entries are contiguous: the first
    /// occurrence gives `min_bound` and the last gives `max_bound`.
    pub fn build_from_sa<SA: SuffixArrayBackend, T: ProteinTextBackend + Sync>(sa: &SA, text: &T, k: usize) -> Self {
        Self::build_kmer_table(sa.len(), |i| sa.get(i) as usize, text.len(), |i| text.get(i), k)
    }

    /// Same as `build_from_sa` but accepts the raw suffix array as a plain slice
    /// and text access via closures — works regardless of whether the mmap feature is active.
    pub fn build_from_raw_sa(sa: &[i64], text_len: usize, get_char: impl Fn(usize) -> u8 + Sync, k: usize) -> Self {
        Self::build_kmer_table(sa.len(), |i| sa[i] as usize, text_len, get_char, k)
    }

    fn build_kmer_table(
        sa_len: usize,
        get_sa: impl Fn(usize) -> usize + Sync,
        text_len: usize,
        get_char: impl Fn(usize) -> u8 + Sync,
        k: usize
    ) -> Self {
        assert!(
            k <= MAX_KMER_K,
            "k={k} exceeds MAX_KMER_K={MAX_KMER_K} (memory cost would be ~{} MB)",
            AMINO_ACID_COUNT.saturating_pow(k as u32).saturating_mul(16) / (1024 * 1024)
        );
        let ascii_array = build_ascii_array();
        let table_size = AMINO_ACID_COUNT.pow(k as u32);

        // Sentinel: (MAX, 0) means "absent". AtomicUsize lets multiple threads update
        // min/max without locks; fetch_min/fetch_max are stable since Rust 1.45.
        let atomic_bounds: Vec<(AtomicUsize, AtomicUsize)> =
            (0..table_size).map(|_| (AtomicUsize::new(usize::MAX), AtomicUsize::new(0))).collect();

        let kmer_index = |suffix_start: usize| -> Option<usize> {
            let mut idx = 0usize;
            for j in 0..k {
                let pos = suffix_start + j;
                if pos >= text_len {
                    return None;
                }
                let char_idx = ascii_array[get_char(pos) as usize];
                if char_idx == 0 {
                    return None;
                }
                idx = idx * AMINO_ACID_COUNT + (char_idx as usize - 1);
            }
            Some(idx)
        };

        (0..sa_len).into_par_iter().for_each(|i| {
            if let Some(idx) = kmer_index(get_sa(i)) {
                // SA is sorted: first occurrence gives min, last gives max.
                // Relaxed ordering is sufficient: we only need the final values after
                // the parallel section, not any inter-thread happens-before guarantees.
                atomic_bounds[idx].0.fetch_min(i, Ordering::Relaxed);
                atomic_bounds[idx].1.fetch_max(i, Ordering::Relaxed);
            }
        });

        let bounds: Vec<(usize, usize)> =
            atomic_bounds.into_iter().map(|(min, max)| (min.into_inner(), max.into_inner())).collect();

        let highest_bound = bounds.iter().filter(|(min, max)| min <= max).map(|(_, max)| *max).max().unwrap_or(0);

        Self { k, ascii_array, bounds, highest_bound }
    }

    /// Maps a byte slice to its flat table index.
    /// Returns `None` if any byte is outside the amino acid alphabet.
    /// The returned index is always `< AMINO_ACID_COUNT^k == self.bounds.len()`.
    #[inline]
    fn bytes_to_kmer_index(&self, kmer: &[u8]) -> Option<usize> {
        let mut idx = 0usize;
        for &c in kmer {
            let char_idx = self.ascii_array[c as usize];
            if char_idx == 0 {
                return None;
            }
            idx = idx * AMINO_ACID_COUNT + (char_idx as usize - 1);
        }
        Some(idx)
    }

    /// The largest SA index any k-mer in this table points at.
    ///
    /// A table is only valid for the suffix array it was built from; see the field docs.
    pub fn highest_bound(&self) -> usize {
        self.highest_bound
    }

    /// Looks up the inclusive `(min_bound, max_bound)` SA range for a k-mer prefix.
    ///
    /// Returns `None` if the k-mer is absent from all proteins.
    /// `kmer` must have exactly `k` bytes; L is treated the same as I.
    #[inline]
    pub fn lookup(&self, kmer: &[u8]) -> Option<(usize, usize)> {
        debug_assert_eq!(kmer.len(), self.k, "kmer length must equal table k");
        let idx = self.bytes_to_kmer_index(kmer)?;
        let &(min, max) = &self.bounds[idx];
        if min > max { None } else { Some((min, max)) }
    }
}

impl WriteBinary for KmerTable {
    fn write_binary<W: Write>(self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        if self.k > u8::MAX as usize {
            return Err(format!("k={} exceeds the maximum serializable value of 255", self.k).into());
        }
        writer.write_all(&[self.k as u8])?;
        writer.write_all(&(AMINO_ACID_COUNT as u64).to_le_bytes())?;
        let mut buf = [0u8; 16];
        for (min, max) in self.bounds {
            buf[..8].copy_from_slice(&(min as u64).to_le_bytes());
            buf[8..].copy_from_slice(&(max as u64).to_le_bytes());
            writer.write_all(&buf)?;
        }
        Ok(())
    }
}

impl ReadBinary for KmerTable {
    fn read_binary<R: BufRead>(reader: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut buf1 = [0u8; 1];
        reader.read_exact(&mut buf1)?;
        let k = buf1[0] as usize;

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let amino_acid_count = u64::from_le_bytes(buf8) as usize;

        if amino_acid_count != AMINO_ACID_COUNT {
            return Err(
                format!("k-mer table: expected amino_acid_count={AMINO_ACID_COUNT}, got {amino_acid_count}").into()
            );
        }

        if k > MAX_KMER_K {
            return Err(format!("k-mer table: k={k} exceeds MAX_KMER_K={MAX_KMER_K}").into());
        }

        let table_size = AMINO_ACID_COUNT.pow(k as u32);
        let mut bounds = Vec::new();
        bounds
            .try_reserve_exact(table_size)
            .map_err(|_| format!("The k-mer table header declares k={k}, whose table cannot be allocated"))?;
        let mut buf16 = [0u8; 16];
        let mut highest_bound = 0usize;
        for _ in 0..table_size {
            reader.read_exact(&mut buf16)?;
            let min = u64::from_le_bytes(buf16[..8].try_into()?) as usize;
            let max = u64::from_le_bytes(buf16[8..].try_into()?) as usize;
            // Absent k-mers carry the `(usize::MAX, 0)` sentinel, so they never raise the maximum.
            if min <= max {
                highest_bound = highest_bound.max(max);
            }
            bounds.push((min, max));
        }

        Ok(Self { k, ascii_array: build_ascii_array(), bounds, highest_bound })
    }
}

/// Only ever the owned route: the table is small relative to the index and has no mmap variant.
impl LoadIndex for KmerTable {
    fn load(path: &std::path::Path) -> Result<Self, Box<dyn Error>> {
        text_compression::load_owned(path)
    }
}

#[cfg(test)]
mod tests {
    use text_compression::InMemoryProteinText;

    use crate::{
        array::{InMemorySA, OriginalSA},
        kmer_table::KmerTable
    };

    fn build_test_table(input: &str, sa_values: Vec<i64>, k: usize) -> KmerTable {
        let text = InMemoryProteinText::from_string(input);
        let sa = InMemorySA::Original(OriginalSA(sa_values, 1));
        KmerTable::build_from_sa(&sa, &text, k)
    }

    /// Three structures from different builds must be refused, not silently mixed.
    ///
    /// Each file is well-formed on its own, so every loader used to succeed and the server started
    /// and reported itself ready — with protein indices resolved against the wrong text. No build
    /// identifier exists in the format to compare, so the check uses two relationships that hold
    /// implicitly and exactly.
    #[test]
    fn structures_from_different_builds_are_rejected() {
        use sa_mappings::proteins::{InMemoryProteins, Protein};
        use text_compression::InMemoryProteinText;

        use crate::{
            array::{InMemorySA, OriginalSA},
            sa_searcher::Searcher,
            suffix_to_protein_index::{InMemorySuffixToProteinMapping, preloaded::BitVecSuffixToProtein}
        };

        let protein = |id: &str| Protein {
            uniprot_id: id.to_string(),
            taxon_id: 1,
            functional_annotations: vec![]
        };

        // A matching set loads: text of 5, sample rate 1, so the SA holds 5 entries.
        let text = InMemoryProteinText::from_string("ACAC$");
        let stp = BitVecSuffixToProtein::new(&text);
        let proteins = InMemoryProteins::new(text, vec![protein("P0")]);
        let sa = InMemorySA::Original(OriginalSA(vec![4, 2, 0, 3, 1], 1));
        assert!(
            Searcher::try_new(sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp)).is_ok(),
            "a matching set must be accepted"
        );

        // A suffix array from a longer text: right shape, wrong build.
        let text = InMemoryProteinText::from_string("ACAC$");
        let stp = BitVecSuffixToProtein::new(&text);
        let proteins = InMemoryProteins::new(text, vec![protein("P0")]);
        let wrong_sa = InMemorySA::Original(OriginalSA(vec![6, 4, 2, 0, 5, 3, 1], 1));
        let err = Searcher::try_new(wrong_sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp))
            .err()
            .expect("a suffix array of the wrong length must be rejected");
        assert!(err.contains("different builds"), "unexpected error: {err}");

        // A mapping built for a longer text, with a matching suffix array.
        let text = InMemoryProteinText::from_string("ACAC$");
        let other = InMemoryProteinText::from_string("ACACACAC$");
        let stp = BitVecSuffixToProtein::new(&other);
        let proteins = InMemoryProteins::new(text, vec![protein("P0")]);
        let sa = InMemorySA::Original(OriginalSA(vec![4, 2, 0, 3, 1], 1));
        let err = Searcher::try_new(sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp))
            .err()
            .expect("a mapping built for another text must be rejected");
        assert!(err.contains("different builds"), "unexpected error: {err}");
    }

    /// A table built from one suffix array must be refused by a searcher holding a different one.
    ///
    /// Nothing in the file format identifies the build, so a rebuilt `sa.bin` with a stale
    /// `kmer_table.bin` used to start cleanly and then fail mid-query — a panic on the preloaded
    /// backend, a fabricated suffix position on the mmap one.
    #[test]
    fn a_table_from_a_different_index_is_rejected() {
        use text_compression::InMemoryProteinText;

        use crate::{
            array::{InMemorySA, OriginalSA},
            sa_searcher::Searcher,
            suffix_to_protein_index::{InMemorySuffixToProteinMapping, preloaded::BitVecSuffixToProtein}
        };

        // A table whose bounds reach index 4, built over a 5-entry array.
        let big_text = InMemoryProteinText::from_string("ACAC$");
        let big_sa = InMemorySA::Original(OriginalSA(vec![4, 2, 0, 3, 1], 1));
        let table = KmerTable::build_from_sa(&big_sa, &big_text, 2);
        assert!(table.highest_bound() > 0, "fixture table should point somewhere");

        // A searcher over a strictly smaller array cannot use it.
        let small_text = InMemoryProteinText::from_string("AC$");
        let small_sa = InMemorySA::Original(OriginalSA(vec![2, 0, 1], 1));
        let stp = BitVecSuffixToProtein::new(&small_text);
        let proteins = sa_mappings::proteins::InMemoryProteins::new(small_text, vec![sa_mappings::proteins::Protein {
            uniprot_id: "P0".to_string(),
            taxon_id: 1,
            functional_annotations: vec![]
        }]);
        let searcher = Searcher::new(small_sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp));

        let err = searcher.try_with_kmer_table(table).err().expect("a mismatched table must be rejected");
        assert!(err.contains("different index"), "unexpected error: {err}");

        // And the matching pair is accepted, so the check is specific.
        let text = InMemoryProteinText::from_string("ACAC$");
        let sa = InMemorySA::Original(OriginalSA(vec![4, 2, 0, 3, 1], 1));
        let good = KmerTable::build_from_sa(&sa, &text, 2);
        let stp = BitVecSuffixToProtein::new(&text);
        let proteins = sa_mappings::proteins::InMemoryProteins::new(text, vec![sa_mappings::proteins::Protein {
            uniprot_id: "P0".to_string(),
            taxon_id: 1,
            functional_annotations: vec![]
        }]);
        let searcher = Searcher::new(sa, proteins, InMemorySuffixToProteinMapping::BitVec(stp));
        assert!(searcher.try_with_kmer_table(good).is_ok(), "the matching table must be accepted");
    }

    #[test]
    fn test_lookup_present() {
        let table = build_test_table("ACAC$", vec![4, 2, 0, 3, 1], 2);
        let result = table.lookup(b"AC");
        assert!(result.is_some(), "AC should be found");
        let (min, max) = result.unwrap();
        assert!(min <= max);
    }

    /// A query byte is a raw `u8`, and nothing on the request path filters the alphabet:
    /// `peptide_search` only uppercases, and `to_uppercase` is Unicode-aware, so a non-ASCII
    /// character arrives here as multi-byte UTF-8 whose every byte is >= 128. The alphabet table
    /// must therefore be indexable by all 256 values — it used to be `[u8; 128]`, which made this
    /// an out-of-bounds index one line *before* the `char_idx == 0` test that rejects such bytes.
    #[test]
    fn test_lookup_rejects_bytes_outside_ascii_without_panicking() {
        let table = build_test_table("ACAC$", vec![4, 2, 0, 3, 1], 2);

        // 'é' (U+00E9) is 0xC3 0xA9 in UTF-8 — what `to_uppercase()` leaves in the query.
        assert!(table.lookup(&[0xC3, 0xA9]).is_none(), "non-ASCII bytes are not in the alphabet");
        // The boundary itself, and the top of the range.
        assert!(table.lookup(&[0x80, 0x80]).is_none());
        assert!(table.lookup(&[0xFF, 0xFF]).is_none());
        // A valid residue followed by one that is out of range.
        assert!(table.lookup(&[b'A', 0xFF]).is_none());
    }

    #[test]
    fn test_lookup_absent() {
        let table = build_test_table("ACAC$", vec![4, 2, 0, 3, 1], 2);
        let result = table.lookup(b"ZZ");
        assert!(result.is_none(), "ZZ should not be found");
    }

    #[test]
    fn test_lookup_separator_returns_none() {
        let table = build_test_table("AC-AC$", vec![5, 2, 0, 3, 4, 1], 2);
        assert!(table.lookup(b"A-").is_none());
    }

    #[test]
    fn test_il_normalization() {
        let table = build_test_table("AIAC$", vec![4, 0, 3, 1, 2], 2);
        let result_i = table.lookup(b"AI");
        let result_l = table.lookup(b"AL");
        assert_eq!(result_i, result_l);
    }

    #[test]
    fn test_roundtrip_serialization() {
        use std::io::{BufReader, BufWriter};

        use crate::{ReadBinary, WriteBinary};

        let table = build_test_table("ACAC$", vec![4, 2, 0, 3, 1], 2);
        let k = table.k;
        let original_lookup = table.lookup(b"AC");

        let mut buf = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buf);
            table.write_binary(&mut writer).unwrap();
        }

        let mut reader = BufReader::new(buf.as_slice());
        let restored = KmerTable::read_binary(&mut reader).unwrap();

        assert_eq!(restored.k, k);
        assert_eq!(restored.lookup(b"AC"), original_lookup);
    }

    /// The table is the one structure with no mmap variant, so its [`LoadIndex`] impl is the only
    /// one no backend-combination test reaches. `sa-server` calls it on every start with
    /// `--kmer-table-file`.
    #[test]
    fn loads_from_a_file() {
        use std::io::Write;

        use crate::{LoadIndex, WriteBinary};

        let table = build_test_table("ACAC$", vec![4, 2, 0, 3, 1], 2);
        let expected = table.lookup(b"AC");

        let mut file = tempfile::NamedTempFile::new().unwrap();
        let mut buf = Vec::new();
        table.write_binary(&mut buf).unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();

        let loaded = KmerTable::load(file.path()).unwrap();

        assert_eq!(loaded.k, 2);
        assert_eq!(loaded.lookup(b"AC"), expected);
    }

    #[test]
    fn loading_a_missing_file_errors() {
        use crate::LoadIndex;

        assert!(KmerTable::load(std::path::Path::new("/nonexistent/kmer_table.bin")).is_err());
    }
}
