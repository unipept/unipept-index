//! Shared test fixtures for the `sa_searcher` submodule tests.
//!
//! # Why these go through files
//!
//! The searcher is generic over the three storage backends, and the `mmap` / `preloaded-*`
//! features pick which ones [`SuffixArray`], `Proteins` and [`SuffixToProteinMapping`] name — nine
//! combinations in all. A fixture that constructed a backend directly could therefore only ever
//! build the preloaded one, which is why these tests used to be gated off under `--features mmap`
//! and never ran in the configuration that ships.
//!
//! So they build the way production does instead: the always-compiled preloaded types write the
//! three files `sa-builder` writes, and [`load_fixture`] reads each one back through whichever
//! loader the active build selects — the same predicates `sa-server`'s `load_by_backend!` uses.
//! Every test below then runs against the backends its build actually uses, and search results
//! must be identical either way, which is the property worth testing.

use std::{
    io::Write,
    ops::{Deref, DerefMut}
};

use sa_mappings::proteins::{InMemoryProteins, Protein, Proteins};
use tempfile::NamedTempFile;
use text_compression::{InMemoryProteinText, WriteBinary};

use crate::{
    KmerTable, SuffixArray,
    array::OriginalSA,
    sa_searcher::Searcher,
    suffix_to_protein_index::{
        BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinMapping
    }
};

/// Reads a fixture file back as `$ty`, mapping it exactly when the active build would.
///
/// The predicates mirror `sa_server::load_by_backend!` one for one; keep them in step. Note the
/// protein one in particular: `proteins.bin` holds both the metadata and the text, so it has to be
/// *mapped* whenever either section is, not merely when the metadata is.
macro_rules! load_fixture {
    ($ty:ty, $file:expr, $mapped_when:meta) => {{
        #[cfg($mapped_when)]
        {
            <$ty as text_compression::ReadBinaryMmap>::read_binary_mmap($file.path()).unwrap()
        }
        #[cfg(not($mapped_when))]
        {
            let mut reader = std::io::BufReader::new(std::fs::File::open($file.path()).unwrap());
            <$ty as text_compression::ReadBinary>::read_binary(&mut reader).unwrap()
        }
    }};
}

/// Serialises one structure to a temporary file. The caller keeps the handle alive for as long as
/// the structure built from it is in use — under `mmap` the mapping borrows this file.
fn write_fixture(value: impl WriteBinary) -> NamedTempFile {
    let mut tmp = NamedTempFile::new().unwrap();
    let mut buf = Vec::new();
    value.write_binary(&mut buf).unwrap();
    tmp.write_all(&buf).unwrap();
    tmp.flush().unwrap();
    tmp
}

/// A fixture searcher, plus the files the active build may still be reading out of.
///
/// Under `mmap` the three structures borrow mappings of the temporary files below, so the handles
/// have to outlive the searcher; a preloaded build has copied everything out and simply drops them
/// at the end of the test. [`Deref`] keeps the wrapper invisible at the call sites.
pub(crate) struct TestSearcher {
    searcher: Searcher<SuffixArray>,
    _backing: Vec<NamedTempFile>
}

impl Deref for TestSearcher {
    type Target = Searcher<SuffixArray>;

    fn deref(&self) -> &Self::Target {
        &self.searcher
    }
}

impl DerefMut for TestSearcher {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.searcher
    }
}

impl TestSearcher {
    /// Forwards to [`Searcher::with_kmer_table`], which consumes the searcher and so cannot be
    /// reached through `Deref`.
    pub(crate) fn with_kmer_table(self, table: KmerTable) -> Self {
        Self {
            searcher: self.searcher.with_kmer_table(table),
            _backing: self._backing
        }
    }
}

/// Which suffix-to-protein representation a fixture searcher should use.
///
/// All three answer identically; tests parameterise over them to check that the three
/// implementations agree, and that search works with whichever the index was built with.
pub(crate) enum Mapping {
    Dense,
    Sparse,
    BitVec
}

/// The one place a fixture searcher is built: write the three structures as `sa-builder` would,
/// then read each back through the active build's loader.
fn build_searcher(text: &str, proteins: Vec<Protein>, sa: Vec<i64>, sparseness: u8, mapping: Mapping) -> TestSearcher {
    let protein_text = InMemoryProteinText::from_string(text);

    let mapping_file = match mapping {
        Mapping::Dense => write_fixture(DenseSuffixToProtein::new(&protein_text)),
        Mapping::Sparse => write_fixture(SparseSuffixToProtein::new(&protein_text)),
        Mapping::BitVec => write_fixture(BitVecSuffixToProtein::new(&protein_text))
    };
    let sa_file = write_fixture(OriginalSA(sa, sparseness));
    let proteins_file = write_fixture(InMemoryProteins::new(protein_text, proteins));

    let searcher = Searcher::new(
        load_fixture!(SuffixArray, sa_file, feature = "mmap"),
        load_fixture!(
            Proteins,
            proteins_file,
            all(feature = "mmap", not(all(feature = "preloaded-text", feature = "preloaded-proteins")))
        ),
        load_fixture!(SuffixToProteinMapping, mapping_file, all(feature = "mmap", not(feature = "preloaded-mapping")))
    );

    TestSearcher {
        searcher,
        _backing: vec![sa_file, proteins_file, mapping_file]
    }
}

/// Example protein set used across the search/batched/retrieval tests.
///
/// Text `"AI-CLACVAA-AC-KCRLY$"` = four proteins (separated by `-`): `AI`, `CLACVAA`,
/// `AC`, `KCRLY`. Each is given a distinct taxon id (10/20/30/40) and uniprot id so
/// retrieval tests can assert which protein a suffix maps to.
const EXAMPLE_TEXT: &str = "AI-CLACVAA-AC-KCRLY$";

/// Four proteins with distinct taxa, matching [`EXAMPLE_TEXT`]'s separators.
fn example_protein_list() -> Vec<Protein> {
    (0..4)
        .map(|i| Protein {
            uniprot_id: format!("P{i}"),
            taxon_id: (i as u32 + 1) * 10,
            functional_annotations: vec![]
        })
        .collect()
}

/// Suffix array of [`EXAMPLE_TEXT`] at sparseness 1.
///
/// Precomputed rather than derived, because several tests assert against specific suffix
/// positions; changing the fixture text means recomputing this.
pub(crate) const EXAMPLE_SA_FULL: [i64; 20] = [19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18];

/// The same text sampled at sparseness 3 — only positions divisible by 3 are indexed.
pub(crate) const EXAMPLE_SA_SPARSE3: [i64; 7] = [9, 0, 3, 12, 15, 6, 18];

/// Builds a searcher over [`EXAMPLE_TEXT`] with an explicit suffix array and mapping.
///
/// Most search tests differ only in these three choices, so they share this rather than
/// re-inlining the three-line construction.
pub(crate) fn example_searcher_with(sa: &[i64], sparseness: u8, mapping: Mapping) -> TestSearcher {
    build_searcher(EXAMPLE_TEXT, example_protein_list(), sa.to_vec(), sparseness, mapping)
}

/// The common case: the full suffix array over the example proteins, with a BitVec mapping.
pub(crate) fn example_searcher() -> TestSearcher {
    example_searcher_with(&EXAMPLE_SA_FULL, 1, Mapping::BitVec)
}

/// A single protein of `n` copies of one residue, terminated — the fixture for anything that
/// needs an SA range larger than a hand-written array can conveniently give.
///
/// `n = 70` clears both defaults that gate the two-pass paths: `prefetch_threshold` (32) and
/// `retrieval_prefetch_distance` (32). Use `'A'` for a residue the I/L rules do not touch, and
/// `'L'` to make `equate_il` matter: `compare` normalises L to I, so searching `"I"` over an
/// all-`L` text matches every position during the bound search and every one of them must then
/// be rejected by validation.
pub(crate) fn repeated_residue_searcher(residue: char, n: usize) -> TestSearcher {
    searcher_over_text(&format!("{}$", residue.to_string().repeat(n)), 1)
}

/// Fixture for the left-extended tryptic search, positions annotated because the test cases
/// depend on them:
///
/// ```text
///   0 M  1 K  2 A  3 P  4 T  5 R  6 V  7 G  8 A  9 K  10 -
///  11 R 12 I 13 Y 14 N 15 K 16 P 17 Q 18 S 19 T  20 -
///  21 P 22 K 23 T 24 R 25 L 26 D 27 E 28 I  29 $
/// ```
///
/// Protein starts at 0, 11, 21; separators at 10, 20; termination at 29. It deliberately contains
/// K/R cut sites (1, 5, 9, 11, 15, 22, 24), a proline-blocked cut (15→16 is K then P), a protein
/// that starts with proline (21, at an *odd* — hence unsampled at sparseness 2 — position), and
/// I/L so `equate_il=false` is exercised.
pub(crate) const TRYPTIC_FIXTURE: &str = "MKAPTRVGAK-RIYNKPQST-PKTRLDEI$";

/// Every pure-amino-acid substring of [`TRYPTIC_FIXTURE`] of length 3..=6, as a peptide corpus.
/// Length 3 is the floor so the corpus stays usable at sparseness 3.
pub(crate) fn tryptic_fixture_peptides() -> Vec<Vec<u8>> {
    let bytes = TRYPTIC_FIXTURE.as_bytes();
    let mut peptides: Vec<Vec<u8>> = Vec::new();
    for len in 3..=6usize {
        for start in 0..bytes.len().saturating_sub(len) {
            let s = &bytes[start..start + len];
            if s.iter().all(|c| c.is_ascii_uppercase()) {
                peptides.push(s.to_vec());
            }
        }
    }
    peptides.sort();
    peptides.dedup();
    peptides
}

/// Builds a `Searcher` over an arbitrary `-`-separated, `$`-terminated protein text, with the
/// suffix array computed by brute force (fine at fixture scale).
///
/// Suffixes are ordered on the L → I normalised text, exactly as the production index is
/// built, so `compare`'s own normalisation is consistent with the array it searches — without
/// that, sparse/L-containing fixtures would exercise a suffix array no real build produces.
///
/// `sparseness` keeps only the SA entries whose text position is a multiple of it, which is
/// the sparse layout `search_matching_suffixes` compensates for with its `skip` loop.
/// Each `-`-separated segment becomes a protein with taxon `10, 20, …`.
pub(crate) fn searcher_over_text(text: &str, sparseness: u8) -> TestSearcher {
    let normalised: Vec<u8> = text.bytes().map(|c| if c == b'L' { b'I' } else { c }).collect();

    let mut positions: Vec<i64> =
        (0..text.len() as i64).filter(|p| (*p as usize).is_multiple_of(sparseness as usize)).collect();
    positions.sort_by(|&a, &b| normalised[a as usize..].cmp(&normalised[b as usize..]));

    let proteins = text
        .trim_end_matches('$')
        .split('-')
        .enumerate()
        .map(|(i, _)| Protein {
            uniprot_id: format!("P{i}"),
            taxon_id: (i as u32 + 1) * 10,
            functional_annotations: vec![]
        })
        .collect();

    build_searcher(text, proteins, positions, sparseness, Mapping::BitVec)
}
