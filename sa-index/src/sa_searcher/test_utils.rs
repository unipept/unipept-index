//! Shared test fixtures for the `sa_searcher` submodule tests.
//!
//! # Why these go through files
//!
//! The searcher is generic over its three storage backends, and both implementations of each are
//! always compiled, so a fixture can build *any* of the sixteen combinations. What it cannot do is
//! construct a mapped backend directly: those types only come into existence by reading a file.
//!
//! So the fixtures build the way production does. The owned types write the three files
//! `sa-builder` writes, and each is read back through [`LoadIndex`], which every backend
//! implements by whichever route it needs. The type parameters alone therefore decide the
//! combination — there is no `#[cfg]` here, and none is mirrored from `sa-server`.
//!
//! Most tests want one combination and take the default, [`PreloadedSearcher`]. The one that wants
//! all sixteen is `super::backend_agreement`, which asserts they answer identically.

use std::{
    io::Write,
    ops::{Deref, DerefMut}
};

use sa_mappings::proteins::{InMemoryProteins, Protein, ProteinsBackend};
use tempfile::NamedTempFile;
use text_compression::{InMemoryProteinText, LoadIndex, WriteBinary};

use crate::{
    KmerTable,
    array::{InMemorySA, OriginalSA, SuffixArrayBackend},
    sa_searcher::{SearchTuning, Searcher},
    suffix_to_protein_index::{
        BitVecSuffixToProtein, DenseSuffixToProtein, InMemorySuffixToProteinMapping, SparseSuffixToProtein,
        SuffixToProteinMappingBackend
    }
};

/// Everything owned — the combination the ordinary tests use, named because three concrete type
/// parameters are unreadable in a signature.
pub(crate) type PreloadedSearcher =
    Searcher<InMemorySA, InMemoryProteins<InMemoryProteinText>, InMemorySuffixToProteinMapping>;

/// Serialises one structure to a temporary file. The caller keeps the handle alive for as long as
/// the structure built from it is in use — a mapped backend borrows this file.
fn write_fixture(value: impl WriteBinary) -> NamedTempFile {
    let mut tmp = NamedTempFile::new().unwrap();
    let mut buf = Vec::new();
    value.write_binary(&mut buf).unwrap();
    tmp.write_all(&buf).unwrap();
    tmp.flush().unwrap();
    tmp
}

/// A fixture searcher, plus the files its backends may still be reading out of.
///
/// A mapped backend borrows a mapping of one of the temporary files below, so the handles have to
/// outlive the searcher; owned backends have copied everything out and simply drop them at the end
/// of the test. [`Deref`] keeps the wrapper invisible at the call sites.
///
/// The parameters default to the owned types, so `TestSearcher` alone means
/// [`PreloadedSearcher`].
pub(crate) struct TestSearcher<
    SA = InMemorySA,
    P = InMemoryProteins<InMemoryProteinText>,
    STPM = InMemorySuffixToProteinMapping
> where
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
{
    searcher: Searcher<SA, P, STPM>,
    _backing: Vec<NamedTempFile>
}

impl<SA, P, STPM> Deref for TestSearcher<SA, P, STPM>
where
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
{
    type Target = Searcher<SA, P, STPM>;

    fn deref(&self) -> &Self::Target {
        &self.searcher
    }
}

impl<SA, P, STPM> DerefMut for TestSearcher<SA, P, STPM>
where
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.searcher
    }
}

impl<SA, P, STPM> TestSearcher<SA, P, STPM>
where
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
{
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
/// then read each back as whichever backend the caller asked for.
pub(crate) fn build_searcher<SA, P, STPM>(
    text: &str,
    proteins: Vec<Protein>,
    sa: Vec<i64>,
    sparseness: u8,
    mapping: Mapping
) -> TestSearcher<SA, P, STPM>
where
    SA: SuffixArrayBackend + LoadIndex,
    P: ProteinsBackend + LoadIndex,
    STPM: SuffixToProteinMappingBackend + LoadIndex
{
    let protein_text = InMemoryProteinText::from_string(text);

    let mapping_file = match mapping {
        Mapping::Dense => write_fixture(DenseSuffixToProtein::new(&protein_text)),
        Mapping::Sparse => write_fixture(SparseSuffixToProtein::new(&protein_text)),
        Mapping::BitVec => write_fixture(BitVecSuffixToProtein::new(&protein_text))
    };
    let sa_file = write_fixture(OriginalSA(sa, sparseness));
    let proteins_file = write_fixture(InMemoryProteins::new(protein_text, proteins));

    let searcher = Searcher::new(
        SA::load(sa_file.path()).unwrap(),
        P::load(proteins_file.path()).unwrap(),
        STPM::load(mapping_file.path()).unwrap()
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
pub(crate) const EXAMPLE_TEXT: &str = "AI-CLACVAA-AC-KCRLY$";

/// Four proteins with distinct taxa, matching [`EXAMPLE_TEXT`]'s separators.
pub(crate) fn example_protein_list() -> Vec<Protein> {
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
/// `n` must clear both defaults that gate the two-pass paths: `prefetch_threshold` and
/// `retrieval_prefetch_distance`, both 32. Callers pass 70. The requirement is asserted rather than
/// left to a comment because a fixture that no longer reaches the two-pass path does not fail — it
/// quietly tests the scalar loop instead and goes on passing, which is the failure mode worth
/// catching. A default raised past 70 is exactly what would trigger it.
///
/// Use `'A'` for a residue the I/L rules do not touch, and `'L'` to make `equate_il` matter:
/// `compare` normalises L to I, so searching `"I"` over an all-`L` text matches every position
/// during the bound search and every one of them must then be rejected by validation.
pub(crate) fn repeated_residue_searcher(residue: char, n: usize) -> TestSearcher {
    let tuning = SearchTuning::default();
    let gate = tuning.prefetch_threshold.max(tuning.retrieval_prefetch_distance);
    assert!(
        n > gate,
        "repeated_residue_searcher: n={n} does not clear the two-pass gate ({gate}); \
         this fixture would silently exercise the scalar path instead"
    );
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
