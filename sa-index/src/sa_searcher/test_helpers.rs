//! Shared test fixtures for the `sa_searcher` submodule tests.

use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
use text_compression::ProteinText;

use crate::{
    array::OriginalSA,
    sa_searcher::Searcher,
    suffix_to_protein_index::{
        BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinMapping,
    },
    SuffixArray,
};

/// Example protein set used across the search/batched/retrieval tests.
///
/// Text `"AI-CLACVAA-AC-KCRLY$"` = four proteins (separated by `-`): `AI`, `CLACVAA`,
/// `AC`, `KCRLY`. Each is given a distinct taxon id (10/20/30/40) and uniprot id so
/// retrieval tests can assert which protein a suffix maps to.
pub(crate) fn get_example_proteins() -> Proteins {
    let text = ProteinText::from_string("AI-CLACVAA-AC-KCRLY$");
    Proteins::new(text, vec![
        Protein { uniprot_id: "P0".to_string(), taxon_id: 10, functional_annotations: vec![] },
        Protein { uniprot_id: "P1".to_string(), taxon_id: 20, functional_annotations: vec![] },
        Protein { uniprot_id: "P2".to_string(), taxon_id: 30, functional_annotations: vec![] },
        Protein { uniprot_id: "P3".to_string(), taxon_id: 40, functional_annotations: vec![] },
    ])
}

/// Suffix array of [`get_example_proteins`]'s text at sparseness 1.
///
/// Precomputed rather than derived, because several tests assert against specific suffix
/// positions; changing the fixture text means recomputing this.
pub(crate) const EXAMPLE_SA_FULL: [i64; 20] =
    [19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18];

/// The same text sampled at sparseness 3 — only positions divisible by 3 are indexed.
pub(crate) const EXAMPLE_SA_SPARSE3: [i64; 7] = [9, 0, 3, 12, 15, 6, 18];

/// Which suffix-to-protein representation a fixture searcher should use.
///
/// All three answer identically; tests parameterise over them to check that the three
/// implementations agree, and that search works with whichever the index was built with.
pub(crate) enum Mapping {
    Dense,
    Sparse,
    BitVec
}

/// Builds a searcher over [`get_example_proteins`] with an explicit suffix array and mapping.
///
/// Most search tests differ only in these three choices, so they share this rather than
/// re-inlining the three-line construction.
pub(crate) fn example_searcher_with(sa: &[i64], sparseness: u8, mapping: Mapping) -> Searcher<SuffixArray> {
    let proteins = get_example_proteins();
    let stp = match mapping {
        Mapping::Dense => SuffixToProteinMapping::Dense(DenseSuffixToProtein::new(proteins.text())),
        Mapping::Sparse => SuffixToProteinMapping::Sparse(SparseSuffixToProtein::new(proteins.text())),
        Mapping::BitVec => SuffixToProteinMapping::BitVec(BitVecSuffixToProtein::new(proteins.text())),
    };
    Searcher::new(SuffixArray::Original(OriginalSA(sa.to_vec(), sparseness)), proteins, stp)
}

/// The common case: the full suffix array over the example proteins, with a BitVec mapping.
pub(crate) fn example_searcher() -> Searcher<SuffixArray> {
    example_searcher_with(&EXAMPLE_SA_FULL, 1, Mapping::BitVec)
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
pub(crate) fn searcher_over_text(text: &str, sparseness: u8) -> Searcher<SuffixArray> {
    let normalised: Vec<u8> =
        text.bytes().map(|c| if c == b'L' { b'I' } else { c }).collect();

    let mut positions: Vec<i64> = (0..text.len() as i64)
        .filter(|p| (*p as usize).is_multiple_of(sparseness as usize))
        .collect();
    positions.sort_by(|&a, &b| normalised[a as usize..].cmp(&normalised[b as usize..]));

    let proteins = Proteins::new(
        ProteinText::from_string(text),
        text.trim_end_matches('$')
            .split('-')
            .enumerate()
            .map(|(i, _)| Protein {
                uniprot_id: format!("P{i}"),
                taxon_id: (i as u32 + 1) * 10,
                functional_annotations: vec![],
            })
            .collect(),
    );
    let stp = BitVecSuffixToProtein::new(proteins.text());
    Searcher::new(
        SuffixArray::Original(OriginalSA(positions, sparseness)),
        proteins,
        SuffixToProteinMapping::BitVec(stp),
    )
}
