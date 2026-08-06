//! Shared test fixtures for the `sa_searcher` submodule tests.

use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
use text_compression::ProteinText;

use crate::{
    array::OriginalSA,
    sa_searcher::Searcher,
    suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping},
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
        .filter(|p| *p as usize % sparseness as usize == 0)
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
