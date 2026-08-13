//! Every storage combination must answer identically. **Tests only** — the module declaration in
//! `super` is `#[cfg(test)]`, so nothing here is compiled into a real build.
//!
//! The searcher is generic over three backends, each with an owned and a mapped implementation,
//! and the protein struct is generic over its text backend on top of that — sixteen combinations
//! in all. Which one a binary uses is a build-time choice made in `sa-server`; which one is
//! *correct* is none, individually, because they are all the same index read three different ways.
//!
//! So this asserts the only property that matters across them: same files in, same answers out.
//! Every other test in this crate runs on the fully-owned combination alone, and this is what
//! entitles them to.
//!
//! It costs sixteen monomorphisations of the whole search path, so the fingerprint below stays
//! deliberately small — one `search_all_matching_suffixes` and one `retrieve_proteins` over a
//! twenty-residue text. It is checking that the backends agree, not re-testing search.

use sa_mappings::proteins::{InMemoryProteins, MmapBackedProteins, ProteinsBackend};
use text_compression::{InMemoryProteinText, LoadIndex, MmapBackedProteinText};

use crate::{
    array::{InMemorySA, MmapBackedSA, SuffixArrayBackend},
    sa_searcher::{
        SearchAllSuffixesResult,
        test_utils::{EXAMPLE_SA_FULL, EXAMPLE_TEXT, Mapping, build_searcher, example_protein_list}
    },
    suffix_to_protein_index::{
        InMemorySuffixToProteinMapping, MmapBackedSuffixToProteinMapping, SuffixToProteinMappingBackend
    }
};

/// One row per (mapping representation, peptide, equate_il, tryptic): the result-variant tag, the
/// sorted matching suffixes, and the `(taxon_id, uniprot_id)` pairs they retrieve.
///
/// Suffixes are sorted because the batched path may emit them in a different order than the scalar
/// one; that ordering is not part of the contract, but the set is.
type Fingerprint = Vec<(&'static str, &'static str, bool, bool, &'static str, Vec<i64>, Vec<(u32, String)>)>;

/// Peptides over [`EXAMPLE_TEXT`] (`"AI-CLACVAA-AC-KCRLY$"`), chosen to hit every branch a backend
/// could get wrong: a hit in several proteins, a single hit, one that only matches with I/L
/// equated, one that only survives without the tryptic filter, and one that matches nothing.
const PEPTIDES: [&[u8]; 5] = [b"AC", b"CLACVAA", b"AL", b"KCRLY", b"WWW"];

/// Searches every peptide against a freshly built searcher of the requested backend combination.
fn fingerprint<SA, P, STPM>() -> Fingerprint
where
    SA: SuffixArrayBackend + LoadIndex,
    P: ProteinsBackend + LoadIndex,
    STPM: SuffixToProteinMappingBackend + LoadIndex
{
    let mut rows = Fingerprint::new();

    // The three suffix-to-protein representations are a runtime choice, not a type parameter, so
    // covering them here is free — it does not multiply the monomorphisations.
    for (mapping_name, mapping) in [("dense", Mapping::Dense), ("sparse", Mapping::Sparse), ("bitvec", Mapping::BitVec)]
    {
        let searcher =
            build_searcher::<SA, P, STPM>(EXAMPLE_TEXT, example_protein_list(), EXAMPLE_SA_FULL.to_vec(), 1, mapping);

        for equate_il in [false, true] {
            for tryptic in [false, true] {
                for peptide in PEPTIDES {
                    let (tag, mut suffixes) = match searcher
                        .search_all_matching_suffixes(&[peptide], usize::MAX, equate_il, tryptic)
                        .remove(0)
                    {
                        SearchAllSuffixesResult::NoMatches => ("none", Vec::new()),
                        SearchAllSuffixesResult::MaxMatches(s) => ("max", s),
                        SearchAllSuffixesResult::SearchResult(s) => ("found", s)
                    };
                    suffixes.sort();

                    let mut proteins: Vec<(u32, String)> = searcher
                        .retrieve_proteins(&suffixes)
                        .iter()
                        .map(|p| (p.taxon_id, p.uniprot_id.to_string()))
                        .collect();
                    proteins.sort();

                    rows.push((
                        mapping_name,
                        std::str::from_utf8(peptide).unwrap(),
                        equate_il,
                        tryptic,
                        tag,
                        suffixes,
                        proteins
                    ));
                }
            }
        }
    }

    rows
}

// The four ways `proteins.bin` can be read. The text and the metadata are independent axes over
// one file, which is why there are four rather than two.
type OwnedProteins = InMemoryProteins<InMemoryProteinText>;
type OwnedMetaMappedText = InMemoryProteins<MmapBackedProteinText>;
type MappedMetaOwnedText = MmapBackedProteins<InMemoryProteinText>;
type MappedProteins = MmapBackedProteins<MmapBackedProteinText>;

#[test]
fn every_backend_combination_returns_identical_results() {
    let expected = fingerprint::<InMemorySA, OwnedProteins, InMemorySuffixToProteinMapping>();

    // An agreement test over a fingerprint that stopped distinguishing anything would pass
    // forever, so check the reference has both kinds of row before comparing against it.
    assert!(
        expected.iter().any(|row| !row.6.is_empty()),
        "the fixture retrieved no proteins — it can no longer tell the backends apart"
    );
    assert!(expected.iter().any(|row| row.4 == "none"), "the fixture has no miss to check");

    // Spelled out rather than generated: sixteen lines that a failure can name, against a macro
    // whose expansion nothing could read.
    let combinations: Vec<(&str, Fingerprint)> = vec![
        (
            "sa=owned  proteins=owned      mapping=mapped",
            fingerprint::<InMemorySA, OwnedProteins, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=owned/map  mapping=owned ",
            fingerprint::<InMemorySA, OwnedMetaMappedText, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=owned/map  mapping=mapped",
            fingerprint::<InMemorySA, OwnedMetaMappedText, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=map/owned  mapping=owned ",
            fingerprint::<InMemorySA, MappedMetaOwnedText, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=map/owned  mapping=mapped",
            fingerprint::<InMemorySA, MappedMetaOwnedText, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=mapped     mapping=owned ",
            fingerprint::<InMemorySA, MappedProteins, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=owned  proteins=mapped     mapping=mapped",
            fingerprint::<InMemorySA, MappedProteins, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=owned      mapping=owned ",
            fingerprint::<MmapBackedSA, OwnedProteins, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=owned      mapping=mapped",
            fingerprint::<MmapBackedSA, OwnedProteins, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=owned/map  mapping=owned ",
            fingerprint::<MmapBackedSA, OwnedMetaMappedText, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=owned/map  mapping=mapped",
            fingerprint::<MmapBackedSA, OwnedMetaMappedText, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=map/owned  mapping=owned ",
            fingerprint::<MmapBackedSA, MappedMetaOwnedText, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=map/owned  mapping=mapped",
            fingerprint::<MmapBackedSA, MappedMetaOwnedText, MmapBackedSuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=mapped     mapping=owned ",
            fingerprint::<MmapBackedSA, MappedProteins, InMemorySuffixToProteinMapping>()
        ),
        (
            "sa=mapped proteins=mapped     mapping=mapped",
            fingerprint::<MmapBackedSA, MappedProteins, MmapBackedSuffixToProteinMapping>()
        ),
    ];

    for (name, rows) in combinations {
        assert_eq!(rows, expected, "{name} disagrees with the fully-owned combination");
    }
}
