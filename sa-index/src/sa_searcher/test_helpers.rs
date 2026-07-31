//! Shared test fixtures for the `sa_searcher` submodule tests.

use sa_mappings::proteins::{Protein, Proteins};
use text_compression::ProteinText;

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
