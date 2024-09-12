use rayon::prelude::*;
use sa_mappings::proteins::Protein;
use serde::Serialize;

use crate::sa_searcher::{SearchAllSuffixesResult, Searcher};

#[derive(Debug, Serialize)]
pub struct SearchResult {
    pub sequence: String,
    pub proteins: Vec<ProteinInfo>,
    pub cutoff_used: bool
}

/// Struct that represents all information known about a certain protein in our database
#[derive(Debug, Serialize)]
pub struct ProteinInfo {
    pub taxon: u32,
    pub uniprot_accession: String,
    pub functional_annotations: String
}

impl From<&Protein> for ProteinInfo {
    fn from(protein: &Protein) -> Self {
        ProteinInfo {
            taxon: protein.taxon_id,
            uniprot_accession: protein.uniprot_id.clone(),
            functional_annotations: protein.get_functional_annotations()
        }
    }
}

/// Searches the `peptide` in the index multithreaded and retrieves the matching proteins
///
/// # Arguments
/// * `searcher` - The Searcher which contains the protein database
/// * `peptide` - The peptide that is being searched in the index
/// * `cutoff` - The maximum amount of matches we want to process from the index
/// * `equate_il` - Boolean indicating if we want to equate I and L during search
/// * `clean_taxa` - Boolean indicating if we want to filter out proteins that are invalid in the
///   taxonomy
///
/// # Returns
///
/// Returns Some if matches are found.
/// The first argument is true if the cutoff is used, otherwise false
/// The second argument is a list of all matching proteins for the peptide
/// Returns None if the peptides does not have any matches, or if the peptide is shorter than the
/// sparseness factor k used in the index
pub fn search_proteins_for_peptide<'a>(
    searcher: &'a Searcher,
    peptide: &str,
    cutoff: usize,
    equate_il: bool
) -> Option<(bool, Vec<&'a Protein>)> {
    let peptide = peptide.trim_end().to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sa.sample_rate() as usize {
        return None;
    }

    let suffix_search = searcher.search_matching_suffixes(peptide.as_bytes(), cutoff, equate_il);
    let (suffixes, cutoff_used) = match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => Some((matched_suffixes, true)),
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => Some((matched_suffixes, false)),
        SearchAllSuffixesResult::NoMatches => None
    }?;

    let proteins = searcher.retrieve_proteins(&suffixes);

    Some((cutoff_used, proteins))
}

pub fn search_peptide(searcher: &Searcher, peptide: &str, cutoff: usize, equate_il: bool) -> Option<SearchResult> {
    let (cutoff_used, proteins) = search_proteins_for_peptide(searcher, peptide, cutoff, equate_il)?;

    Some(SearchResult {
        sequence: peptide.to_string(),
        proteins: proteins.iter().map(|&protein| protein.into()).collect(),
        cutoff_used
    })
}

/// Searches the list of `peptides` in the index and retrieves all related information about the
/// found proteins This does NOT perform any of the analyses
///
/// # Arguments
/// * `searcher` - The Searcher which contains the protein database
/// * `peptides` - List of peptides we want to search in the index
/// * `cutoff` - The maximum amount of matches we want to process from the index
/// * `equate_il` - Boolean indicating if we want to equate I and L during search
/// * `clean_taxa` - Boolean indicating if we want to filter out proteins that are invalid in the
///   taxonomy
///
/// # Returns
///
/// Returns an `OutputData<SearchOnlyResult>` object with the search results for the peptides
pub fn search_all_peptides(
    searcher: &Searcher,
    peptides: &Vec<String>,
    cutoff: usize,
    equate_il: bool
) -> Vec<SearchResult> {
    peptides
        .par_iter()
        .filter_map(|peptide| search_peptide(searcher, peptide, cutoff, equate_il))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_json_eq(generated_json: &str, expected_json: &str) {
        assert_eq!(
            generated_json.parse::<serde_json::Value>().unwrap(),
            expected_json.parse::<serde_json::Value>().unwrap(),
        );
    }

    #[test]
    fn test_serialize_protein_info() {
        let protein_info = ProteinInfo {
            taxon: 1,
            uniprot_accession: "P12345".to_string(),
            functional_annotations: "GO:0001234;GO:0005678".to_string()
        };

        let generated_json = serde_json::to_string(&protein_info).unwrap();
        let expected_json =
            "{\"taxon\":1,\"uniprot_accession\":\"P12345\",\"functional_annotations\":\"GO:0001234;GO:0005678\"}";

        assert_json_eq(&generated_json, expected_json);
    }

    #[test]
    fn test_serialize_search_result() {
        let search_result = SearchResult {
            sequence: "MSKIAALLPSV".to_string(),
            proteins: vec![],
            cutoff_used: true
        };

        let generated_json = serde_json::to_string(&search_result).unwrap();
        let expected_json = "{\"sequence\":\"MSKIAALLPSV\",\"proteins\":[],\"cutoff_used\":true}";

        assert_json_eq(&generated_json, expected_json);
    }
}
