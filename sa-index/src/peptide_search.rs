use rayon::prelude::*;
use sa_mappings::{
    functionality::FunctionalAggregation,
    proteins::Protein
};
use serde::Serialize;

use crate::sa_searcher::{
    SearchAllSuffixesResult,
    Searcher
};

/// Struct representing a collection of `SearchResultWithAnalysis` or `SearchOnlyResult` results
#[derive(Debug, Serialize)]
pub struct OutputData<T: Serialize> {
    pub result: Vec<T>
}

/// Struct representing the search result of the `sequence` in the index, including the analyses
#[derive(Debug, Serialize)]
pub struct SearchResultWithAnalysis {
    pub sequence: String,
    pub taxa: Vec<u32>,
    pub uniprot_accession_numbers: Vec<String>,
    pub fa: Option<FunctionalAggregation>,
    pub cutoff_used: bool
}

/// Struct representing the search result of the `sequence` in the index (without the analyses)
#[derive(Debug, Serialize)]
pub struct SearchOnlyResult {
    pub sequence:    String,
    pub proteins:    Vec<ProteinInfo>,
    pub cutoff_used: bool
}

/// Struct that represents all information known about a certain protein in our database
#[derive(Debug, Serialize)]
pub struct ProteinInfo {
    pub taxon:                  u32,
    pub uniprot_accession:      String,
    pub functional_annotations: Vec<String>
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
    let peptide = peptide.strip_suffix('\n').unwrap_or(peptide).to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sa.sample_rate() as usize {
        return None;
    }

    let suffix_search =
        searcher.search_matching_suffixes(peptide.as_bytes(), cutoff, equate_il);
    let mut cutoff_used = false;
    let suffixes = match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => {
            cutoff_used = true;
            matched_suffixes
        }
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => matched_suffixes,
        SearchAllSuffixesResult::NoMatches => {
            eprintln!("No matches found for peptide: {}", peptide);
            return None;
        }
    };

    let proteins = searcher.retrieve_proteins(&suffixes);

    Some((cutoff_used, proteins))
}

/// Searches the `peptide` in the index multithreaded and retrieves the protein information from the
/// database This does NOT perform any of the analyses, it only retrieves the functional and
/// taxonomic annotations
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
/// Returns Some(SearchOnlyResult) if the peptide has matches
/// Returns None if the peptides does not have any matches, or if the peptide is shorter than the
/// sparseness factor k used in the index
pub fn search_peptide_retrieve_annotations(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equate_il: bool
) -> Option<SearchOnlyResult> {
    let (cutoff_used, proteins) =
        search_proteins_for_peptide(searcher, peptide, cutoff, equate_il)?;

    let annotations = searcher.get_all_functional_annotations(&proteins);

    let mut protein_info: Vec<ProteinInfo> = vec![];
    for (&protein, annotations) in proteins.iter().zip(annotations) {
        protein_info.push(ProteinInfo {
            taxon:                  protein.taxon_id as u32,
            uniprot_accession:      protein.uniprot_id.clone(),
            functional_annotations: annotations
        })
    }

    Some(SearchOnlyResult {
        sequence: peptide.to_string(),
        proteins: protein_info,
        cutoff_used
    })
}

/// Searches the `peptide` in the index multithreaded and performs the taxonomic and functional
/// analyses
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
/// Returns Some(SearchResultWithAnalysis) if the peptide has matches
/// Returns None if the peptides does not have any matches, or if the peptide is shorter than the
/// sparseness factor k used in the index
pub fn analyse_peptide(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equate_il: bool
) -> Option<SearchResultWithAnalysis> {
    let (cutoff_used, proteins) =
        search_proteins_for_peptide(searcher, peptide, cutoff, equate_il)?;

    let mut uniprot_accession_numbers = vec![];
    let mut taxa = vec![];

    for protein in &proteins {
        taxa.push(protein.taxon_id as u32);
        uniprot_accession_numbers.push(protein.uniprot_id.clone());
    }

    let fa = searcher.retrieve_function(&proteins);
    // output the result
    Some(SearchResultWithAnalysis {
        sequence: peptide.to_string(),
        cutoff_used,
        uniprot_accession_numbers,
        taxa,
        fa
    })
}

/// Searches the list of `peptides` in the index multithreaded and performs the functional and
/// taxonomic analyses
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
/// Returns an `OutputData<SearchResultWithAnalysis>` object with the search and analyses results
/// for the peptides
pub fn analyse_all_peptides(
    searcher: &Searcher,
    peptides: &Vec<String>,
    cutoff: usize,
    equate_il: bool
) -> OutputData<SearchResultWithAnalysis> {
    let res: Vec<SearchResultWithAnalysis> = peptides
        .par_iter()
        // calculate the results and remove None's
        .filter_map(|peptide| analyse_peptide(searcher, peptide, cutoff, equate_il))
        .collect();

    OutputData {
        result: res
    }
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
) -> OutputData<SearchOnlyResult> {
    let res: Vec<SearchOnlyResult> = peptides
        .par_iter()
        // calculate the results and remove None's
        .filter_map(|peptide| {
            search_peptide_retrieve_annotations(
                searcher,
                peptide,
                cutoff,
                equate_il
            )
        })
        .collect();

    OutputData {
        result: res
    }
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
    fn test_serialize_output_data() {
        let output_data = OutputData {
            result: vec![1, 2, 3]
        };

        let generated_json = serde_json::to_string(&output_data).unwrap();
        let expected_json = "{\"result\":[1,2,3]}";

        assert_json_eq(&generated_json, expected_json);
    }

    #[test]
    fn test_serialize_search_result_with_analysis() {
        let search_result = SearchResultWithAnalysis {
            sequence: "MSKIAALLPSV".to_string(),
            taxa: vec![1, 2, 3],
            uniprot_accession_numbers: vec!["P12345".to_string(), "P23456".to_string()],
            fa: None,
            cutoff_used: true
        };

        let generated_json = serde_json::to_string(&search_result).unwrap();
        let expected_json = "{\"sequence\":\"MSKIAALLPSV\",\"taxa\":[1,2,3],\"uniprot_accession_numbers\":[\"P12345\",\"P23456\"],\"fa\":null,\"cutoff_used\":true}";

        assert_json_eq(&generated_json, expected_json);
    }

    #[test]
    fn test_serialize_protein_info() {
        let protein_info = ProteinInfo {
            taxon:                  1,
            uniprot_accession:      "P12345".to_string(),
            functional_annotations: vec!["GO:0001234".to_string(), "GO:0005678".to_string()]
        };

        let generated_json = serde_json::to_string(&protein_info).unwrap();
        let expected_json = "{\"taxon\":1,\"uniprot_accession\":\"P12345\",\"functional_annotations\":[\"GO:0001234\",\"GO:0005678\"]}";

        assert_json_eq(&generated_json, expected_json);
    }

    #[test]
    fn test_serialize_search_only_result() {
        let search_result = SearchOnlyResult {
            sequence:    "MSKIAALLPSV".to_string(),
            proteins:    vec![],
            cutoff_used: true
        };

        let generated_json = serde_json::to_string(&search_result).unwrap();
        let expected_json = "{\"sequence\":\"MSKIAALLPSV\",\"proteins\":[],\"cutoff_used\":true}";

        assert_json_eq(&generated_json, expected_json);
    }
}
