//! The public entry point: peptides in, proteins out.
//!
//! Wraps the searcher with the parts a caller actually wants — normalising the query, dropping
//! peptides too short to be searchable, resolving suffixes to proteins, and decoding annotations.
//! `sa-server` calls [`search_all_peptides`] and serialises the result directly.

use rayon::prelude::*;
use sa_mappings::proteins::{ProteinRef, ProteinsBackend};
use serde::Serialize;

use crate::{
    array::SuffixArrayBackend,
    sa_searcher::{DEFAULT_MLP_BATCH, SearchAllSuffixesResult, Searcher},
    suffix_to_protein_index::SuffixToProteinMappingBackend
};

/// Everything found for one peptide. Serialised straight to JSON by the server.
#[derive(Debug, Serialize)]
pub struct SearchResult {
    /// The peptide as the caller wrote it, before normalisation.
    pub sequence: String,
    /// Every protein containing the peptide.
    pub proteins: Vec<ProteinInfo>,
    /// Whether the match cutoff was hit, meaning `proteins` is a truncated sample rather than the
    /// complete set.
    pub cutoff_used: bool
}

/// One matching protein, with its annotations already decoded.
#[derive(Debug, Serialize)]
pub struct ProteinInfo {
    /// NCBI taxon id.
    pub taxon: u32,
    /// UniProt accession, e.g. `P12345`.
    pub uniprot_accession: String,
    /// Functional annotations in their `GO:`/`EC:`/`IPR:` text form.
    pub functional_annotations: String
}

impl From<ProteinRef<'_>> for ProteinInfo {
    fn from(protein: ProteinRef<'_>) -> Self {
        ProteinInfo {
            taxon: protein.taxon_id,
            uniprot_accession: protein.uniprot_id.to_string(),
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
/// * `tryptic` - Boolean indicating if we only want tryptic matches.
///
/// # Returns
///
/// Returns Some if matches are found.
/// The first argument is true if the cutoff is used, otherwise false
/// The second argument is a list of all matching proteins for the peptide
/// Returns None if the peptides does not have any matches, or if the peptide is shorter than the
/// sparseness factor k used in the index
pub fn search_proteins_for_peptide<
    'a,
    SA: SuffixArrayBackend,
    P: ProteinsBackend,
    STPM: SuffixToProteinMappingBackend
>(
    searcher: &'a Searcher<SA, P, STPM>,
    peptide: &str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<(bool, Vec<ProteinRef<'a>>)> {
    let peptide = peptide.trim_end().to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sa.sample_rate() as usize {
        return None;
    }

    let suffix_search = searcher.search_matching_suffixes(peptide.as_bytes(), cutoff, equate_il, tryptic);
    let (suffixes, cutoff_used) = match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => Some((matched_suffixes, true)),
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => Some((matched_suffixes, false)),
        SearchAllSuffixesResult::NoMatches => None
    }?;

    let proteins = searcher.retrieve_proteins(&suffixes);

    Some((cutoff_used, proteins))
}

/// Searches one peptide and packages the result.
///
/// The single-peptide path. `sa-server` uses [`search_all_peptides`] instead, which batches the
/// suffix searches across the whole request; this one is kept for callers with a single query.
///
/// Returns `None` when the peptide has no matches or is shorter than the sparseness factor.
pub fn search_peptide<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptide: &str,
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Option<SearchResult> {
    let (cutoff_used, proteins) = search_proteins_for_peptide(searcher, peptide, cutoff, equate_il, tryptic)?;

    Some(SearchResult {
        sequence: peptide.to_string(),
        proteins: proteins.into_iter().map(|protein| protein.into()).collect(),
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
/// * `tryptic` - Boolean indicating if we only want tryptic matches.
///
/// # Returns
///
/// Returns an `OutputData<SearchOnlyResult>` object with the search results for the peptides
pub fn search_all_peptides<SA: SuffixArrayBackend, P: ProteinsBackend, STPM: SuffixToProteinMappingBackend>(
    searcher: &Searcher<SA, P, STPM>,
    peptides: &[String],
    cutoff: usize,
    equate_il: bool,
    tryptic: bool
) -> Vec<SearchResult> {
    let sample_rate = searcher.sa.sample_rate() as usize;

    // Normalise once and keep only searchable peptides (anything shorter than the sample rate
    // cannot be searched), remembering each peptide's original index so the returned `sequence`
    // stays verbatim — matching the previous per-peptide behaviour.
    let prepared: Vec<(usize, String)> = peptides
        .iter()
        .enumerate()
        .map(|(index, peptide)| (index, peptide.trim_end().to_uppercase()))
        .filter(|(_, normalised)| normalised.len() >= sample_rate)
        .collect();

    let byte_peptides: Vec<&[u8]> = prepared.iter().map(|(_, peptide)| peptide.as_bytes()).collect();

    // One batched (MLP) search over the whole list — the same code path the benchmark measures
    // and the only place the batch size is chosen.
    let suffix_results =
        searcher.search_all_matching_suffixes(&byte_peptides, cutoff, equate_il, tryptic, DEFAULT_MLP_BATCH);

    // Retrieve the proteins for each hit and build the result, dropping peptides with no matches
    // (preserves the previous filter_map semantics and result ordering).
    //
    // Retrieval is per peptide by design: a cross-query batched variant was built and measured
    // (run3) and moved throughput by a median of +1.7%, never clearing the noise floor. See the
    // module comment on `sa_searcher::orchestrate` for the full result.
    prepared
        .par_iter()
        .zip(suffix_results.par_iter())
        .filter_map(|((original_index, _), suffix_result)| {
            let (suffixes, cutoff_used) = match suffix_result {
                SearchAllSuffixesResult::MaxMatches(matched) => (matched, true),
                SearchAllSuffixesResult::SearchResult(matched) => (matched, false),
                SearchAllSuffixesResult::NoMatches => return None
            };

            let proteins = searcher.retrieve_proteins(suffixes);

            Some(SearchResult {
                sequence: peptides[*original_index].to_string(),
                proteins: proteins.into_iter().map(|protein| protein.into()).collect(),
                cutoff_used
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(not(feature = "mmap"))]
    #[test]
    fn test_search_all_peptides_matches_scalar_reference() {
        use sa_mappings::proteins::{Protein, Proteins, ProteinsBackend as _};
        use text_compression::ProteinText;

        use crate::{
            SuffixArray,
            array::OriginalSA,
            suffix_to_protein_index::{BitVecSuffixToProtein, SuffixToProteinMapping}
        };

        // Example DB "AI-CLACVAA-AC-KCRLY$", sample rate 1 (so single characters are searchable).
        let text = ProteinText::from_string("AI-CLACVAA-AC-KCRLY$");
        let proteins = Proteins::new(text, vec![
            Protein {
                uniprot_id: "P0".to_string(),
                taxon_id: 10,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: "P1".to_string(),
                taxon_id: 20,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: "P2".to_string(),
                taxon_id: 30,
                functional_annotations: vec![]
            },
            Protein {
                uniprot_id: "P3".to_string(),
                taxon_id: 40,
                functional_annotations: vec![]
            },
        ]);
        let sa = SuffixArray::Original(OriginalSA(
            vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18],
            1
        ));
        let stp = BitVecSuffixToProtein::new(proteins.text());
        let searcher = Searcher::new(sa, proteins, SuffixToProteinMapping::BitVec(stp));

        // Mixed case (normalisation), whitespace + empty (trim/length filter), a no-match ("ZZZ").
        let peptides: Vec<String> =
            ["A", "ai", "CLA", "kcrly", "VAA", "ZZZ", "", "  ", "AC"].iter().map(|s| s.to_string()).collect();

        let as_json = |v: &[SearchResult]| serde_json::to_string(v).unwrap();
        for equate_il in [true, false] {
            for tryptic in [true, false] {
                // Batched production path vs the scalar per-peptide reference (the old behaviour).
                let got = search_all_peptides(&searcher, &peptides, 1000, equate_il, tryptic);
                let reference: Vec<SearchResult> =
                    peptides.iter().filter_map(|p| search_peptide(&searcher, p, 1000, equate_il, tryptic)).collect();
                assert_eq!(as_json(&got), as_json(&reference), "il={} tryptic={}", equate_il, tryptic);
            }
        }
    }

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
