use std::{cmp::min, ops::Deref};

use sa_mappings::proteins::{Protein, Proteins, SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use text_compression::ProteinTextSlice;

use crate::{
    sa_searcher::BoundSearch::{Maximum, Minimum},
    suffix_to_protein_index::{DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinIndex},
    Nullable, SuffixArray
};

/// Enum indicating if we are searching for the minimum, or maximum bound in the suffix array
#[derive(Clone, Copy, PartialEq)]
enum BoundSearch {
    Minimum,
    Maximum
}

/// Enum representing the minimum and maximum bound of the found matches in the suffix array
#[derive(PartialEq, Debug)]
pub enum BoundSearchResult {
    NoMatches,
    SearchResult((usize, usize))
}

/// Enum representing the matching suffixes after searching a peptide in the suffix array
/// Both the MaxMatches and SearchResult indicate found suffixes, but MaxMatches is used when the
/// cutoff is reached.
#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    NoMatches,
    MaxMatches(Vec<i64>),
    SearchResult(Vec<i64>)
}

/// Custom implementation of partialEq for SearchAllSuffixesResult
/// We consider 2 SearchAllSuffixesResult equal if they exist of the same key, and the Vec contains
/// the same values, but the order can be different
impl PartialEq for SearchAllSuffixesResult {
    fn eq(&self, other: &Self) -> bool {
        /// Returns true if `arr1` and `arr2` contains the same elements, the order of the elements
        /// is ignored
        ///
        /// # Arguments
        /// * `arr1` - The first array used in the comparison
        /// * `arr2` - The second array used in the comparison
        ///
        /// # Returns
        ///
        /// Returns true if arr1 and arr2 contains the same elements, the order of the elements is
        /// ignored
        fn array_eq_unordered(arr1: &[i64], arr2: &[i64]) -> bool {
            let mut arr1_copy = arr1.to_owned();
            let mut arr2_copy = arr2.to_owned();

            arr1_copy.sort();
            arr2_copy.sort();

            arr1_copy == arr2_copy
        }

        match (self, other) {
            (SearchAllSuffixesResult::MaxMatches(arr1), SearchAllSuffixesResult::MaxMatches(arr2)) => {
                array_eq_unordered(arr1, arr2)
            }
            (SearchAllSuffixesResult::SearchResult(arr1), SearchAllSuffixesResult::SearchResult(arr2)) => {
                array_eq_unordered(arr1, arr2)
            }
            (SearchAllSuffixesResult::NoMatches, SearchAllSuffixesResult::NoMatches) => true,
            _ => false
        }
    }
}

pub struct SparseSearcher(Searcher);

impl SparseSearcher {
    pub fn new(sa: SuffixArray, proteins: Proteins) -> Self {
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));
        Self(searcher)
    }
}

impl Deref for SparseSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub struct DenseSearcher(Searcher);

impl DenseSearcher {
    pub fn new(sa: SuffixArray, proteins: Proteins) -> Self {
        let suffix_index_to_protein = DenseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));
        Self(searcher)
    }
}

impl Deref for DenseSearcher {
    type Target = Searcher;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Struct that contains all the elements needed to search a peptide in the suffix array
/// This struct also contains all the functions used for search
///
/// # Arguments
/// * `sa` - The sparse suffix array representing the protein database
/// * `sparseness_factor` - The sparseness factor used by the suffix array
/// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix is
///   part of
/// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the
///   taxonomic analysis provided by Unipept
/// * `function_aggregator` - Object used to retrieve the functional annotations and to calculate
///   the functional analysis provided by Unipept
pub struct Searcher {
    pub sa: SuffixArray,
    pub proteins: Proteins,
    pub suffix_index_to_protein: Box<dyn SuffixToProteinIndex>
}

impl Searcher {
    /// Creates a new Searcher object
    ///
    /// # Arguments
    /// * `sa` - The sparse suffix array representing the protein database
    /// * `sparseness_factor` - The sparseness factor used by the suffix array
    /// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix
    ///   is part of
    /// * `proteins` - List of all the proteins where the suffix array is build on
    /// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the
    ///   taxonomic analysis provided by Unipept
    /// * `function_aggregator` - Object used to retrieve the functional annotations and to
    ///   calculate the functional analysis provided by Unipept
    ///
    /// # Returns
    ///
    /// Returns a new Searcher object
    pub fn new(sa: SuffixArray, proteins: Proteins, suffix_index_to_protein: Box<dyn SuffixToProteinIndex>) -> Self {
        Self { sa, proteins, suffix_index_to_protein }
    }

    /// Compares the `search_string` to the `suffix`
    /// During search this function performs extra logic since the suffix array is build with I ==
    /// L, while ` self.proteins.input_string` is the original text where I != L
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide being searched in the suffix array
    /// * `suffix` - The current suffix from the suffix array we are comparing with in the binary
    ///   search
    /// * `skip` - How many characters we can skip in the comparison because we already know these
    ///   match
    /// * `bound` - Indicates if we are searching for the min of max bound
    ///
    /// # Returns
    ///
    /// The first argument is true if `bound` == `Minimum` and `search_string` <= `suffix` or if
    /// `bound` == `Maximum` and `search_string` >= `suffix` The second argument indicates how
    /// far the `suffix` and `search_string` matched
    fn compare(&self, search_string: &[u8], suffix: i64, skip: usize, bound: BoundSearch) -> (bool, usize) {
        let mut index_in_suffix = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;

        // Depending on if we are searching for the min of max bound our condition is different
        let condition_check = match bound {
            Minimum => |a: u8, b: u8| a < b,
            Maximum => |a: u8, b: u8| a > b
        };

        // match as long as possible
        while index_in_search_string < search_string.len()
            && index_in_suffix < self.proteins.text.len()
            && (search_string[index_in_search_string] == self.proteins.text.get(index_in_suffix)
                || (search_string[index_in_search_string] == b'L' && self.proteins.text.get(index_in_suffix) == b'I')
                || (search_string[index_in_search_string] == b'I' && self.proteins.text.get(index_in_suffix) == b'L'))
        {
            index_in_suffix += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty
        // search string should not be found)
        if !search_string.is_empty() {
            if index_in_search_string == search_string.len() {
                is_cond_or_equal = true
            } else if index_in_suffix < self.proteins.text.len() {
                // in our index every L was replaced by a I, so we need to replace them if we want
                // to search in the right direction
                let peptide_char = if search_string[index_in_search_string] == b'L' {
                    b'I'
                } else {
                    search_string[index_in_search_string]
                };

                let protein_char = if self.proteins.text.get(index_in_suffix) == b'L' {
                    b'I'
                } else {
                    self.proteins.text.get(index_in_suffix)
                };

                is_cond_or_equal = condition_check(peptide_char, protein_char);
            }
        }

        (is_cond_or_equal, index_in_search_string)
    }

    /// Searches for the minimum or maximum bound for a string in the suffix array
    ///
    /// # Arguments
    /// * `bound` - Indicates if we are searching the minimum or maximum bound
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// The first argument is true if a match was found
    /// The second argument indicates the index of the minimum or maximum bound for the match
    /// (depending on `bound`)
    fn binary_search_bound(&self, bound: BoundSearch, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last
        // iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(center), skip, bound);

            found |= lcp_center == search_string.len();

            // update the left and right bound, depending on if we are searching the min or max
            // bound
            if retval && bound == Minimum || !retval && bound == Maximum {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }

        // handle edge case to search at index 0
        if right == 1 && left == 0 {
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(0), min(lcp_left, lcp_right), bound);

            found |= lcp_center == search_string.len();

            if bound == Minimum && retval {
                right = 0;
            }
        }

        match bound {
            Minimum => (found, right),
            Maximum => (found, left)
        }
    }

    /// Searches for the minimum and maximum bound for a string in the suffix array
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// Returns the minimum and maximum bound of all matches in the suffix array, or `NoMatches` if
    /// no matches were found
    pub fn search_bounds(&self, search_string: &[u8]) -> BoundSearchResult {
        let (found_min, min_bound) = self.binary_search_bound(Minimum, search_string);

        if !found_min {
            return BoundSearchResult::NoMatches;
        }

        let (_, max_bound) = self.binary_search_bound(Maximum, search_string);

        BoundSearchResult::SearchResult((min_bound, max_bound + 1))
    }

    /// Searches for the suffixes matching a search string
    /// During search I and L can be equated
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    /// * `max_matches` - The maximum amount of matches processed, if more matches are found we
    ///   don't process them
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    /// * `tryptic` - Boolean indicating if we only want tryptic matches.
    ///
    /// # Returns
    ///
    /// Returns all the matching suffixes
    #[inline]
    pub fn search_matching_suffixes(
        &self,
        search_string: &[u8],
        max_matches: usize,
        equate_il: bool,
        tryptic: bool
    ) -> SearchAllSuffixesResult {
        let mut matching_suffixes: Vec<i64> = vec![];
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        let mut skip: usize = 0;
        while skip < self.sa.sample_rate() as usize {
            let mut il_locations_start = 0;
            while il_locations_start < il_locations.len() && il_locations[il_locations_start] < skip {
                il_locations_start += 1;
            }
            let il_locations_current_suffix = &il_locations[il_locations_start..];
            let current_search_string_prefix = &search_string[..skip];
            let current_search_string_suffix = &search_string[skip..];
            let search_bound_result = self.search_bounds(&search_string[skip..]);
            // if the shorter part is matched, see if what goes before the matched suffix matches
            // the unmatched part of the prefix
            if let BoundSearchResult::SearchResult((min_bound, max_bound)) = search_bound_result {
                // try all the partially matched suffixes and store the matching suffixes in an
                // array (stop when our max number of matches is reached)
                let mut sa_index = min_bound;
                while sa_index < max_bound {
                    let suffix = self.sa.get(sa_index) as usize;

                    if suffix >= skip {
                        let match_start = suffix - skip;
                        let match_end = suffix + search_string.len() - skip;

                        // filter away matches where I was wrongfully equalized to L, and check the
                        // unmatched prefix when I and L equalized, we only need to
                        // check the prefix, not the whole match, when the prefix is 0, we don't need to
                        // check at all
                        if (skip == 0
                            || Self::check_prefix(
                                current_search_string_prefix,
                                ProteinTextSlice::new(&self.proteins.text, match_start, suffix),
                                equate_il
                            ))
                            && Self::check_suffix(
                                skip,
                                il_locations_current_suffix,
                                current_search_string_suffix,
                                ProteinTextSlice::new(&self.proteins.text, suffix, match_end),
                                equate_il
                            )
                            && (!tryptic
                                || ((self.check_start_of_protein(match_start) || self.check_tryptic_cut(match_start))
                                    && (self.check_end_of_protein(match_end) || self.check_tryptic_cut(match_end))))
                        {
                            matching_suffixes.push((suffix - skip) as i64);

                            // return if max number of matches is reached
                            if matching_suffixes.len() >= max_matches {
                                return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                            }
                        }
                    }

                    sa_index += 1;
                }
            }
            skip += 1;
        }

        if matching_suffixes.is_empty() {
            SearchAllSuffixesResult::NoMatches
        } else {
            SearchAllSuffixesResult::SearchResult(matching_suffixes)
        }
    }

    /// Check if a cut is the start of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the start of a protein.
    #[inline]
    fn check_start_of_protein(&self, cut_index: usize) -> bool {
        cut_index == 0 || self.proteins.text.get(cut_index - 1) == SEPARATION_CHARACTER
    }

    /// Check if a cut is the end of a protein.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is at the end of a protein.
    #[inline]
    fn check_end_of_protein(&self, cut_index: usize) -> bool {
        self.proteins.text.get(cut_index) == TERMINATION_CHARACTER
            || self.proteins.text.get(cut_index) == SEPARATION_CHARACTER
    }

    /// Check if a cut is a tryptic cut, so check if the amino acid preceding the cut is K or R and the amino acid at the cut is not P.
    ///
    /// # Arguments
    /// * `cut_index` - The index of the cut in the text of proteins.
    ///
    /// # Returns
    ///
    /// Returns true if the cut is a tryptic cut.
    #[inline]
    fn check_tryptic_cut(&self, cut_index: usize) -> bool {
        (self.proteins.text.get(cut_index - 1) == b'K' || self.proteins.text.get(cut_index - 1) == b'R')
            && self.proteins.text.get(cut_index) != b'P'
    }

    /// Returns true of the prefixes are the same
    /// if `equate_il` is set to true, L and I are considered the same
    ///
    /// # Arguments
    /// * `search_string_prefix` - The unchecked prefix of the string/peptide that is searched
    /// * `index_prefix` - The unchecked prefix from the protein from the suffix array
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string_prefix` and `index_prefix` are considered the same, otherwise
    /// false
    #[inline]
    fn check_prefix(search_string_prefix: &[u8], index_prefix: ProteinTextSlice, equate_il: bool) -> bool {
        index_prefix.equals_slice(search_string_prefix, equate_il)
    }

    /// Returns true of the search_string and index_string are equal
    /// This is automatically true if `equate_il` is set to true, since there matched during
    /// search where I = L If `equate_il` is set to false, we need to check if the I and
    /// L locations have the same character
    ///
    /// # Arguments
    /// * `skip` - The used skip factor during the search iteration
    /// * `il_locations` - The locations of the I's and L's in the **original** peptide
    /// * `search_string` - The peptide that is being searched, but already with the skipped prefix
    ///   removed from it
    /// * `index_string` - The suffix that search_string matches with when I and L were equalized
    ///   during search
    /// * `equate_il` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string` and `index_string` are considered the same, otherwise false
    fn check_suffix(
        skip: usize,
        il_locations: &[usize],
        search_string: &[u8],
        text_slice: ProteinTextSlice,
        equate_il: bool
    ) -> bool {
        if equate_il { true } else { text_slice.check_il_locations(skip, il_locations, search_string) }
    }

    /// Returns all the proteins that correspond with the provided suffixes
    ///
    /// # Arguments
    /// * `suffixes` - List of suffix indices
    ///
    /// # Returns
    ///
    /// Returns the proteins that every suffix is a part of
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &Vec<i64>) -> Vec<&Protein> {
        let mut res = vec![];
        for &suffix in suffixes {
            let protein_index = self.suffix_index_to_protein.suffix_to_protein(suffix);
            if !protein_index.is_null() {
                res.push(&self.proteins[protein_index as usize]);
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use sa_mappings::proteins::{Protein, Proteins};
    use text_compression::ProteinText;

    use crate::{
        sa_searcher::{BoundSearchResult, SearchAllSuffixesResult, Searcher},
        suffix_to_protein_index::SparseSuffixToProtein,
        SuffixArray
    };

    #[test]
    fn test_partial_eq_search_all_suffixes_result() {
        let search_all_suffixes_result_1 = SearchAllSuffixesResult::SearchResult(vec![1, 2, 3]);
        let search_all_suffixes_result_2 = SearchAllSuffixesResult::SearchResult(vec![3, 2, 1]);
        let search_all_suffixes_result_3 = SearchAllSuffixesResult::SearchResult(vec![1, 2, 4]);
        let search_all_suffixes_result_4 = SearchAllSuffixesResult::MaxMatches(vec![1, 2, 3]);
        let search_all_suffixes_result_5 = SearchAllSuffixesResult::MaxMatches(vec![3, 2, 1]);
        let search_all_suffixes_result_6 = SearchAllSuffixesResult::MaxMatches(vec![1, 2, 4]);
        let search_all_suffixes_result_7 = SearchAllSuffixesResult::NoMatches;
        let search_all_suffixes_result_8 = SearchAllSuffixesResult::NoMatches;

        assert_eq!(search_all_suffixes_result_1, search_all_suffixes_result_2);
        assert_ne!(search_all_suffixes_result_1, search_all_suffixes_result_3);
        assert_eq!(search_all_suffixes_result_4, search_all_suffixes_result_5);
        assert_ne!(search_all_suffixes_result_4, search_all_suffixes_result_6);
        assert_eq!(search_all_suffixes_result_7, search_all_suffixes_result_8);
        assert_ne!(search_all_suffixes_result_1, search_all_suffixes_result_7);
        assert_ne!(search_all_suffixes_result_4, search_all_suffixes_result_7);
    }

    fn get_example_proteins() -> Proteins {
        let input_string = "AI-CLACVAA-AC-KCRLY$";
        let text = ProteinText::from_string(input_string);

        Proteins {
            text,
            proteins: vec![
                Protein {
                    uniprot_id: String::new(),
                    taxon_id: 0,
                    functional_annotations: vec![]
                },
                Protein {
                    uniprot_id: String::new(),
                    taxon_id: 0,
                    functional_annotations: vec![]
                },
                Protein {
                    uniprot_id: String::new(),
                    taxon_id: 0,
                    functional_annotations: vec![]
                },
                Protein {
                    uniprot_id: String::new(),
                    taxon_id: 0,
                    functional_annotations: vec![]
                },
            ]
        }
    }

    #[test]
    fn test_search_simple() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1);

        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

        // search bounds 'A'
        let bounds_res = searcher.search_bounds(&[b'A']);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((4, 9)));

        // search bounds '$'
        let bounds_res = searcher.search_bounds(&[b'$']);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((0, 1)));

        // search bounds 'AC'
        let bounds_res = searcher.search_bounds(&[b'A', b'C']);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((6, 8)));
    }

    #[test]
    fn test_search_sparse() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![9, 0, 3, 12, 15, 6, 18], 3);

        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

        // search suffix 'VAA'
        let found_suffixes = searcher.search_matching_suffixes(&[b'V', b'A', b'A'], usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![7]));

        // search suffix 'AC'
        let found_suffixes = searcher.search_matching_suffixes(&[b'A', b'C'], usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![5, 11]));
    }

    #[test]
    fn test_il_equality() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 4, 17, 14, 16, 7, 18], 1);

        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

        let bounds_res = searcher.search_bounds(&[b'I']);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((13, 16)));

        // search bounds 'RIZ' with equal I and L
        let bounds_res = searcher.search_bounds(&[b'R', b'I', b'Y']);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult((17, 18)));
    }

    #[test]
    fn test_il_equality_sparse() {
        let proteins = get_example_proteins();
        let sa = SuffixArray::Original(vec![9, 0, 3, 12, 15, 6, 18], 3);

        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sa, proteins, Box::new(suffix_index_to_protein));

        // search bounds 'RIZ' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(&[b'R', b'I', b'Y'], usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![16]));

        // search bounds 'RIZ' without equal I and L
        let found_suffixes = searcher.search_matching_suffixes(&[b'R', b'I', b'Y'], usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::NoMatches);
    }

    // test edge case where an I or L is the first index in the sparse SA.
    #[test]
    fn test_l_first_index_in_sa() {
        let input_string = "LMPYY$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![0, 2, 4], 2);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(&[b'I', b'M'], usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0]));
    }

    #[test]
    fn test_il_missing_matches() {
        let input_string = "AAILLL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![6, 0, 1, 5, 4, 3, 2], 1);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        let found_suffixes = searcher.search_matching_suffixes(&[b'I'], usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![2, 3, 4, 5]));
    }

    #[test]
    fn test_il_duplication() {
        let input_string = "IIIILL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![6, 5, 4, 3, 2, 1, 0], 1);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        let found_suffixes = searcher.search_matching_suffixes(&[b'I', b'I'], usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4]));
    }

    #[test]
    fn test_il_suffix_check() {
        let input_string = "IIIILL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![6, 4, 2, 0], 2);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        // search all places where II is in the string IIIILL, but with a sparse SA
        // this way we check if filtering the suffixes works as expected
        let found_suffixes = searcher.search_matching_suffixes(&[b'I', b'I'], usize::MAX, false, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2]));
    }

    #[test]
    fn test_il_duplication2() {
        let input_string = "IILLLL$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![6, 5, 4, 3, 2, 1, 0], 1);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(&[b'I', b'I'], usize::MAX, true, false);
        assert_eq!(found_suffixes, SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4]));
    }

    #[test]
    fn test_tryptic_search() {
        let input_string = "PAA-AAKPKAPAA$";
        let text = ProteinText::from_string(input_string);

        let proteins = Proteins {
            text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                taxon_id: 0,
                functional_annotations: vec![]
            }]
        };

        let sparse_sa = SuffixArray::Original(vec![13, 3, 12, 11, 1, 4, 2, 5, 9, 8, 6, 10, 0, 7], 1);
        let suffix_index_to_protein = SparseSuffixToProtein::new(&proteins.text);
        let searcher = Searcher::new(sparse_sa, proteins, Box::new(suffix_index_to_protein));

        let found_suffixes_1 = searcher.search_matching_suffixes(&[b'P', b'A', b'A'], usize::MAX, false, true);
        assert_eq!(found_suffixes_1, SearchAllSuffixesResult::SearchResult(vec![0]));

        let found_suffixes_2 = searcher.search_matching_suffixes(&[b'A', b'P', b'A', b'A'], usize::MAX, false, true);
        assert_eq!(found_suffixes_2, SearchAllSuffixesResult::SearchResult(vec![9]));
    }
}
