//! What a search returns.
//!
//! Two enums, both of which distinguish "found nothing" from "found something" as a variant rather
//! than through an empty collection, because the callers act differently on the two: a peptide with
//! no matches is dropped before retrieval instead of being looked up.

/// Enum representing the minimum and maximum bound of the found matches in the suffix array
#[derive(PartialEq, Debug)]
pub enum BoundSearchResult {
    /// No suffix in the array has the search string as a prefix.
    NoMatches,
    /// The half-open SA index range `[min, max)` whose suffixes share that prefix.
    SearchResult((usize, usize))
}

/// Enum representing the matching suffixes after searching a peptide in the suffix array
///
/// Both `MaxMatches` and `SearchResult` indicate found suffixes. The distinction is exact and
/// user-visible — it becomes `cutoff_used` in the response, which a caller reads to decide whether
/// a protein list is exhaustive:
///
/// * `SearchResult` — every match is here. The set is complete.
/// * `MaxMatches` — there were **strictly more** than `max_matches` matches, and this is a sample
///   of exactly `max_matches` of them.
///
/// A set of exactly `max_matches` matches is therefore `SearchResult`, not `MaxMatches`: it is
/// complete, and nothing was dropped. Every producer decides this by collecting one match *past*
/// the cutoff — reaching `max_matches + 1` is what proves the set was truncated — and then hands
/// the accumulator to `SearchAllSuffixesResult::truncated`, which drops that extra element.
#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    /// Nothing matched.
    NoMatches,
    /// More than `max_matches` suffixes matched; this is a sample of exactly `max_matches`.
    MaxMatches(Vec<i64>),
    /// Every matching text position, complete.
    SearchResult(Vec<i64>)
}

impl SearchAllSuffixesResult {
    /// Builds the truncated-sample result from an accumulator that ran one past the cutoff.
    ///
    /// Shared by the scalar and batched searchers rather than duplicated, because the two must
    /// return identical results for identical input and this is exactly the kind of one-element
    /// detail that drifts apart.
    pub(crate) fn truncated(mut suffixes: Vec<i64>, max_matches: usize) -> Self {
        debug_assert!(
            suffixes.len() > max_matches,
            "truncated() called with {} matches, which does not exceed the cutoff of {max_matches}",
            suffixes.len()
        );
        suffixes.truncate(max_matches);
        Self::MaxMatches(suffixes)
    }
}

/// Custom implementation of partialEq for SearchAllSuffixesResult
/// We consider 2 SearchAllSuffixesResult equal if they exist of the same key, and the Vec contains
/// the same values, but the order can be different
impl PartialEq for SearchAllSuffixesResult {
    fn eq(&self, other: &Self) -> bool {
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

#[cfg(test)]
mod tests {
    use super::SearchAllSuffixesResult;

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
}
