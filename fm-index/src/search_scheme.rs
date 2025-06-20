//! Search scheme representation for approximate pattern matching.
//!
//! This module provides structures to define and iterate over search
//! strategies used in approximate string matching. Each search strategy
//! consists of a sequence of parts (pattern segments), and constraints
//! on mismatches for each segment. This is particularly useful for bidirectional
//! FM-index searching where flexible matching strategies are required.

use std::error::Error;
use std::fs;
use std::path::Path;

/// A single search configuration within a search scheme.
///
/// Represents a search over a pattern divided into parts, specifying:
/// - the order of part traversal,
/// - minimum and maximum allowed mismatches per part.
pub struct Search {
    pub order: Vec<u8>,
    pub min_mismatches: Vec<u8>,
    pub max_mismatches: Vec<u8>,
}

/// Iterator over the parts in a `Search`.
///
/// Each item contains:
/// - the part index,
/// - the minimum mismatches allowed for that part,
/// - the maximum mismatches allowed for that part.
pub struct SearchIter<'a> {
    search: &'a Search,
    index: usize,
}

impl<'a> Iterator for SearchIter<'a> {
    type Item = (u8, u8, u8);

    /// Returns the next `(part, min_mismatches, max_mismatches)` tuple.
    ///
    /// The iterator ends after all parts in the search have been visited.
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.search.order.len() {
            let i = self.index;
            self.index += 1;
            Some((
                self.search.order[i],
                self.search.min_mismatches[i],
                self.search.max_mismatches[i],
            ))
        } else {
            None
        }
    }
}

impl Search {

    /// Constructs a new `Search` instance.
    ///
    /// A `Search` defines the order in which pattern segments are searched, along with
    /// lower and upper bounds on the number of allowed mismatches for each segment.
    ///
    /// # Arguments
    ///
    /// * `order` - A vector specifying the order in which pattern parts are searched.
    ///   Each value is an index into the partitioned pattern.
    /// * `min_mismatches` - A vector of the same length as `order`, specifying the minimum
    ///   number of mismatches allowed for each corresponding pattern segment.
    /// * `max_mismatches` - A vector of the same length as `order`, specifying the maximum
    ///   number of mismatches allowed for each corresponding pattern segment.
    pub fn new(order: Vec<u8>, min_mismatches: Vec<u8>, max_mismatches: Vec<u8>) -> Search {
        Search { order, min_mismatches, max_mismatches }
    }

    /// Returns an iterator over the search parts.
    pub fn iter(&self) -> SearchIter {
        SearchIter {
            search: self,
            index: 0,
        }
    }

    /// Returns whether the direction of search at the given index is leftward.
    ///
    /// Compares current part's position to the previous part in the order.
    pub fn get_direction_left(&self, idx: usize) -> bool {
        if idx == 0 {
            return self.order[1] < self.order[0];
        }

        self.order[idx] < self.order[idx-1]
    }

    /// Gets the maximum number of mismatches allowed at the specified part index.
    pub fn get_upperbound(&self, idx: usize) -> u8 {
        self.max_mismatches[idx]
    }

    /// Gets the minimum number of mismatches required at the specified part index.
    pub fn get_lowerbound(&self, idx: usize) -> u8 {
        self.min_mismatches[idx]
    }

    /// Gets the pattern part index from the search order at the given step.
    pub fn get_part(&self, idx: usize) -> u8 {
        self.order[idx]
    }

    /// Returns the number of parts in the search.
    pub fn len(&self) -> usize {
        self.order.len()
    }

    /// Validates the internal consistency of the search definition.
    ///
    /// Ensures:
    /// - All fields have the same length.
    /// - max ≥ min mismatches for each step.
    /// - mismatch bounds do not decrease.
    /// - parts are contiguous (adjacent to the searched part).
    ///
    /// # Errors
    /// Returns an error if any consistency check fails.
    pub fn validate(&self) -> Result<(), Box<dyn Error>> {
        let len = self.order.len();
        if self.min_mismatches.len() != len || self.max_mismatches.len() != len {
            return Err("Length mismatch in Search fields".into());
        }

        let mut seen = (0, 0);
        let mut prev_min = 0;
        let mut prev_max = 0;

        for (i, &pos) in self.order.iter().enumerate() {
            // Check bounds
            let min = self.min_mismatches[i];
            let max = self.max_mismatches[i];

            if max < min {
                return Err(format!("max_mismatches[{}] < min_mismatches[{}]", i, i).into());
            }

            if i > 0 {
                if min < prev_min {
                    return Err(format!("min_mismatches decreased at step {}", i).into());
                }
                if max < prev_max {
                    return Err(format!("max_mismatches decreased at step {}", i).into());
                }
            }

            // Check contiguity
            if i == 0 {
                seen = (pos, pos);
            } else {
                let (begin, end) = seen;
                if begin == pos + 1 && pos < begin {
                    seen = (pos, end);
                } else if end + 1 == pos && pos < len as u8 {
                    seen = (begin, pos);
                } else {
                    return Err(format!("Part at pos {} does not border seen parts", pos).into());
                }
            }

            prev_min = min;
            prev_max = max;
        }

        Ok(())
    }
}

/// A search scheme, i.e., a set of `Search` configurations.
///
/// These are typically loaded from a file that encodes each search line
/// with three brace-enclosed lists:
/// `{order} {min_mismatches} {max_mismatches}`
pub struct SearchScheme {
    searches: Vec<Search>
}

impl SearchScheme {

    /// Creates a new `SearchScheme` from a list of individual `Search` objects.
    ///
    /// A `SearchScheme` defines a strategy for approximate string matching, where each `Search`
    /// specifies partitioning and mismatch bounds for searching a pattern.
    ///
    /// # Arguments
    ///
    /// * `searches` - A vector of `Search` instances, each representing one traversal strategy
    ///   through a partitioned pattern with allowed mismatches.
    ///
    /// # Returns
    ///
    /// A `SearchScheme` containing the provided `Search` objects.
    pub fn new(searches: Vec<Search>) -> SearchScheme {
        SearchScheme { searches }
    }

    /// Loads a `SearchScheme` from a file path.
    ///
    /// # Format
    /// Each line should contain three brace-enclosed lists:
    /// `{order} {min} {max}`
    ///
    /// # Errors
    /// Returns an error if the file cannot be read or is malformed.
    pub fn from_file(path: &Path) -> Result<Self, Box<dyn Error>> {
        let content = fs::read_to_string(path)?;
        let mut searches = Vec::new();

        for line in content.lines() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(format!("Invalid line format: {}", line).into());
            }

            let order = Self::parse_braced_numbers(parts[0])?;
            let min_mismatches = Self::parse_braced_numbers(parts[1])?;
            let max_mismatches = Self::parse_braced_numbers(parts[2])?;
            
            searches.push(Search {
                order,
                min_mismatches,
                max_mismatches
            });
        }

        Ok(SearchScheme { searches })
    }

    /// Validates all `Search` entries in the scheme.
    pub fn validate(&self) -> Result<(), Box<dyn Error>> {
        for (i, search) in self.searches.iter().enumerate() {
            search.validate().map_err(|e| format!("Search {} invalid: {}", i, e))?;
        }

        Ok(())
    }

    /// Parses a comma-separated list of numbers within braces.
    fn parse_braced_numbers(s: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let s = s.trim();
        if !(s.starts_with('{') && s.ends_with('}')) {
            return Err(format!("Invalid format for list: {}", s).into());
        }

        let numbers = s[1..s.len() - 1]
            .split(',')
            .map(|num| num.trim().parse::<u8>())
            .collect::<Result<Vec<_>, _>>()?;

        Ok(numbers)
    }

    /// Returns the number of parts in the search scheme (from the first search).
    pub fn get_parts_amount(&self) -> u8 {
        self.searches[0].order.len() as u8
    }

    /// Retrieves the search at the specified index.
    pub fn get_search(&self, index: usize) -> &Search {
        &self.searches[index]
    }

}

impl<'a> IntoIterator for &'a SearchScheme {
    type Item = &'a Search;
    type IntoIter = std::slice::Iter<'a, Search>;

    /// Creates an iterator over all `Search` instances in the `SearchScheme`.
    ///
    /// This allows iterating directly over a reference to a `SearchScheme`:
    fn into_iter(self) -> Self::IntoIter {
        self.searches.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_iter_yields_correct_tuples() {
        let search = Search {
            order: vec![0, 1, 2],
            min_mismatches: vec![0, 1, 1],
            max_mismatches: vec![1, 2, 3],
        };

        let results: Vec<_> = search.iter().collect();
        assert_eq!(results, vec![
            (0, 0, 1),
            (1, 1, 2),
            (2, 1, 3),
        ]);
    }

    #[test]
    fn test_get_direction_left() {
        let search = Search {
            order: vec![2, 1, 0],
            min_mismatches: vec![0, 0, 0],
            max_mismatches: vec![0, 0, 0],
        };

        assert_eq!(search.get_direction_left(0), true); // 1 < 2
        assert_eq!(search.get_direction_left(1), true); // 0 < 1
    }

    #[test]
    fn test_upper_lower_bound_access() {
        let search = Search {
            order: vec![0, 1],
            min_mismatches: vec![1, 2],
            max_mismatches: vec![3, 4],
        };

        assert_eq!(search.get_upperbound(0), 3);
        assert_eq!(search.get_lowerbound(1), 2);
    }

    #[test]
    fn test_get_part_and_len() {
        let search = Search {
            order: vec![1, 2, 3],
            min_mismatches: vec![0, 0, 0],
            max_mismatches: vec![0, 0, 0],
        };

        assert_eq!(search.get_part(1), 2);
        assert_eq!(search.len(), 3);
    }

    #[test]
    fn test_validate_valid_search() {
        let search = Search {
            order: vec![1, 2, 0],
            min_mismatches: vec![0, 1, 2],
            max_mismatches: vec![1, 2, 3],
        };

        assert!(search.validate().is_ok());
    }

    #[test]
    fn test_validate_invalid_length_mismatch() {
        let search = Search {
            order: vec![0, 1],
            min_mismatches: vec![0],
            max_mismatches: vec![1, 2],
        };

        let err = search.validate().unwrap_err().to_string();
        assert!(err.contains("Length mismatch"));
    }

    #[test]
    fn test_validate_invalid_bounds_order() {
        let search = Search {
            order: vec![0, 1, 2],
            min_mismatches: vec![0, 2, 1], // decreases
            max_mismatches: vec![1, 3, 2], // decreases
        };

        let err = search.validate().unwrap_err().to_string();
        assert!(err.contains("min_mismatches decreased") || err.contains("max_mismatches decreased"));
    }

    #[test]
    fn test_validate_non_contiguous_parts() {
        let search = Search {
            order: vec![0, 2, 4],
            min_mismatches: vec![0, 1, 2],
            max_mismatches: vec![1, 2, 3],
        };

        let err = search.validate().unwrap_err().to_string();
        assert!(err.contains("does not border"));
    }

    #[test]
    fn test_parse_braced_numbers_valid() {
        let input = "{1,2,3}";
        let parsed = SearchScheme::parse_braced_numbers(input).unwrap();
        assert_eq!(parsed, vec![1, 2, 3]);
    }

    #[test]
    fn test_parse_braced_numbers_invalid_format() {
        let input = "1,2,3";
        let err = SearchScheme::parse_braced_numbers(input).unwrap_err().to_string();
        assert!(err.contains("Invalid format"));
    }

    #[test]
    fn test_scheme_validate_all_valid() {
        let search1 = Search {
            order: vec![0, 1],
            min_mismatches: vec![0, 1],
            max_mismatches: vec![1, 2],
        };
        let search2 = Search {
            order: vec![1, 0],
            min_mismatches: vec![0, 1],
            max_mismatches: vec![1, 2],
        };

        let scheme = SearchScheme {
            searches: vec![search1, search2],
        };

        assert!(scheme.validate().is_ok());
    }

    #[test]
    fn test_scheme_validate_with_invalid_search() {
        let search = Search {
            order: vec![0, 2],
            min_mismatches: vec![0, 1],
            max_mismatches: vec![1, 2],
        };

        let scheme = SearchScheme {
            searches: vec![search],
        };

        let err = scheme.validate().unwrap_err().to_string();
        assert!(err.contains("invalid"));
    }
}