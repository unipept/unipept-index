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

