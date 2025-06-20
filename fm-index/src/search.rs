//! Approximate pattern matching using FM-index and search schemes.
//!
//! This module implements parallelized approximate matching for multiple patterns
//! using an FM-index and flexible search schemes. It supports banded dynamic
//! programming to handle mismatches and provides mechanisms to locate all valid
//! matches in the indexed text.

use crate::fm_index::{FMIndex, FMIndexRange};
use crate::search_scheme::{SearchScheme, Search};
use crate::search_pattern::SearchPattern;
use crate::banded_matrix::BandedMatrix;
use std::collections::HashSet;
use std::error::Error;
use rayon::prelude::*;


/// Represents an intermediate occurrence in the FM-index search process.
/// 
/// - `range`: The range in the FM-index corresponding to the current match.
/// - `mismatches`: Number of mismatches encountered so far.
/// - `match_length`: Length of the currently matched segment.
pub struct FMOcc {
    pub range: FMIndexRange,
    pub mismatches: u8,
    pub match_length: usize
}


/// Represents a partial match extension during the search traversal.
/// 
/// - `range`: Current FM-index range for the partial match.
/// - `depth`: Current depth of the search.
/// - `c`: Character used to extend the match.
pub struct FMMatchToExplore {
    pub range: FMIndexRange,
    pub depth: usize,
    pub c: u8
}

/// Represents a completed match in the text.
/// 
/// - `start_position`: The starting position of the match in the original text.
/// - `length`: The length of the matched substring.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct FMMatch {
    pub start_position: usize,
    pub length: usize
}

/// Searches for multiple patterns in parallel using the provided FM-index and search scheme.
///
/// # Arguments
/// * `fm_index` - Reference to an FMIndex instance.
/// * `patterns` - A vector of patterns to search for (each as a `Vec<u8>`).
/// * `search_scheme` - A search scheme guiding the approximate matching.
///
/// # Returns
/// A `Result` containing a set of `FMMatch` if successful, or an error if any occurs.
///
/// # Parallelism
/// Utilizes Rayon to perform each pattern search in parallel.
pub fn search_multiple(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, search_scheme: SearchScheme) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {

    let matches: HashSet<FMMatch> = patterns
        .into_par_iter() // Parallel iterator from rayon
        .map(|pattern| approximate_search(fm_index, pattern, &search_scheme))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect();

    Ok(matches)
}

/// Performs approximate pattern matching for a single pattern using a given search scheme.
///
/// # Arguments
/// * `fm_index` - Reference to an FMIndex instance.
/// * `pattern` - The pattern to search for (as a `Vec<u8>`).
/// * `search_scheme` - A reference to a search scheme defining the allowed mismatches.
///
/// # Returns
/// A `Result` containing a set of `FMMatch` instances representing matches.
pub fn approximate_search(fm_index: &FMIndex, mut pattern: Vec<u8>, search_scheme: &SearchScheme) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {

    let mut matches: HashSet<FMMatch> = HashSet::new();

    translate_l_to_i(&mut pattern);
    let pattern = fm_index.map_pattern(&pattern);
    let pattern = SearchPattern::new(pattern, search_scheme.get_parts_amount() as usize)?;

    for search in search_scheme.into_iter() {

        let range = FMIndexRange { begin: 0, end: fm_index.len(), begin_rev: 0, end_rev: fm_index.len() };
        let start_occ = FMOcc { range, mismatches: 0, match_length: 0 };

        let mut occs: Vec<FMOcc> = Vec::new();

        approximate_search_rec(&fm_index, search, start_occ, &pattern, 0, &mut occs);
        
        for occ in occs {
            for pos in occ.range.begin..occ.range.end {
                let _ = matches.insert(FMMatch { start_position: fm_index.locate(pos), length: occ.match_length });
            }
        }
    }

    Ok(matches)

}

/// Recursively performs approximate matching for a segment of the search scheme.
///
/// This function extends matches while keeping track of mismatches using a banded matrix.
///
/// # Arguments
/// * `fm_index` - Reference to the FM-index.
/// * `search` - The current search object defining direction and bounds.
/// * `occ` - The current partial match state.
/// * `pattern` - The structured search pattern split by parts.
/// * `idx_in_search` - Current index in the search scheme.
/// * `occs` - Accumulator vector for completed occurrences.
fn approximate_search_rec(fm_index: &FMIndex, search: &Search, occ: FMOcc, pattern: &SearchPattern, idx_in_search: u8, occs: &mut Vec<FMOcc>) {
    
    let idx = idx_in_search as usize;
    let direction = search.get_direction_left(idx);
    let part_idx = search.get_part(idx);

    let part_size = pattern.get_part_len(part_idx);
    let width = search.get_upperbound(idx) - occ.mismatches;
    let mut bandedmatrix = BandedMatrix::new(part_size, width, occ.mismatches);

    let mut stack: Vec<FMMatchToExplore> = Vec::with_capacity(fm_index.get_alphabet_size() as usize * pattern.len());
    extend_match(fm_index, &occ.range, 0, &mut stack, direction);

    while !stack.is_empty() {

        let current_pos = stack.pop().unwrap();

        let part: &Vec<u8> = &pattern.get_part(search.get_part(idx), direction);
        let mismatches = bandedmatrix.update_matrix_row(part, current_pos.depth, current_pos.c);
        if mismatches <= search.get_upperbound(idx) {

            extend_match(fm_index, &current_pos.range, current_pos.depth, &mut stack, direction);

            if bandedmatrix.is_final_column(current_pos.depth) {

                let distance: u8 = bandedmatrix.get_value_in_final_column(current_pos.depth);

                if distance <= search.get_upperbound(idx) && distance >= search.get_lowerbound(idx) {

                    let matched_length: usize = occ.match_length + current_pos.depth;
                    let next_occ = FMOcc { range: current_pos.range, mismatches: distance, match_length: matched_length };

                    if idx == search.len() - 1 {
                        occs.push(next_occ);
                    } else {
                        approximate_search_rec(fm_index, search, next_occ, pattern, idx_in_search + 1, occs);
                    }
                }
            }
        }
    }
}

/// Pushes all valid single-character extensions of the current match onto the stack.
///
/// # Arguments
/// * `fm_index` - Reference to the FM-index.
/// * `range` - Current FM-index range to extend.
/// * `depth` - Current search depth.
/// * `stack` - Stack of matches to explore further.
/// * `direction_left` - Whether to extend leftward or rightward.
fn extend_match(fm_index: &FMIndex, range: &FMIndexRange, depth: usize, stack: &mut Vec<FMMatchToExplore>, direction_left: bool) {
    let alphabet_size = fm_index.get_alphabet_size();
    for c in 0..alphabet_size {
        let new_range = if direction_left { 
            fm_index.left_extension(c, range.clone())
        } else {
            fm_index.right_extension(c, range.clone())
        };

        if ! new_range.empty() {
            stack.push(FMMatchToExplore { range: new_range.clone(), depth: depth + 1, c: c });
        }
    }
}

/// Translates all 'L' characters in a given sequence to 'I'.
///
/// # Arguments
/// * `text` - A mutable slice of bytes representing the sequence to normalize.
pub fn translate_l_to_i(text: &mut [u8]) {
    for character in text.iter_mut() {
        if *character == b'L' {
            *character = b'I'
        }
    }
}