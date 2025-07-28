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
use std::hash::Hash;
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


pub fn search_multiple_casanovo(fm_index: &FMIndex, patterns: Vec<Vec<u8>>)  -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {
    let all_patterns = transform_patterns(patterns);

    let matches: HashSet<FMMatch> = search_multiple_exact_grouped(fm_index, all_patterns)?
        .into_iter()
        .flatten()
        .collect();

    Ok(matches)
}

fn transform_patterns(patterns: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut result = Vec::with_capacity(patterns.len() * 3); 

    for pattern in patterns {
        let mut transformed_first = pattern.clone();
        let mut transformed_last = pattern.clone();

        // Invert the first 3 elements (if available)
        transformed_first[0] = pattern[2];
        transformed_first[2] = pattern[0];

        let length = pattern.len();
        transformed_last[length-1] = pattern[length-3];
        transformed_last[length-3] = pattern[length-1];

        // Push original and transformed versions
        result.push(pattern);
        result.push(transformed_first);
        result.push(transformed_last);
    }

    result
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
pub fn search_multiple(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, search_scheme: &SearchScheme) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {

    let matches: HashSet<FMMatch> = search_multiple_grouped(fm_index, patterns, search_scheme)?
        .into_iter()
        .flatten()
        .collect();

    Ok(matches)
}

/// Searches for multiple patterns in parallel using the provided FM-index and search scheme, only search unique matches in the text.
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
pub fn search_multiple_unique(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, search_scheme: &SearchScheme) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {

    let matches: HashSet<FMMatch> = search_multiple_grouped_unique(fm_index, patterns, search_scheme)?
        .into_iter()
        .flatten()
        .collect();

    Ok(matches)
}

/// Searches for multiple patterns in parallel using the provided FM-index and search scheme, only retrieve unique matches.
///
/// # Arguments
/// * `fm_index` - Reference to an FMIndex instance.
/// * `patterns` - A vector of patterns to search for (each as a `Vec<u8>`).
/// * `search_scheme` - A search scheme guiding the approximate matching.
///
/// # Returns
/// A `Result` containing a vector of sets of `FMMatch` if successful, or an error if any occurs.
/// Each set contains the matches for one pattern.
///
/// # Parallelism
/// Utilizes Rayon to perform each pattern search in parallel.
/// use std::collections::HashSet;
pub fn search_multiple_grouped_unique(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, search_scheme: &SearchScheme) -> Result<Vec<HashSet<FMMatch>>, Box<dyn Error + Send + Sync>> {

    let matches: Vec<HashSet<FMMatch>> = patterns
        .into_par_iter() // Parallel iterator from rayon
        .map(|pattern| approximate_search(fm_index, pattern, search_scheme))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(|occs| {
            let mut matches: HashSet<FMMatch> = HashSet::new();
            for occ in occs {
                let _ = matches.insert(FMMatch { start_position: fm_index.locate(occ.range.begin), length: occ.match_length });
            }
            return matches;
        })
        .collect();

    Ok(matches)
}


/// Searches for multiple patterns in parallel using the provided FM-index and search scheme.
///
/// # Arguments
/// * `fm_index` - Reference to an FMIndex instance.
/// * `patterns` - A vector of patterns to search for (each as a `Vec<u8>`).
/// * `search_scheme` - A search scheme guiding the approximate matching.
///
/// # Returns
/// A `Result` containing a vector of sets of `FMMatch` if successful, or an error if any occurs.
/// Each set contains the matches for one pattern.
///
/// # Parallelism
/// Utilizes Rayon to perform each pattern search in parallel.
pub fn search_multiple_grouped(fm_index: &FMIndex, patterns: Vec<Vec<u8>>, search_scheme: &SearchScheme) -> Result<Vec<HashSet<FMMatch>>, Box<dyn Error + Send + Sync>> {
    let matches: Vec<HashSet<FMMatch>> = patterns
        .into_par_iter() // Parallel iterator from rayon
        .map(|pattern| approximate_search(fm_index, pattern, search_scheme))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(|occs| {
            let mut matches: HashSet<FMMatch> = HashSet::new();
            for occ in occs {
                for pos in occ.range.begin..occ.range.end {
                    let _ = matches.insert(FMMatch { start_position: fm_index.locate(pos), length: occ.match_length });
                }
            }
            return matches;
        })
        .collect();

    Ok(matches)
}

/// Searches for multiple patterns using exact matching on the FM-index.
/// 
/// # Arguments
/// * `fm_index` - Reference to the FMIndex instance.
/// * `patterns` - Vector of patterns to search for (each as a `Vec<u8>`).
/// 
/// # Returns
/// A `Result` containing a set of positions (`FMMatch`) in the original text where exact matches start.
/// 
/// # Parallelism
/// Uses Rayon to perform all pattern searches in parallel for increased throughput.
pub fn search_multiple_exact(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {
    let matches: HashSet<FMMatch> = search_multiple_exact_grouped(fm_index, patterns)?
        .into_iter()
        .flatten()
        .collect();

    Ok(matches)
}

/// Searches for multiple patterns using exact matching on the FM-index.
/// 
/// # Arguments
/// * `fm_index` - Reference to the FMIndex instance.
/// * `patterns` - Vector of patterns to search for (each as a `Vec<u8>`).
/// 
/// # Returns
/// A `Result` containing a vector of sets of positions (`FMMatch`) in the original text where exact matches start.
/// Each set contains the matches for one pattern.
/// 
/// # Parallelism
/// Uses Rayon to perform all pattern searches in parallel for increased throughput.
pub fn search_multiple_exact_grouped(fm_index: &FMIndex, patterns: Vec<Vec<u8>>) -> Result<Vec<HashSet<FMMatch>>, Box<dyn Error + Send + Sync>> {
    let matches: Vec<HashSet<FMMatch>> = patterns
        .into_par_iter() // Parallel iterator from rayon
        .map(|pattern| exact_search(fm_index, pattern))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(matches)
}


/// Searches for exact matches of a single pattern in the FM-index.
/// 
/// # Arguments
/// * `fm_index` - Reference to the FMIndex instance.
/// * `pattern` - The pattern to search (as a `Vec<u8>`).
/// 
/// # Returns
/// A `Result` containing a set of positions (`FMMatch`) where the exact match was found in the original text.
pub fn exact_search(fm_index: &FMIndex, mut pattern: Vec<u8>) -> Result<HashSet<FMMatch>, Box<dyn Error + Send + Sync>> {

    let mut matches: HashSet<FMMatch> = HashSet::new();

    translate_l_to_i(&mut pattern);
    let pattern = fm_index.map_pattern(&pattern);

    let range = FMIndexRange { begin: 0, end: fm_index.len(), begin_rev: 0, end_rev: 0 };
    let new_range = fm_index.match_exact(range, &pattern);

    for pos in new_range.begin..new_range.end {
        let _ = matches.insert(FMMatch { start_position: fm_index.locate(pos), length: pattern.len() });
    }

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
/// A `Result` containing a vector of `FMOcc` instances representing matches.
pub fn approximate_search(fm_index: &FMIndex, mut pattern: Vec<u8>, search_scheme: &SearchScheme) -> Result<Vec<FMOcc>, Box<dyn Error + Send + Sync>> {

    let mut matches: Vec<FMOcc> = Vec::new();

    translate_l_to_i(&mut pattern);
    let pattern = fm_index.map_pattern(&pattern);
    let pattern = SearchPattern::new(pattern, search_scheme.get_parts_amount() as usize)?;

    for search in search_scheme.into_iter() {

        let range = FMIndexRange { begin: 0, end: fm_index.len(), begin_rev: 0, end_rev: fm_index.len() };
        let start_occ = FMOcc { range, mismatches: 0, match_length: 0 };

        let mut occs: Vec<FMOcc> = Vec::new();

        approximate_search_rec(&fm_index, search, start_occ, &pattern, 0, &mut occs);
        
        matches.extend(occs.into_iter());
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm_index::FMIndex;
    use crate::search_scheme::{SearchScheme, Search};
    use qwt::QWT256;
    use succinct::{BitVecPush, BitVector};
    use succinct::rank::Rank9;

    fn get_fmindex() -> FMIndex {
        // Create a minimal mock of a BWT over "BANANA$" -> "ANNB$AA" with the mapped alphabet
        let bwt_data = vec![1, 3, 3, 2, 0, 1, 1];
        let bwt = QWT256::from(bwt_data.clone());
        // Create a BWT for "ANANAB$" -> "BNN$AAA" with mapped alphabet
        let rev_bwt = vec![2, 3, 3, 0, 1, 1, 1];
        let bwt_rev = QWT256::from(rev_bwt.clone());

        let counts = vec![0, 1, 4, 5]; // simplistic, only for test
        let ssa = vec![6, 3, 0]; // some sampled suffix array

        // Mark all suffixes as sampled
        let mut bv = BitVector::new();
        bv.push_bit(true);
        bv.push_bit(false);
        bv.push_bit(true);
        bv.push_bit(false);
        bv.push_bit(true);
        bv.push_bit(false);
        bv.push_bit(false);
        let ssa_occs = Rank9::new(bv);

        let mut char_to_id = vec![0u8; 256];
        for (i, c) in b"$ABN".iter().enumerate() {
            char_to_id[*c as usize] = i as u8;
        }

        FMIndex::new(bwt, bwt_rev, counts, ssa, ssa_occs, char_to_id)
    }

    fn get_search_scheme() -> SearchScheme {
        let search_1 = Search::new( vec![1, 0], vec![0, 0], vec![0, 1]);
        let search_2 = Search::new( vec![0, 1], vec![0, 0], vec![0, 1]);
        
        SearchScheme::new(vec![search_1, search_2])
    }

    #[test]
    fn test_translate_l_to_i() {
        let mut seq = b"ALLL".to_vec();
        translate_l_to_i(&mut seq);
        assert_eq!(&seq, b"AIII");
    }

    #[test]
    fn test_approximate_search_basic() {
        let fm_index = get_fmindex();
        let search_scheme = get_search_scheme();

        let pattern = b"NA".to_vec();

        let result = approximate_search(&fm_index, pattern, &search_scheme);

        assert!(result.is_ok());
        let matches = result.unwrap();

        // Expect at least one match since pattern is 'NA'
        assert!(!matches.is_empty());

        // Matches contain positions within the text length
        for m in matches {
            assert!(fm_index.locate(m.range.begin) < fm_index.len());
            assert!(m.match_length >= 1);
        }
    }

    #[test]
    fn test_search_multiple_parallel() {
        let fm_index = get_fmindex();
        let search_scheme = get_search_scheme();

        let pattern_1 = b"AB".to_vec();
        let pattern_2 = b"NA".to_vec();
        let patterns = vec![pattern_1, pattern_2];

        let result = search_multiple(&fm_index, patterns, &search_scheme);

        assert!(result.is_ok());
        let matches = result.unwrap();

        assert!(matches.len() >= 1);

        for m in matches {
            assert!(m.start_position < fm_index.len());
            assert!(m.length > 0);
        }
    }

    #[test]
    fn test_exact_search_single() {
        let fm_index = get_fmindex();

        let pattern = b"NA".to_vec(); // should appear in the mock text "BANANA$"

        let result = exact_search(&fm_index, pattern);
        assert!(result.is_ok());

        let matches = result.unwrap();

        assert!(!matches.is_empty());
        for pos in matches {
            assert!(pos.start_position < fm_index.len());
        }
    }

    #[test]
    fn test_search_multiple_exact_parallel() {
        let fm_index = get_fmindex();

        let patterns = vec![
            b"NA".to_vec(),
            b"BA".to_vec(),
            b"AN".to_vec()
        ];

        let result = search_multiple_exact(&fm_index, patterns);
        assert!(result.is_ok());

        let matches = result.unwrap();
        assert!(matches.iter().all(|&pos| pos.start_position < fm_index.len()));
    }

    #[test]
    fn test_search_multiple_exact_empty() {
        let fm_index = get_fmindex();
        let patterns = vec![];

        let result = search_multiple_exact(&fm_index, patterns);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }
}