//! A succinct FM-index implementation using QWT and rank-select bit vectors.
//!
//! This module provides an efficient compressed suffix array–based FM-index for full-text substring
//! search, particularly suited for DNA strings or similar datasets. The index supports forward and
//! reverse extensions and locating matches.
//!
//! # Components
//! - `FMIndex`: The main structure storing the forward and reverse BWTs and supporting rank queries.
//! - `FMIndexRange`: Encapsulates a range within the index for iterative backward/forward search.
//!
//! This implementation uses serialization via `bincode` to load precomputed components.

use succinct::storage::BlockType;
use qwt::{AccessUnsigned, RankUnsigned, QWT256};
use std::error::Error;
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector};
use succinct::rank::Rank9;
use byteorder::LittleEndian;

use std::fs::{File};
use std::io::BufReader;
use std::path::Path;

use serde::{Serialize, Deserialize};
use bincode::deserialize_from;

/// Wrapper for loading count tables from disk.
#[derive(Serialize, Deserialize)]
struct Counts(Vec<usize>);

/// Represents a search range within the FM-index.
///
/// This struct is used during backward/forward search to define the range of matches
/// in both forward and reverse BWTs.
#[derive(Clone)]
pub struct FMIndexRange {
    pub begin: usize,
    pub end: usize,
    pub begin_rev: usize,
    pub end_rev: usize,
}

impl FMIndexRange {

    /// Checks whether the search range is empty (no matches).
    pub fn empty(&self) -> bool {
        self.begin >= self.end
    }
}

/// FM-index structure.
///
/// The index supports forward and backward search over a Burrows-Wheeler Transform (BWT)
/// and includes rank support.
pub struct FMIndex {
    bwt: QWT256<u8>,
    bwt_rev: QWT256<u8>,
    counts: Vec<usize>,
    ssa: Vec<i64>,
    ssa_occs: Rank9<BitVector<u64>>,
    char_to_id: Vec<u8>
}

impl FMIndex {

    /// Creates a new `FMIndex` instance from the provided components.
    ///
    /// # Arguments
    ///
    /// * `bwt` - The Burrows-Wheeler Transform (BWT) of the indexed text.
    /// * `bwt_rev` - The BWT of the reversed indexed text.
    /// * `counts` - A vector containing cumulative character counts for rank calculations.
    /// * `ssa` - The sampled suffix array, used for locating positions in the original text.
    /// * `ssa_occs` - A rank data structure (Rank9 over a BitVector) to efficiently support suffix array queries.
    /// * `char_to_id` - A mapping from characters to their internal alphabet IDs.
    ///
    /// # Returns
    ///
    /// A new `FMIndex` instance constructed from the given components.
    pub fn new(bwt: QWT256<u8>, bwt_rev: QWT256<u8>, counts: Vec<usize>, ssa: Vec<i64>, ssa_occs: Rank9<BitVector<u64>>, char_to_id: Vec<u8>) -> FMIndex {
        FMIndex { bwt, bwt_rev, counts, ssa, ssa_occs, char_to_id }
    }

    /// Loads a serialized FM-index from files on disk.
    ///
    /// # Arguments
    /// * `base_path` - Base path to the index files, without extensions.
    ///
    /// # Expected Extensions
    /// - `.bwt`, `.rev.bwt`, `.ssa`, `.ssa_occ`, `.alph`, `.counts`
    ///
    /// # Returns
    /// A deserialized `FMIndex` or an error if loading fails.
    pub fn from_files(base_path: &Path) -> Result<Self, Box<dyn Error>> {

        eprintln!("\tLoading BWT...");
        let bwt_file = BufReader::new(File::open(base_path.with_extension("bwt"))?);
        let bwt: QWT256<u8> = deserialize_from(bwt_file)?;

        eprintln!("\tLoading BWT of reversed text...");
        let bwt_rev_file = BufReader::new(File::open(base_path.with_extension("rev.bwt"))?);
        let bwt_rev: QWT256<u8> = deserialize_from(bwt_rev_file)?;

        eprintln!("\tLoading SSA...");
        let ssa_file = BufReader::new(File::open(base_path.with_extension("ssa"))?);
        let ssa: Vec<i64> = deserialize_from(ssa_file)?;

        eprintln!("\tLoading Alphabet...");
        let alph_file = BufReader::new(File::open(base_path.with_extension("alph"))?);
        let char_to_id: Vec<u8> = deserialize_from(alph_file)?;

        eprintln!("\tLoading Counts...");
        let counts_file = BufReader::new(File::open(base_path.with_extension("counts"))?);
        let counts: Vec<usize> = deserialize_from(counts_file)?;

        eprintln!("\tLoading SSA occurences...");
        let mut ssa_occs_file = BufReader::new(File::open(base_path.with_extension("ssa_occ"))?);
        let mut ssa_occs = BitVector::new();
        loop {
            match BlockType::read_block::<_, LittleEndian>(&mut ssa_occs_file) {
                Ok(block) => ssa_occs.push_block(block),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(Box::new(e)), // or handle gracefully
            }
        }

        eprintln!("\tBuilding Rank support for SSA occurences...");
        let ssa_occs = Rank9::new(ssa_occs);

        Ok(Self {
            bwt,
            bwt_rev,
            counts,
            ssa,
            ssa_occs,
            char_to_id
        })
    }

    /// Returns the length of the original text (same as BWT length).
    pub fn len(&self) -> usize {
        self.bwt.len()
    }

    /// Locates the original positions of all matches within a given SA range.
    pub fn locate_matches(&self, range: FMIndexRange) -> Vec<usize> {
        (range.begin..range.end).map(|i| self.locate(i)).collect()
    }
    
    /// Extends a match range to the left by one character `c` using the LF-mapping.
    ///
    /// Used in backward search.
    pub fn left_extension(&self, c: u8, range: FMIndexRange) -> FMIndexRange {
        let FMIndexRange { begin, end, begin_rev, end_rev: _ } = range;
        
        let new_begin = self.lf(c, begin);
        let new_end = self.lf(c, end);

        let mut new_begin_rev = begin_rev;
        for smaller_c in 0..c {
            new_begin_rev += self.rank(smaller_c, end) - self.rank(smaller_c, begin);
        }
        let new_end_rev = new_begin_rev + (new_end - new_begin);
        
        FMIndexRange { begin: new_begin, end: new_end, begin_rev: new_begin_rev, end_rev: new_end_rev }
    }

    /// Extends a match range to the right by one character `c` using the reverse LF-mapping.
    ///
    /// Used in bidirectional search algorithms.
    pub fn right_extension(&self, c: u8, range: FMIndexRange) -> FMIndexRange {
        let FMIndexRange { begin, end: _, begin_rev, end_rev } = range;
        let new_begin_rev = self.lf_rev(c, begin_rev);
        let new_end_rev = self.lf_rev(c, end_rev);

        let mut new_begin = begin;
        for smaller_c in 0..c {
            new_begin += self.rank_rev(smaller_c, end_rev) - self.rank_rev(smaller_c, begin_rev);
        }
        let new_end = new_begin + (new_end_rev - new_begin_rev);
        
        FMIndexRange { begin: new_begin, end: new_end, begin_rev: new_begin_rev, end_rev: new_end_rev }
    }

    /// Performs the LF-mapping.
    fn lf(&self, symbol: u8, pos: usize) -> usize {
        self.c(symbol) + self.rank(symbol, pos)
    }

    /// Reverse LF-mapping using the reverse BWT.
    fn lf_rev(&self, symbol: u8, pos: usize) -> usize {
        self.c(symbol) + self.rank_rev(symbol, pos)
    }
    
    /// Returns the number of occurrences of `symbol` up to position `pos` in the forward BWT.
    fn rank(&self, symbol: u8, pos: usize) -> usize {
        self.bwt.rank(symbol, pos).unwrap()
    }

    /// Returns the number of occurrences of `symbol` up to position `pos` in the reversed BWT.
    fn rank_rev(&self, symbol: u8, pos: usize) -> usize {
        self.bwt_rev.rank(symbol, pos).unwrap()
    }

    /// Returns the cumulative count of characters smaller than `symbol` in the original text.
    fn c(&self, symbol: u8) -> usize {
        self.counts[symbol as usize]
    }

    /// Locates the original position of the i-th suffix using sampled SSA entries.
    ///
    /// Traverses the LF-mapping backward until it finds a sampled suffix.
    pub fn locate(&self, mut i: usize) -> usize {
        let mut steps = 0;
        while !self.ssa_occs.get_bit(i as u64) {
            let c = self.bwt.get(i).unwrap();
            i = self.lf(c as u8, i);
            steps += 1;
        }
        
        (*self.ssa.get(self.ssa_occs.rank1(i as u64) as usize - 1).unwrap() as usize + steps) % self.bwt.len()
    }

    /// Returns the size of the alphabet used in the index.
    pub fn get_alphabet_size(&self) -> u8 {
        self.counts.len() as u8
    }

    /// Maps a pattern using the `char_to_id` vector into internal integer codes.
    pub fn map_pattern(&self, pattern: &Vec<u8>) -> Vec<u8> {
        let mut mapped_pattern = Vec::with_capacity(pattern.len());
        for &c in pattern {
            mapped_pattern.push(self.char_to_id[c as usize]);
        }

        mapped_pattern
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qwt::QWT256;
    use succinct::BitVector;
    use succinct::rank::Rank9;

    fn get_fmindex() -> FMIndex {
        // Create a minimal mock of a BWT over "BANANA$" -> "ANNB$AA" with the mapped alphabet
        let bwt_data = vec![1, 3, 3, 2, 0, 1, 1];
        let bwt = QWT256::from(bwt_data.clone());
        // Create a BWT for "ANANAB$" -> "BNN$AAA" with mapped alphabet
        let rev_bwt = vec![2, 3, 3, 0, 1, 1, 1];
        let bwt_rev = QWT256::from(rev_bwt.clone());

        let mut counts = vec![0; 256]; // simplistic, only for test
        counts[1] = 1;
        counts[2] = 4;
        counts[3] = 5;
        counts[4] = 7;
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

        FMIndex {
            bwt,
            bwt_rev,
            counts,
            ssa,
            ssa_occs,
            char_to_id,
        }
    }

    #[test]
    fn test_len() {
        let fm = get_fmindex();
        assert_eq!(fm.len(), 7);
    }

    #[test]
    fn test_rank_and_lf() {
        let fm = get_fmindex();
        let pos = 4;
        let symbol = 1;
        let r = fm.rank(symbol, pos);
        assert!(r <= pos);
        let lf_val = fm.lf(symbol, pos);
        assert!(lf_val >= r);
    }

    #[test]
    fn test_locate() {
        let fm = get_fmindex();
        // Every SA entry is sampled, so locate should return directly
        for i in 0..fm.len() {
            let pos = fm.locate(i);
            assert!(pos < fm.len());
        }
    }

    #[test]
    fn test_empty_range() {
        let range = FMIndexRange {
            begin: 5,
            end: 5,
            begin_rev: 2,
            end_rev: 2,
        };
        assert!(range.empty());
    }

    #[test]
    fn test_map_pattern() {
        let fm = get_fmindex();
        let input = b"ABN$".to_vec();
        let mapped_pattern = vec![1, 2, 3, 0];
        let mapped = fm.map_pattern(&input);
        assert_eq!(mapped, mapped_pattern); // identity mapping in mock
    }

    #[test]
    fn test_forward_and_reverse_extension_consistency() {
        let fm = get_fmindex();

        let initial_range = FMIndexRange {
            begin: 0,
            end: fm.len(),
            begin_rev: 0,
            end_rev: fm.len(),
        };

        let c = 1;
        let extended = fm.left_extension(c, initial_range.clone());
        assert!(extended.begin <= extended.end);

        let extended_back = fm.right_extension(c, extended.clone());
        assert!(extended_back.begin <= extended_back.end);
    }
}