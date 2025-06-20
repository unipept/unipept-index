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
use qwt::{AccessUnsigned, RankUnsigned, SpaceUsage, QWT256};
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