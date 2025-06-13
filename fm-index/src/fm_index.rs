use succinct::storage::BlockType;
use sucds::char_sequences::WaveletMatrix;
use sucds::int_vectors::CompactVector;
use sucds::bit_vectors::Rank9Sel;
use std::error::Error;
use succinct::{BitRankSupport, BitVec, BitVecPush, BitVector};
use succinct::rank::Rank9;
use byteorder::LittleEndian;

use std::fs::{File};
use std::io::{BufReader, Read};
use std::path::Path;

use serde::{Serialize, Deserialize};
use bincode::deserialize_from;

#[derive(Serialize, Deserialize)]
struct Counts(Vec<usize>);

#[derive(Clone)]
pub struct FMIndexRange {
    pub begin: usize,
    pub end: usize,
    pub begin_rev: usize,
    pub end_rev: usize,
}

impl FMIndexRange {

    pub fn empty(&self) -> bool {
        self.begin >= self.end
    }
}

pub struct FMIndex {
    bwt: WaveletMatrix<Rank9Sel>,
    bwt_rev: WaveletMatrix<Rank9Sel>,
    counts: Vec<usize>,
    ssa: Vec<i64>,
    ssa_occs: Rank9<BitVector<u64>>,
    char_to_id: Vec<u8>
}

impl FMIndex {

    fn build_wavelet_matrix(bwt: Vec<u8>) -> Result<WaveletMatrix<Rank9Sel>, Box<dyn Error>> {
        let mut seq: CompactVector = CompactVector::with_capacity(bwt.len(), 8)?;
        seq.extend(bwt.into_iter().map(|e| e as usize))?;
        WaveletMatrix::<Rank9Sel>::new(seq).map_err(|_| "Could not create Wavelet Matrix".into())
    }

    pub fn from_files(base_path: &Path) -> Result<Self, Box<dyn Error>> {
        // Load BWT
        let mut bwt = Vec::new();
        BufReader::new(File::open(base_path.with_extension("bwt"))?)
            .read_to_end(&mut bwt)?;
        let bwt = Self::build_wavelet_matrix(bwt)?;

        // Load BWT of reverse text
        let mut bwt_rev = Vec::new();
        BufReader::new(File::open(base_path.with_extension("rev.bwt"))?)
            .read_to_end(&mut bwt_rev)?;
        let bwt_rev = Self::build_wavelet_matrix(bwt_rev)?;

        // Load SSA
        let ssa_file = BufReader::new(File::open(base_path.with_extension("ssa"))?);
        let ssa: Vec<i64> = deserialize_from(ssa_file)?;

        // Load char_to_id mapping
        let alph_file = BufReader::new(File::open(base_path.with_extension("alph"))?);
        let char_to_id: Vec<u8> = deserialize_from(alph_file)?;

        // Load counts
        let counts_file = BufReader::new(File::open(base_path.with_extension("counts"))?);
        let counts: Vec<usize> = deserialize_from(counts_file)?;

        // Load SSA occurences
        let mut ssa_occs_file = BufReader::new(File::open(base_path.with_extension("ssa_occ"))?);
        let mut ssa_occs = BitVector::new();
        loop {
            match BlockType::read_block::<_, LittleEndian>(&mut ssa_occs_file) {
                Ok(block) => ssa_occs.push_block(block),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(Box::new(e)), // or handle gracefully
            }
        }

        // Add rank/select support:
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

    pub fn len(&self) -> usize {
        self.bwt.len()
    }

    pub fn locate_matches(&self, range: FMIndexRange) -> Vec<usize> {
        (range.begin..range.end).map(|i| self.locate(i)).collect()
    }
    
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

    fn lf(&self, symbol: u8, pos: usize) -> usize {
        self.c(symbol) + self.rank(symbol, pos)
    }

    fn lf_rev(&self, symbol: u8, pos: usize) -> usize {
        self.c(symbol) + self.rank_rev(symbol, pos)
    }
    
    fn rank(&self, symbol: u8, pos: usize) -> usize {
        self.bwt.rank(pos, symbol as usize).unwrap()
    }

    fn rank_rev(&self, symbol: u8, pos: usize) -> usize {
        self.bwt_rev.rank(pos, symbol as usize).unwrap()
    }

    fn c(&self, symbol: u8) -> usize {
        self.counts[symbol as usize]
    }

    pub fn locate(&self, mut i: usize) -> usize {
        let mut steps = 0;
        while !self.ssa_occs.get_bit(i as u64) {
            let c = self.bwt.access(i).unwrap();
            i = self.c(c as u8) + self.bwt.rank(i, c).unwrap();
            steps += 1;
        }
        
        (*self.ssa.get(self.ssa_occs.rank1(i as u64) as usize - 1).unwrap() as usize + steps) % self.bwt.len()
    }

    pub fn get_alphabet_size(&self) -> u8 {
        self.counts.len() as u8
    }

    pub fn map_pattern(&self, pattern: &Vec<u8>) -> Vec<u8> {
        let mut mapped_pattern = Vec::with_capacity(pattern.len());
        for &c in pattern {
            mapped_pattern.push(self.char_to_id[c as usize]);
        }

        mapped_pattern
    }
}