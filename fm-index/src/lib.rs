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

pub struct FMIndex {
    bwt_wavelet: WaveletMatrix<Rank9Sel>,
    bwt: Vec<u8>,
    counts: Vec<usize>,
    ssa: Vec<i64>,
    ssa_occs: Rank9<BitVector<u64>>
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
        let bwt_wavelet = Self::build_wavelet_matrix(bwt.clone())?;

        // Load SSA
        let ssa_file = BufReader::new(File::open(base_path.with_extension("ssa"))?);
        let ssa: Vec<i64> = deserialize_from(ssa_file)?;

        // Load counts
        let counts_file = BufReader::new(File::open(base_path.with_extension("counts"))?);
        let counts = deserialize_from(counts_file)?;

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
            bwt_wavelet,
            bwt,
            counts,
            ssa,
            ssa_occs
        })
    }
    
    pub fn extend(&self, c: u8, sp: usize, ep: usize) -> (usize, usize) {
        let new_sp = self.c(c) + self.rank(c, sp).unwrap();
        let new_ep = self.c(c) + self.rank(c, ep).unwrap();
        (new_sp, new_ep)
    }

    
    pub fn rank(&self, symbol: u8, pos: usize) -> Option<usize> {
        self.bwt_wavelet.rank(pos, symbol as usize)
    }

    pub fn c(&self, symbol: u8) -> usize {
        self.counts[symbol as usize]
    }

    pub fn locate(&self, mut i: usize) -> usize {
        let mut steps = 0;
        while !self.ssa_occs.get_bit(i as u64) {
            let c = self.bwt_wavelet.access(i).unwrap();
            i = self.c(c as u8) + self.bwt_wavelet.rank(i, c).unwrap();
            steps += 1;
        }
        
        (*self.ssa.get(self.ssa_occs.rank1(i as u64) as usize - 1).unwrap() as usize + steps) % self.bwt.len()
    }
}