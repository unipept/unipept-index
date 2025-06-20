//! FM-Index Construction
//!
//! This module defines the logic for building a bidirectional FM-index from a given
//! concatenated sequence file (e.g., protein sequences). The FM-index is stored in
//! multiple serialized components including the BWT, suffix array samples, and
//! occurrence bit vectors for fast pattern searching.
//!
//! Steps include:
//! - Character normalization (`L` → `I`)
//! - Text transformation and alphabet mapping
//! - Suffix array construction (using `libsais64_rs`)
//! - Burrows-Wheeler Transform (BWT) computation
//! - Suffix array sampling and serialization
//! - Bidirectional FM-index support (reversed text)
//!
//! This utility is intended to be run as a CLI application, controlled by the [`Arguments`] struct.

use sa_builder::{translate_l_to_i};
use succinct::storage::BlockType;
use std::error::Error;
use std::path::Path;
use std::io::BufWriter;
use std::fs::{File};
use succinct::{BitVec, BitVecMut, BitVector};
use byteorder::LittleEndian;
use std::mem::drop;
use qwt::QWT256;

use bincode::serialize_into;

use clap::Parser;

/// Command-line arguments for the FM-index builder.
///
/// This struct defines the parameters accepted from the command line, such as input
/// and output file paths, and the sparseness factor for suffix array sampling.
#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used. All the proteins are expected to be
    /// concatenated.
    #[arg(short, long)]
    pub database_file: String,
    /// Output location where to store the FM-index
    #[arg(short, long)]
    pub output: String,
    /// The sparseness_factor used on the suffix array samples
    #[arg(short, long, default_value_t = 8)]
    pub sparseness_factor: u8
}

/// Build a FM-index from the given text
///
/// # Arguments
/// * `text` - The text on which we want to build the FM-index
/// * `sparseness_factor` - The sparseness factor used on the suffix array
///
/// # Returns
///
/// Returns the constructed FM-index
///
/// # Errors
///
/// The errors that occurred during the building of the suffix array itself
pub fn build_fm(
    mut text: Vec<u8>,
    sparseness_factor: u8,
    output_path: &Path
) -> Result<(), Box<dyn Error>> {
    eprintln!("Building bidirectional FM-index");

    eprintln!("\tTranslating L to I in text...");
    translate_l_to_i(&mut text);

    eprintln!("\tTransforming text and building alphabet file");
    let (char_to_id, alph_size) = transform_text(&mut text);
    let alph_file = BufWriter::new(File::create(output_path.with_extension("alph"))?);
    let _ = serialize_into(alph_file, &char_to_id).unwrap();
    eprintln!("\t\tWritten {}", output_path.with_extension("alph").to_str().unwrap());

    eprintln!("\tBuilding counts table");
    let counts = build_counts(&text, alph_size);
    let counts_file = BufWriter::new(File::create(output_path.with_extension("counts"))?);
    let _ = serialize_into(counts_file, &counts).unwrap();
    eprintln!("\t\tWritten {}", output_path.with_extension("counts").to_str().unwrap());

    eprintln!("\tBuilding suffix array");
    let mut sa: Vec<i64> = libsais64_rs::sais64(text.clone(), 1)?;
    
    eprintln!("\tBuilding BWT...");
    let bwt = build_bwt(&text, &sa);
    let bwt_file = BufWriter::new(File::create(output_path.with_extension("bwt"))?);
    let _ = serialize_into(bwt_file, &bwt).unwrap();
    eprintln!("\t\tWritten {}", output_path.with_extension("bwt").to_str().unwrap());
    drop(bwt);

    eprintln!("\tSampling suffix array");
    let ssa_occs = sample_sa(&mut sa, sparseness_factor);
    let sa_file: BufWriter<File> = BufWriter::new(File::create(output_path.with_extension("ssa"))?);
    let _ = serialize_into(sa_file, &sa);
    eprintln!("\t\tWritten {}", output_path.with_extension("ssa").to_str().unwrap());
    let mut ssa_occs_file: BufWriter<File> = BufWriter::new(File::create(output_path.with_extension("ssa_occ"))?);
    for i in 0..ssa_occs.block_len() {
        let _ = ssa_occs.get_block(i).write_block::<_, LittleEndian>(&mut ssa_occs_file);
    }
    eprintln!("\t\tWritten {}", output_path.with_extension("ssa_occ").to_str().unwrap());
    drop(sa);
    drop(ssa_occs);

    eprintln!("\tReversing text...");
    text.reverse();

    eprintln!("\tBuilding suffix array of reversed text...");
    let sa_rev: Vec<i64> = libsais64_rs::sais64(text.clone(), 1)?;

    eprintln!("\tBuilding BWT of reversed text...");
    let bwt_rev = build_bwt(&text, &sa_rev);
    let bwt_rev_file = BufWriter::new(File::create(output_path.with_extension("rev.bwt"))?);
    let _ = serialize_into(bwt_rev_file, &bwt_rev).unwrap();
    eprintln!("\t\tWritten {}", output_path.with_extension("rev.bwt").to_str().unwrap());

    Ok(())
}

/// Sample the suffix array with the given sparseness factor
///
/// # Arguments
/// * `sa` - The suffix array that we want to sample
/// * `sparseness_factor` - The sparseness factor used for sampling
///
/// # Returns
///
/// The sampled suffix array
pub fn sample_sa(sa: &mut Vec<i64>, sparseness_factor: u8) -> BitVector {
    let mut ssa_occs = BitVector::with_fill(sa.len() as u64, false);

    let mut current_sampled_index = 0;
    for i in 0..sa.len() {
        let current_sa_val = sa[i];
        if current_sa_val % sparseness_factor as i64 == 0 {
            ssa_occs.set_bit(i as u64, true);
            sa[current_sampled_index] = current_sa_val;
            current_sampled_index += 1;
        }
    }

    // make shorter
    sa.resize(current_sampled_index, 0);

    ssa_occs
}

/// Constructs the Burrows-Wheeler Transform (BWT) from a suffix array.
///
/// # Arguments
/// * `text` - Original input text.
/// * `sa` - Suffix array built from the text.
///
/// # Returns
/// A `QWT256<u8>` object representing the BWT.
fn build_bwt(text: &Vec<u8>, sa: &Vec<i64>) -> QWT256<u8> {
    let bwt: Vec<u8> = sa.iter()
        .map(|&i| if i == 0 { text[text.len() - 1] } else { text[i as usize - 1] })
        .collect();

    QWT256::from(bwt)
}

/// Builds the count table used for FM-index rank support.
///
/// # Arguments
/// * `text` - Transformed input text.
/// * `alph_size` - Size of the reduced alphabet.
///
/// # Returns
/// A vector where `counts[i]` is the number of symbols in the text less than `i`.
fn build_counts(text: &Vec<u8>, alph_size: u8) -> Vec<usize> {
    let mut freqs = [0usize; 256];
    for &c in text {
        freqs[c as usize] += 1;
    }

    let mut counts = Vec::with_capacity(alph_size as usize);
    let mut total = 0;
    for i in 0..alph_size {
        counts.push(total);
        total += freqs[i as usize];
    }

    counts
}

/// Transforms the input text into a reduced integer alphabet and remaps characters.
///
/// # Arguments
/// * `text` - The input text to be transformed (modified in-place).
///
/// # Returns
/// A tuple of (char_to_id, alphabet_size), where `char_to_id` maps original characters
/// to their compressed form.
fn transform_text(text: &mut Vec<u8>) -> (Vec<u8>, u8) {
    let mut occs = [false; 256];
    for &c in text.iter() {
        occs[c as usize] = true;
    }

    let mut char_to_id: Vec<u8> = vec![0;256];
    let mut id = 0;
    for i in 0..256 {
        if occs[i] {
            char_to_id[i] = id;
            id += 1;
        }
    }

    for c in text.iter_mut() {
        *c = char_to_id[*c as usize];
    }

    (char_to_id, id)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use qwt::AccessUnsigned;

    #[test]
    fn test_transform_text_basic() {
        let mut text = b"BANANA$".to_vec();
        let (char_to_id, alph_size) = transform_text(&mut text);

        assert_eq!(alph_size, 4); // B, A, N, $
        for &c in &text {
            assert!(c < alph_size);
        }

        // Check that mapping is consistent
        let mut seen = HashMap::new();
        let mapped_chars = &[char_to_id[b'B' as usize], char_to_id[b'A' as usize], char_to_id[b'N' as usize], char_to_id[b'$' as usize]];
        for (orig, &mapped) in b"BAN$".iter().zip(mapped_chars) {
            seen.insert(*orig, mapped);
        }

        for (i, &c) in b"BANANA".iter().enumerate() {
            assert_eq!(text[i], seen[&c]);
        }
    }

    #[test]
    fn test_build_counts_correctness() {
        let text = vec![2, 0, 1, 1, 0, 2]; // Already transformed
        let counts = build_counts(&text, 3);

        assert_eq!(counts, vec![0, 2, 4]);
    }

    #[test]
    fn test_sample_sa_correctness() {
        let mut sa = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let sparseness_factor = 3;

        let sampled = sample_sa(&mut sa, sparseness_factor);

        // Should retain only values divisible by 3
        assert_eq!(sa, vec![0, 3, 6, 9]);

        for i in 0..10 {
            let expected = i % 3 == 0;
            assert_eq!(sampled.get_bit(i as u64), expected);
        }
    }

    #[test]
    fn test_build_bwt_basic() {
        let text = b"BANANA$".to_vec();
        let sa = vec![6, 5, 3, 1, 0, 4, 2]; // Suffix array for "BANANA$"
        let bwt: QWT256<u8> = build_bwt(&text, &sa);
        let bwt_vec: Vec<u8> = (0..bwt.len()).map(|i| bwt.get(i).unwrap()).collect();

        assert_eq!(bwt_vec, vec![b'A', b'N', b'N', b'B', b'$', b'A', b'A']);
    }
}