use sa_builder::{translate_l_to_i};
use succinct::storage::BlockType;
use std::error::Error;
use std::path::Path;
use std::io::{BufWriter, Write};
use std::fs::{File};
use succinct::{BitVec, BitVecMut, BitVector};
use byteorder::LittleEndian;
use log::info;

use bincode::serialize_into;

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
    info!("------------Building bidirectional FM-index----------------");
    env_logger::init();

    // translate all L's to a I
    info!("Translating L to I in text...");
    translate_l_to_i(&mut text);

    info!("Building counts table");
    let counts: Vec<usize> = build_counts(&text);
    let counts_file = BufWriter::new(File::create(output_path.with_extension("counts"))?);
    let _ = serialize_into(counts_file, &counts).unwrap();
    info!("\tWritten {}", output_path.with_extension("counts").to_str().unwrap());

    info!("Building suffix array");
    let mut sa: Vec<i64> = libsais64_rs::sais64(text.clone(), 1)?;
    
    info!("Building BWT...");
    let bwt: Vec<u8> = build_bwt(&text, &sa);
    let bwt_file = BufWriter::new(File::create(output_path.with_extension("bwt"))?);
    bwt_file.into_inner()?.write_all(&bwt)?;
    info!("\tWritten {}", output_path.with_extension("bwt").to_str().unwrap());

    info!("Sampling suffix array");
    let ssa_occs = sample_sa(&mut sa, sparseness_factor);
    let sa_file: BufWriter<File> = BufWriter::new(File::create(output_path.with_extension("ssa"))?);
    let _ = serialize_into(sa_file, &sa);
    info!("\tWritten {}", output_path.with_extension("ssa").to_str().unwrap());
    let mut ssa_occs_file: BufWriter<File> = BufWriter::new(File::create(output_path.with_extension("ssa_occ"))?);
    for i in 0..ssa_occs.block_len() {
        let _ = ssa_occs.get_block(i).write_block::<_, LittleEndian>(&mut ssa_occs_file);
    }
    info!("\tWritten {}", output_path.with_extension("ssa_occ").to_str().unwrap());

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

fn build_bwt(text: &Vec<u8>, sa: &Vec<i64>) -> Vec<u8> {
    sa.iter()
        .map(|&i| if i == 0 { text[text.len() - 1] } else { text[i as usize - 1] })
        .collect()
}

fn build_counts(text: &Vec<u8>) -> Vec<usize> {
    let mut freqs = [0usize; 256];
    for &c in text {
        freqs[c as usize] += 1;
    }

    let mut counts = vec![0; 256];
    let mut total = 0;
    for i in 0..256 {
        counts[i] = total;
        total += freqs[i];
    }

    counts
}
