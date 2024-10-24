use std::error::Error;
use clap::{Parser, ValueEnum};

/// Build a (sparse, compressed) suffix array from the given text
#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be
    /// concatenated using a hashtag `#`.
    #[arg(short, long)]
    pub database_file: String,
    /// Output location where to store the suffix array
    #[arg(short, long)]
    pub output: String,
    /// The sparseness_factor used on the suffix array (default value 1, which means every value in
    /// the SA is used)
    #[arg(short, long, default_value_t = 1)]
    pub sparseness_factor: u8,
    /// The algorithm used to construct the suffix array (default value LibSais)
    #[arg(short('a'), long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    pub construction_algorithm: SAConstructionAlgorithm,
    /// If the suffix array should be compressed (default value true)
    #[arg(short, long, default_value_t = false)]
    pub compress_sa: bool
}

/// Enum representing the two possible algorithms to construct the suffix array
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SAConstructionAlgorithm {
    LibDivSufSort,
    LibSais
}

/// Build a sparse suffix array from the given text
///
/// # Arguments
/// * `text` - The text on which we want to build the suffix array
/// * `construction_algorithm` - The algorithm used during construction
/// * `sparseness_factor` - The sparseness factor used on the suffix array
///
/// # Returns
///
/// Returns the constructed (sparse) suffix array
///
/// # Errors
///
/// The errors that occurred during the building of the suffix array itself
pub fn build_ssa(
    text: &mut Vec<u8>,
    construction_algorithm: &SAConstructionAlgorithm,
    sparseness_factor: u8
) -> Result<Vec<i64>, Box<dyn Error>> {
    // translate all L's to a I
    translate_l_to_i(text);

    // Build the suffix array using the selected algorithm
    let mut sa = match construction_algorithm {
        SAConstructionAlgorithm::LibSais => libsais64(text, sparseness_factor)?,
        SAConstructionAlgorithm::LibDivSufSort => libdivsufsort_rs::divsufsort64(text).ok_or("Building suffix array failed")?
    };

    // make the SA sparse and decrease the vector size if we have sampling (sampling_rate > 1)
    if *construction_algorithm == SAConstructionAlgorithm::LibDivSufSort {
        sample_sa(&mut sa, sparseness_factor);
    }

    Ok(sa)
}

// Max sparseness for libsais because it creates a bucket for each element of the alphabet (2 ^ (sparseness * bits_per_char) buckets).
const MAX_SPARSENESS: usize = 5;
fn libsais64(text: &Vec<u8>, sparseness_factor: u8) -> Result<Vec<i64>, &str> {
    let sparseness_factor = sparseness_factor as usize;

    // set libsais_sparseness to highest sparseness factor fitting in 32-bit value and sparseness factor divisible by libsais sparseness
    // max 28 out of 32 bits used, because a bucket is created for every element of the alfabet 8 * 2^28).
    let mut libsais_sparseness = MAX_SPARSENESS;
    while sparseness_factor % libsais_sparseness != 0 {
        libsais_sparseness -= 1;
    }
    let sample_rate = sparseness_factor / libsais_sparseness;
    eprintln!("\tSparseness factor: {}", sparseness_factor);
    eprintln!("\tLibsais sparseness factor: {}", libsais_sparseness);
    eprintln!("\tSample rate: {}", sample_rate);

    let mut sa = libsais64_rs::sais64(text, libsais_sparseness)?;

    if sample_rate > 1 {
        sample_sa(&mut sa, sample_rate as u8);
    }

    Ok(sa)
}

/// Translate all L's to I's in the given text
///
/// # Arguments
/// * `text` - The text in which we want to translate the L's to I's
///
/// # Returns
///
/// The text with all L's translated to I's
fn translate_l_to_i(text: &mut [u8]) {
    for character in text.iter_mut() {
        if *character == b'L' {
            *character = b'I'
        }
    }
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
fn sample_sa(sa: &mut Vec<i64>, sparseness_factor: u8) {
    if sparseness_factor <= 1 {
        return;
    }

    let mut current_sampled_index = 0;
    for i in 0..sa.len() {
        let current_sa_val = sa[i];
        if current_sa_val % sparseness_factor as i64 == 0 {
            sa[current_sampled_index] = current_sa_val;
            current_sampled_index += 1;
        }
    }

    // make shorter
    sa.resize(current_sampled_index, 0);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arguments() {
        let args = Arguments::parse_from(&[
            "sa-builder",
            "--database-file",
            "database.fa",
            "--output",
            "output.fa",
            "--sparseness-factor",
            "2",
            "--construction-algorithm",
            "lib-div-suf-sort",
            "--compress-sa"
        ]);

        assert_eq!(args.database_file, "database.fa");
        assert_eq!(args.output, "output.fa");
        assert_eq!(args.sparseness_factor, 2);
        assert_eq!(args.construction_algorithm, SAConstructionAlgorithm::LibDivSufSort);
        assert_eq!(args.compress_sa, true);
    }

    #[test]
    fn test_sa_construction_algorithm() {
        assert_eq!(
            SAConstructionAlgorithm::from_str("lib-div-suf-sort", false),
            Ok(SAConstructionAlgorithm::LibDivSufSort)
        );
        assert_eq!(SAConstructionAlgorithm::from_str("lib-sais", false), Ok(SAConstructionAlgorithm::LibSais));
    }

    #[test]
    fn test_build_ssa_libsais() {
        let mut text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibSais, 1).unwrap();
        assert_eq!(sa, vec![11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]);
    }

    #[test]
    fn test_build_ssa_libsais_empty() {
        let mut text = b"".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibSais, 1).unwrap();
        assert_eq!(sa, vec![]);
    }

    #[test]
    fn test_build_ssa_libsais_sparse() {
        let mut text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibSais, 2).unwrap();
        assert_eq!(sa, vec![10, 0, 8, 4, 6, 2]);
    }

    #[test]
    fn test_build_ssa_libdivsufsort() {
        let mut text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibDivSufSort, 1).unwrap();
        assert_eq!(sa, vec![11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]);
    }

    #[test]
    fn test_build_ssa_libdivsufsort_empty() {
        let mut text = b"".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibDivSufSort, 1).unwrap();
        assert_eq!(sa, vec![]);
    }

    #[test]
    fn test_build_ssa_libdivsufsort_sparse() {
        let mut text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(&mut text, &SAConstructionAlgorithm::LibDivSufSort, 2).unwrap();
        assert_eq!(sa, vec![10, 0, 8, 4, 6, 2]);
    }

    #[test]
    fn test_translate_l_to_i() {
        let mut text = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ$-".to_vec();
        translate_l_to_i(&mut text);
        assert_eq!(text, b"ABCDEFGHIJKIMNOPQRSTUVWXYZ$-".to_vec());
    }

    #[test]
    fn test_sample_sa_1() {
        let mut sa = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        sample_sa(&mut sa, 1);
        assert_eq!(sa, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_sample_sa_2() {
        let mut sa = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        sample_sa(&mut sa, 2);
        assert_eq!(sa, vec![0, 2, 4, 6, 8]);
    }
}
