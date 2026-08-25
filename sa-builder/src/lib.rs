//! Builds the on-disk index that `sa-server` reads: the suffix array, the protein store, the
//! suffix-to-protein mapping and an optional k-mer bounds table.
//!
//! This is the only writer in the workspace. Every format the other crates read is produced here,
//! through the `WriteBinary` implementations documented in `binary-traits` — so a format question
//! is answered at its writer, not here.
//!
//! One run produces one index. The four files are only usable as a set: the suffix array indexes
//! positions in the text stored in `proteins.bin`, and the mapping resolves those same positions
//! to entries in it, so mixing files from different builds yields wrong answers rather than
//! errors. The binary therefore writes every section to a temporary sibling and renames them all
//! only once the last one has succeeded; see `main.rs`.
//!
//! This crate holds the pieces worth testing on their own — the command line in [`Arguments`] and
//! suffix array construction in [`build_ssa`]. The file writing lives in the binary.
#![warn(missing_docs)]

use std::error::Error;

use clap::{Parser, ValueEnum};
use sa_index::kmer_table::MAX_KMER_K;

/// Build a (sparse, compressed) suffix array from the given text
#[derive(Parser, Debug)]
pub struct Arguments {
    /// Tab-separated file with one protein per row: `uniprot_id`, `taxon_id`, `sequence` and
    /// `annotations`. The sequences are concatenated into a single text internally, separated by
    /// `-` and terminated by `$`; neither character may appear in a sequence.
    #[arg(short, long)]
    pub database_file: String,
    /// Output location where to store the suffix array
    #[arg(long)]
    pub output_sa: String,
    /// Output location where to store the proteins binary
    #[arg(long)]
    pub output_proteins: String,
    /// Output location where to store the suffix-to-protein mapping binary
    #[arg(long)]
    pub output_mapping: String,
    /// The sparseness factor used on the suffix array: only every n-th text position is indexed
    /// (default 1, which indexes every position).
    ///
    /// Larger values shrink the suffix array at the cost of search work. Peptides shorter than n
    /// cannot be found at all.
    // Rejected below 1 rather than left to the arithmetic, because every step downstream accepts
    // zero and produces a well-formed index that matches nothing. `0.is_multiple_of(5)` is true,
    // so `libsais64` keeps `libsais_sparseness == MAX_SPARSENESS` and derives a `sample_rate` of
    // `0 / 5 == 0`, which skips the second sampling pass; the zero then reaches the file header
    // unchanged through `dump_suffix_array`. The array is built at a stride of 5 while the header
    // claims 0, and a claimed stride of zero means no peptide is ever long enough to search. The
    // build reports success throughout.
    #[arg(short, long, default_value_t = 1, value_parser = clap::value_parser!(u8).range(1..))]
    pub sparseness_factor: u8,
    /// The algorithm used to construct the suffix array (default value LibSais)
    #[arg(short('a'), long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    pub construction_algorithm: SAConstructionAlgorithm,
    /// If the suffix array should be compressed (default value false)
    #[arg(short, long, default_value_t = false)]
    pub compress_sa: bool,
    /// The style of suffix-to-protein mapping to build (default value BitVec)
    #[arg(long, value_enum, default_value_t = SuffixToProteinMappingStyle::BitVec)]
    pub mapping_style: SuffixToProteinMappingStyle,
    /// Output location where to store the k-mer bounds table (optional).
    /// When set, a k-mer lookup table is built and written to this path.
    #[arg(long)]
    pub output_kmer_table: Option<String>,
    /// The k-mer size used when building the k-mer bounds table (default 5, maximum 7).
    ///
    /// The table is dense at 24^k entries of 16 bytes, so its size depends only on k and not on
    /// the database: k=5 is 127 MB (0.12 GiB) and k=6 is 3.06 GB (2.85 GiB), a 24x step for one
    /// more level of the probe chain.
    ///
    /// 5 is the default because it is the size that pays in both regimes. With the index fully
    /// resident, the 6-mer's edge over it sits inside the noise floor on most length regimes and
    /// reaches only +4.1% on large peptides, which does not buy 2.9 GB. Raise it to 6 only under
    /// a memory ceiling, where the table's value is working-set size rather than probe count: a
    /// 5-mer narrows the search to ~7 SA pages per query and a 6-mer to ~1, which is +18.4%
    /// against no table where the 5-mer manages +3.2%.
    // The upper bound is enforced here rather than left to the `k <= MAX_KMER_K` assertion inside
    // `KmerTable::build_kmer_table`, which does not run until the suffix array has already been
    // built — hours into a full-database build. The bound itself is taken from `sa_index` so the
    // two cannot drift. The lower bound matters for the same reason `--sparseness-factor` has one:
    // k=0 indexes every suffix into a single bucket, producing a table that narrows nothing and
    // reports no error.
    #[arg(long, default_value_t = 5,
          value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=MAX_KMER_K as u64))]
    pub kmer_size: usize
}

/// The library used to construct the suffix array.
///
/// Both produce the same array; they differ in how the sparseness factor is applied. `libsais`
/// can sample during construction and so never materialises the full array, which is what makes
/// a sparse build of a UniProt-sized database feasible.
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SAConstructionAlgorithm {
    /// `libdivsufsort`. Always builds the dense array first and samples it afterwards, so peak
    /// memory is that of the dense array whatever the sparseness factor.
    LibDivSufSort,
    /// `libsais`. Samples during construction where the factor allows it. The default.
    LibSais
}

/// Which suffix-to-protein mapping is written to `--output-mapping`.
///
/// All three answer the same question — which protein contains this text position — and the
/// style is recorded in the file, so the server picks its reader from what it finds rather than
/// from configuration. The choice is purely a space against lookup-cost trade made here, and
/// `sa_index::suffix_to_protein_index` documents the three representations in full.
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SuffixToProteinMappingStyle {
    /// One `u32` per text position: a single load per lookup, at 4 bytes per residue — 0.84 GB
    /// over a 209 M-position text, and ~256 GB at full UniProt scale, which exceeds the whole
    /// rest of the index. For small databases only.
    Dense,
    /// The start position of each protein, binary-searched. Smallest, at O(log m) dependent
    /// loads per lookup for m proteins, each likely a cache miss.
    Sparse,
    /// A bit per text position marking the separators and the terminator, with a rank structure
    /// over it. Near-dense speed at ~1.25 bits per position. The default.
    BitVec
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
    mut text: Vec<u8>,
    construction_algorithm: &SAConstructionAlgorithm,
    sparseness_factor: u8
) -> Result<Vec<i64>, Box<dyn Error>> {
    // translate all L's to a I
    translate_l_to_i(&mut text);

    // Build the suffix array using the selected algorithm
    let mut sa = match construction_algorithm {
        SAConstructionAlgorithm::LibSais => libsais64(text, sparseness_factor)?,
        SAConstructionAlgorithm::LibDivSufSort => {
            libdivsufsort_rs::divsufsort64(&text).ok_or("Building suffix array failed")?
        }
    };

    // make the SA sparse and decrease the vector size if we have sampling (sampling_rate > 1)
    if *construction_algorithm == SAConstructionAlgorithm::LibDivSufSort {
        sample_sa(&mut sa, sparseness_factor);
    }

    Ok(sa)
}

/// The largest number of text characters `libsais` may fold into a single symbol.
///
/// `libsais` allocates one bucket per symbol, and a symbol spanning `s` characters of the 5-bit
/// protein alphabet needs `2^(5*s)` of them. At `s = 5` that is `2^25` buckets, and every further
/// step multiplies it by 32.
const MAX_SPARSENESS: usize = 5;

/// Builds the suffix array with `libsais`, splitting the sparseness factor between the sampling
/// `libsais` performs itself and a second pass over the result.
fn libsais64(text: Vec<u8>, sparseness_factor: u8) -> Result<Vec<i64>, &'static str> {
    let sparseness_factor = sparseness_factor as usize;

    // libsais can only sample at a factor it folds into its symbols, so take the largest such
    // factor that divides the requested sparseness and leave the remainder to `sample_sa`. The
    // walk down from MAX_SPARSENESS always terminates because 1 divides everything.
    let mut libsais_sparseness = MAX_SPARSENESS;
    while !sparseness_factor.is_multiple_of(libsais_sparseness) {
        libsais_sparseness -= 1;
    }
    let sample_rate = sparseness_factor / libsais_sparseness;

    let mut sa = libsais64_rs::sais64(text, libsais_sparseness)?;

    if sample_rate > 1 {
        sample_sa(&mut sa, sample_rate as u8);
    }

    Ok(sa)
}

/// Translate all L's to I's in the given text, in place
///
/// Leucine and isoleucine have the same mass, so a mass-spectrometry search cannot tell them
/// apart. Folding them together here means the suffix array never contains an `L` at all.
///
/// # Arguments
/// * `text` - The text in which we want to translate the L's to I's
fn translate_l_to_i(text: &mut [u8]) {
    for character in text.iter_mut() {
        if *character == b'L' {
            *character = b'I'
        }
    }
}

/// Sample the suffix array with the given sparseness factor, in place
///
/// Keeps only the entries whose text position is a multiple of the factor, compacting them to the
/// front and truncating `sa` to what survived.
///
/// # Arguments
/// * `sa` - The suffix array that we want to sample
/// * `sparseness_factor` - The sparseness factor used for sampling
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
        let args = Arguments::parse_from([
            "sa-builder",
            "--database-file",
            "database.fa",
            "--output-sa",
            "output.fa",
            "--output-proteins",
            "output.proteins",
            "--sparseness-factor",
            "2",
            "--construction-algorithm",
            "lib-div-suf-sort",
            "--compress-sa",
            "--output-mapping",
            "output.mapping",
            "--mapping-style",
            "dense",
            "--output-kmer-table",
            "output.kmer",
            "--kmer-size",
            "5"
        ]);

        assert_eq!(args.database_file, "database.fa");
        assert_eq!(args.output_sa, "output.fa");
        assert_eq!(args.output_proteins, "output.proteins");
        assert_eq!(args.sparseness_factor, 2);
        assert_eq!(args.construction_algorithm, SAConstructionAlgorithm::LibDivSufSort);
        assert!(args.compress_sa);
        assert_eq!(args.output_mapping, "output.mapping".to_string());
        assert_eq!(args.mapping_style, SuffixToProteinMappingStyle::Dense);
        assert_eq!(args.output_kmer_table, Some("output.kmer".to_string()));
        assert_eq!(args.kmer_size, 5);
    }

    /// The k-mer table is optional; omitting its path must leave it unbuilt.
    #[test]
    fn test_arguments_without_kmer_table() {
        let args = Arguments::parse_from([
            "sa-builder",
            "--database-file",
            "database.fa",
            "--output-sa",
            "output.fa",
            "--output-proteins",
            "output.proteins",
            "--output-mapping",
            "output.mapping"
        ]);

        assert_eq!(args.output_kmer_table, None);
        assert_eq!(args.sparseness_factor, 1, "default sparseness");
        assert_eq!(args.kmer_size, 5, "default k-mer size");
        assert!(!args.compress_sa, "compression is opt-in");
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
        let text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibSais, 1).unwrap();
        assert_eq!(sa, vec![11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]);
    }

    #[test]
    fn test_build_ssa_libsais_empty() {
        let text = b"".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibSais, 1).unwrap();
        assert_eq!(sa, Vec::<i64>::new());
    }

    #[test]
    fn test_build_ssa_libsais_sparse() {
        let text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibSais, 2).unwrap();
        assert_eq!(sa, vec![10, 0, 8, 4, 6, 2]);
    }

    #[test]
    fn test_build_ssa_libdivsufsort() {
        let text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibDivSufSort, 1).unwrap();
        assert_eq!(sa, vec![11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2]);
    }

    #[test]
    fn test_build_ssa_libdivsufsort_empty() {
        let text = b"".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibDivSufSort, 1).unwrap();
        assert_eq!(sa, Vec::<i64>::new());
    }

    #[test]
    fn test_build_ssa_libdivsufsort_sparse() {
        let text = b"ABRACADABRA$".to_vec();
        let sa = build_ssa(text, &SAConstructionAlgorithm::LibDivSufSort, 2).unwrap();
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
