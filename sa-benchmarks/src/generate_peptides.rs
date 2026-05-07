use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rand::Rng;
use sa_server::load_proteins_file;
use text_compression::ProteinTextBackend as _;

/// Characters that delimit protein boundaries in the text.
/// b'-' is SEPARATION_CHARACTER, b'$' is TERMINATION_CHARACTER.
const BOUNDARY_CHARS: [u8; 2] = [b'-', b'$'];

/// Sample real protein subsequences from proteins.bin for use as benchmark peptides.
///
/// The output file will contain `--amount` peptides (one per line), each drawn
/// from an actual protein in the index so that all peptides are guaranteed hits.
#[derive(Parser, Debug)]
#[command(about = "Sample real protein subsequences from proteins.bin for use as benchmark peptides")]
struct Args {
    /// Folder containing proteins.bin
    #[arg(short, long)]
    index_dir: PathBuf,

    /// Output file (one peptide per line)
    #[arg(short, long)]
    output_file: PathBuf,

    /// Number of peptides to sample
    #[arg(long, default_value_t = 10_000)]
    amount: usize,

    /// Minimum peptide length
    #[arg(long, default_value_t = 5)]
    min_len: usize,

    /// Maximum peptide length
    #[arg(long, default_value_t = 50)]
    max_len: usize,

}

/// Scans `text` once and returns all protein runs as `(start, len)` pairs,
/// keeping only runs whose length is at least `min_len`.
fn collect_protein_runs(
    text: &text_compression::ProteinText,
    min_len: usize,
) -> Vec<(usize, usize)> {
    let total = text.len();
    let mut runs = Vec::new();
    let mut run_start = 0usize;

    for i in 0..total {
        let c = text.get(i);
        if BOUNDARY_CHARS.contains(&c) {
            let run_len = i - run_start;
            if run_len >= min_len {
                runs.push((run_start, run_len));
            }
            run_start = i + 1;
        }
    }

    // Handle any trailing run before the terminator (if the text doesn't end with '$')
    let trailing = total - run_start;
    if trailing >= min_len {
        runs.push((run_start, trailing));
    }

    runs
}

fn sample_peptides(
    text: &text_compression::ProteinText,
    amount: usize,
    min_len: usize,
    max_len: usize,
) -> Vec<String> {
    let runs = collect_protein_runs(text, min_len);

    // Prefix sums over the number of valid start positions per run.
    // A run of length L has (L - min_len + 1) valid start positions.
    let prefix_sums: Vec<usize> = runs
        .iter()
        .scan(0usize, |acc, (_, len)| {
            *acc += len - min_len + 1;
            Some(*acc)
        })
        .collect();

    let total_valid_starts = *prefix_sums.last().unwrap_or(&0);
    assert!(
        total_valid_starts > 0,
        "No valid start positions found — the index may be too small for min_len = {}",
        min_len
    );

    let mut rng = rand::thread_rng();
    let mut peptides = Vec::with_capacity(amount);

    for _ in 0..amount {
        // Pick a uniform random position among all valid starts.
        let pos = rng.gen_range(0..total_valid_starts);

        // Map pos → (run_index, offset_within_run) via binary search on prefix sums.
        let run_idx = prefix_sums.partition_point(|&s| s <= pos);
        let run_base = if run_idx == 0 { 0 } else { prefix_sums[run_idx - 1] };
        let offset = pos - run_base;

        let (run_start, run_len) = runs[run_idx];
        let start = run_start + offset;
        let available = run_len - offset;
        let len = rng.gen_range(min_len..=max_len.min(available));

        let peptide: String = (0..len)
            .map(|i| text.get(start + i) as char)
            .collect();
        peptides.push(peptide);
    }

    peptides
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    if args.min_len > args.max_len {
        return Err("--min-len must be <= --max-len".into());
    }
    if args.min_len == 0 {
        return Err("--min-len must be >= 1".into());
    }

    let proteins_path = args.index_dir.join("proteins.bin");
    eprintln!("Loading proteins from {}...", proteins_path.display());
    let proteins = load_proteins_file(proteins_path.to_str().unwrap())?;

    let text = proteins.text();
    eprintln!("  Text length: {} characters", text.len());
    eprintln!("Sampling {} peptides (length {}-{})...", args.amount, args.min_len, args.max_len);

    let peptides = sample_peptides(text, args.amount, args.min_len, args.max_len);

    let mut out = BufWriter::new(File::create(&args.output_file)?);
    for p in &peptides {
        writeln!(out, "{}", p)?;
    }

    eprintln!("Wrote {} peptides to {}", peptides.len(), args.output_file.display());
    Ok(())
}
