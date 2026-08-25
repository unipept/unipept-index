use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    path::PathBuf,
    time::{SystemTime, SystemTimeError, UNIX_EPOCH}
};

use binary_traits::WriteBinary;
use clap::Parser;
use sa_builder::{Arguments, SuffixToProteinMappingStyle, build_ssa};
use sa_index::{
    KmerTable,
    array::{dump_compressed_suffix_array, dump_suffix_array},
    suffix_to_protein_index::{BitVecSuffixToProtein, DenseSuffixToProtein, SparseSuffixToProtein}
};
use sa_mappings::proteins::{InMemoryProteins, ProteinsBackend as _};
use text_compression::ProteinTextBackend as _;

fn main() {
    let Arguments {
        database_file,
        output_sa,
        output_proteins,
        sparseness_factor,
        construction_algorithm,
        compress_sa,
        output_mapping,
        mapping_style,
        output_kmer_table,
        kmer_size
    } = Arguments::parse();

    let proteins = timed("loading the proteins", || {
        InMemoryProteins::load_from_tsv(&database_file).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()))
    });
    let bits_per_value = (proteins.text().len() as f64).log2().ceil() as usize;

    // Not the most efficient way to get the text, but still acceptable
    let text: Vec<u8> = proteins.text().iter().collect();
    let sa = timed("building the suffix array", || {
        build_ssa(text, &construction_algorithm, sparseness_factor)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()))
    });
    eprintln!("\tAmount of items: {}", sa.len());

    // Build the k-mer table while the SA Vec is still in memory (before it is consumed below).
    let kmer_table: Option<KmerTable> = output_kmer_table.as_ref().map(|_| {
        timed(&format!("building k-mer table (k={})", kmer_size), || {
            KmerTable::build_from_raw_sa(&sa, proteins.text().len(), |i| proteins.text().get(i), kmer_size)
        })
    });

    // Every section is written to a temporary sibling and renamed only once all of them have
    // succeeded, so a failed build leaves the previous index untouched. See `commit_all`.
    let mut pending: Vec<PendingRename> = Vec::new();
    let mut wrote: Vec<String> = Vec::new();

    // open the output file
    let mut sa_file =
        Output::open(&output_sa, 100 * 1024 * 1024).unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("dumping the suffix array", || {
        eprintln!("\tAmount of bits per item: {}", if compress_sa { bits_per_value } else { 64 });
        if compress_sa {
            if let Err(err) = dump_compressed_suffix_array(sa, sparseness_factor, bits_per_value, sa_file.writer()) {
                eprint_and_exit(err.to_string().as_str());
            };
        } else {
            if let Err(err) = dump_suffix_array(sa, sparseness_factor, sa_file.writer()) {
                eprint_and_exit(err.to_string().as_str());
            }
        }
    });
    pending.push(sa_file.seal().unwrap_or_else(|err| eprint_and_exit(&err)));
    wrote.push(output_sa);

    let mut mapping_file = Output::open(&output_mapping, 100 * 1024 * 1024)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("writing suffix-to-protein mapping binary", || {
        let text_len = proteins.text().len();
        // The three arms differ only in the concrete mapping type; each writes its own type byte,
        // which is what `read_binary` later dispatches on. A generic helper would need the closure
        // to be nameable, so the repetition is left explicit.
        let char_at = |i: usize| proteins.text().get(i);
        let result = match &mapping_style {
            SuffixToProteinMappingStyle::Dense => {
                DenseSuffixToProtein::from_text_parts(text_len, char_at).write_binary(mapping_file.writer())
            }
            SuffixToProteinMappingStyle::Sparse => {
                SparseSuffixToProtein::from_text_parts(text_len, char_at).write_binary(mapping_file.writer())
            }
            SuffixToProteinMappingStyle::BitVec => {
                BitVecSuffixToProtein::from_text_parts(text_len, char_at).write_binary(mapping_file.writer())
            }
        };
        result.unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    });
    pending.push(mapping_file.seal().unwrap_or_else(|err| eprint_and_exit(&err)));
    wrote.push(output_mapping);

    let mut proteins_file = Output::open(&output_proteins, 100 * 1024 * 1024)
        .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

    timed("writing proteins binary", || {
        proteins
            .write_binary(proteins_file.writer())
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
    });
    pending.push(proteins_file.seal().unwrap_or_else(|err| eprint_and_exit(&err)));
    wrote.push(output_proteins);

    if let (Some(table), Some(kmer_table_path)) = (kmer_table, output_kmer_table) {
        let mut kmer_file = Output::open(&kmer_table_path, 64 * 1024 * 1024)
            .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));

        timed("writing k-mer table", || {
            table
                .write_binary(kmer_file.writer())
                .unwrap_or_else(|err| eprint_and_exit(err.to_string().as_str()));
        });
        pending.push(kmer_file.seal().unwrap_or_else(|err| eprint_and_exit(&err)));
        wrote.push(kmer_table_path);
    }

    // Nothing has entered the index directory until this line.
    commit_all(pending).unwrap_or_else(|err| eprint_and_exit(&err));

    eprintln!();
    eprintln!("Wrote:");
    for path in wrote {
        eprintln!("\t{path}");
    }
}

/// An index file that only replaces its target once it is complete.
///
/// Writes land in a sibling `.tmp` file and are renamed into place by [`Output::finish`]. Two
/// separate problems that fixes:
///
/// * The builder used to `truncate(true)` the real path when it *opened* it, before writing a
///   byte. A build that then failed — a malformed input row, a full disk, a kill — had already
///   destroyed the previous index.
/// * The `proteins.bin` size guard fires after `sa.bin` and the mapping are written, so a database
///   too large for the offset format left the directory holding a mix of two builds.
///
/// A rename within a directory is atomic, so a reader sees either the whole old file or the whole
/// new one. A failed build leaves `.tmp` files behind on purpose: that is the cost of not touching
/// the index that still works.
struct Output {
    writer: BufWriter<File>,
    temporary: PathBuf,
    destination: PathBuf
}

impl Output {
    /// Opens `destination`'s temporary sibling, creating or truncating only that.
    fn open(destination: &str, buffer_size: usize) -> std::io::Result<Self> {
        let destination = PathBuf::from(destination);
        let mut temporary = destination.clone().into_os_string();
        temporary.push(".tmp");
        let temporary = PathBuf::from(temporary);

        let file = OpenOptions::new().create(true).write(true).truncate(true).open(&temporary)?;

        Ok(Self {
            writer: BufWriter::with_capacity(buffer_size, file),
            temporary,
            destination
        })
    }

    /// Borrows the writer for the duration of one `write_binary` call.
    fn writer(&mut self) -> &mut BufWriter<File> {
        &mut self.writer
    }

    /// Flushes and fsyncs the temporary file, returning the rename still owed.
    ///
    /// `BufWriter` flushes on drop but **discards the error**, so relying on that produced a
    /// silently truncated index that then loaded without complaint — the failure `binary_traits`
    /// requires `write_all` to prevent, reinstated one layer down. With buffers of 100 MB that is
    /// most of a section. `sync_all` follows because a successful flush only means the bytes
    /// reached the kernel, and an index written just before a crash is worth having on disk.
    ///
    /// The rename is deliberately *not* done here. Sealing frees the buffer as soon as a section is
    /// written, but nothing enters the index directory until every section has succeeded — see
    /// [`commit_all`].
    fn seal(self) -> Result<PendingRename, String> {
        let Self { mut writer, temporary, destination } = self;
        let name = destination.display();

        writer.flush().map_err(|err| format!("Could not write {name}: {err}"))?;
        writer
            .into_inner()
            .map_err(|err| format!("Could not write {name}: {}", err.error()))?
            .sync_all()
            .map_err(|err| format!("Could not flush {name} to disk: {err}"))?;

        Ok(PendingRename { temporary, destination })
    }
}

/// A written-and-synced temporary file waiting to be renamed over its destination.
struct PendingRename {
    temporary: PathBuf,
    destination: PathBuf
}

/// Moves every completed section into place, once all of them have been written.
///
/// The index is only usable as a set: `sa.bin` indexes positions in `proteins.bin`, and the mapping
/// resolves those positions to entries in it. Renaming each section as it finished would leave the
/// directory holding a *mix* of two builds whenever a later section failed — which is exactly what
/// the `proteins.bin` size guard does, since it fires after the suffix array and the mapping are
/// already written. Deferring every rename to here means a failed build changes nothing at all, and
/// the previous index keeps working.
///
/// The renames themselves are individually atomic but not atomic as a group; a crash *between* them
/// can still leave a mix. Making that impossible needs a directory swap, which is a bigger change
/// than this is worth — the window is microseconds of `rename` calls rather than the minutes a
/// build takes.
fn commit_all(pending: Vec<PendingRename>) -> Result<(), String> {
    for PendingRename { temporary, destination } in pending {
        std::fs::rename(&temporary, &destination).map_err(|err| {
            format!("Could not move {} into place at {}: {err}", temporary.display(), destination.display())
        })?;
    }
    Ok(())
}

fn eprint_and_exit(err: &str) -> ! {
    eprintln!("{}", err);
    std::process::exit(1);
}

fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}

fn timed<T, F: FnOnce() -> T>(msg: &str, f: F) -> T {
    eprintln!();
    eprintln!("📋 Started {}...", msg);
    let start = get_time_ms().unwrap();
    let result = f();
    eprintln!("✅ Successfully {} in {} seconds!", msg, (get_time_ms().unwrap() - start) / 1000.0);
    result
}
