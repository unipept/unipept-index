//! The shared TSV fixture, used by both backends' test modules.
//!
//! Both build their index from the same rows, so their assertions are directly comparable — which
//! is what lets `matches_the_preloaded_backend_field_for_field` in `mmap` write a file with one
//! backend and read it back with the other. Each test picks how many rows it wants.

use std::{fs::File, io::Write, path::PathBuf};

use tempfile::TempDir;

/// `(uniprot_id, taxon_id, sequence, annotations)`, in the order they appear in the file.
pub(crate) const TEST_PROTEINS: [(&str, u32, &str, &str); 4] = [
    ("P12345", 1, "MLPGLALLLLAAWTARALEV", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P54321", 2, "PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P67890", 6, "KWDSDPSGTKTCIDT", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P13579", 17, "KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH", "GO:0009279;IPR:IPR016364;IPR:IPR008816")
];

/// Writes `proteins` as a UniProt TSV into `tmp_dir` and returns its path.
pub(crate) fn write_database_file(tmp_dir: &TempDir, proteins: &[(&str, u32, &str, &str)]) -> PathBuf {
    let path = tmp_dir.path().join("database.tsv");
    let mut f = File::create(&path).unwrap();
    for (uid, taxon, sequence, annotations) in proteins {
        writeln!(f, "{uid}\t{taxon}\t{sequence}\t{annotations}").unwrap();
    }
    path
}
