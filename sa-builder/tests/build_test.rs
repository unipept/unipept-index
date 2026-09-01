use std::{io::Write, process::Command};

use sa_index::{Nullable, ReadBinaryMmap, SuffixArray, suffix_to_protein_index::legacy::SuffixToProteinMapping};
use sa_mappings::proteins::Proteins;

/// Four proteins used as test input, matching the fixture in sa-mappings unit tests.
/// Text layout (L→I translation happens internally during SA construction):
///   pos   0-19  : MLPGLALLLLAAWTARALEV  (protein 0, taxon 1)
///   pos  20     : separator '-'
///   pos  21-50  : PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG  (protein 1, taxon 2)
///   pos  51     : separator '-'
///   pos  52-66  : KWDSDPSGTKTCIDT  (protein 2, taxon 6)
///   pos  67     : separator '-'
///   pos  68-112 : KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH  (protein 3, taxon 17)
///   pos 113     : terminator '$'
///   total length: 114
const PROTEINS: &[(&str, u32, &str, &str)] = &[
    ("P12345", 1, "MLPGLALLLLAAWTARALEV", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P54321", 2, "PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P67890", 6, "KWDSDPSGTKTCIDT", "GO:0009279;IPR:IPR016364;IPR:IPR008816"),
    ("P13579", 17, "KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH", "GO:0009279;IPR:IPR016364;IPR:IPR008816")
];

const TEXT_LENGTH: usize = 114; // sum of sequence lengths + 3 separators + 1 terminator

fn write_tsv(path: &std::path::Path) {
    let mut f = std::fs::File::create(path).unwrap();
    for (uid, taxon, seq, fa) in PROTEINS {
        writeln!(f, "{}\t{}\t{}\t{}", uid, taxon, seq, fa).unwrap();
    }
}

#[test]
fn test_build_creates_all_output_files() {
    let dir = tempfile::tempdir().unwrap();

    let db = dir.path().join("proteins.tsv");
    let out_sa = dir.path().join("sa.bin");
    let out_proteins = dir.path().join("proteins.bin");
    let out_mapping = dir.path().join("mapping.bin");

    write_tsv(&db);

    let status = Command::new(env!("CARGO_BIN_EXE_sa-builder"))
        .args([
            "--database-file",
            db.to_str().unwrap(),
            "--output-sa",
            out_sa.to_str().unwrap(),
            "--output-proteins",
            out_proteins.to_str().unwrap(),
            "--output-mapping",
            out_mapping.to_str().unwrap()
        ])
        .status()
        .expect("failed to run sa-builder");

    assert!(status.success(), "sa-builder exited with non-zero status");
    assert!(out_sa.exists(), "suffix array file not created");
    assert!(out_proteins.exists(), "proteins binary not created");
    assert!(out_mapping.exists(), "mapping binary not created");
}

#[test]
fn test_suffix_array_output() {
    let dir = tempfile::tempdir().unwrap();

    let db = dir.path().join("proteins.tsv");
    let out_sa = dir.path().join("sa.bin");
    let out_proteins = dir.path().join("proteins.bin");
    let out_mapping = dir.path().join("mapping.bin");

    write_tsv(&db);

    Command::new(env!("CARGO_BIN_EXE_sa-builder"))
        .args([
            "--database-file",
            db.to_str().unwrap(),
            "--output-sa",
            out_sa.to_str().unwrap(),
            "--output-proteins",
            out_proteins.to_str().unwrap(),
            "--output-mapping",
            out_mapping.to_str().unwrap()
        ])
        .status()
        .unwrap();

    let sa = SuffixArray::read_binary_mmap(&out_sa).unwrap();

    assert_eq!(sa.sample_rate(), 1, "expected sparseness factor 1");
    assert_eq!(sa.bits_per_value(), 64, "expected uncompressed (64 bits per value)");
    assert_eq!(sa.len(), TEXT_LENGTH, "SA length must equal text length");

    // Every value must be a valid text position.
    for i in 0..sa.len() {
        let v = sa.get(i);
        assert!(v >= 0 && v < TEXT_LENGTH as i64, "SA value {} out of range at index {}", v, i);
    }
}

#[test]
fn test_proteins_output() {
    let dir = tempfile::tempdir().unwrap();

    let db = dir.path().join("proteins.tsv");
    let out_sa = dir.path().join("sa.bin");
    let out_proteins = dir.path().join("proteins.bin");
    let out_mapping = dir.path().join("mapping.bin");

    write_tsv(&db);

    Command::new(env!("CARGO_BIN_EXE_sa-builder"))
        .args([
            "--database-file",
            db.to_str().unwrap(),
            "--output-sa",
            out_sa.to_str().unwrap(),
            "--output-proteins",
            out_proteins.to_str().unwrap(),
            "--output-mapping",
            out_mapping.to_str().unwrap()
        ])
        .status()
        .unwrap();

    let proteins = Proteins::read_binary_mmap(&out_proteins).unwrap();

    assert_eq!(proteins.len(), PROTEINS.len(), "protein count mismatch");

    for (i, (uid, taxon, _, _)) in PROTEINS.iter().enumerate() {
        assert_eq!(proteins.get(i).taxon_id, *taxon, "taxon mismatch for protein {}", i);
        assert_eq!(proteins.get(i).uniprot_id, *uid, "uniprot id mismatch for protein {}", i);
    }

    assert_eq!(proteins.text().len(), TEXT_LENGTH, "protein text length mismatch");
}

#[test]
fn test_mapping_output() {
    let dir = tempfile::tempdir().unwrap();

    let db = dir.path().join("proteins.tsv");
    let out_sa = dir.path().join("sa.bin");
    let out_proteins = dir.path().join("proteins.bin");
    let out_mapping = dir.path().join("mapping.bin");

    write_tsv(&db);

    Command::new(env!("CARGO_BIN_EXE_sa-builder"))
        .args([
            "--database-file",
            db.to_str().unwrap(),
            "--output-sa",
            out_sa.to_str().unwrap(),
            "--output-proteins",
            out_proteins.to_str().unwrap(),
            "--output-mapping",
            out_mapping.to_str().unwrap()
        ])
        .status()
        .unwrap();

    let mapping = SuffixToProteinMapping::read_binary_mmap(&out_mapping).unwrap();
    let idx = &*mapping.0;

    // Positions inside each protein resolve to the correct protein index.
    assert_eq!(idx.suffix_to_protein(0), 0, "position 0 should map to protein 0");
    assert_eq!(idx.suffix_to_protein(21), 1, "position 21 should map to protein 1");
    assert_eq!(idx.suffix_to_protein(52), 2, "position 52 should map to protein 2");
    assert_eq!(idx.suffix_to_protein(68), 3, "position 68 should map to protein 3");

    // Separator and terminator positions resolve to NULL.
    assert!(idx.suffix_to_protein(20).is_null(), "position 20 (separator) should be NULL");
    assert!(idx.suffix_to_protein(51).is_null(), "position 51 (separator) should be NULL");
    assert!(idx.suffix_to_protein(67).is_null(), "position 67 (separator) should be NULL");
    assert!(idx.suffix_to_protein(113).is_null(), "position 113 (terminator) should be NULL");
}
