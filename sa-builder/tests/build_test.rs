//! End-to-end tests: run the real `sa-builder` binary and inspect what it writes.
//!
//! The builder writes with the owned types in every configuration, and both readers of every
//! structure are always compiled, so each assertion below reads the output back through *both* —
//! which is the property that matters for a builder: what it writes has to be readable either way.
//! These used to be gated off behind `--features mmap` and only ever checked the mapped reader.

use std::{io::Write, path::Path, process::Command};

use binary_traits::LoadIndex;
use protein_metadata::{InMemoryProteins, MmapBackedProteins, ProteinsBackend};
use sa_index::{
    Nullable,
    array::{InMemorySA, MmapBackedSA, SuffixArrayBackend},
    suffix_to_protein_index::{
        InMemorySuffixToProteinMapping, MmapBackedSuffixToProteinMapping, SuffixToProteinMappingBackend
    }
};
use tempfile::TempDir;
use text_compression::{InMemoryProteinText, MmapBackedProteinText, ProteinTextBackend as _};

/// Four proteins used as test input, matching the fixture in protein-metadata unit tests.
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

/// The three files one builder run produces. The [`TempDir`] is returned with them because
/// dropping it deletes all three.
struct BuiltIndex {
    _dir: TempDir,
    sa: std::path::PathBuf,
    proteins: std::path::PathBuf,
    mapping: std::path::PathBuf
}

/// Writes a TSV of [`PROTEINS`] and runs the real builder binary over it.
fn build_index() -> BuiltIndex {
    let dir = tempfile::tempdir().unwrap();

    let db = dir.path().join("proteins.tsv");
    let sa = dir.path().join("sa.bin");
    let proteins = dir.path().join("proteins.bin");
    let mapping = dir.path().join("mapping.bin");

    let mut f = std::fs::File::create(&db).unwrap();
    for (uid, taxon, seq, fa) in PROTEINS {
        writeln!(f, "{}\t{}\t{}\t{}", uid, taxon, seq, fa).unwrap();
    }
    drop(f);

    let status = Command::new(env!("CARGO_BIN_EXE_sa-builder"))
        .args([
            "--database-file",
            db.to_str().unwrap(),
            "--output-sa",
            sa.to_str().unwrap(),
            "--output-proteins",
            proteins.to_str().unwrap(),
            "--output-mapping",
            mapping.to_str().unwrap()
        ])
        .status()
        .expect("failed to run sa-builder");
    assert!(status.success(), "sa-builder exited with non-zero status");

    BuiltIndex { _dir: dir, sa, proteins, mapping }
}

#[test]
fn test_build_creates_all_output_files() {
    let index = build_index();

    assert!(index.sa.exists(), "suffix array file not created");
    assert!(index.proteins.exists(), "proteins binary not created");
    assert!(index.mapping.exists(), "mapping binary not created");
}

#[test]
fn test_suffix_array_output() {
    let index = build_index();

    fn assert_sa(sa: &impl SuffixArrayBackend) {
        assert_eq!(sa.sample_rate(), 1, "expected sparseness factor 1");
        assert_eq!(sa.bits_per_value(), 64, "expected uncompressed (64 bits per value)");
        assert_eq!(sa.len(), TEXT_LENGTH, "SA length must equal text length");

        // Every value must be a valid text position.
        for i in 0..sa.len() {
            let v = sa.get(i);
            assert!(v >= 0 && v < TEXT_LENGTH as i64, "SA value {} out of range at index {}", v, i);
        }
    }

    assert_sa(&InMemorySA::load(&index.sa).unwrap());
    assert_sa(&MmapBackedSA::load(&index.sa).unwrap());
}

#[test]
fn test_proteins_output() {
    let index = build_index();

    fn assert_proteins(proteins: &impl ProteinsBackend) {
        assert_eq!(proteins.len(), PROTEINS.len(), "protein count mismatch");

        for (i, (uid, taxon, _, _)) in PROTEINS.iter().enumerate() {
            assert_eq!(proteins.get(i).taxon_id, *taxon, "taxon mismatch for protein {}", i);
            assert_eq!(proteins.get(i).uniprot_id, *uid, "uniprot id mismatch for protein {}", i);
        }

        assert_eq!(proteins.text().len(), TEXT_LENGTH, "protein text length mismatch");
    }

    // All four pairings of the two independent axes — `proteins.bin` holds the metadata and the
    // text, and each may be owned or mapped.
    fn check<P: ProteinsBackend + LoadIndex>(path: &Path) {
        assert_proteins(&P::load(path).unwrap());
    }
    check::<InMemoryProteins<InMemoryProteinText>>(&index.proteins);
    check::<InMemoryProteins<MmapBackedProteinText>>(&index.proteins);
    check::<MmapBackedProteins<InMemoryProteinText>>(&index.proteins);
    check::<MmapBackedProteins<MmapBackedProteinText>>(&index.proteins);
}

#[test]
fn test_mapping_output() {
    let index = build_index();

    // The mapping enum implements the backend trait directly — no inner box to unwrap.
    fn assert_mapping(idx: &impl SuffixToProteinMappingBackend) {
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

    assert_mapping(&InMemorySuffixToProteinMapping::load(&index.mapping).unwrap());
    assert_mapping(&MmapBackedSuffixToProteinMapping::load(&index.mapping).unwrap());
}
