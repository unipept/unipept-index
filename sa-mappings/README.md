# Suffix Array Mappings

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=sa-mappings&logo=codecov)
![Static Badge](https://img.shields.io/badge/doc-rustdoc-blue)

A suffix array search returns positions in the concatenated protein text, and the suffix-to-protein
mapping turns those positions into protein indices. The `sa-mappings` library is the last step: it
turns an index into the metadata a search result actually reports — the UniProt accession, the NCBI
taxon id, and the functional annotations.

It also owns the concatenated text itself. `InMemoryProteins::load_from_tsv` reads a UniProt TSV
(`uniprot_id`, `taxon_id`, `sequence`, `annotations`), upper-cases every sequence and joins them
with `-` between and `$` at the end, which is the layout the suffix array is built over.

## Two backends

Both are always compiled. A caller picks a backend by naming its type:

* `InMemoryProteins` — accessions and encoded annotations held in owned memory. Always available,
  because it also owns the `WriteBinary` implementation that `sa-builder` uses to produce
  `proteins.bin` for *both* backends to read.
* `MmapBackedProteins` — the same fields decoded straight out of a memory mapping of that file,
  which is what keeps a multi-gigabyte protein table servable in bounded RSS.

Both implement `proteins::ProteinsBackend` and both hand out a borrowed `ProteinRef`, so reading
code is written once and never copies per result. The on-disk format of `proteins.bin` is
documented at its writer, on the `impl WriteBinary for InMemoryProteins` block in
`src/proteins/preloaded.rs`.

## Two axes, not one

`proteins.bin` holds two things: the concatenated protein text and the metadata table. They are
the hottest and the biggest structures in the index respectively, so the best storage for one is
not the best storage for the other — and both structs above are generic over their text backend so
the two can be chosen independently.

| pairing | metadata | text |
|---|---|---|
| `InMemoryProteins<InMemoryProteinText>` | owned | owned |
| `MmapBackedProteins<MmapBackedProteinText>` | mapped | mapped |
| `MmapBackedProteins<InMemoryProteinText>` | mapped | owned |
| `InMemoryProteins<MmapBackedProteinText>` | owned | mapped |

All four load from the same file, and every reader that needs the mapping lives in
`src/proteins/mmap.rs`, sharing one header parser so the combinations cannot drift apart on how
they bound it. The third row is the interesting one: it keeps the multi-gigabyte metadata table
mapped while the ~190 MB text that search reads once per character compared sits in owned RAM.

## Example

```rust
use sa_mappings::proteins::{InMemoryProteins, ProteinsBackend};

fn main() {
    let proteins = InMemoryProteins::load_from_tsv("database.tsv").unwrap();

    let protein = proteins.get(0);

    // "P12345"
    println!("{}", protein.uniprot_id);

    // 1
    println!("{}", protein.taxon_id);

    // Annotations are stored encoded and decoded on demand, once per reported result:
    // "GO:0009279;IPR:IPR016364;IPR:IPR008816"
    println!("{}", protein.get_functional_annotations());
}
```
