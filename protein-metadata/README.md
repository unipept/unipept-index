# protein-metadata

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

Protein metadata — accessions, taxon ids and functional annotations — addressed by index. A suffix
array search returns positions in the concatenated protein text; the suffix-to-protein mapping
turns those positions into protein indices; this crate is the last step, turning an index into the
metadata a search result actually reports.

It also owns the loader for the text itself. `InMemoryProteins::load_from_tsv` reads a UniProt TSV
(`uniprot_id`, `taxon_id`, `sequence`, `annotations`), upper-cases every sequence and joins them
with `-` between and `$` at the end, which is the layout the suffix array is built over. The packed
text lives in [`protein-text`](../protein-text/README.md); the two share one `proteins.bin`.

## Two backends

Both are always compiled. A caller picks one by naming its type:

| type | where the bytes live |
|---|---|
| `InMemoryProteins` | accessions and encoded annotations in owned memory. Always available, because it also owns the `WriteBinary` implementation `sa-builder` uses to produce `proteins.bin` for *both* backends |
| `MmapBackedProteins` | the same fields decoded straight out of a mapping of that file, which is what keeps a multi-gigabyte protein table servable in bounded RSS |

Both implement `ProteinsBackend` and both hand out a borrowed `ProteinRef`, so reading code is
written once and never copies per result. The on-disk format of `proteins.bin` is documented at its
writer, on the `impl WriteBinary for InMemoryProteins` block in `src/preloaded.rs`.

## Two axes, not one

`proteins.bin` holds two things: the concatenated protein text and the metadata table. The text is
the hottest structure in the index and the metadata table the one that grows most when preloaded —
roughly tripling, since it becomes owned strings rather than bytes in a file — so the best storage
for one is not the best storage for the other — and both structs above are generic over their text backend so the two can be chosen
independently.

| pairing | metadata | text |
|---|---|---|
| `InMemoryProteins<InMemoryProteinText>` | owned | owned |
| `MmapBackedProteins<MmapBackedProteinText>` | mapped | mapped |
| `MmapBackedProteins<InMemoryProteinText>` | mapped | owned |
| `InMemoryProteins<MmapBackedProteinText>` | owned | mapped |

All four load from the same file, and every reader that needs the mapping lives in `src/mmap.rs`,
sharing one header parser so the combinations cannot drift apart on how they bound it. The third
row is the interesting one: it keeps the multi-gigabyte metadata table mapped while the text that
search reads once per character compared sits in owned RAM.

## Example

```rust
use protein_metadata::{InMemoryProteins, ProteinsBackend};

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

`get` is unchecked in both backends and neither guarantees a panic on an out-of-range index —
`InMemoryProteins` indexes a `Vec` and panics, `MmapBackedProteins` may decode an entry-sized window
of whatever follows the table and return fabricated metadata. Callers bound the index themselves.

## Where it sits

Depends on `binary-traits`, `protein-text`, `fa-compression` and `memory-hints`. Used by `sa-index`,
`sa-builder` and `sa-server`.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p protein-metadata --open`
