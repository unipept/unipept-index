# sa-builder

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

Builds the on-disk index that [`sa-server`](../sa-server/README.md) reads: the suffix array, the
protein store, the suffix-to-protein mapping and an optional k-mer bounds table.

**This is the only writer in the workspace.** Every format the other crates read is produced here,
through the `WriteBinary` implementations documented in
[`binary-traits`](../binary-traits/README.md) — so a format question is answered at its writer, not
here. It names only the preloaded types, since there is one writer per structure whichever backend
reads it later.

## Usage

```bash
cargo build --release -p sa-builder

./target/release/sa-builder \
  --database-file    proteins.tsv \
  --output-sa        sa.bin \
  --output-proteins  proteins.bin \
  --output-mapping   mapping.bin \
  --sparseness-factor 3 \
  --compress-sa
```

Input rows are `uniprot_id<TAB>taxon_id<TAB>sequence<TAB>annotations`. Sequences are concatenated
into a single text, separated by `-` and terminated by `$`; neither character may appear in a
sequence.

| option | default | what it does |
|---|---|---|
| `--sparseness-factor n` | 1 | index only every n-th text position. Shrinks the array at the cost of search work — **peptides shorter than n cannot be searched at all** |
| `--compress-sa` | off | pack entries at the minimum width the text length needs rather than 64 bits, roughly halving the file |
| `--construction-algorithm` | `lib-sais` | `lib-sais` samples during construction; `lib-div-suf-sort` always builds the dense array first, so peak memory is that of the dense array whatever the sparseness |
| `--mapping-style` | `bit-vec` | `bit-vec`, `sparse` or `dense`. The style is recorded in the file, so the server picks its reader from what it finds |
| `--output-kmer-table` | none | also build a k-mer bounds table |
| `--kmer-size` | 5 | k for that table, maximum 7 |

## One run produces one index

The four files are only usable as a set: the suffix array indexes positions in the text stored in
`proteins.bin`, and the mapping resolves those same positions to entries in it. Mixing files from
different builds yields **wrong answers rather than errors** — `Searcher::try_new` is what refuses
to start on a mismatched set, and only once the files meet.

The binary therefore writes every section to a temporary sibling and renames them all only once the
last one has succeeded, so an interrupted build cannot leave a half-set behind.

## The k-mer table

Dense at `24^k` entries of 16 bytes, so its size depends only on `k` and never on the database:
k=5 is 127 MB and k=6 is 3.06 GB, a 24x step for one more level of the probe chain. That makes `k`
the whole memory decision.

5 is the default because it is the size that pays in both regimes. With the index fully resident,
the 6-mer's edge over it sits inside the noise floor on most length regimes and does not buy the
extra 2.9 GB. Raise it to 6 only under a memory ceiling, where the table's value is working-set
size rather than probe count — a 5-mer narrows the search to ~7 suffix-array pages per query and a
6-mer to ~1, which is +18.4% against no table where the 5-mer manages +3.2%.

The table can only be built during a full index build, since it is derived from the suffix array
before that is written out.

## Where it sits

Depends on `binary-traits`, `libsais64-rs`, `libdivsufsort-rs`, `protein-text`, `protein-metadata`
and `sa-index`. Nothing depends on it. Building it needs `libsais64-rs`' toolchain — see the
[workspace README](../README.md).

The pieces worth testing on their own live in `src/lib.rs` — the command line in `Arguments` and
suffix-array construction in `build_ssa`. File writing lives in the binary.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p sa-builder --open`
