# Unipept Index

![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&logo=Codecov)

The Unipept index, written entirely in `Rust`: given a peptide, find every protein that contains
it. This repository is a Cargo workspace of several crates that depend on each other.

## The one thing to know first: there are two builds

Every storage structure in the index has two implementations — one holding **owned memory**, one
borrowing a **memory mapping** — and which you get is a *compile-time* choice made by the `mmap`
feature. There is no runtime switch.

**No crate declares `default = [...]`, so a plain `cargo build --release` gives you the preloaded
backend.** The production server is built with `--features mmap`.

| | preloaded (default) | `--features mmap` |
|---|---|---|
| Where the index lives | resident in the process | a memory mapping; the kernel decides what stays resident |
| Memory | pays the full index size in RSS | bounded, at the cost of page faults |
| Startup | slow: everything is read up front | fast, but the first queries pay the faults unless warmed |
| Use when | the index comfortably fits in RAM | the index is large relative to available memory |

A running server reports which one it is compiled as:

```
Storage backend: mmap
```

## Installation

> [!NOTE]
> To build and use the Unipept Index, you need to have Rust installed. If you don't have Rust
> installed, you can get it from [rust-lang.org](https://www.rust-lang.org/).

```bash
git clone https://github.com/unipept/unipept-index.git
cd unipept-index

cargo build --release                    # preloaded backend
cargo build --release --features mmap    # memory-mapped backend
```

## Building an index

`sa-builder` turns a protein TSV into the three files the server needs. Rows are
`uniprot_id<TAB>taxon_id<TAB>sequence<TAB>annotations`.

```bash
./target/release/sa-builder \
  --database-file    proteins.tsv \
  --output-sa        sa.bin \
  --output-proteins  proteins.bin \
  --output-mapping   mapping.bin \
  --sparseness-factor 3 \
  --compress-sa
```

All three outputs are mandatory and all three must come from the same run — the suffix array
indexes positions in the text stored in `proteins.bin`, and `mapping.bin` maps those same
positions to entries in it. Mixing files from different builds produces wrong answers rather than
errors.

Two options worth understanding:

* `--sparseness-factor n` indexes only every n-th text position, shrinking the suffix array at the
  cost of search work. **Peptides shorter than `n` cannot be searched at all.**
* `--compress-sa` packs entries at the minimum width the text length needs rather than 64 bits,
  roughly halving the file.

### Optional: a k-mer bounds table

```bash
  --output-kmer-table kmer_table.bin --kmer-size 5
```

This precomputes the suffix-array range for every k-mer, so a search starts from those bounds
instead of binary-searching the whole array. It is an accelerator only — results are identical
with and without it — and it helps most on short peptides.

Note it is a *dense* table: at `--kmer-size 5` it is roughly 127 MB regardless of database size,
because it has one entry per possible k-mer. It can only be built during a full index build,
since it is derived from the suffix array before that is written out.

## Running the server

```bash
./target/release/sa-server \
  --database-file proteins.bin \
  --index-file    sa.bin \
  --mapping-file  mapping.bin \
  --kmer-table-file kmer_table.bin \   # optional
  --address 0.0.0.0:3000               # optional, this is the default
```

Then POST to `/search`:

```bash
curl -X POST http://localhost:3000/search \
  -H 'Content-Type: application/json' \
  -d '{"peptides": ["MLPGLALLLL"], "equate_il": false, "tryptic": false, "cutoff": 10000}'
```

| field | default | meaning |
|---|---|---|
| `peptides` | — | the peptides to search |
| `equate_il` | `false` | treat I and L as interchangeable |
| `tryptic` | `false` | return only matches at tryptic boundaries |
| `cutoff` | `10000` | stop after this many matches; the response flags `cutoff_used` |

## Testing

Because the backend is a compile-time choice, so is what gets tested:

```bash
cargo test                  # preloaded backend
cargo test --all-features   # memory-mapped backend, plus the metrics instrumentation
```

Both are run by CI. Lints:

```bash
cargo clippy --all-targets --all-features -- -D warnings
cargo +nightly fmt --all --check    # nightly is required, see below
```

`cargo fmt` **must** run on nightly: `.rustfmt.toml` enables `unstable_features` along with
`imports_granularity`, `group_imports` and `normalize_comments`, all of which stable rustfmt
silently ignores.

## The crates

| crate | what it does |
|---|---|
| `sa-server` | HTTP server; the deployed artifact |
| `sa-builder` | builds an index from a protein TSV |
| `sa-index` | the search itself: suffix array, k-mer table, suffix→protein mapping |
| `sa-mappings` | protein metadata — accessions, taxa, annotations |
| `text-compression` | the concatenated protein text, packed at 5 bits per residue |
| `bitarray` | dense arrays of fixed-width values |
| `fa-compression` | encoding for functional annotations |
| `binary-traits` | the read/write traits every on-disk structure implements |
| `prefetch` | one portable software-prefetch hint |
| `libsais64-rs` | bindings to the suffix-array construction library |
| `sa-benchmarks` | measurement harness (see below) |

`sa-benchmarks` is **excluded from `default-members`**, so ordinary builds skip it and its extra
dependencies. Build it explicitly:

```bash
cargo build --release -p sa-benchmarks
```

Beware that `cargo build --workspace` *overrides* `default-members` and will include it again.

## Further reading

* [`docs/design/`](docs/design/) — why the index is shaped the way it is: the optimization levers
  that were measured, and the ones that were measured and rejected. Read
  `mlp-batching-design.md` first.
* [`scripts/`](scripts/) — dev tooling. `run_index.sh` boots a server and replays a peptide file;
  `baseline.sh` sweeps the full configuration matrix and diffs two captures, which is how a
  refactor is shown not to have changed any answer.
* Crate-level `cargo doc` for `sa-index` explains the search pipeline and the feature axes in
  detail.
