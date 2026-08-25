# Unipept Index

![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&logo=Codecov)

The Unipept index, written entirely in `Rust`: given a peptide, find every protein that contains
it. This repository is a Cargo workspace of several crates that depend on each other.

## The one thing to know first: storage is a build-time choice

Every storage structure in the index has two implementations — one holding **owned memory**, one
borrowing a **memory mapping**. Both are always compiled, and the library crates are generic over
which they are handed; **the choice is made by `sa-server`, at compile time, through features on
that crate alone** (`sa-server/src/backends.rs`). There is no runtime switch, and no dispatch
anywhere in the search path.

**`sa-server` declares no `default = [...]`, so a plain `cargo build --release -p sa-server` gives
you the preloaded backend.** Build with `--features mmap` for the memory-mapped one.

| | preloaded (default) | `--features mmap` |
|---|---|---|
| Where the index lives | resident in the process | a memory mapping; the kernel decides what stays resident |
| Memory | pays the full index size in RSS | bounded, at the cost of page faults |
| Startup | slow: everything is read up front | fast, but the first queries pay the faults unless warmed |
| Use when | the index comfortably fits in RAM | the index is large relative to available memory |

### Mixing the two

The choice is per structure, not per build. `mmap` maps all four; each `preloaded-*` feature then
pulls **one** of them back into owned memory:

| feature | structure it un-maps |
|---|---|
| `preloaded-text` | the concatenated protein text (~190 MB at UniProt scale) |
| `preloaded-proteins` | the protein metadata table (accessions + annotations) |
| `preloaded-mapping` | the suffix-to-protein mapping |

The suffix array always follows `mmap`; there is no `preloaded-sa`. A `preloaded-*` feature has no
effect without `mmap`, since everything is preloaded already.

This exists because the structures are not alike: the text is read once per character compared
while the metadata table is the largest thing in the index and read once per reported result. So
`--features mmap,preloaded-text` keeps the index mapped but the hottest structure resident — worth
measuring on your own data before choosing a default.

**If you preload anything, preload the text as well as the metadata.**
`--features mmap,preloaded-text,preloaded-proteins` is the fastest configuration measured that
still leaves the 160 GB suffix array mapped: on the full UniProt index it beats
`mmap,preloaded-proteins` alone in all 16 default cells, and ties the fully preloaded build.
Preloading only the metadata is the configuration to avoid — it carries the highest resident
footprint of any arm (250 GB against 242 GB fully preloaded) and still leaves search mapped, which
is where the remaining gap lives.

**Any of them is a bet on the index fitting.** Preloaded memory is non-evictable, so under pressure
it displaces the page cache the mapped structures live in rather than being reclaimed. Measured
against a falling ceiling, plain `mmap` is ahead of `mmap,preloaded-proteins` at every ceiling that
binds, and at 78 GB the preloaded arms lose almost all their throughput to the fault rate. Preload
when residency is guaranteed; do not where the ceiling might move. The numbers are in the "When the
index does not fit in RAM" section of `sa-index/src/lib.rs`.

These nine feature combinations are the ones the server exposes; the types themselves compose into
sixteen, and `sa-index`'s `backend_agreement` test builds every one of them and asserts they return
identical results.

A running server reports what it was compiled with:

```
Storage backends: sa=mmap text=preloaded proteins=mmap mapping=mmap
```

## Running an mmap build when the index does not fit in RAM

Two settings matter far more than the storage flags above, and neither is on by default. Measured
on the full 223 GB index with cgroup ceilings; see the `sa-index` crate docs for the method.

**1. Raise the thread count.** A page fault blocks the thread that takes it, so with rayon at the
core count every faulting thread idles a core. Raising it does not reduce faults — it overlaps
them:

| RAM available | default threads | tuned | gain |
|---|---|---|---|
| more than the index | 35,710 qps | 35,046 (48 threads) | **−1.9%** |
| 75% of the index | 15,739 qps | 26,071 (48 threads) | **+65.6%** |
| 50% of the index | 10,561 qps | 19,654 (96 threads) | **+86.1%** |

```bash
RAYON_NUM_THREADS=96 ./target/release/sa-server ...
```

The gain scales with how much the index overflows RAM, and it *costs* up to 10% when everything
fits — so set it for constrained deployments and leave it alone otherwise.

**2. Build a 6-mer k-mer table, not a 5-mer.** Under pressure the 6-mer is +18.4% and takes 27.9%
fewer major faults than no table; a 5-mer is +3.2%, barely better than nothing. It costs 3.06 GB
against 127 MB, which is why the resident-case measurement rejected it.

```bash
./target/release/sa-builder ... --output-kmer-table kmer-tables/6mer_table.bin --kmer-size 6
```

With both applied, a box holding 75% of the index runs at ~74% of its unconstrained throughput
instead of ~36%.

## Installation

> [!NOTE]
> To build and use the Unipept Index, you need to have Rust installed. If you don't have Rust
> installed, you can get it from [rust-lang.org](https://www.rust-lang.org/).

```bash
git clone https://github.com/unipept/unipept-index.git
cd unipept-index

cargo build --release -p sa-server                    # preloaded backend
cargo build --release -p sa-server --features mmap    # memory-mapped backend

# ... or mix: map everything except the protein text
cargo build --release -p sa-server --features mmap,preloaded-text
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
  --output-kmer-table kmer_table.bin        # --kmer-size defaults to 6
```

This precomputes the suffix-array range for every k-mer, so a search starts from those bounds
instead of binary-searching the whole array. It is an accelerator only — results are identical
with and without it — and it helps most on short peptides.

Note it is a *dense* table: it has one entry per possible k-mer, so its size depends only on `k`
and not on the database. That makes `k` the whole memory decision — roughly 127 MB at
`--kmer-size 5` against 2.85 GB at the default `6`. The default assumes the mapped deployment,
where the larger table pays for itself in avoided page faults; see "Choosing a storage backend"
above. Pass `--kmer-size 5` if the index is guaranteed resident. The table can only be built during
a full index build, since it is derived from the suffix array before that is written out.

## Running the server

`sa-server` exists to exercise the index over HTTP — for testing, and for serving an index directly
when that is what you want. **It is not how the search reaches Unipept:** the API site builds its
own endpoints on top of `sa-index`'s public functions (`peptide_search::search_all_peptides_json`
and friends) rather than proxying this binary.

That matters for anyone adding limits or validation. `sa-server` is not the untrusted boundary, so
hardening it protects nothing that runs in production; the caller that exposes these functions to
the network is the place for that. See the note on resource limits in `sa-index::peptide_search`.

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

Both backends of every structure are always compiled, so one run covers them all — the searcher's
fixtures build each combination by naming its types, and `backend_agreement` asserts the sixteen
give identical answers:

```bash
cargo test                                  # everything, both backends
cargo test -p sa-index --features measure   # the one remaining feature
```

CI additionally checks that all nine of `sa-server`'s feature combinations typecheck. Lints:

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
| `sa-server` | HTTP server for testing and for serving an index directly; see the note below |
| `sa-builder` | builds an index from a protein TSV |
| `sa-index` | the search itself: suffix array, k-mer table, suffix→protein mapping |
| `sa-mappings` | protein metadata — accessions, taxa, annotations |
| `text-compression` | the concatenated protein text, packed at 5 bits per residue |
| `bitarray` | dense arrays of fixed-width values |
| `fa-compression` | encoding for functional annotations |
| `binary-traits` | the read/write/load traits every on-disk structure is written and read through |
| `prefetch` | one portable software-prefetch hint |
| `libsais64-rs` | bindings to the suffix-array construction library |
| `sa-benchmarks` | measurement harness (see below) |

`sa-benchmarks` is **excluded from `default-members`**, so ordinary builds skip it and its extra
dependencies. Build it explicitly:

```bash
cargo build --release -p sa-benchmarks
```

Beware that `cargo build --workspace` *overrides* `default-members` and will include it again.

## Benchmarking

Everything needed to measure the index — the harness, the driver, the suite definitions and the
machine profiles — lives in `sa-benchmarks/`. See [its README](sa-benchmarks/README.md).

```bash
cp sa-benchmarks/profiles/example.toml sa-benchmarks/profiles/local.toml   # once per machine
./sa-benchmarks/run.sh defaults      # the regression gate — run after any change to the search path
sudo ./sa-benchmarks/run.sh all      # every suite in one session, into one report.md
```

Five suites: `defaults` (throughput at production defaults), `detail` (where the time goes inside a
search), `startup` (what each storage configuration costs before the first query), `ram` (scaling as
the RAM ceiling falls) and `threads` (whether oversubscription pays). `ram` and `threads` need root
for cgroup ceilings and `drop_caches`; `run.sh all` skips them without it and says so in the report.

Every run writes a self-contained `report.html` — sidebar, row filter, sortable columns — next to
its results, plus `report.md` and, for `all`, a `report.json` a later run can use as its baseline.

Measurement code stays out of what ships. The hot-path instrumentation is behind `sa-index`'s
`measure` feature, which compiles to nothing when off and which CI proves `sa-server` and
`sa-builder` never enable; run-level measurement lives in `sa-benchmarks`, which never ships at all.

## Further reading

* Crate-level `cargo doc` for `sa-index` explains the search pipeline and the feature axes in
  detail, including why the two backends exist and which decisions were measured and rejected.
