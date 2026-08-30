# Unipept Index

![Test](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github&label=test)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&logo=Codecov)

The Unipept index, written entirely in `Rust`: given a peptide, find every protein that contains
it. This repository is a Cargo workspace of several crates that depend on each other.

Every performance figure below comes from one session on the full UniProt index (223 GB, 12-core
Xeon Silver 4410Y, 295 GB RAM), run at `660befd7ee`. The harness that produced it is in
[`sa-benchmarks/`](sa-benchmarks/README.md).

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
| Memory | pays the full index size in RSS, none of it evictable | bounded, at the cost of page faults |
| Startup | slow: 466 s to read the 223 GB index up front | 3 s to map it, then 316 s of optional warmup |
| Use when | the index comfortably fits in RAM | the index is large relative to available memory |

### Mixing the two

The choice is per structure, not per build. `mmap` maps all four; each `preloaded-*` feature then
pulls **one** of them back into owned memory:

| feature | structure it un-maps |
|---|---|
| `preloaded-text` | the concatenated protein text (~43 GB on the full UniProt index) |
| `preloaded-proteins` | the protein metadata table — accessions and annotations, ~8 GB on disk and ~24 GB once resident |
| `preloaded-mapping` | the suffix-to-protein mapping (~11 GB) |

The suffix array always follows `mmap`; there is no `preloaded-sa`. A `preloaded-*` feature has no
effect without `mmap`, since everything is preloaded already.

This exists because the structures are not alike: the text is read once per character compared,
while the metadata table is read once per reported result and triples in size when it moves to the
heap. So `--features mmap,preloaded-text` keeps the index mapped but the hottest structure
resident — worth measuring on your own data before choosing a default.

**If you preload anything, preload the text as well as the metadata.**
`--features mmap,preloaded-text,preloaded-proteins`, and the same plus `preloaded-mapping`, are the
two fastest configurations measured. They cannot be told apart from each other, and together they
are ahead of every other arm by more than the noise floor in 6 of the 16 throughput cells — by 12%
at the median. Both beat `mmap,preloaded-proteins` alone in **all 16** cells, and both match the
fully preloaded build while leaving the 160 GB suffix array mapped, pinning roughly 78 GB of
non-evictable memory where the fully preloaded build pins 242 GB.

Preloading only the metadata is the configuration to avoid. It closes the retrieval gap and none of
the search one, which is the larger half, and it carries the highest resident footprint of any arm
(250 GB, against 242 GB for the fully preloaded build).

**Any of them is a bet on the index fitting.** Preloaded memory is non-evictable, so under pressure
it displaces the page cache the mapped structures live in rather than being reclaimed. Swept
against a falling cgroup ceiling, plain `mmap` is level with or ahead of every preloading arm at
every ceiling that binds, and at 78 GB — roughly a third of the index — the preloaded arms do not
degrade, they collapse:

| ceiling | `mmap` | `mmap,preloaded-proteins` | `mmap,preloaded-text,preloaded-proteins` |
|---|---|---|---|
| none | 35,725 qps | 43,603 | 52,198 |
| 167 GB (75%) | 14,748 | 13,832 | 15,161 |
| 112 GB (50%) | 10,679 | 9,582 | 10,212 |
| 78 GB (35%) | **7,411** | 292 | 169 |

The unconstrained `mmap` figure here (35,725 qps) and the one in the thread-count table below
(38,081) are the same build measured by two different suites in the same session. The 6.6% between
them is inside both cells' resolution floors — compare within a table, not across.

At 78 GB the preloaded-text arm takes 55x `mmap`'s major faults, and adding `preloaded-mapping` on
top does not fit at all. Preload when residency is guaranteed; do not where the ceiling might move.
The method, and the alternatives that were measured and rejected, are in the "When the index does
not fit in RAM" section of `sa-index/src/lib.rs`.

These nine feature combinations are the ones the server exposes; the types themselves compose into
sixteen, and `sa-index`'s `every_backend_combination_returns_identical_results` test builds every
one of them and asserts they return identical results.

A running server reports what it was compiled with:

```
Storage backends: sa=mmap text=preloaded proteins=mmap mapping=mmap
```

## Running an mmap build when the index does not fit in RAM

Two settings matter far more than the storage flags above, and neither is on by default.

**1. Raise the thread count.** A page fault blocks the thread that takes it, so with rayon at the
core count every faulting thread idles a core. Raising it does not reduce faults — it overlaps
them. Measured on the `mmap` arm, with a 6-mer table attached:

| RAM available | default threads | tuned | change |
|---|---|---|---|
| more than the index | 38,081 qps | 34,264 (48 threads) | **−10.0%** |
| 75% of the index (167 GB) | 15,518 qps | 25,734 (96 threads) | **+65.8%** |
| 50% of the index (112 GB) | 10,503 qps | 20,473 (96 threads) | **+94.9%** |

```bash
RAYON_NUM_THREADS=96 ./target/release/sa-server ...
```

The gain scales with how much the index overflows RAM, and it *costs* up to 12% when everything
fits — so set it for constrained deployments and leave it alone otherwise. The preloading arms gain
more still under a ceiling (`mmap,preloaded-text,preloaded-proteins` is +102.7% at 167 GB), but
they are also the arms that collapse first when the ceiling keeps falling.

**2. Raise the k-mer table to 6, from the default 5.** Under a ceiling the 6-mer is +18.4% and
takes 27.9% fewer major faults than no table; a 5-mer is +3.2%, barely better than nothing. The
difference is working-set size rather than probe count — a 6-mer narrows the search to about one
suffix-array page per query where a 5-mer leaves seven. It costs 3.06 GB against 127 MB, which is
why the resident-case measurement rejected it and why 5 is the default; this is the one regime
where the 24x is worth paying.

```bash
./target/release/sa-builder ... --output-kmer-table kmer-tables/6mer_table.bin --kmer-size 6
```

With both applied, a box holding 75% of the index runs at ~68% of its unconstrained throughput
instead of ~41%.

## Installation

> [!NOTE]
> Building needs Rust — get it from [rust-lang.org](https://www.rust-lang.org/) — plus `git`,
> `cmake`, `make`, a C compiler and `libclang`. `libsais64-rs` builds the suffix-array construction
> library from source: its build script clones
> [`unipept/libsais-packed`](https://github.com/unipept/libsais-packed), compiles it with CMake and
> generates bindings with `bindgen`, so the first build needs network access.

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
errors; `Searcher::try_new` is what refuses to start on a mismatched set.

Two options worth understanding:

* `--sparseness-factor n` indexes only every n-th text position, shrinking the suffix array at the
  cost of search work. **Peptides shorter than `n` cannot be searched at all.**
* `--compress-sa` packs entries at the minimum width the text length needs rather than 64 bits,
  roughly halving the file.

### Optional: a k-mer bounds table

```bash
  --output-kmer-table kmer_table.bin        # --kmer-size defaults to 5
```

This precomputes the suffix-array range for every k-mer, so a search starts from those bounds
instead of binary-searching the whole array. It is an accelerator only — results are identical
with and without it.

It helps most on **long** peptides, which is the opposite of what the mechanism suggests: attaching
a table is worth +14% to +29% on 35-50-residue queries and nothing that clears the noise floor on
5-9-residue ones. A short peptide matches an enormous suffix-array range, so its cost is the scan
over that range rather than the binary-search descent the table removes.

Note it is a *dense* table: it has one entry per possible k-mer, so its size depends only on `k`
and not on the database. That makes `k` the whole memory decision — roughly 127 MB at the default
`--kmer-size 5` against 3.06 GB at `6`, and `7` is the maximum the builder accepts. The default is
5 because the 6-mer's extra 2.9 GB does not pay for itself with the index resident: nothing
separates the two sizes above the noise floor in any resident cell. Pass `--kmer-size 6` only under
a memory ceiling, where it does; see "Running an mmap build" above. The table can only be built
during a full index build, since it is derived from the suffix array before that is written out.

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

The request body is capped at 5 MB, which is the only limit the server imposes.

## Testing

Both backends of every structure are always compiled, so one run covers them all — the searcher's
fixtures build each combination by naming its types, and
`sa_index::sa_searcher::tests::every_backend_combination_returns_identical_results` asserts the
sixteen give identical answers:

```bash
cargo test   # everything, both backends
```

CI additionally builds and runs `sa-server`'s tests under all nine of its feature combinations, so
each one reports the backends it actually selected. Lints:

```bash
cargo clippy --all-targets --all-features -- -D warnings
cargo +nightly fmt --all --check    # nightly is required, see below
```

`cargo fmt` **must** run on nightly: `.rustfmt.toml` enables `unstable_features` along with
`imports_granularity`, `group_imports` and `normalize_comments`, all of which stable rustfmt
silently ignores. The nightly is pinned in `rust-toolchain.toml`; every other CI job overrides it
with stable, which is what the workspace ships from.

## The crates

| crate | what it does |
|---|---|
| [`sa-server`](sa-server/README.md) | HTTP server for testing and for serving an index directly; the only place a storage feature is read |
| [`sa-builder`](sa-builder/README.md) | builds an index from a protein TSV — the only writer in the workspace |
| [`sa-index`](sa-index/README.md) | the search itself: suffix array, k-mer table, suffix→protein mapping |
| [`protein-metadata`](protein-metadata/README.md) | protein metadata — accessions, taxa, annotations |
| [`protein-text`](protein-text/README.md) | the concatenated protein text, packed at 5 bits per residue |
| [`bitarray`](bitarray/README.md) | dense arrays of fixed-width values |
| [`fa-compression`](fa-compression/README.md) | encoding for functional annotations |
| [`binary-traits`](binary-traits/README.md) | the read/write/load traits every on-disk structure is written and read through |
| [`memory-hints`](memory-hints/README.md) | prefetch, transparent-huge-page and page-warmup hints to the memory subsystem |
| [`libsais64-rs`](libsais64-rs/README.md) | bindings to the suffix-array construction library |
| [`sa-benchmarks`](sa-benchmarks/README.md) | measurement harness (see below) |

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
sudo ./sa-benchmarks/run.sh all      # every suite in one session, into one report
```

Six suites: `defaults` (throughput at production defaults), `kmer` (what each k-mer table buys),
`stream` (how throughput depends on the number of peptides in one call), `startup` (what each
storage configuration costs before the first query), `ram` (scaling as the RAM ceiling falls) and
`threads` (whether oversubscription pays). The last three need root, for `drop_caches` and cgroup
ceilings. `ram` and `threads` are optional and are reported as skipped without it; `startup` is
not, so a rootless `run.sh all` aborts in the preflight rather than producing a partial report.

Every run writes a self-contained `report.html` — sidebar, row filter, sortable columns — next to
its results, plus `report.md` and, for `all`, a `report.json` a later run can use as its baseline.

Measurement code stays out of what ships. Nothing in the search path reads a clock or bumps a
counter — the `measure` feature that used to gate that instrumentation is gone, and its findings
are recorded in `sa-index`'s crate docs. Run-level measurement lives in `sa-benchmarks`, which
never ships at all.

## Further reading

* Crate-level `cargo doc` for `sa-index` explains the search pipeline and the feature axes in
  detail, including why the two backends exist and which decisions were measured and rejected.
* Each crate's own README, linked in the table above, covers what that crate is for and how it
  fits; the rustdoc under it carries the arguments.
