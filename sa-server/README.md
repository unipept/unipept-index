# sa-server

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

An HTTP server over a prebuilt index: one route, `POST /search`, taking a list of peptides and
returning the proteins that contain them.

**This is not how search reaches Unipept.** The API site calls `sa_index::peptide_search` directly
rather than proxying this binary, so this server is a testing and direct-serving tool, not the
untrusted boundary. Hardening it protects nothing that runs in production; the limits belong in the
caller that exposes those functions to the network. See the note on resource limits in
[`sa-index`](../sa-index/README.md).

## Running

```bash
cargo build --release -p sa-server

./target/release/sa-server \
  --database-file proteins.bin \
  --index-file    sa.bin \
  --mapping-file  mapping.bin \
  --kmer-table-file kmer_table.bin \   # optional
  --address 0.0.0.0:3000               # optional, this is the default
```

All four files come from one `sa-builder` run and only mean anything as a set. Each is loaded and
validated on its own, so a mismatched set gets past loading; `Searcher::try_new` is what refuses to
start on one.

```bash
curl -X POST http://localhost:3000/search \
  -H 'Content-Type: application/json' \
  -d '{"peptides": ["MLPGLALLLL"], "equate_il": false, "tryptic": false, "cutoff": 10000}'
```

| field | default | meaning |
|---|---|---|
| `peptides` | — | the peptides to search, answered independently and in request order |
| `equate_il` | `false` | treat I and L as interchangeable |
| `tryptic` | `false` | return only matches at tryptic boundaries |
| `cutoff` | `10000` | stop after this many matches; the result then flags `cutoff_used` |

Every field but `peptides` is optional, so an older caller sending only a peptide list keeps getting
the same answers. The request body is capped at 5 MB — the only limit the server imposes.

## This crate owns the storage decision

`src/backends.rs` is the **only** place in the workspace where a storage feature is read. The
libraries compile both backends of every structure unconditionally and are generic over which one
they are handed; the features exist here, at the top of the dependency graph, purely to name one
concrete type per structure.

| feature | effect |
|---|---|
| *(none)* | everything preloaded into owned memory. **This is what a plain build gives you** — no crate declares `default = [...]` |
| `mmap` | map all four structures |
| `preloaded-text` | with `mmap`: pull the concatenated protein text back into owned memory |
| `preloaded-proteins` | with `mmap`: pull the protein metadata table back |
| `preloaded-mapping` | with `mmap`: pull the suffix-to-protein mapping back |

Nine configurations in all. The suffix array follows `mmap` and has no override — there is no
`preloaded-sa` — and a `preloaded-*` feature without `mmap` is a no-op, since everything is
preloaded already. Cargo features are additive and cannot be negated by a dependent crate, so they
only ever *remove* mapping, never add it.

Which combination to build is a measured question; the [workspace README](../README.md) has the
numbers. A running server reports what it was compiled with, because nothing else at runtime
reveals it:

```
Storage backends: sa=mmap text=preloaded proteins=mmap mapping=mmap
```

CI builds and runs this crate's tests under all nine combinations, so each one is checked to report
the backends it actually selected.

## Why the response is built by hand

`search_all_peptides_json` hands back the body already serialised, one chunk per peptide, built on
the rayon workers that did the search. Serialising in the handler instead would put the whole answer
— hundreds of megabytes on a large non-tryptic request — through a single-threaded
`serde_json::to_vec`, which is the largest serial stretch a request has. It also could not borrow:
a `SearchResult` holds references into the `Searcher`, so it never satisfies the `'static` bound a
response body needs.

## Where it sits

Depends on `axum`, `tokio`, `clap`, `binary-traits`, `sa-index`, `protein-text` and
`protein-metadata`. `sa-benchmarks` forwards its storage features here and nowhere else, which is
what keeps the harness and the loaders it calls from ever resolving to different backends.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p sa-server --open`
