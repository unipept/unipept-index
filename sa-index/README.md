# sa-index

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

The search itself: given a peptide, find every protein that contains it. This is the crate the
Unipept API builds its endpoints on — `peptide_search::search_all_peptides_json` and friends —
rather than proxying [`sa-server`](../sa-server/README.md).

## The pipeline

1. `sa_searcher` binary-searches the suffix array for the range of suffixes sharing the peptide as
   a prefix, optionally starting from bounds looked up in a `kmer_table`;
2. it validates each candidate in that range against the text, since a sparse suffix array only
   indexes every n-th position, and I/L equating and tryptic filtering add further conditions;
3. `suffix_to_protein_index` maps surviving text positions to protein indices;
4. `peptide_search` turns those into results with accessions and annotations.

| module | what it holds |
|---|---|
| `array` | the suffix array, plain and compressed |
| `kmer_table` | precomputed suffix-array bounds per k-mer. An accelerator only; results are identical with and without it |
| `sa_searcher` | binary search, candidate validation, tryptic filtering, retrieval |
| `suffix_to_protein_index` | text position → protein index, in three representations |
| `peptide_search` | the public entry point, including JSON serialisation |

`suffix_to_protein_index` is chosen at build time (`sa-builder --mapping-style`) and recorded in the
file: **BitVec** (a bit per position plus rank; near-dense speed at ~1.25 bits per position, and the
default), **Sparse** (protein starts, binary-searched; smallest, O(log m) dependent loads) or
**Dense** (one `u32` per position; one load per lookup, but 4 bytes per residue — ~300 GB at
full-UniProt scale, larger than the suffix array itself, so small databases only).

## Storage: two backends per structure, and no opinion about which

Every storage structure has two implementations, one owned and one borrowing a memory mapping.
Both are always compiled, the searcher is generic over all three of them, and **nothing here names
a concrete one**.

| structure | owned | mapped |
|---|---|---|
| suffix array | `array::InMemorySA` | `array::MmapBackedSA` |
| protein text | `protein_text::InMemoryProteinText` | `protein_text::MmapBackedProteinText` |
| protein metadata | `protein_metadata::InMemoryProteins<T>` | `protein_metadata::MmapBackedProteins<T>` |
| suffix→protein | `suffix_to_protein_index::InMemorySuffixToProteinMapping` | `…::MmapBackedSuffixToProteinMapping` |

The choice is made once per build by the binary — `sa-server`'s `backends` module is the only place
in the workspace a storage feature is read. Selection is by type, so there is no runtime branch and
no dispatch anywhere in the search path. Sixteen combinations are constructible;
`sa_searcher::tests::every_backend_combination_returns_identical_results` builds every one and
asserts they answer identically.

## Resource limits are the caller's, and there are none here

Nothing in `peptide_search` bounds the work a request can ask for, and two things multiply:

* **The peptide list.** Every entry costs an independent search plus per-hit retrieval, and short
  peptides are the expensive ones — a single residue matches a large fraction of the index.
* **`cutoff`.** It is *not* an upper bound on work in the direction that matters. It caps a result
  set only while the match range is larger than it; once the range is smaller, the whole range is
  collected regardless. A very large `cutoff` means "return everything", not "return at most this
  many".

A caller exposing these functions to untrusted input owns that. The alphabet and length filters in
this crate are correctness filters, not limits.

## Tuning is compile-time, and there is nothing left to sweep

The searcher's performance parameters — the cross-query MLP batch, the two-pass validation batch and
both prefetch distances — are constants in `src/sa_searcher/tuning.rs`. They were runtime fields
until the full-database run at `660befd7ee` could not separate three of them from noise anywhere
(`validate_batch`: 0 of 40 contexts cleared their own floor; `prefetch_threshold` ×
`retrieval_prefetch_distance`: 0 of 80 pairs), and the fourth had no value that won everywhere.
`tuning.rs` records which sweep retired each one and argues the case before you spend time
restoring a runtime path.

## No instrumentation

Nothing in the search path reads a clock or bumps a counter. It used to, behind a `measure` feature
whose atomics perturbed the very numbers they produced (~2% at `mlp_batch=1`); the feature is gone
and the findings it settled are recorded in the crate docs. Run-level measurement lives in
[`sa-benchmarks`](../sa-benchmarks/README.md), which never ships.

The crate docs are the long form of all of this: why the code is written the way it is, what LTO
was measured at and rejected for, and what changes when the index does not fit in RAM.

## Where it sits

Depends on `binary-traits`, `bitarray`, `protein-text`, `protein-metadata`, `fa-compression` and
`memory-hints`. Used by `sa-builder`, `sa-server`, and by the Unipept API directly.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p sa-index --open`
