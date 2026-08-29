# The c4f0f30 baseline

This is a git worktree of **c4f0f30** (`Mappings.bin renamed to mapping.bin in releases`) — the
commit that introduced the **initial, pre-optimisation mmap/preloaded split** — carrying a port of
`feature/preloaded-sa-improvements`'s benchmark harness.

It is the middle point of a three-commit progression:

| | commits from fbc9328 | storage configurations |
|---|---|---|
| `fbc9328` | 0 | one: everything owned, protein store and mapping rebuilt from the TSV at startup |
| **`c4f0f30`** | **32** | **two: `mmap` / `preloaded`, one runtime bool covering all three structures** |
| branch tip | 206 | five: per-structure selection at compile time |

Companion to the fbc9328 port on `baseline/fbc9328-harness`; the two share the staging mechanism and
the branch-side suite files, which are byte-identical in both trees. Measuring this commit separates
**having the split** from **the optimisation work done on top of it** — the branch's throughput gain
against fbc9328 is otherwise one undifferentiated number.

## Two arms, and they are all-or-nothing

At this commit the storage choice is a **runtime bool** — `sa_server::load_*_file(path, mmap)` — and
one bool covers the suffix array, the proteins and the mapping together. So:

* there are exactly two configurations, `mmap` and `preloaded`;
* the four `*_BACKEND` constants in a record always agree with each other;
* **`pprot` / `ptext` / `pmap` have no counterpart here and never will.** Mixing per structure is
  itself something the branch added.

The bool is fed from `cfg!(feature = "mmap")` rather than a CLI flag. That is not a change to what
is measured — the same two configurations, selected the same way at the same point — but the driver
picks arms by cargo feature and refuses two arms that compile to the same bytes, and this gives it
two binaries that genuinely differ. `stage.sh` asserts they do before staging.

## The index is read unchanged

`sa.bin`, `proteins.bin` and `mapping.bin` all load directly from the branch's index directory, in
both arms. Nothing is rebuilt, nothing is converted, and unlike fbc9328 there is **no
`proteins.tsv` dependency and no startup memory problem** — the prebuilt files already exist at this
commit. Verified by answering 2,000 peptides identically to the branch, in both arms and at both
`equate_il` settings.

## What the port changed, and what it did not

The record schema is **identical** to the branch's — `SCHEMA_VERSION`, every field, the JSONL
layout — verified by diffing the key sets of a record from each tree: no field differs. The port was
done by deletion where possible, so anything not listed here is byte-identical to the branch's.

| the branch has | c4f0f30 has | what the port does |
|---|---|---|
| `ActiveSearcher`, four `*_BACKEND` constants | `Searcher::new(sa, proteins, mapping)`, one runtime bool | alias + four constants derived from `cfg!(feature = "mmap")` |
| `KmerTable` and the k-mer axis | nothing | axis removed, `kmer_k` pinned to 0, a cell asking for a table is **rejected** |
| `search_all_matching_suffixes_batched` | `search_matching_suffixes` per peptide | the flat `par_iter` this commit's own `search_all_peptides` uses |
| `frame_chunks` / `json_chunk`, borrowed `ProteinInfo` | nothing; owned `ProteinInfo` | reimplemented with identical framing |
| `touch_all_pages` via `memory_hints` | nothing, and private mmap handles | accessor-driven sweep — see below |

### The page sweep is the one caveat

The branch walks the raw `Mmap` one byte per page under `MADV_SEQUENTIAL`. Neither `memory-hints`
nor access to the mmap handles exists here, so the port sweeps through the public accessors instead
— one `get` per page-worth of entries. It faults in the same pages and populates the same PTEs, so
**the warmup does its job and `baseline` is unaffected**. But:

* **`GB swept` is comparable** with the branch's (computed analytically from the same entry counts).
* **`GB/s` is NOT.** An accessor call per page is far more work than a byte load per page, and there
  is no `MADV_SEQUENTIAL`. Read this commit's sweep rate against its own other arm, never against
  the branch's. What the column exists for — telling a cold sweep from one handed a warm page cache
  — still works within a commit.

The preloaded arm sweeps nothing and reports zero, exactly as the branch's does.

## Runbook

One session, all five arms, interleaved. That matters: `baseline` is palindrome-ordered so
between-invocation drift is measured rather than absorbed into an arm difference, and separate runs
cannot be interleaved.

```bash
BRANCH=/path/to/unipept-index                 # feature/preloaded-sa-improvements
OLD9328=/path/to/unipept-index-fbc9328
OLDC4=/path/to/unipept-index-c4f0f30          # this tree
SESSION=/data/bench/progression-$(date +%Y%m%d)

# 1. install the suite definitions into the branch checkout (4 untracked files).
#    Identical in both older trees — install from either.
bash "$OLDC4/sa-benchmarks/branch-side/install.sh" "$BRANCH"

# 2. THE GATES. Both older trees must answer identically to the branch before any timing counts.
bash "$OLD9328/sa-benchmarks/check_answers_cross.sh" "$BRANCH" 2000
bash "$OLDC4/sa-benchmarks/check_answers_cross.sh"   "$BRANCH" 2000   # checks BOTH arms

# 3. stage all three old arms into the session's bin/ (once; covers both suites)
bash "$OLD9328/sa-benchmarks/stage.sh" "$SESSION"
bash "$OLDC4/sa-benchmarks/stage.sh"   "$SESSION"

# 4. run
cd "$BRANCH"
./sa-benchmarks/run.sh baseline --check --out "$SESSION"       # expect 5 arms, 20 grid cells
./sa-benchmarks/run.sh baseline --out "$SESSION"
sudo ./sa-benchmarks/run.sh baseline_startup --out "$SESSION"

# 5. re-render without re-running
./sa-benchmarks/run.sh baseline --out "$SESSION" --report-only

# afterwards
bash "$OLDC4/sa-benchmarks/branch-side/install.sh" "$BRANCH" --uninstall
```

### Cost

`baseline` is 16 cells x 100 reps x 2 slots per arm. Five arms is 2.5x a two-arm session. If wall
clock binds, drop `c4f0f30-mmap` (the mapped path is the one the branch changed least) before
dropping `c4f0f30-preloaded`, which is the arm that isolates the optimisation work.

`baseline_startup` is `sequential`, so adding arms to it later is purely additive and it resumes
into an existing session. `baseline` is palindrome and does **not** — its slot letters are
positional, so more arms renames the reverse-half labels and orphans the old files. Re-run that one
into a fresh session.

### What the server needs

* **Root**, for `baseline_startup` only.
* **A default toolchain that can build c4f0f30.** No `rust-toolchain.toml` at this commit, so it
  takes the box default; `stage.sh` and the gate both honour `CARGO="rustup run <version> cargo"`.
* Nothing else. No `proteins.tsv`, no cmake, no vendored `libsais` — `sa-server` here does not
  depend on `sa-builder`.
