# The fbc9328 baseline

This is a git worktree of **fbc9328** (`Merge pull request #31 from
unipept/feature/suffix-to-protein-optimization`) — the commit before
`feature/preloaded-sa-improvements` — carrying a port of that branch's benchmark harness. It exists
for one reason: fbc9328 predates the harness entirely, so the only way to say what the branch
changed is to measure the commit before it with the same instrument.

Nothing in the branch's repository is modified by any of this. Two untracked files are added there
(`sa-benchmarks/suites/baseline.toml` and `sa-benchmarks/bench/suites/baseline.py`); delete them and
the branch is exactly as it was.

## Why no index rebuild is needed

Both trees run off **one index directory**, read-only:

| file | branch reads | fbc9328 reads |
|---|---|---|
| `sa.bin` | yes | yes — the on-disk header is unchanged (`[bits_per_value][sparseness][count][data]`), so the same file loads in both |
| `proteins.tsv` | no | **yes** — fbc9328 has no `proteins.bin`; it parses the database TSV the index was built from and reconstructs the protein store at startup. The branch's index directory still ships this file. |
| `proteins.bin` | yes | no — a branch artefact |
| `mapping.bin` | yes | **type byte only** — fbc9328 rebuilds the mapping from the text, and reads this file solely to learn which KIND to build, so both trees walk the same mapping (bitvec, on the current index) |
| peptide files | yes | yes — plain text |

## The one axis that cannot be compared: the k-mer table

The k-mer bound table arrived **with** the branch. fbc9328 has no cell that can sit opposite a
`kmer = 5` one, so the comparison runs at `kmer_k = 0` on both sides — which is why this is the
`baseline` suite and not `defaults` (that one pins `kmer = [5]`).

The harness **rejects** a grid cell asking for a table rather than running it untabled, so this
cannot drift by accident. For what the table buys on top of the delta, read `run.sh kmer` on the
branch, which already sweeps `kmer = [0, 5, 6]`.

## What the delta contains

Everything the branch did to the search path, together: the batched search orchestrator and its MLP
interleaving, the borrowed `ProteinInfo` and per-peptide JSON serialisation, retrieval and
suffix-to-protein work, and the tuning constants.

**Compare the `fbc9328` column against `preloaded`.** fbc9328 holds every structure in owned memory,
so that is its counterpart; against `mmap` it is answering two questions in one column.

Two `startup` fields are **not** comparable across the trees, and are recorded anyway because "what
did each pay before its first query" is a real question:

* `load_proteins_ms` — fbc9328 parses a ~293 MB TSV, uppercases and concatenates every sequence,
  fa-compresses every annotation and bit-packs the text; the branch reads a prebuilt file.
* `load_mapping_ms` — fbc9328 builds the mapping; the branch reads it.

The `warmup_*_bytes` columns are zero here because nothing is mapped. That is the same zero the
branch's `preloaded` arm reports, for the same reason.

Throughput, the search/retrieval phase split, `response_*` and the fault counters **are**
comparable.

## The two suites

| suite | mode | needs root | what it answers |
|---|---|---|---|
| `baseline` | matrix | no | steady-state throughput across 16 cells — four length regimes × both search options, plus the response phase |
| `baseline_startup` | single | **yes** | what each configuration pays before its first query, on a cold page cache |

`baseline_startup` is where the branch's index-build change shows up: fbc9328 has no `proteins.bin`
or `mapping.bin` and reconstructs both from the database TSV at startup, so **read its `total`
column first and its phase split second** — `sa` is the same file read the same way in both trees,
but `proteins` and `mapping` are a parse and a build on one side and a file read on the other. The
suite's own notes say this at length; it is the one place in the comparison where a column invites
the wrong reading.

## The runbook

The whole comparison is **one session with three arms**, not two runs stitched together. That
matters: `palindrome` ordering interleaves the arms so between-invocation drift is measured rather
than absorbed into an arm difference, and two separate runs cannot be interleaved.

The driver builds every declared arm from the tree it runs in, and it cannot build this one. What it
*can* do is skip an arm whose binary is already in `<session>/bin` with a matching `.features`
manifest — its resume path. `stage.sh` puts the ported binary there.

```bash
# 0. Paths. Nothing below is machine-specific; both scripts locate this tree from their own path.
BRANCH=/path/to/unipept-index                  # the feature/preloaded-sa-improvements checkout
OLD=/path/to/unipept-index-fbc9328             # this tree
SESSION=/data/bench/baseline-$(date +%Y%m%d)

# 1. Install the two suite definitions into the branch checkout (four untracked files).
bash "$OLD/sa-benchmarks/branch-side/install.sh" "$BRANCH"

# 2. THE GATE. Both trees must answer identically before any timing is worth reading.
#    Runs each tree's own sa-server against one index and diffs /search peptide by peptide.
#    Index and peptide paths come from the BRANCH tree's machine profile.
bash "$OLD/sa-benchmarks/check_answers_cross.sh" "$BRANCH" 2000

# 3. Stage the fbc9328 harness into the session's bin/.
bash "$OLD/sa-benchmarks/stage.sh" "$SESSION"

# 4. Run all three arms, from the BRANCH tree. Two suites, one session, one bin/.
cd "$BRANCH"
./sa-benchmarks/run.sh baseline --check --out "$SESSION"           # confirm the box is fit to run
./sa-benchmarks/run.sh baseline --out "$SESSION"                   # throughput
sudo ./sa-benchmarks/run.sh baseline_startup --out "$SESSION"      # time to first query (NEEDS ROOT)

# 5. Re-render either one without re-running, as often as you like.
./sa-benchmarks/run.sh baseline --out "$SESSION" --report-only
./sa-benchmarks/run.sh baseline_startup --out "$SESSION" --report-only
```

`baseline_startup` **needs root**, exactly as the branch's own `startup` suite does, because it
drops the page cache before every arm. Without that the first arm's RSS evicts the cache it filled,
the second faults the index in from the device and the third sweeps what the second left resident —
which is how two arms doing identical work come out an order of magnitude apart. It is a cold-boot
measurement; a restart on a warm box is much faster and is not what it measures.

Both suites share `<session>/bin`, so step 3 is done once and covers both.

Step 2 must be re-run whenever this tree's harness changes — the driver skips a staged arm without
checking how old it is, and a stale binary writes perfectly well-formed records.

### Getting this onto a server

This tree is a branch in the same repository, so it travels with the repository:

```bash
# on the machine that has it
git push origin baseline/fbc9328-harness

# on the server, from the branch checkout
git fetch origin baseline/fbc9328-harness
git worktree add ../unipept-index-fbc9328 baseline/fbc9328-harness
```

Everything the comparison needs is in that branch — the ported harness, both suite definitions, the
two scripts, and the branch-side files `install.sh` copies across. There is no second thing to copy.

### What the server needs

* **No machine profile in this tree.** The run is driven from the branch checkout and uses ITS
  profile for every arm, the staged one included. That is what makes the cells line up:
  `peptide_source` is identical by construction rather than by two profiles agreeing. Only a
  standalone debug run of this tree alone needs `profiles/local.toml`, and `example.toml` is the
  template for it.
* **Root**, for `baseline_startup` only. `baseline` runs unprivileged.
* **`proteins.tsv` in the index directory.** fbc9328 has no `proteins.bin`; the gate checks for it
  up front, because it is the one file the branch no longer needs at runtime and so the one an index
  directory can plausibly be missing.
* **cmake, make and network on first build**, for the answer gate only. It builds fbc9328's
  `sa-server`, which pulls in `sa-builder` -> `libsais64-rs`, whose build script clones and builds
  `libsais-packed`. Any machine that already builds this repository has this. The harness itself does
  not depend on it.
* **A default toolchain that can build fbc9328.** This commit has no `rust-toolchain.toml`, so it
  takes the box default. Set `CARGO="rustup run <version> cargo"` for both scripts if that fails.

## What was changed in the port, and what was not

`src/main.rs` is the branch's harness with the pieces fbc9328 cannot support removed. The record
schema is **identical** — `SCHEMA_VERSION`, every field, the JSONL layout — which is the whole point:
the driver joins cells on `(peptide_source, equate_il, tryptic, kmer_k, amount_of_peptides)` and
splits them by the `arm` dim, so these records land beside the branch's own arms in one table with
no translation. Verified by diffing the key sets of a record from each tree: no field differs.

| the branch has | fbc9328 has | what the port does |
|---|---|---|
| `sa_server::{ActiveSearcher, load_*_file, *_BACKEND}` | no `sa-server` library at all | inlined; the four `*_BACKEND` constants pinned to `"preloaded"`, which is what this tree does |
| `KmerTable` and the whole k-mer axis | nothing | axis removed, `kmer_k` pinned to 0, a cell asking for a table is rejected |
| `search_all_matching_suffixes_batched` | `search_matching_suffixes` per peptide | the flat `par_iter` this commit's own `search_all_peptides` uses — measuring fbc9328 with the branch's batching would attribute the branch's win to something else |
| `frame_chunks` / `json_chunk` | nothing | reimplemented with identical framing, against this commit's owned `SearchResult` |
| four storage backend traits | one build | one arm; the page sweep reports zero because nothing is mapped |
| `proteins.bin` + `mapping.bin` | `proteins.tsv` | parsed and rebuilt at startup, as this commit's server does |

Everything else — the grid parser, the rep loop, the percentile and band statistics, the fault and
RSS counters, the JSONL emission — is byte-identical to the branch's, because the port was done by
deletion rather than rewriting.

**If the branch's schema moves, this moves with it, or the baseline stops being readable.**
