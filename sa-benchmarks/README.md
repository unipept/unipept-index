# sa-benchmarks

Everything needed to measure the suffix-array index lives in this directory: the Rust harness, the
Python driver, the suite definitions, and the machine profiles.

The crate is a workspace member but is **excluded from `default-members`**, so a plain `cargo build`
or `cargo test` skips it and its extra dependencies. Build it explicitly with `-p sa-benchmarks`.

```
./sa-benchmarks/run.sh defaults          # the regression gate — run this after any search change
./sa-benchmarks/run.sh ram --dry-run     # plan a sweep without touching the index
sudo ./sa-benchmarks/run.sh all          # every suite, into one report.md
```

## The suites

| Suite | Answers | Needs root |
|---|---|---|
| `defaults` | What this version does at production defaults, across small/medium/large peptides, on both storage backends. The regression gate. | no |
| `detail` | Where the time goes inside a search: binary search against range scan, candidate acceptance rate, and the MLP batch curve. Instrumented, so its throughput is perturbed. | no |
| `startup` | What each of the nine storage configurations costs before it can answer the first query. | only for `--cold` |
| `ram` | How the storage arms scale as the RAM ceiling falls, and where they cross over. | yes |
| `threads` | Whether thread oversubscription pays, by how much, and what it costs when RAM is ample. | yes |
| `all` | Every suite in one session, into one `report.md` + `report.json`. | partially — see below |

`run.sh all` runs each suite in order into one session directory, sharing one built binary per arm.
Suites it cannot run on this machine are **skipped and reported as skipped**, never silently
omitted — so a report from a laptop shows filled `defaults`/`detail`/`startup` sections and an
explicit `not run — needs root` for `ram` and `threads`.

## The report

Every run — a single suite or the whole session — writes **`report.html`** next to its results, plus
`report.md` for pasting into a PR. `all` additionally writes `report.json`, which is what a later run
consumes as its `--baseline`.

The page is one self-contained file with no external references, so it survives being `scp`'d off
the benchmark server and opened from disk:

- a **sidebar** built from the report's own headings, with a status dot per suite, so a skipped
  `ram` is visible before you scroll;
- a **filter box** that hides non-matching rows across every table at once — how a 34-row grid
  becomes the four rows in question;
- **click-to-sort** columns, parsing the numbers out of `1,234,567` and `+4.2%`, with cells that
  aren't numbers (`VOID`, `did not fit`) sorting to the bottom rather than interleaving;
- **marked cells** for the outcomes that must not be skimmed past: `VOID`, `REGRESSION`, did-not-fit,
  and signed deltas.

The renderer never changes what a suite said — it classifies cells by matching the strings the
suites already emit, so there is still one analysis behind the terminal, the markdown and the page.

## Setting up a machine

```bash
cp sa-benchmarks/profiles/example.toml sa-benchmarks/profiles/local.toml
$EDITOR sa-benchmarks/profiles/local.toml     # index dir, peptide files, k-mer tables, scratch
```

Suites name peptide files and k-mer tables (`mixed`, `small`, `k6`), never paths, so the same suite
definition runs anywhere a profile exists. Every path is validated at startup, because a sweep that
discovers a missing bucket four hours in has wasted four hours.

The length-bucketed peptide files come from `bucket_peptides.sh`; the peptide files themselves come
from the `generate-peptides` binary, which samples real subsequences out of `proteins.bin` so every
query is a guaranteed hit.

`profiles/local.toml` is gitignored. Results never go in the repo — they land under the profile's
`scratch`.

## Reading the output

Every suite prints its tables with the noise floor beside the deltas, and follows them with the
prose that says how to read them. Four statistics recur, and none is decoration:

* **band** — half a cell's p10..p90 spread. How steady that one cell was.
* **slots** — the gap between a cell's own two invocations under a palindrome ordering. This is the
  floor on what the experiment can resolve: a delta smaller than it is *no answer*, not a small
  effect.
* **drift** — a cell's first quarter of reps against its last. A capped cell climbs out of whatever
  the page sweep left in the cache; large drift means it never reached steady state.
* **void** — a capped cell whose RSS landed above its ceiling was never constrained, and is
  discarded rather than read.

The measured run-to-run noise floor on the full database is **3.9%**. Deltas below it are noise.

## Adding a suite

1. Write `suites/<name>.toml`: `[[arms]]` (feature sets to build) and `[axes]` (values to sweep).
   Axis names have defined effects and an unknown one is an error — see `bench/config.py`.
2. Optionally write `bench/suites/<name>.py` with an `analyse(report, suite, records, out_dir)`
   function. Without one, a generic table of every cell is printed.
3. Add it to `suites/all.toml` if it belongs in the master run.

`analyse` appends to a `Report` and never prints, which is what lets the terminal and `report.md`
show exactly the same analysis.

```bash
PYTHONPATH=sa-benchmarks python3 -m bench.selftest
```

runs every analysis against fabricated records with known shapes. `ram` and `threads` need root,
cgroup v2 and Linux, so on a development machine this is the only thing that exercises them — and
they are where the crossover detection, the void-cell rule and the fault-flatness warning live. It
asserts that a cell which was OOM-killed, one whose cap did not bind, and one that never ran each
appear in the output and stay distinguishable from one another, and that all three renderings
survive: the page is well-formed, self-contained, and marks those cells.

## Measurement code and shipping binaries

`metrics` (on `sa-index`) is the only gate for hot-path instrumentation, and nothing that ships
enables it: with the feature off the counters are zero-sized and every write compiles away. CI
resolves `sa-server`'s and `sa-builder`'s feature graphs on every push to prove it stays that way.

Measurement that is a property of a whole run rather than of the hot path — load timings, page-fault
counts — lives in this crate instead, which never ships at all.

## Correctness gate

```bash
bash sa-benchmarks/check_answers.sh
```

Starts `sa-server` once per storage configuration and checks that all nine return byte-identical
answers. Run it before any storage comparison: a fast arm is worthless if it answers differently,
and no suite would notice.

## Full-database runbook

On the benchmark server, from the repo root:

```bash
git pull && git rev-parse --short HEAD
$EDITOR sa-benchmarks/profiles/local.toml      # check the dated index directory name for this box

bash sa-benchmarks/check_answers.sh            # correctness first
./sa-benchmarks/run.sh all --dry-run           # eyeball the plan and the cell count
sudo ./sa-benchmarks/run.sh all                # -> <scratch>/<commit>-<ts>/report.md
```

Then, after the next change:

```bash
sudo ./sa-benchmarks/run.sh all --baseline <scratch>/<previous-commit>-<ts>
```

Anything can be run on its own, and everything is resumable at cell granularity — an interrupted
overnight sweep restarts where it stopped:

```bash
./sa-benchmarks/run.sh defaults
./sa-benchmarks/run.sh startup --cold          # first-boot rather than warm-restart load times
sudo ./sa-benchmarks/run.sh ram
sudo ./sa-benchmarks/run.sh threads
```

Two things to check before starting a long run: the tree is clean (a dirty tree makes the numbers
unattributable, and the driver says so), and nothing else is running on the box (a co-tenant job
invalidates every comparison — the driver warns when the load average is high).

## Plotting

```bash
python3 sa-benchmarks/bench/plot.py --stat median <results>/*.jsonl
```

Needs matplotlib. Nothing else here has third-party dependencies, so the driver runs on a benchmark
server without a virtualenv.
