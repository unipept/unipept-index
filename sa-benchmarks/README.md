# sa-benchmarks

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

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

Listed in the order `all` runs them — cheapest and root-free first, so a session that is going to
fail on privileges has already produced the throughput numbers by the time it gets there.

| Suite | Answers | Needs root |
|---|---|---|
| `defaults` | What this version does at production defaults, across the length regimes and every storage arm, varying only `equate_il` and `tryptic`. The regression gate. | no |
| `kmer` | What each k-mer table buys against attaching none, per length regime, with the table's resident cost beside the win. | no |
| `stream` | How throughput depends on the number of peptides in one call — the only suite that varies the query count, which every other one holds fixed. | no |
| `startup` | What each storage arm costs before it can answer the first query, from cold. | yes, and **not optional** |
| `ram` | How the storage arms scale as the RAM ceiling falls, and where they cross over. | yes |
| `threads` | Whether thread oversubscription pays, by how much, and what it costs when RAM is ample. | yes |
| `all` | Every suite in one session, into one `report.md` + `report.json`. | partially — see below |

The three root-requiring suites are not equivalent. `ram` and `threads` are listed as `optional` in
`suites/all.toml`, so without root they are **skipped and reported as skipped**. `startup` is not:
its arms share a page cache and it has to drop that cache between them, so a rootless run of it
would measure the cache rather than the arm. A rootless `run.sh all` therefore **aborts in the
preflight** rather than producing a partial report — run the root-free suites by name instead.

Every suite is in `all`, so one session produces one report. A suite's answer only means something
against what the same box did on the same commit at the shipped settings, and splitting a session
compares two machines' moods rather than two configurations. Any suite still runs on its own when
that is all you need.

Five storage arms run in every suite, and they are a **ladder** rather than independent choices —
each one preloads everything the one before it does, plus one more structure: **mmap** (everything
mapped), **pprot** (+ the protein metadata, which is what retrieval walks), **ptext** (+ the protein
text, which search reads once per residue compared), **pmap** (+ the suffix-to-protein mapping) and
**preloaded** (everything, the suffix array included).

Read them as rungs, not as rivals: the difference between two adjacent arms IS what that one
structure is worth. `ptext` and `pmap` tie at the top and `pprot` is a stopping point to walk past,
not a deployment target — it closes the retrieval gap and none of the search one.

`tryptic` is swept in every suite; `equate_il` only in `defaults`. Both are what the caller asked
for rather than tuning, and both change the work — but not in the same way, which is what decides
where each is swept. From the instrumented cells of `defaults` at 660befd7ee:

* **`tryptic` changes the SHAPE of a request.** Acceptance collapses from 6-92% of the candidates
  examined to 0.1-0.8%, so retrieval all but disappears and a request becomes almost pure search.
  It does *not* reduce the candidates examined — those are level on `mixed` and `small`, and 1.7x
  and 7x **higher** on `large` and `medium`, because almost nothing is accepted and the cutoff
  never stops the scan early.
* **`equate_il` changes the VOLUME.** It cuts candidates examined roughly 2-3x (191M to 81M on
  `small`, 1.35M to 0.44M on `large`) while leaving acceptance where it was, so the phase mix a
  knob is measured against barely moves.

A knob's answer therefore has to hold across tryptic before it is an answer about the default,
which is why every suite crosses it; `defaults` is where the caller-visible cost of `equate_il` is
priced instead. Crossing `equate_il` everywhere cost 45% of the sweep and changed a resolved
verdict in 2 of 23 contexts.

`defaults` answers one question at high precision; `kmer` and `stream` each answer one more, at
whatever precision fits. That split is deliberate: a regression gate wants few cells resolved
tightly enough that a 3% move is visible, and a sweep wants many cells at whatever precision the
budget allows. One grid cannot be both, and when `defaults` tried to be, it was the gate that lost
— a 34-row matrix re-measuring the k-mer and batch questions on every commit.

`run.sh all` runs each suite in order into one session directory, sharing one built binary per arm.
Optional suites it cannot run on this machine are **skipped and reported as skipped**, never
silently omitted — so a rootless session that got past the preflight shows filled
`defaults`/`kmer`/`stream` sections and an explicit `not run — needs root` for `ram` and `threads`.

## The report

Every run — a single suite or the whole session — writes **`report.html`** next to its results, plus
`report.md` for pasting into a PR. `all` additionally writes `report.json`: which commit, which
suites ran, and what each cost.

`--baseline` names the session **directory**, not that file. Each suite reads its own jsonl out of
it and diffs cell against cell, so the comparison is per cell and needs the records rather than the
summary.

The page is one self-contained file with no external references, so it survives being `scp`'d off
the benchmark server and opened from disk.

**Every suite opens with its answer.** A row of stat tiles states what it found — the best value
against the shipped one, the gain, and how many of its contexts cleared their own floor — followed
by one sentence reading it: a value that wins everywhere past its floor is a candidate, one that
wins somewhere is a deployment choice, one that wins nowhere changes nothing. Then the figure, then
the grid it came from, folded.

**One figure per suite, not one per length regime.** The regimes are panels of a single
small-multiple grid sharing one scale and one legend. They used to be four sibling sections, each
scaled to its own maximum — 176 of the page's 204 figures, and the one comparison the split invites
was the one it made impossible. Knob curves facet the same way, at (regime x context), so no panel
carries more than the storage arms themselves; the old single-axes version drew twenty-four lines,
three times what the palette can tell apart.

**Knob planes are still rendered, though no current suite emits one.** A two-knob plane carries a
switch choosing what a cell says: `throughput` on a ramp scaled to that plane, or `vs shipped`, the
signed difference against the shipped pair, neutral inside the noise floor. The second is the
finding and the first is what makes it readable, so it opens on throughput. Because the ramp is
scaled per plane, each one prints its own spread beside the floor: a plane that varies by less than
its noise has no shape, whatever the colours suggest. That machinery outlived the five knob suites
it was built for — see [Re-tuning](#re-tuning) — and is what a restored one would draw into.

**Only `defaults` times the whole request.** A request has three phases — `search`,
`retrieval`, then `response` (unpacking each hit's annotations and writing the JSON).
Phase 3 is measured only where a sweep sets `response = true` in its `[[sweep]]` block, which is
`defaults.toml` and nowhere else, because paying for it on every cell of every
knob sweep would cost more than the knobs being measured. So the `time split` in the other suites
has two segments rather than three, and every throughput number outside `defaults` describes
search plus retrieval only. What fraction of a real request that is — 19% to 96% on the
full-database run at `660befd7ee`, depending almost entirely on how big the answer is — is the
`what a request actually costs` section of `defaults`, and it is the number to read every other
suite's gains through.

**The phase split is a share, not a duration.** Composition is a ratio, and the length regimes
differ by two orders of magnitude in absolute time, which left the short ones a few pixels tall
exactly where the split mattered. Every column is normalised to itself and the absolute
milliseconds stay in the hover; the magnitudes are the `search time` and `retrieval time` variants
of the same switch.

Other things the page does:

- a **sidebar** built from the report's own headings, with a status dot per suite, so a skipped
  `ram` is visible before you scroll;
- a **filter box** that hides non-matching rows across every table at once, plus per-table chips —
  how a ninety-row sweep becomes the four rows in question;
- a **search-mode control** in the header, since `tryptic` is swept in every suite and half the page
  is therefore a tryptic figure beside its non-tryptic twin. It hides the figures of the mode you
  are not reading and moves the `tryptic` chip on every table that has one, so the tables narrow
  with the pictures. It opens on non-tryptic and remembers the pick. A figure that carries tryptic
  as an axis of its own — the `defaults` grid — is not a duplicate of anything and is never hidden;
- **click-to-sort** columns, parsing the numbers out of `1,234,567` and `+4.2%`, with cells that
  aren't numbers (`VOID`, `did not fit`) sorting to the bottom rather than interleaving;
- **marked cells** for the outcomes that must not be skimmed past: `VOID`, `REGRESSION`, did-not-fit,
  and signed deltas;
- **hover glossaries** on the column headings, which is where the definitions of `band`, `slots`,
  `drift` and `floor` live rather than in a paragraph repeated under every suite.

The renderer never changes what a suite said — it classifies cells by matching the strings the
suites already emit, so there is still one analysis behind the terminal, the markdown and the page.

### Colour

The storage arms take **one hue each** (`--arm-1..5` in `bench/html.py`): `mmap` blue, `pprot`
orange, `preloaded` green, `ptext` purple, `pmap` amber. The hue is keyed by arm NAME in
`charts.ARM_HUE`, not by position, so inserting an arm into the residency ladder cannot silently
recolour the ones after it — an arm keeps its colour across releases and whatever a filter leaves on
screen. Everything else categorical takes `--s1..s5` in fixed order. `bench/selftest.py` fails if any
legend shows two series in one colour.

They were a single-hue ordinal ramp, because the arms genuinely are an ordinal axis — how much is
resident. It did not survive contact: reading a ramp means judging which of two blues is darker, and
`mmap` against `pprot` is the one comparison this report exists to make. Separate hues separate at a
glance, at facet size, and in a screenshot; the residency ORDER still lives in `charts.ARM_ORDER`,
where it can be read exactly instead of estimated from a tint.

The hues for the two newest arms are taken from the data-viz reference palette rather than
invented — no ordering of five hand-picked ones cleared the gates.

Before substituting any of them, re-run the data-viz validator — the arms are a categorical palette,
so they face the categorical gates (lightness band, chroma floor, CVD separation, normal-vision
floor, contrast). Validate them **in draw order**, `charts.ARM_ORDER`, not in the order the suite
file declares them: bars and lines are judged on the ADJACENT pairlist, so which hue sits next to
which is the question, and the same five hues pass or fail depending on it. The current five pass in
both modes (worst adjacent CVD ΔE 9.1 light / 8.4 dark against an ≥8 target; normal-vision 22.9 /
19.8 against ≥15). Three of them sit under 3:1 on the light surface and the documented relief
applies — every chart here sits beside the table holding the same numbers, and every legend is
labelled.

### Re-rendering without re-running

Changing how the report LOOKS never needs the sweep again — `--report-only` rebuilds `report.html`
and `report.md` from a finished session's jsonl:

```bash
./sa-benchmarks/run.sh all --report-only --out <scratch>/<commit>-<ts>
```

## Setting up a machine

```bash
cp sa-benchmarks/profiles/example.toml sa-benchmarks/profiles/local.toml
$EDITOR sa-benchmarks/profiles/local.toml     # index dir, peptide files, k-mer tables, scratch
```

Suites name peptide files and k-mer tables (`mixed`, `small`, `k6`), never paths, so the same suite
definition runs anywhere a profile exists. Every path is validated at startup, because a sweep that
discovers a missing bucket four hours in has wasted four hours.

Every peptide file comes from the `generate-peptides` binary, which samples real subsequences out
of `proteins.bin` so every query is a guaranteed hit. The length band is an argument, so the
bucketed files are four invocations rather than a post-processing step — and each one gets the
count it was asked for, instead of however many happened to fall in the band:

```bash
cargo build --release -p sa-benchmarks --bin generate-peptides
GP=target/release/generate-peptides
IDX=<index_dir>                                 # the same one the profile names

$GP -i $IDX -o $IDX/../peptides/peptides_5_50.txt --amount 1000000               # mixed
$GP -i $IDX -o $IDX/../peptides/small.txt        --amount 200000 --min-len 5  --max-len 9
$GP -i $IDX -o $IDX/../peptides/medium.txt       --amount 200000 --min-len 15 --max-len 25
$GP -i $IDX -o $IDX/../peptides/large.txt        --amount 200000 --min-len 35 --max-len 50
```

`mixed` needs `runs * amount` lines because the single-mode suites read it sequentially — a million
covers `runs = 100, amount = 10000`. The three buckets are only read by the matrix suites, which
re-read the same prefix per cell, so they need the largest single `amount` and no more.

`profiles/local.toml` is gitignored. Results never go in the repo — they land under the profile's
`scratch`.

## Before anything runs: the preflight

Every run — one suite or `all` — starts by printing what it is about to do and whether the box can
do it. Nothing is built until that block is on screen, and a `FAIL` in it aborts the session.

```
== preflight — 6 suites ==

  ok    host         AMD EPYC 7502P · 64 cores · Linux 6.1.0 x86_64
  warn  tree         DIRTY — these numbers are not attributable to that commit
  warn  load         14.02 / 11.80 / 9.61 (1 / 5 / 15 min) on 64 cores — a co-tenant job ...
  ok    privileges   running as root
  ok    cgroups      v2 memory controller available · systemd-run present
  ok    index        184.20 GB in /mnt/data/uniprot-2026-01/suffix-array
  ok    scratch      412.6 GB free at /mnt/data/tmp/bench

  suite         processes   grid cells        queries   status
  defaults              4           64     12,800,000   run
       queries (needs/has): mixed 10,000/1,000,000 · small 10,000/138,062 · ...
  ...
  ram                  24           24     24,000,000   skip   needs root (cgroup ceilings ...)
  total                28         1006    184,335,500   4 suite(s) to run

  ok to run, with the warnings above
```

Two things it exists for. **Nothing that can be asked up front is asked late**: under `all` the
suites that need root run last, so without this a session could spend the evening on `defaults` and
only then discover it was never going to be able to run `ram`. The peptide supply, the cgroup probe,
the toolchain, the free space and every profile path are checked once, before the first build, for
every suite in the session. **The size is stated before it is paid for**: "28 processes, 1,006 grid
cells, 184 million queries" is the difference between a coffee break and an overnight run.

A `warn` prints and the run continues — a dirty tree or a busy box invalidates a comparison, but
whether to run anyway is the operator's call. A `FAIL` (a short peptide file, a missing k-mer bucket,
a non-optional suite this machine cannot run) stops the session with nothing run. The
`status` column says `run`, `skip` (optional and blocked — it will be reported as skipped) or `FAIL`,
and a resumed session says how many of each suite's cells already have results.

Three ways to ask without running:

```bash
./sa-benchmarks/run.sh all --check       # the checks and the plan table, nothing else
./sa-benchmarks/run.sh ram --dry-run     # the same, plus the per-cell plan and each cell's command
./sa-benchmarks/run.sh all --check && nohup sudo ./sa-benchmarks/run.sh all &
```

The toolchain check is asked the way `bench/build.py` invokes cargo — as the invoking user under a
login shell — because `sudo` replaces PATH with `secure_path` and a rustup toolchain lives in that
user's `~/.cargo/bin`. It never fails the run either way: refusing to start a session that would
have built fine is worse than the seconds cargo takes to say so itself, and the builds all happen
before the first cell.

`--check` builds nothing, writes nothing — not even a session directory — and exits **1** when the
run could not start, so it gates an overnight sweep from a shell script. It reports a fresh session,
so nothing counts as already complete; pass `--out <dir>` to ask about a session being resumed.
`--dry-run` never gates (a plan is usually read on a machine other than the one that will run it).

## While it runs: the progress bar

```
[############..................]  41.3%  ·  kmer 2/3 in 6.1m  ·  elapsed 38.2m  ·  eta ~54.1m
```

One line, re-drawn in place at the bottom of the terminal, spanning the whole session rather than
the current suite, with the current suite's own clock beside the session's. Each suite's bar is left
on screen when it ends, so the scrollback holds one line per suite saying what it cost. Progress is measured in **timed queries, not cells**: a matrix suite's process
sweeps a whole grid and a tryptic cell at a fifth the query count is a fifth of the cost, so counting
cells would make the bar lurch and the ETA lie. The ETA extrapolates from the queries this session
has already completed and the seconds they took — it is empty until the first cell finishes, and
cells skipped as already-complete are kept out of that average, since they cost no time. A suite that
is skipped or aborts leaves the total, so the bar still reaches 100%.

When the output is not a terminal (`nohup`, CI, a pipe to `tee`) it degrades to one `progress:` line
per cell with the same numbers, rather than filling the log with escape codes.

### Where the time went

Each suite's wall clock is reported four ways: live in the bar, as `-- kmer: ok in 21.4 min` when it
finishes, in the report's **Suites** table (with each suite's share of the session), and in
`report.json`. Under `all` it is also written to `<session>/timings.json` **after every suite**, so a
session that is interrupted — or one still running — can still be asked where its evening went:

```bash
cat "$SESSION/timings.json"
```

A suite's clock includes the arm builds it paid for, and every later suite reuses those binaries, so
the first suite in the order carries the build cost for all of them. That is the honest number for
"what does adding this suite to the session cost" only for suites after the first.

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

`slots` is the one worth insisting on, because the alternative to it is a systematic error that
reads exactly like a result. Every suite runs one process per arm — in matrix mode a single
invocation emits *every cell* of that arm — so the reps inside it are not independent samples of
the arm. They are one sample, repeated. Whatever state that process happened to start in shifts all
of its cells together, in the same direction, and a table of per-cell bands will report that as a
resolved effect in a dozen cells at once. The palindrome ordering exists to bound it: run each arm
twice and let the gap between its two invocations set the floor. A suite that ran only one
invocation per arm reports `slots` as `-`, and its arm-vs-arm deltas are worth less than they look.

The same failure has a page-cache form, which is why `startup` drops the cache before every
configuration and needs root to do it. Suites run their arms in a fixed order over one shared page
cache, so an arm that runs after another has already faulted the index in is measuring the cache
rather than itself. Reading `GB/s` beside `warmup` is what catches it: two arms sweeping the same
structure do the same work, and if one of them is an order of magnitude faster per byte, it was
handed a warm cache and the pair is not comparable.

**Throughput here is search plus retrieval, and that is not a whole request.** Production then turns
every hit into the `ProteinInfo` it returns — an fa-compression decode of the annotations plus a
`String` for the accession — and serialises the result to JSON. `defaults` times that, in its "what
a request actually costs" section, and reports the share. A knob that buys 20% of the measured part
buys less than 20% to a caller, and that ratio is the first thing to read before any verdict in the
rest of the report.

Schema v13 changed how the phase is measured: it now runs in production's shape, with the decode
parallel across peptides. v12 timed the decode serially while timing search and retrieval in
parallel, which overstated the response share by up to the core count — **any share quoted from a
v12 session, including the "~12%" that used to sit here, is wrong and needs re-measuring.**

## Re-tuning

There is nothing to sweep. The searcher's performance parameters — the cross-query MLP batch, the
two-pass validation batch, and both prefetch distances — are compile-time constants in
`sa-index/src/sa_searcher/tuning.rs`, and that file records which sweep retired each one and what it
found. Five suites lived here (`mlp`, `validate`, `prefetch`, `mlp_validate`, `combos`), one per
knob plus their crosses; they went when the knobs did.

Re-opening one means restoring the runtime path first: give the searcher a field, thread it to the
constant's use, and teach this harness to set it. That is deliberately more work than editing a
number — a parameter no measurement could move is one that costs more to carry than to rebuild — and
`tuning.rs` argues the case before you spend the time.

What the harness still sweeps is everything that is not a searcher setting: the workload
(`defaults`, `stream`), the index build (`kmer`), and the machine (`startup`, `ram`, `threads`).

Read each one's **resolution** table before its results. These suites interleave a reference cell
through every process, so drift is measured and removed rather than folded into the numbers, and
the reference's leftover scatter is the floor on what that process could resolve at all. A residual
wider than the effects below it means the run answered nothing, which on a shared box is the common
outcome rather than the rare one.

## Adding a suite

1. Write `suites/<name>.toml`: `[[arms]]` (feature sets to build) and `[axes]` (values to sweep),
   or `[[sweep]]` blocks for `mode = "matrix"`. Axis names have defined effects and an unknown one
   is an error — see `bench/config.py`.
2. Optionally write `bench/suites/<name>.py` with an `analyse(report, suite, records, out_dir)`
   function. Without one, a generic table of every cell is printed.
3. Add it to `suites/all.toml` if it belongs in the master run.

### Suites that sweep many cells

A suite in `mode = "matrix"` loads the index once per process and sweeps in-process, which at
full-database scale is the difference between a sweep and a day of index loads. Every such suite is
built from `[[sweep]]` blocks (see `bench/grid.py`) — the harness carries no grid of its own, so a
sweep is described in exactly one place. A narrow suite is one block (`defaults`, `kmer`), a
wider one is several, and the rule that keeps a wider one affordable is that blocks
**add** rather than multiply:

* one block varying one coordinate against **one** background — the curve;
* one block at **one** point against every context — what each context costs.

That is a few dozen cells where a cross-product is thousands, and it covers strictly more. Widening
a block is a line of TOML; narrowing a cross-product after the run is not, because the run is
already over.

Two things such a suite gets for free, both of which make lower per-cell rep counts safe:

* **`base_every`** re-measures the reference configuration through the process. Drift then becomes a
  measured series: every cell is rescaled by the reference interpolated to its own position, and the
  reference's leftover scatter becomes that process's resolution floor, reported in that suite's
  **resolution** table. Set it on any suite whose process runs more than a handful of cells —
  `kmer` and `stream` both do, and on a shared box their processes routinely
  move by more than the effects they are measuring. Not `defaults`: its two arms are separate
  processes, so its exposure is drift *between* them, which an in-process cadence cannot address.
* **`runs_target_band`** stops a cell once its own spread is tight enough to read. A fixed rep count
  spends the noisiest cell's budget on every cell.

Only `threads` and `ceiling_gb` may be `[axes]` in matrix mode — a rayon pool is built once per
process and a cgroup scope wraps one, so neither can change while an index stays loaded. Everything
else belongs in a block, where it costs a cell rather than an index load. Pinning threads needs no
root; only a memory ceiling does.

`analyse` appends to a `Report` and never prints, which is what lets the terminal and `report.md`
show exactly the same analysis.

```bash
PYTHONPATH=sa-benchmarks python3 -m bench.selftest
```

runs the `defaults`, `kmer`, `ram` and `threads` analyses against fabricated records with known
shapes. `stream` and `startup` have no fixtures yet, so their `analyse` is not covered — worth
knowing for `startup`, which is not optional and runs first on the benchmark server. `ram` and
`threads` need root,
cgroup v2 and Linux, so on a development machine this is the only thing that exercises them — and
they are where the crossover detection, the void-cell rule and the fault-flatness warning live. It
asserts that a cell which was OOM-killed, one whose cap did not bind, and one that never ran each
appear in the output and stay distinguishable from one another, and that all three renderings
survive: the page is well-formed, self-contained, and marks those cells.

The grid expander is checked separately: that two blocks describing the same measurement
collapse while the same measurement at two precisions does not, that cells land in the process their
block named, and that a misspelled context key — or a suite still carrying a `tune` table — is an
error rather than a silently different sweep.

## Measurement code and shipping binaries

There is no hot-path instrumentation. `sa-index` used to carry a `measure` feature that swapped the
searcher's zero-sized counters for real atomics; it perturbed throughput by ~2% wherever it was on,
so only one arm of one suite ever used it, and it was removed once the two questions it answered —
the binary-search/range-scan split, and tryptic's candidate acceptance rate — were settled. Both
findings are recorded above and in `sa-index`'s crate docs.

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
sudo ./sa-benchmarks/run.sh startup            # cold by default; see below for why it has to be
sudo ./sa-benchmarks/run.sh ram
sudo ./sa-benchmarks/run.sh threads
```

Two things to check before starting a long run: the tree is clean (a dirty tree makes the numbers
unattributable, and the driver says so), and nothing else is running on the box (a co-tenant job
invalidates every comparison — the driver warns when the load average is high).

Every suite renders its own charts into `report.html`; there is nothing extra to run and no
third-party dependency, so the driver works on a benchmark server without a virtualenv.

## Where it sits

Depends on `sa-index`, `sa-server` and `protein-text`, plus `clap`, `serde`, `rand`, `sysinfo` and
`rayon`. Its storage features forward to `sa-server` and nowhere else, which is what keeps this
crate and the loaders it calls from ever resolving to different backends. Nothing depends on it,
and it is excluded from `default-members`, so it never ships.

The Python driver has no third-party dependencies — it needs python >= 3.11 for `tomllib` and
nothing else, so it runs on a benchmark server without a virtualenv.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p sa-benchmarks --open`
