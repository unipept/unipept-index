"""Running cells: one harness invocation each, under the right ceiling, in the right order.

The rules encoded here are the ones that separate a measurement from a pile of numbers, and each of
them was a comment in one or more of the scripts this replaces:

* **Drop the page cache before every capped cell.** cgroup v2 charges a page-cache page to the
  cgroup that FIRST faults it in, so a page left resident by the previous cell is reused without
  being charged again and the ceiling silently does not apply.
* **`MemorySwapMax=0`.** The failure mode at the floor should be a clean OOM, not swap thrash, which
  would measure the swap device instead of the thing under test.
* **An OOM is an answer.** A cell killed under its ceiling records "this arm cannot run here" and the
  sweep continues; it is a fact about the arm, not a failed run.
* **Everything is resumable.** A cell with results, or with a did-not-fit marker, is skipped, so an
  interrupted overnight sweep restarts where it stopped.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import grid
from .config import Cell, Suite
from .profile import Profile
from .rig import RigError, check_peptide_supply, drop_caches, is_root

#: How a SIGKILL — under a cgroup ceiling, the OOM killer — reaches this driver.
#:
#: Two spellings, because two things report it. `subprocess` returns the negated signal number for
#: a child killed by one, so that is what the runner compares against; `137` is the shell's
#: `128 + signal` convention, which is what the bash scripts this package replaces wrote into their
#: markers and therefore what a marker from one of those sessions still holds. Both are accepted
#: wherever a status is judged, so a marker written before this package is read the way it was
#: meant. `OOM_EXIT` is the value written into new markers, in the shell spelling, so a marker
#: stays readable by anything that predates this driver.
OOM_EXIT = 137
OOM_STATUSES = (OOM_EXIT, -signal.SIGKILL)

#: Warmup peptide count used when a suite asks for the page-touch sweep ("all") but the arm has
#: nothing mapped to touch. See `_warmup_for`.
PRELOADED_FALLBACK_WARMUP = 1000


@dataclass
class CellResult:
    cell: Cell
    status: str  # "ok" | "skipped" | "did-not-fit" | "failed"
    seconds: float = 0.0
    detail: str = ""


class Runner:
    """Executes a suite's cells into `out_dir`."""

    def __init__(
        self,
        suite: Suite,
        profile: Profile,
        binaries: dict[str, Path],
        out_dir: Path,
        echo=print,
        progress=None,
    ) -> None:
        self.suite = suite
        self.profile = profile
        self.binaries = binaries
        self.out_dir = out_dir
        self.echo = echo
        #: Optional `bench.progress.Progress`; the runner is the only place that knows when a cell
        #: has finished and what it was worth, so it is what advances the bar.
        self.progress = progress
        # Not created here: `--check` and `--dry-run` build a Runner purely to expand and weigh the
        # cell list, and a question about a run should not leave a session directory behind.
        self._grid: dict[tuple, list[dict]] | None = None

    # -- planning

    def settings(self, cell: Cell) -> dict[str, Any]:
        """Suite defaults with this cell's axis values laid over them."""
        return {**self.suite.defaults, **cell.axes}

    @property
    def grid(self) -> dict[tuple, list[dict]]:
        """`(arm, threads, ceiling_gb) -> [grid cell]` for a `[[sweep]]` suite; empty otherwise.

        Expanded from the suite file alone, so it is the same before and after the arms are built —
        which is what makes a dry run's cell count exact rather than an upper bound.
        """
        if self._grid is None:
            self._grid = self._expand_grid()
        return self._grid

    def _expand_grid(self) -> dict[tuple, list[dict]]:
        if not self.suite.sweeps:
            return {}
        return grid.expand(
            self.suite.sweeps,
            suite_axes=self.suite.axes,
            suite_files=self.suite.defaults.get("files", []),
            suite_defaults=self.suite.defaults,
        )

    def grid_for(self, cell: Cell) -> list[dict]:
        """This cell's in-process grid: the blocks that named its arm, thread count and ceiling."""
        key = (cell.arm.name, cell.axes.get("threads", "default"), cell.axes.get("ceiling_gb", 0))
        return self.grid.get(key, [])

    def applicable(self, cells: list[Cell]) -> list[Cell]:
        """Drops the processes no `[[sweep]]` block asked for.

        `[axes]` is a product over every arm, but a block may name only some of them — the thread
        ladder is swept on the mapped arm and only sampled on the preloaded one, because a preloaded
        index load is the most expensive thing in the suite and there are no faults there to
        overlap. Those combinations are absences by design, not gaps, so they are dropped here
        rather than each becoming a process that pays an index load to run nothing.
        """
        if not self.suite.sweeps:
            return cells
        kept = [cell for cell in cells if self.grid_for(cell)]
        for cell in cells:
            if not self.grid_for(cell):
                self.echo(f"  no sweep covers {cell.label} — not run")
        return kept

    def weight(self, cell: Cell) -> int:
        """Timed queries this cell will run — its share of the session's cost.

        What the progress bar apportions itself over, and what a plan totals. Not the cell count:
        cells differ by orders of magnitude, so counting them would describe a different run.
        """
        settings = self.settings(cell)
        runs, amount = int(settings.get("runs", 100)), int(settings.get("amount", 10_000))
        if self.suite.sweeps:
            return grid.query_count(self.grid_for(cell), runs, amount)
        return runs * amount

    def requirements(self, cells: list[Cell]) -> dict[str, int]:
        """Peptide-file name -> lines these cells will consume from it.

        Separate from `check_supply` so the same arithmetic can be reported before a run starts
        and enforced when it does, rather than the two drifting apart.
        """
        needed: dict[str, int] = {}
        for cell in cells:
            settings = self.settings(cell)
            if self.suite.mode == "matrix":
                # Every grid cell re-reads the same prefix of the file, so the requirement is the
                # largest `amount` any cell asks for — not the sum, and not `runs * amount`.
                for entry in self.grid_for(cell):
                    name = entry["file"]
                    wanted = int(entry.get("amount", settings.get("amount", 10_000)))
                    needed[name] = max(needed.get(name, 0), wanted)
                continue
            name = settings["peptides"]
            required = int(settings.get("runs", 100)) * int(settings.get("amount", 10_000))
            needed[name] = max(needed.get(name, 0), required)
        return needed

    def check_supply(self, cells: list[Cell]) -> None:
        """Enforces the query supply for every distinct (peptide file, runs x amount) the plan needs.

        The harness consumes lines sequentially, so a short file does not fail — it shortens the
        later reps, and the run completes having measured something else. `bench.preflight` reports
        this before the session starts; this is where it is enforced, so a cell can never run
        against a file that cannot feed it however the runner was reached.
        """
        for name, required in self.requirements(cells).items():
            check_peptide_supply(self.profile.peptide_file(name), required)

    # -- execution

    def completed(self, cell: Cell) -> bool:
        """True when `run_cell` would skip this cell: it already has results, or a did-not-fit marker.

        Resuming is cell-granular, so this is also what the plan and the progress bar have to weigh
        a resumed session by — counting work that will not be done again would predict a session
        several times longer than it is.
        """
        jsonl = self.out_dir / f"{cell.label}.jsonl"
        return (jsonl.exists() and jsonl.stat().st_size > 0) or _is_oom_marker(
            self.out_dir / f"{cell.label}.oom"
        )

    def run(self, cells: list[Cell]) -> list[CellResult]:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        # Weighed once, before the first cell: a cell that finishes writes its own jsonl, and asking
        # afterwards whether it was already complete would answer yes for everything.
        weights = {cell.label: 0 if self.completed(cell) else self.weight(cell) for cell in cells}
        if self.progress:
            self.progress.begin_suite(self.suite.name, sum(weights.values()), len(cells))

        results = []
        for cell in cells:
            result = self.run_cell(cell)
            results.append(result)
            if self.progress:
                self.progress.cell_done(
                    weights[cell.label], ran=result.status != "skipped", seconds=result.seconds
                )
        return results

    def run_cell(self, cell: Cell) -> CellResult:
        jsonl = self.out_dir / f"{cell.label}.jsonl"
        marker = self.out_dir / f"{cell.label}.oom"
        if jsonl.exists() and jsonl.stat().st_size > 0:
            self.echo(f"  skip {cell.label} (results exist)")
            return CellResult(cell, "skipped", detail="results exist")
        if _is_oom_marker(marker):
            self.echo(f"  skip {cell.label} (recorded as did-not-fit)")
            return CellResult(cell, "did-not-fit", detail="recorded earlier")
        if marker.exists():
            # A marker left by an older driver, which wrote one for every non-zero exit. It records
            # a crash, not a ceiling, so the cell is retried rather than being skipped forever as a
            # measurement nobody made.
            self.echo(f"  retrying {cell.label} (its marker records exit {_marker_exit(marker)}, not an OOM)")
            marker.unlink()

        if self.suite.drop_caches and is_root():
            drop_caches()

        command = self._command(cell)
        log = self.out_dir / f"{cell.label}.log"
        self.echo(f"[{time.strftime('%H:%M:%S')}] {cell.label}")

        started = time.monotonic()
        with log.open("w") as handle:
            completed = subprocess.run(
                command, stdout=handle, stderr=subprocess.STDOUT, env=self._env(cell), check=False
            )
        elapsed = time.monotonic() - started

        if completed.returncode == 0:
            return CellResult(cell, "ok", elapsed)

        if completed.returncode in OOM_STATUSES:
            # ONLY an OOM writes a marker, and the distinction is load-bearing twice over. A marker
            # is a RESULT — `records.unfit_cells` turns it into "this arm cannot run at this
            # ceiling" in the `ram` and `threads` tables — and it is also what `completed()` reads,
            # so a cell that has one is never attempted again. Writing one for every non-zero exit
            # meant a panic or a bad path was reported as a fact about the arm's memory behaviour
            # and then permanently skipped on resume, silently, with no warning the second time.
            #
            # The marker carries the cell's dims, not just its exit code: a cell that produced no
            # records still has to appear in the report, and recovering that from the file name
            # would reintroduce the parsing the dims envelope exists to remove.
            # `OOM_EXIT`, not `completed.returncode`: the two spellings mean the same thing, and
            # writing the shell's is what keeps a marker readable by anything that predates this
            # driver. The status is normalised on the way in, not on the way out, so a marker never
            # carries the platform detail that a SIGKILL reaches Python as a negative number.
            marker.write_text(
                json.dumps({"exit": OOM_EXIT, "dims": cell.dims, "label": cell.label}) + "\n"
            )
            self.echo(f"  -> killed under its ceiling; recorded as did-not-fit (see {log.name})")
            return CellResult(cell, "did-not-fit", elapsed, "OOM-killed under its ceiling")

        # Anything else is a failure, not an answer. Nothing is written, so the cell is retried on
        # the next run rather than being frozen into the session as a measurement nobody made.
        self.echo(f"  -> exit {completed.returncode} (see {log})")
        return CellResult(cell, "failed", elapsed, f"exit {completed.returncode}, see {log}")

    # -- command construction

    def _command(self, cell: Cell) -> list[str]:
        settings = self.settings(cell)
        command = self._wrapper(settings)
        command.append(str(self.binaries[cell.arm.name]))
        command += ["--index-dir", str(self.profile.index_dir)]
        command += ["--output", str(self.out_dir), "--label", cell.label]
        command += ["--suite", self.suite.name]
        for key, value in cell.dims.items():
            command += ["--dim", f"{key}={value}"]
        command += ["--runs", str(settings.get("runs", 100))]
        command += ["--amount-of-peptides", str(settings.get("amount", 10_000))]
        command += ["--max-matches", str(settings.get("max_matches", 10_000))]

        if self.suite.mode == "matrix":
            command += self._sweep_args(cell, settings)
        else:
            command += ["--peptide-file", str(self.profile.peptide_file(settings["peptides"]))]
            command += self._single_args(cell, settings)

        if settings.get("no_theoretical_memory"):
            command.append("--no-theoretical-memory")
        return command

    def _sweep_args(self, cell: Cell, settings: dict[str, Any]) -> list[str]:
        """A `[[sweep]]` suite: the driver expands the grid and hands the harness the cell list.

        The grid file is written beside the results rather than to a temporary path, because it is
        the record of what this process was asked to run — a report that cannot be traced back to a
        cell list is a report that cannot be reproduced.
        """
        cells = self.grid_for(cell)
        if not cells:
            raise RigError(
                f"no [[sweep]] block covers {cell.label} — every cell of a matrix suite must have "
                f"work to do, or the process pays an index load for nothing"
            )

        # Suites name peptide files by profile key ("mixed"); the harness knows them by the stem of
        # the path it was handed ("peptides_5_50"). Translated here, at the boundary between the
        # two, rather than by teaching either side about the other's naming. The key travels along
        # as `bucket` and is what the records carry, so a report reads in the suite's vocabulary on
        # every machine — the stem is one profile's filename, not a property of the measurement.
        files = sorted({entry["file"] for entry in cells})
        paths = {name: self.profile.peptide_file(name) for name in files}
        cells = [{**entry, "bucket": entry["file"], "file": paths[entry["file"]].stem} for entry in cells]
        # In `grids/`, not beside the results: `records.load_dir` globs `*.jsonl` over the results
        # directory, and a grid file is JSONL too — left there it would be read back as a few
        # hundred records with no measurements in them.
        grid_path = grid.write(cells, self.out_dir / "grids" / f"{cell.label}.jsonl")

        args = ["--grid-file", str(grid_path)]
        args += ["--matrix", "--matrix-files", ",".join(str(path) for path in paths.values())]
        # Matrix mode always page-sweeps; this is what additionally asks it for the pipeline half.
        # It used to be single-mode only, which meant the two modes warmed differently and the
        # difference landed on the arms unevenly — a preloaded structure lives in anonymous memory
        # that no page sweep reaches, so only real queries warm it. See `run_matrix` in main.rs.
        warmup = _warmup_for(cell, settings.get("warmup"))
        if warmup:
            args += ["--warmup", warmup]
        for k in sorted({entry["kmer_k"] for entry in cells} - {0}):
            args += self._kmer_file_arg(k)
        if settings.get("runs_target_band"):
            args += ["--runs-target-band", str(settings["runs_target_band"])]
            args += ["--min-runs", str(settings.get("min_runs", 5))]
        return args

    def _kmer_file_arg(self, k: int) -> list[str]:
        """`--kmer-file k=<path>`, or nothing when the table has to be built in-process.

        Same trade as `_kmer_args`: a table the profile does not have is built rather than dropped,
        because dropping it would quietly measure a different index while building it only costs
        startup time.
        """
        table = self.profile.kmer_table(f"k{k}")
        if table:
            return ["--kmer-file", f"{k}={table}"]
        self.echo(
            f"  note: profile '{self.profile.name}' has no pre-built 'k{k}' table, so it is built "
            f"in-process instead (same table, paid at startup)"
        )
        return []

    def _single_args(self, cell: Cell, settings: dict[str, Any]) -> list[str]:
        args: list[str] = []
        warmup = _warmup_for(cell, settings.get("warmup"))
        if warmup:
            args += ["--warmup", warmup]
        args += self._kmer_args(settings)
        if "equate_il" in settings:
            args += ["--equate-il", _bool(settings["equate_il"])]
        if "tryptic" in settings:
            args += ["--tryptic", _bool(settings["tryptic"])]
        return args

    def _kmer_args(self, settings: dict[str, Any]) -> list[str]:
        """Attaches the k-mer table this suite asked for, loading it or building it.

        A table the profile does not have is built in-process instead. Running without it
        would quietly measure a different index — the table removes most of the binary search's
        probes — whereas building it only costs startup time. The warning says which happened,
        because it makes the `startup` suite's `kmer` column mean something different.
        """
        k = int(settings.get("kmer") or 0)
        if k == 0:
            return []
        table = self.profile.kmer_table(f"k{k}")
        if table:
            return ["--kmer-table-file", str(table)]
        self.echo(
            f"  note: profile '{self.profile.name}' has no pre-built 'k{k}' table, so a k={k} "
            f"table is built in-process instead (same table, paid at startup)"
        )
        return ["--build-kmer-table", str(k)]

    def _wrapper(self, settings: dict[str, Any]) -> list[str]:
        """`systemd-run --scope` prefix imposing this cell's memory ceiling.

        Only a ceiling needs the scope. Pinning threads is setting an environment variable, which
        `_env` does directly — routing that through systemd-run too would have made every thread
        sweep need root, cgroup v2 and a working `--setenv`, none of which the thread count itself
        requires. When a cell has both, the scope is already there and carries the variable across
        it (`rig.setenv_reaches_child` is the probe that this actually happens).
        """
        ceiling = float(settings.get("ceiling_gb", 0) or 0)
        if ceiling <= 0:
            return []

        wrapper = ["systemd-run", "--scope", "--quiet"]
        threads = settings.get("threads", "default")
        if str(threads) != "default":
            wrapper.append(f"--setenv=RAYON_NUM_THREADS={threads}")
        # Swap off, so the floor is a clean OOM rather than a measurement of the swap device.
        wrapper += ["-p", f"MemoryMax={ceiling:g}G", "-p", "MemorySwapMax=0"]
        wrapper.append("--")
        return wrapper

    def _env(self, cell: Cell) -> dict[str, str]:
        """The child's environment: this process's, plus the cell's thread count when it pins one.

        Set unconditionally, even when `_wrapper` has already passed it through `--setenv`, so the
        two paths cannot disagree about what a cell ran at.
        """
        env = dict(os.environ)
        threads = self.settings(cell).get("threads", "default")
        if str(threads) != "default":
            env["RAYON_NUM_THREADS"] = str(threads)
        else:
            # An inherited value would silently pin every "default" cell to whatever the shell had.
            env.pop("RAYON_NUM_THREADS", None)
        return env

    # -- dry run

    def plan(self, cells: list[Cell]) -> list[str]:
        """The cell list and the exact command for each, without touching the index."""
        lines = [
            f"suite      : {self.suite.name} ({self.suite.mode} mode, {self.suite.ordering} ordering)",
            f"arms       : " + ", ".join(f"{arm.name}[{arm.feature_string or 'default'}]" for arm in self.suite.arms),
            f"processes  : {len(cells)}" if self.suite.sweeps else f"cells      : {len(cells)}",
            f"results    : {self.out_dir}",
            "",
        ]
        total_queries = 0
        for cell in cells:
            settings = self.settings(cell)
            runs, amount = int(settings.get("runs", 100)), int(settings.get("amount", 10_000))
            if self.suite.sweeps:
                # A process, not a measurement: what costs time is the grid inside it, and the
                # index load it pays once for all of them.
                inner = self.grid_for(cell)
                queries = grid.query_count(inner, runs, amount)
                lines.append(f"  {cell.describe()}   ({len(inner)} grid cells, {queries:,} queries)")
            else:
                queries = runs * amount
                lines.append(f"  {cell.describe()}   ({queries:,} queries)")
            total_queries += queries

        if self.suite.sweeps:
            grid_cells = sum(len(self.grid_for(cell)) for cell in cells)
            lines += [
                "",
                f"total      : {len(cells)} index loads, {grid_cells} grid cells, "
                f"{total_queries:,} timed queries",
            ]
        return lines


def _marker_exit(marker: Path) -> int | None:
    """The exit status a did-not-fit marker recorded, or None when it cannot be read.

    A marker predating the `exit` field is read as an OOM, which is what it meant when it was
    written: at that point the runner only reached the marker path from an OOM branch.
    """
    try:
        return int(json.loads(marker.read_text()).get("exit", OOM_EXIT))
    except (json.JSONDecodeError, OSError, TypeError, ValueError):
        return None


def _is_oom_marker(marker: Path) -> bool:
    """True when this marker is a RESULT — a cell killed under its ceiling — rather than a crash.

    The distinction decides both whether the cell is ever retried and whether the report says "this
    arm cannot run at this ceiling" about it, so an unreadable marker is treated as an OOM: it was
    written by a run that only wrote them for OOMs, and re-running a cell on the strength of a
    truncated file would discard a real answer. `records.unfit_cells` applies the same rule.
    """
    if not marker.exists():
        return False
    recorded = _marker_exit(marker)
    return recorded is None or recorded in OOM_STATUSES


def _warmup_for(cell: Cell, warmup: Any) -> str | None:
    """Translates a suite's warmup setting for this arm.

    "all" touches every page of every mapped region. A fully preloaded arm has nothing mapped to
    touch, so "all" there is meaningless — it degrades to the pipeline warmup, which is what the
    preloaded arm needs anyway (CPU caches and TLB, not the page cache).
    """
    if warmup in (None, "", 0):
        return None
    warmup = str(warmup)
    if "mmap" in cell.arm.features:
        return warmup
    if warmup.startswith("all:"):
        return warmup.split(":", 1)[1]
    if warmup == "all":
        return str(PRELOADED_FALLBACK_WARMUP)
    return warmup


def _bool(value: Any) -> str:
    return "true" if str(value).lower() in ("true", "1", "yes", "on") else "false"

