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
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import TUNE_PREFIX, Cell, Suite
from .profile import Profile
from .rig import check_peptide_supply, drop_caches, is_root

#: Exit status of a process killed by SIGKILL, which under a cgroup ceiling means the OOM killer.
OOM_EXIT = 137

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
    ) -> None:
        self.suite = suite
        self.profile = profile
        self.binaries = binaries
        self.out_dir = out_dir
        self.echo = echo
        out_dir.mkdir(parents=True, exist_ok=True)

    # -- planning

    def settings(self, cell: Cell) -> dict[str, Any]:
        """Suite defaults with this cell's axis values laid over them."""
        return {**self.suite.defaults, **cell.axes}

    def preflight(self, cells: list[Cell]) -> None:
        """Checks the query supply for every distinct (peptide file, runs x amount) the plan needs.

        Up front, because the harness consumes lines sequentially: a short file shortens the later
        reps rather than failing, so the run completes and measures something else.
        """
        needed: dict[str, int] = {}
        for cell in cells:
            settings = self.settings(cell)
            if self.suite.mode == "matrix":
                # Matrix mode re-reads the first `amount` lines for every cell, so it needs `amount`
                # lines, not `runs * amount`.
                for name in self.suite.matrix.get("files", []):
                    needed[name] = max(needed.get(name, 0), int(settings.get("amount", 10_000)))
                continue
            name = settings["peptides"]
            required = int(settings.get("runs", 100)) * int(settings.get("amount", 10_000))
            needed[name] = max(needed.get(name, 0), required)

        for name, required in needed.items():
            check_peptide_supply(self.profile.peptide_file(name), required)

    # -- execution

    def run(self, cells: list[Cell]) -> list[CellResult]:
        results = []
        for cell in cells:
            results.append(self.run_cell(cell))
        return results

    def run_cell(self, cell: Cell) -> CellResult:
        jsonl = self.out_dir / f"{cell.label}.jsonl"
        marker = self.out_dir / f"{cell.label}.oom"
        if jsonl.exists() and jsonl.stat().st_size > 0:
            self.echo(f"  skip {cell.label} (results exist)")
            return CellResult(cell, "skipped", detail="results exist")
        if marker.exists():
            self.echo(f"  skip {cell.label} (recorded as did-not-fit)")
            return CellResult(cell, "did-not-fit", detail="recorded earlier")

        if self.suite.drop_caches and is_root():
            drop_caches()

        command = self._command(cell)
        log = self.out_dir / f"{cell.label}.log"
        self.echo(f"[{time.strftime('%H:%M:%S')}] {cell.label}")

        started = time.monotonic()
        with log.open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
        elapsed = time.monotonic() - started

        if completed.returncode == 0:
            return CellResult(cell, "ok", elapsed)

        # The marker carries the cell's dims, not just its exit code: a cell that produced no
        # records still has to appear in the report as "this arm cannot run at this ceiling", and
        # recovering that from the file name would reintroduce the parsing the dims envelope
        # exists to remove.
        marker.write_text(
            json.dumps({"exit": completed.returncode, "dims": cell.dims, "label": cell.label}) + "\n"
        )
        if completed.returncode == OOM_EXIT:
            self.echo(f"  -> killed under its ceiling; recorded as did-not-fit (see {log.name})")
            return CellResult(cell, "did-not-fit", elapsed, "OOM-killed under its ceiling")
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
            command += self._matrix_args(settings)
        else:
            command += ["--peptide-file", str(self.profile.peptide_file(settings["peptides"]))]
            command += self._single_args(cell, settings)

        if settings.get("no_theoretical_memory"):
            command.append("--no-theoretical-memory")
        return command

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
        # Every `tune.<field>` setting, passed straight through. This driver never learns the field
        # names: the harness validates them against the real `SearchTuning`, so a knob added there
        # is sweepable the same day without a change here.
        for key, value in sorted(settings.items()):
            if key.startswith(TUNE_PREFIX):
                args += ["--tune", f"{key[len(TUNE_PREFIX):]}={_tune_value(value)}"]
        return args

    def _kmer_args(self, settings: dict[str, Any]) -> list[str]:
        """Attaches the k-mer table this suite asked for, loading it or building it.

        A named table the profile does not have is built in-process instead. Running without it
        would quietly measure a different index — the table removes most of the binary search's
        probes — whereas building it only costs startup time. The warning says which happened,
        because it makes the `startup` suite's `kmer` column mean something different.
        """
        name = _table_name(settings.get("kmer_table"))
        if name is None:
            return []
        table = self.profile.kmer_table(name)
        if table:
            return ["--kmer-table-file", str(table)]
        k = _kmer_k(name)
        self.echo(
            f"  note: profile '{self.profile.name}' has no pre-built '{name}' table, so a k={k} "
            f"table is built in-process instead (same table, paid at startup)"
        )
        return ["--build-kmer-table", str(k)]

    def _matrix_args(self, settings: dict[str, Any]) -> list[str]:
        matrix = self.suite.matrix
        files = [str(self.profile.peptide_file(name)) for name in matrix.get("files", [])]
        args = ["--matrix", "--matrix-files", ",".join(files)]
        if "batches" in matrix:
            args += ["--matrix-batches", ",".join(str(batch) for batch in matrix["batches"])]
        for key, flag in (("k5", "--kmer5-file"), ("k6", "--kmer6-file")):
            if key in self.profile.kmer_tables:
                args += [flag, str(self.profile.kmer_tables[key])]
        if matrix.get("kmer6"):
            args.append("--matrix-kmer6")
        return args

    def _wrapper(self, settings: dict[str, Any]) -> list[str]:
        """`systemd-run --scope` prefix imposing this cell's memory ceiling and thread count."""
        ceiling = float(settings.get("ceiling_gb", 0) or 0)
        threads = settings.get("threads", "default")
        constrained = ceiling > 0
        pinned = str(threads) != "default"
        if not (constrained or pinned):
            return []

        wrapper = ["systemd-run", "--scope", "--quiet"]
        if pinned:
            wrapper.append(f"--setenv=RAYON_NUM_THREADS={threads}")
        if constrained:
            # Swap off, so the floor is a clean OOM rather than a measurement of the swap device.
            wrapper += ["-p", f"MemoryMax={ceiling:g}G", "-p", "MemorySwapMax=0"]
        wrapper.append("--")
        return wrapper

    # -- dry run

    def plan(self, cells: list[Cell]) -> list[str]:
        """The cell list and the exact command for each, without touching the index."""
        lines = [
            f"suite      : {self.suite.name} ({self.suite.mode} mode, {self.suite.ordering} ordering)",
            f"arms       : " + ", ".join(f"{arm.name}[{arm.feature_string or 'default'}]" for arm in self.suite.arms),
            f"cells      : {len(cells)}",
            f"results    : {self.out_dir}",
            "",
        ]
        for cell in cells:
            settings = self.settings(cell)
            reps = int(settings.get("runs", 100)) * int(settings.get("amount", 10_000))
            lines.append(f"  {cell.describe()}   ({reps:,} queries)")
        return lines


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


def _table_name(value: Any) -> str | None:
    """`kmer_table = "none"` (or absent) means run with no table attached."""
    if value in (None, "", "none"):
        return None
    return str(value)


def _kmer_k(name: str) -> int:
    """The k in a table name. Names are `k<N>` precisely so a missing file can still be rebuilt."""
    if not (name.startswith("k") and name[1:].isdigit()):
        raise ValueError(
            f"k-mer table name '{name}' must be 'k<N>' (e.g. k5, k6) so the k is recoverable when "
            f"the profile has no pre-built file"
        )
    return int(name[1:])


def _bool(value: Any) -> str:
    return "true" if str(value).lower() in ("true", "1", "yes", "on") else "false"


def _tune_value(value: Any) -> str:
    """TOML booleans reach us as Python bools, which stringify as `True` — Rust wants `true`."""
    return _bool(value) if isinstance(value, bool) else str(value)
