"""Expanding `[[sweep]]` blocks into the cell list one process executes.

A suite is a list of blocks, each varying one thing against a fixed background, rather than one
grid crossing everything with everything. Their costs then add instead of multiplying: a block that
turns out to matter gets widened, where a cross-product cannot be narrowed after the fact, because
the run is already over.

Everything a block does not name is inherited: contexts from the suite's `[axes]`, precision from
its `[defaults]`. That inheritance is what keeps a block down to the two or three lines that say
what it is *for*.

## What a block can vary

The contexts in `CONTEXT_KEYS` — the workload, the build, the machine — and nothing else. Blocks
used to carry a `[sweep.tune]` table and a `strategy` (`ofat`, `pairs`, `full`) describing how to
walk it, because the searcher had runtime performance knobs. It does not any more: those became
compile-time constants once no sweep could separate their values from noise, so every block is now
the one shape that used to be called `base` — one measurement per context.

Deduplication is on the fully-resolved cell, not on bookkeeping about which block emitted it. Two
blocks that happen to describe the same measurement therefore collapse into one, which is what lets
a block be added without auditing every other one for overlap.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any, Iterable

#: Context keys a block may set. Each is a list, and the block runs the product of those it names.
#: `arms` and `threads` and `ceiling_gb` are *process* coordinates — they decide which process a
#: cell belongs to rather than travelling in the grid file — because a rayon pool is built once per
#: process and a cgroup scope wraps one.
PROCESS_KEYS = ("arms", "threads", "ceiling_gb")
CELL_KEYS = ("files", "kmer", "equate_il", "tryptic", "amounts", "cutoffs")
CONTEXT_KEYS = PROCESS_KEYS + CELL_KEYS

#: Everything else a block may set, beyond its contexts.
SETTING_KEYS = ("name", "runs", "amount", "base_every", "response")


class GridError(Exception):
    """A `[[sweep]]` block is malformed, or asks for something that cannot be expanded."""


# ---------------------------------------------------------------------------
# Contexts
# ---------------------------------------------------------------------------


def contexts(block: dict, suite_axes: dict[str, list], suite_files: list[str]) -> list[dict[str, Any]]:
    """The (process coordinate, cell coordinate) points this block measures at.

    A key the block does not name falls back to the suite's `[axes]`, then to a single sensible
    value. Inheriting rather than requiring each block to restate every context is what keeps a
    block readable as the one thing it varies.
    """
    values = {
        "arms": block.get("arms"),
        "threads": block.get("threads", suite_axes.get("threads")) or ["default"],
        "ceiling_gb": block.get("ceiling_gb", suite_axes.get("ceiling_gb")) or [0],
        "files": block.get("files", suite_files) or [],
        "kmer": block.get("kmer", [5]),
        "equate_il": block.get("equate_il", [True]),
        "tryptic": block.get("tryptic", [False]),
        # `amounts` is the AXIS; the scalar `amount` setting is the value when nothing sweeps it.
        # Query count is a coordinate, not a precision dial: two cells that ran different stream
        # lengths are not comparable, so a suite that varies it is measuring it.
        "amounts": block.get("amounts", [None]),
        # The match cutoff. `None` leaves the invocation's `--max-matches` alone. This one changes
        # the ANSWER, not just the time, so it is a coordinate rather than a setting.
        "cutoffs": block.get("cutoffs", [None]),
    }
    if not values["arms"]:
        raise GridError(f"sweep '{block.get('name', '?')}': no arms — name them, or the block runs nothing")
    if not values["files"]:
        raise GridError(
            f"sweep '{block.get('name', '?')}': no peptide files — set `files` on the block or on the suite"
        )

    names = list(values)
    return [dict(zip(names, point)) for point in product(*(values[name] for name in names))]


# ---------------------------------------------------------------------------
# Expansion
# ---------------------------------------------------------------------------


def expand(
    sweeps: list[dict],
    *,
    suite_axes: dict[str, list] | None = None,
    suite_files: list[str] | None = None,
    suite_defaults: dict[str, Any] | None = None,
) -> dict[tuple, list[dict]]:
    """Every block, expanded and grouped by the process that must run it.

    Returns `(arm, threads, ceiling_gb) -> [grid cell]`, each cell in the JSON shape the harness's
    `--grid-file` reads. Cells keep the order their blocks were declared in, so a suite controls
    what runs early and what runs late — which matters, because a process that drifts drifts in
    cell order.
    """
    suite_axes = suite_axes or {}
    suite_files = suite_files or []
    suite_defaults = suite_defaults or {}

    processes: dict[tuple, list[dict]] = {}
    seen: dict[tuple, set[tuple]] = {}

    for block in sweeps:
        _validate(block)
        runs = block.get("runs", suite_defaults.get("runs"))
        amount = block.get("amount", suite_defaults.get("amount"))

        for context in contexts(block, suite_axes, suite_files):
            process = (context["arms"], context["threads"], context["ceiling_gb"])
            cells = processes.setdefault(process, [])
            known = seen.setdefault(process, set())
            cell = _cell(block, context, runs, amount)
            # Identity is the measurement, not the block that asked for it. `sweep` and
            # `grid_slot` are excluded so two blocks describing the same cell collapse rather
            # than running it twice at (possibly) different precision.
            identity = _identity(cell)
            if identity in known:
                continue
            known.add(identity)
            cells.append(cell)

    for process, cells in processes.items():
        processes[process] = _with_drift_cadence(cells, sweeps)
    return processes


def _cell(block: dict, context: dict, runs, amount) -> dict:
    amount = context["amounts"] if context.get("amounts") is not None else amount
    cell = {
        "file": context["files"],
        "kmer_k": int(context["kmer"]),
        "equate_il": bool(context["equate_il"]),
        "tryptic": bool(context["tryptic"]),
        "sweep": block.get("name", ""),
        "grid_slot": "a",
    }
    if context.get("cutoffs") is not None:
        cell["max_matches"] = int(context["cutoffs"])
    if block.get("response"):
        cell["response"] = True
    if runs is not None:
        cell["runs"] = int(runs)
    if amount is not None:
        cell["amount"] = int(amount)
    return cell


def _identity(cell: dict) -> tuple:
    return (
        cell["file"],
        cell["kmer_k"],
        cell["equate_il"],
        cell["tryptic"],
        cell.get("runs"),
        cell.get("amount"),
        cell.get("max_matches"),
    )


def _validate(block: dict) -> None:
    unknown = set(block) - set(CONTEXT_KEYS) - set(SETTING_KEYS)
    if unknown:
        retired = {"tune", "strategy", "pairs"} & unknown
        hint = (
            f"\n  {', '.join(sorted(retired))} went with `SearchTuning`: the searcher's performance "
            f"parameters are compile-time constants, so there is nothing left to sweep."
            if retired
            else ""
        )
        raise GridError(
            f"sweep '{block.get('name', '?')}': unknown key(s) {', '.join(sorted(unknown))}. A key "
            f"with no defined effect would narrow nothing while still reading as though it did.\n"
            f"  contexts: {', '.join(CONTEXT_KEYS)}\n"
            f"  settings: {', '.join(SETTING_KEYS)}{hint}"
        )
    for key in CONTEXT_KEYS:
        if key in block and not isinstance(block[key], list):
            raise GridError(f"sweep '{block.get('name', '?')}': '{key}' must be a list, got {block[key]!r}")


# ---------------------------------------------------------------------------
# Drift
# ---------------------------------------------------------------------------


def _with_drift_cadence(cells: list[dict], sweeps: list[dict]) -> list[dict]:
    """Interleaves a reference cell every `base_every` cells, in slots "z0", "z1", ...

    A process running a hundred cells takes long enough that machine drift across it can exceed the
    knob effects it is measuring, and a single reference at the start cannot tell the two apart. Run
    the reference repeatedly and drift becomes a measured series: the analysis subtracts its trend
    and reads what is left of it as this process's resolution floor.

    Placed here rather than by the caller because it has to happen after dedup — the whole point is
    that these repeats are NOT deduplicated away, which is what the distinct slots encode.
    """
    cadence = next((block.get("base_every") for block in sweeps if block.get("base_every")), None)
    if not cadence or not cells:
        return cells

    reference = dict(cells[0])
    reference["sweep"] = "drift"

    out: list[dict] = []
    marks = 0
    for index, cell in enumerate(cells):
        if index % cadence == 0:
            out.append({**reference, "grid_slot": f"z{marks}"})
            marks += 1
        out.append(cell)
    # Closing mark, so the last stretch of cells is bracketed rather than extrapolated past.
    out.append({**reference, "grid_slot": f"z{marks}"})
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write(cells: Iterable[dict], path: Path) -> Path:
    """Writes a grid file. One JSON object per line, in run order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(cell, sort_keys=True) + "\n" for cell in cells))
    return path


def query_count(cells: Iterable[dict], runs: int, amount: int) -> int:
    """Timed queries these cells will run, which is what the wall clock tracks.

    Not the cell count: once cells may differ in size, a tryptic cell at a fifth the query count is
    a fifth of the cost, and counting configurations hides exactly that difference.
    """
    return sum(cell.get("runs", runs) * cell.get("amount", amount) for cell in cells)
