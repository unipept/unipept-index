"""Expanding `[[sweep]]` blocks into the cell list one process executes.

A sweep of N knobs across M contexts has two possible shapes, and the choice between them decides
whether the run takes an hour or a day.

**Multiplicative** — every knob value crossed with every context — asks "does every knob's optimum
depend on every context?". Four knobs at their measured value counts, three k-mer tables, four
thread counts, both search options and both backends is 3,456 cells, nearly all of them
re-measuring a curve that was flat the first three times.

**Additive** — several blocks, each varying one thing against a fixed background — asks the three
questions anyone actually has, and adds their costs instead of multiplying them:

    what is each knob's curve?          every knob,  ONE context     ->  ~36 cells
    what does each context cost?        ONE tuning point, every one  ->  ~96 cells
    does a knob's optimum MOVE?         only the pairs with a mechanism behind them

so a suite is a list of blocks rather than one grid. A block that turns out to matter gets widened;
a cross-product cannot be narrowed after the fact, because the run is already over.

Everything a block does not name is inherited: contexts from the suite's `[axes]`, precision from
its `[defaults]`, tuning from the knob values this binary reports as shipped. That inheritance is
what keeps a block down to the two or three lines that say what it is *for*.

## Strategies

    base    the shipped defaults, one point. For blocks whose subject is the context.
    ofat    one knob moved at a time, every other knob at its shipped value. `1 + sum(len-1)`
            points, so knobs add rather than multiply.
    pairs   ofat, plus the full 2-D product of each named knob pair. A plane already contains both
            ofat lines through the default point, so the lines it subsumes are dropped rather than
            re-run — a pair costs |A|x|B| instead of |A|x|B| + (|A|-1) + (|B|-1).
    full    the complete cross-product of every knob in the block. Honest for two or three values
            of two knobs, ruinous beyond that.

Deduplication is on the fully-resolved cell, not on bookkeeping about which strategy emitted it.
Two blocks that happen to describe the same measurement therefore collapse into one, whichever
routes led there — which is what lets a block be added without auditing every other one for
overlap.
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

#: Everything else a block may set, beyond its contexts and its `tune` table.
SETTING_KEYS = ("name", "strategy", "runs", "amount", "base_every", "pairs", "response")

STRATEGIES = ("base", "ofat", "pairs", "full")


class GridError(Exception):
    """A `[[sweep]]` block is malformed, or asks for something that cannot be expanded."""


# ---------------------------------------------------------------------------
# Tuning points
# ---------------------------------------------------------------------------


def tuning_points(block: dict, defaults: dict[str, Any]) -> list[dict[str, Any]]:
    """The `SearchTuning` overrides this block measures, as a list of `{field: value}` dicts.

    Every point is complete — it names every knob, at its swept or shipped value — so two points
    are equal exactly when they describe the same measurement, and dedup needs no knowledge of how
    either was produced.
    """
    knobs = block.get("tune") or {}
    strategy = block.get("strategy", "base" if not knobs else "ofat")
    if strategy not in STRATEGIES:
        raise GridError(
            f"sweep '{block.get('name', '?')}': strategy must be one of {STRATEGIES}, got '{strategy}'"
        )
    if strategy != "base" and not knobs:
        raise GridError(f"sweep '{block.get('name', '?')}': strategy '{strategy}' needs a [sweep.tune] table")

    base = dict(defaults)
    if strategy == "base":
        return [base]

    if strategy == "full":
        names = sorted(knobs)
        return [{**base, **dict(zip(names, values))} for values in product(*(knobs[name] for name in names))]

    # ofat / pairs. Planes first, so the lines they subsume are already in `seen` and get skipped
    # rather than being emitted and then deduplicated away — the difference is visible in the cell
    # ORDER, and cell order is what the drift cadence interleaves against.
    points: list[dict[str, Any]] = [base]
    seen = {_key(base)}

    for pair in block.get("pairs") or []:
        axes = pair.get("axes") if isinstance(pair, dict) else pair
        if not isinstance(axes, list) or len(axes) != 2:
            raise GridError(
                f"sweep '{block.get('name', '?')}': a pair must name exactly two knobs, got {axes!r}. "
                f"A plane is two-dimensional; three knobs at once is `strategy = \"full\"`."
            )
        for axis in axes:
            if axis not in knobs:
                raise GridError(
                    f"sweep '{block.get('name', '?')}': pair axis '{axis}' has no values in "
                    f"[sweep.tune] (has: {', '.join(sorted(knobs)) or 'nothing'})"
                )
        for left, right in product(knobs[axes[0]], knobs[axes[1]]):
            point = {**base, axes[0]: left, axes[1]: right}
            if _key(point) not in seen:
                seen.add(_key(point))
                points.append(point)

    for name in sorted(knobs):
        for value in knobs[name]:
            point = {**base, name: value}
            if _key(point) not in seen:
                seen.add(_key(point))
                points.append(point)

    return points


def _key(point: dict[str, Any]) -> tuple:
    return tuple(sorted(point.items()))


# ---------------------------------------------------------------------------
# Contexts
# ---------------------------------------------------------------------------


def contexts(block: dict, suite_axes: dict[str, list], suite_files: list[str]) -> list[dict[str, Any]]:
    """The (process coordinate, cell coordinate) points this block runs its tuning points at.

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
        # The match cutoff. `None` leaves the invocation's `--max-matches` alone. Unlike a tuning
        # knob this one changes the ANSWER, so it is a coordinate rather than a setting.
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
    defaults: dict[str, Any],
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
        points = tuning_points(block, defaults)
        runs = block.get("runs", suite_defaults.get("runs"))
        amount = block.get("amount", suite_defaults.get("amount"))

        for context in contexts(block, suite_axes, suite_files):
            process = (context["arms"], context["threads"], context["ceiling_gb"])
            cells = processes.setdefault(process, [])
            known = seen.setdefault(process, set())
            for point in points:
                cell = _cell(block, context, point, runs, amount)
                # Identity is the measurement, not the block that asked for it. `sweep` and
                # `grid_slot` are excluded so two blocks describing the same cell collapse rather
                # than running it twice at (possibly) different precision.
                identity = _identity(cell)
                if identity in known:
                    continue
                known.add(identity)
                cells.append(cell)

    for process, cells in processes.items():
        processes[process] = _with_drift_cadence(cells, sweeps, defaults)
    return processes


def _cell(block: dict, context: dict, point: dict, runs, amount) -> dict:
    amount = context["amounts"] if context.get("amounts") is not None else amount
    cell = {
        "file": context["files"],
        "kmer_k": int(context["kmer"]),
        "equate_il": bool(context["equate_il"]),
        "tryptic": bool(context["tryptic"]),
        "tune": dict(point),
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
        _key(cell["tune"]),
        cell.get("runs"),
        cell.get("amount"),
        cell.get("max_matches"),
    )


def _validate(block: dict) -> None:
    unknown = set(block) - set(CONTEXT_KEYS) - set(SETTING_KEYS) - {"tune"}
    if unknown:
        raise GridError(
            f"sweep '{block.get('name', '?')}': unknown key(s) {', '.join(sorted(unknown))}. A key "
            f"with no defined effect would narrow nothing while still reading as though it did.\n"
            f"  contexts: {', '.join(CONTEXT_KEYS)}\n"
            f"  settings: {', '.join(SETTING_KEYS)}, tune"
        )
    for key in CONTEXT_KEYS:
        if key in block and not isinstance(block[key], list):
            raise GridError(f"sweep '{block.get('name', '?')}': '{key}' must be a list, got {block[key]!r}")


# ---------------------------------------------------------------------------
# Drift
# ---------------------------------------------------------------------------


def _with_drift_cadence(cells: list[dict], sweeps: list[dict], defaults: dict[str, Any]) -> list[dict]:
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
    reference["tune"] = dict(defaults)
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
    a fifth of the cost, and counting configurations hides exactly the thing that was tuned.
    """
    return sum(cell.get("runs", runs) * cell.get("amount", amount) for cell in cells)
