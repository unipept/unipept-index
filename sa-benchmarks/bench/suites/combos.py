"""The three accelerators crossed, and the one question that needs the cross.

`kmer`, `mlp` and `validate` each report one knob's curve against a fixed background. Stacking their
three winners into a configuration is only valid if the knobs are separable, and that is what this
suite tests: does the best tuple move across contexts, or is there one that wins everywhere?

The output is therefore a table of winners, not a wall of heatmaps. A full cross of three axes over
four search-option pairs and three backends is thirty-six planes, which is not a reading — so the
planes are shown only for the contexts where something actually cleared its floor, and the rest are
counted rather than drawn.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record, delta_pct
from ..report import Report, Table, band, pct, qps
from .shared import (
    by_cell,
    fmt_tune,
    floor_of,
    held_and_swept,
    kmer_label,
    knob_planes,
    resolution_table,
    tips_for,
)

#: The three axes, in the order they read as a tuple.
AXES = ("kmer_k", "mlp_batch", "validate_batch")

#: How many planes are worth drawing. Past a handful they stop being read.
MAX_PLANES = 6


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    # Built into a scratch report so the aggregate — which is this suite's whole point — can be
    # stated above the table it comes from rather than in a paragraph below it.
    body = Report()
    winners, decided, rows = _winners(body, loaded, suite)
    _verdict_tiles(report, winners, decided, rows)

    report.heading("summary", level=3)
    held_and_swept(report, loaded)
    resolution_table(report, loaded)
    report.extend(body)
    if not winners:
        return

    report.heading("the planes", level=3)
    _planes(report, loaded, winners)

    if suite.notes:
        report.note(suite.notes)


def _verdict_tiles(report: Report, winners: list[tuple], decided: set[tuple], rows: int) -> None:
    """Whether the three accelerators' individual optima are simultaneously reachable.

    That is the only reason to pay for a multiplicative sweep, so it is the only thing the tiles say.
    A stable winner licenses stacking the three single-knob suites' answers; a winner that moves
    means the tuple has to be chosen together, and the single-knob suites are each one slice.
    """
    if not rows:
        return
    if not winners:
        status, value, reading = "flat", "none", (
            f"**No context resolved.** Across {rows} of them, no tuple beat the shipped one by more "
            f"than that context's own floor — which is itself the separability result: crossing the "
            f"three knobs found nothing their individual suites would have missed."
        )
    elif len(decided) == 1:
        status, value = "good", _tuple_label(next(iter(decided)))
        reading = (
            f"**One tuple wins wherever anything resolves**, in {len(winners)} of {rows} contexts. "
            f"A stable winner across k-mer tables, search options and backends is a candidate "
            f"configuration rather than a deployment one."
        )
    else:
        status, value = "warn", f"{len(decided)} tuples"
        reading = (
            f"**The winning tuple MOVES**: {len(decided)} different tuples win across "
            f"{len(winners)} of {rows} contexts. The three accelerators' individual optima are not "
            f"simultaneously reachable, so the tuple is a deployment choice, not a default."
        )
    report.verdict(
        [
            ("best tuple", value, "/".join(AXES), ""),
            ("contexts resolved", f"{len(winners)} / {rows}", "cleared the floor", status),
        ],
        reading,
    )


def _winners(report: Report, loaded: list[Record], suite: Suite) -> tuple[list[tuple], set[tuple], int]:
    """One row per context: the best tuple, and whether it beats the shipped one past the floor.

    Returns the contexts that resolved (so the plane section can draw those and skip the rest), the
    distinct winning tuples, and how many contexts were measured at all.
    """
    keys = ("peptide_source", "equate_il", "tryptic", "amount_of_peptides", *AXES)
    cells = by_cell(loaded, keys, correct_drift=True)
    defaults = next((r.config.get("tuning_defaults") for r in loaded if r.config), {}) or {}
    arms = [arm.name for arm in suite.arms]

    # The shipped tuple, with the k-mer table the other suites treat as production. Read off the
    # binary for the two knobs; the table is not a `SearchTuning` field, so it is named here.
    shipped_kmer = 5 if any(key[keys.index("kmer_k")] == 5 for key in cells) else None
    shipped = (shipped_kmer, defaults.get("mlp_batch"), defaults.get("validate_batch"))

    # One column per swept coordinate, so each is its own filter — a joined `context` string is a
    # lookup key, not something a chip can narrow.
    headers = ["file", "equate_il", "tryptic", "arm", "shipped tuple", "best tuple", "gain", "floor", "reading"]
    table = Table(
        headers=headers,
        aligns=["<", "<", "<", "<", ">", ">", ">", ">", "<"],
        chips=["file", "equate_il", "tryptic", "arm"],
        tips=tips_for(headers),
    )
    resolved_contexts: list[tuple] = []
    decided: set[tuple] = set()
    rows = 0

    grouped: dict[tuple, dict[str, dict]] = {}
    for key, per_arm in cells.items():
        context = (key[0], key[1], key[2], key[3])
        tuple_key = tuple(key[keys.index(axis)] for axis in AXES)
        for arm, cell in per_arm.items():
            grouped.setdefault((context, arm), {})[tuple_key] = cell

    for (context, arm), points in sorted(grouped.items(), key=str):
        if arm not in arms:
            continue
        usable = {values: cell for values, cell in points.items() if cell["p50"]}
        reference = usable.get(shipped)
        if not usable or not reference:
            continue
        best = max(usable, key=lambda values: usable[values]["p50"])
        difference = delta_pct(usable[best]["p50"], reference["p50"])
        floor = floor_of(usable[best], reference)
        clears = abs(difference) > floor
        if clears:
            resolved_contexts.append((context, arm, best, difference))
            decided.add(best)
        rows += 1
        table.row(
            context[0],
            # Bare booleans, not `context_label`: that helper names the coordinate because it writes
            # into a joined context string where nothing else would, and here the column is already
            # called `equate_il`. `tryptic=true` under a header reading `tryptic` says it twice, and
            # it is not one of the spellings the page paints as a pill.
            fmt_tune(context[1]),
            fmt_tune(context[2]),
            arm,
            _tuple_label(shipped),
            _tuple_label(best),
            pct(difference),
            band(floor),
            "the shipped tuple is the peak"
            if best == shipped
            else (f"{_tuple_label(best)} wins" if clears else "inside the floor — no better tuple shown"),
        )
    report.table(table)

    # The reading of the whole table is the verdict row above it — see `_verdict_tiles`. What stays
    # here is only what the COLUMNS mean, which the tiles have no room for.
    report.para(
        f"Tuples are {'/'.join(AXES)}. `gain` is the best tuple against the shipped one in the same "
        "context, and it only means something past the floor beside it."
    )
    return resolved_contexts, decided, rows


def _planes(report: Report, loaded: list[Record], winners: list[tuple]) -> None:
    """The `mlp_batch` x `validate_batch` planes, drawn only where something resolved.

    A full cross over every context is dozens of heatmaps, and a reader who has to scroll past
    thirty of them to find the interesting one has not been given a picture, only a pile.
    """
    if not winners:
        report.para(
            "No plane is drawn: nothing resolved, so every one of them would be a grid of neutral "
            "cells. The winners table above is the whole result."
        )
        return

    # Ranked by gain, but budgeted per search mode rather than over the pile. The page filters
    # figures by `tryptic`, and one global top-six can be entirely tryptic — which would leave the
    # reader looking at an empty section under the other half of the filter, with the prose below it
    # still describing planes. Each mode gets its own share of the budget, so neither view is blank.
    ranked = sorted(winners, key=lambda entry: -abs(entry[3]))
    budget = max(1, MAX_PLANES // 2)
    taken: dict[object, int] = {}
    interesting = []
    for entry in ranked:
        mode = entry[0][2]
        if taken.get(mode, 0) >= budget:
            continue
        taken[mode] = taken.get(mode, 0) + 1
        interesting.append(entry)
    # A run that swept only one mode should still get the whole budget rather than half of it.
    if len(taken) == 1:
        interesting = ranked[:MAX_PLANES]
    wanted = {(context, arm) for context, arm, _, _ in interesting}
    kept = [
        record
        for record in loaded
        if (
            (
                record.config.get("peptide_source"),
                record.config.get("equate_il"),
                record.config.get("tryptic"),
                record.config.get("amount_of_peptides"),
            ),
            record.dims.get("arm"),
        )
        in wanted
        or record.config.get("sweep") == "drift"
    ]
    knob_planes(report, kept)
    omitted = len(winners) - len(interesting)
    report.para(
        f"Drawn for the {len(interesting)} context(s) with the largest resolved gain"
        + (f"; {omitted} other resolved context(s) are in the table above but not plotted." if omitted else ".")
        + " A cell inside its floor is painted the neutral midpoint, so a plane that looks empty is "
        "one where the two batches did nothing this run can separate."
    )


def _tuple_label(values: tuple) -> str:
    kmer, mlp, validate = values
    return f"{kmer_label(kmer)}/{fmt_tune(mlp)}/{fmt_tune(validate)}"
