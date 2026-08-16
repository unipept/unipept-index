"""Both prefetch distances as one plane, plus the marginal curve of each.

The plane is the primary reading here, not a supplement: these two knobs consume the same hardware
at two points in the same search, so their curves are only separable if they do not interact, and
that is the thing being tested.

The report is built to make a NULL result readable. A suite whose honest answer is "neither knob
pays, delete both" has to say so in a way that survives being skimmed, which means the flatness has
to be stated against a floor rather than left as a grid of small percentages.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record, delta_pct
from ..report import Report, Table, band, pct
from .shared import by_cell, cell_band, floor_of, held_and_swept, knob_planes, resolution_table, tips_for

KNOBS = ("prefetch_threshold", "retrieval_prefetch_distance")


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    table, totals = _verdict(loaded, suite)
    _verdict_tiles(report, totals)

    report.heading("summary", level=3)
    held_and_swept(report, loaded)
    resolution_table(report, loaded)

    report.table(table)
    report.para(
        f"`shipped pair` and `best pair` are {KNOBS[0]}/{KNOBS[1]}. `cells` counts the pairs that "
        "produced a measurement, and the reading is about all of them at once: a plane where "
        "nothing clears its floor has not found a small effect, it has found no effect, and these "
        "two fields could then leave `SearchTuning` entirely."
    )

    report.heading("the plane", level=3)
    knob_planes(report, loaded)
    report.para(
        "Each cell is that pair of distances against the shipped pair. A cell inside its floor is "
        "painted the neutral midpoint rather than a faint tint, so a plane that looks empty IS the "
        "result: neither distance changed anything this run can resolve."
    )

    if suite.notes:
        report.note(suite.notes)


def _verdict(loaded: list[Record], suite: Suite) -> tuple[Table, list[tuple]]:
    """Does anything on the plane beat the shipped pair by more than its floor?

    One row per backend, and the whole point is the last column. A grid of sixteen numbers between
    -3% and +3% reads as sixteen results; the same grid summarised against its floor reads as the
    one result it is.

    Returns the table and `(best pair, gain, floor, resolved, total)` per arm, so the aggregate can
    lead the suite instead of closing it — a null result in particular has to survive being skimmed.
    """
    keys = ("peptide_source", *KNOBS)
    cells = by_cell(loaded, keys, correct_drift=True)
    defaults = next((record.config.get("tuning_defaults") for record in loaded if record.config), {}) or {}
    shipped = tuple(defaults.get(knob) for knob in KNOBS)
    arms = [arm.name for arm in suite.arms]

    table = Table(
        headers=["arm", "cells", "shipped pair", "best pair", "gain", "floor", "reading"],
        aligns=["<", ">", ">", ">", ">", ">", "<"],
        tips=tips_for(["arm", "gain", "floor"]),
    )
    totals: list[tuple] = []
    for arm in arms:
        points = {key[1:]: per_arm[arm] for key, per_arm in cells.items() if arm in per_arm}
        reference = points.get(shipped)
        usable = {pair: cell for pair, cell in points.items() if cell["p50"]}
        if not usable or not reference:
            continue
        best = max(usable, key=lambda pair: usable[pair]["p50"])
        difference = delta_pct(usable[best]["p50"], reference["p50"])
        floor = floor_of(usable[best], reference)
        resolved = sum(
            1
            for cell in usable.values()
            if abs(delta_pct(cell["p50"], reference["p50"])) > floor_of(cell, reference)
        )
        totals.append((best, difference, floor, resolved, len(usable), shipped))
        table.row(
            arm,
            len(usable),
            "/".join(str(value) for value in shipped),
            "/".join(str(value) for value in best),
            pct(difference),
            band(floor),
            _reading(difference, floor, resolved, len(usable)),
        )
    return table, totals


def _verdict_tiles(report: Report, totals: list[tuple]) -> None:
    """The plane's answer across every backend at once.

    This suite is the one most likely to conclude NOTHING, and a null result is the hardest kind to
    publish legibly: sixteen cells between -3% and +3% look like sixteen findings. Stated as
    "0 of 48 pairs cleared their floor" it reads as what it is — grounds for deleting two fields
    from `SearchTuning` rather than for tuning them.
    """
    if not totals:
        return
    resolved = sum(entry[3] for entry in totals)
    measured = sum(entry[4] for entry in totals)
    best = max(totals, key=lambda entry: entry[1])
    shipped = "/".join(str(value) for value in totals[0][5])

    if resolved == 0:
        status, reading = "flat", (
            f"**FLAT** — not one of the {measured} pairs measured, on any backend, beats the shipped "
            f"pair by more than its own noise floor. This is not a small effect that needs a longer "
            f"run to see; it is grounds for dropping both fields from `SearchTuning`."
        )
    else:
        status, reading = "good", (
            f"{resolved} of {measured} pairs clear their floor. Read the plane below before "
            f"changing anything: whether the ridge is straight decides if these two knobs can be "
            f"tuned separately at all."
        )
    report.verdict(
        [
            ("best pair", "/".join(str(value) for value in best[0]), f"shipped {shipped}", ""),
            ("best gain", f"{best[1]:+.1f}%", f"floor ±{best[2]:.1f}%", ""),
            ("pairs resolved", f"{resolved} / {measured}", "cleared the floor", status),
        ],
        reading,
    )


def _reading(difference: float, floor: float, resolved: int, total: int) -> str:
    if difference != difference:
        return "no data"
    if resolved == 0:
        return f"FLAT — 0 of {total} pairs clear the floor; grounds for dropping both knobs"
    if abs(difference) <= floor:
        return f"the best pair is inside its floor, though {resolved} of {total} cells are not"
    return f"the best pair wins by more than its floor; {resolved} of {total} pairs resolve"
