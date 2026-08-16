"""RAM scaling: one block per ceiling, and where the arms cross over.

The delta between arms is printed with the floor it has to clear — the wider of the two cells' own
slot spreads. Under a palindrome ordering each arm runs twice per ceiling, and the gap between an
arm's own two invocations is the honest limit on what that ceiling can resolve.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, by_residency, lines
from ..config import Suite
from ..records import Record, Summary, delta_pct, group, noise_floor, summarise, unfit_cells, verdict
from ..report import Report, Table, band, caveats, count, gb, pct, qps
from .shared import tips_for


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = _cells(loaded)
    arms = [arm.name for arm in suite.arms]
    ceilings = _ceilings(suite, cells)
    unfit = _unfit_arms(out_dir)

    # The per-ceiling table is built into a scratch report first: the crossover falls out of that
    # loop, and the crossover is this suite's answer — so it has to be known before anything is
    # emitted, and printed above the table it came from.
    body = Report()
    body.heading("per ceiling", level=3, folded=True)
    table = Table(
        headers=["ceiling", "arm", "n", "qps", "band", "slots", "drift", "majflt/rep", "RSS GB"],
        aligns=["<", "<", ">", ">", ">", ">", ">", ">", ">"],
        tips=tips_for(["arm", "qps", "band", "drift"]),
    )
    crossover: list[tuple[str, float, float, str]] = []

    for ceiling in ceilings:
        name = "none" if ceiling == 0 else f"{ceiling:g}G"
        present: dict[str, Summary] = {}
        for arm in arms:
            summary = cells.get((ceiling, arm))
            if summary is None:
                table.row(name, arm, "-", "did not fit" if arm in unfit.get(ceiling, set()) else "(not run)", *[""] * 5)
                continue
            present[arm] = summary
            table.row(
                name,
                arm,
                summary.n,
                "VOID" if summary.void_reason else qps(summary.qps),
                band(summary.band),
                band(summary.slot_spread),
                pct(summary.drift),
                count(summary.major_faults),
                gb(summary.rss_gb),
            )

        if len(arms) == 2 and all(arm in present and present[arm].usable for arm in arms):
            base, other = present[arms[0]], present[arms[1]]
            difference = delta_pct(other.qps, base.qps)
            floor = noise_floor(base, other)
            call = verdict(difference, floor, better=arms[1], worse=arms[0])
            crossover.append((name, difference, floor, call))
            table.row("", f"-> {arms[1]} vs {arms[0]}", "", pct(difference), "", band(floor), "", "", call)

    body.table(table, raw=True)

    if crossover:
        body.heading("crossover", level=3)
        readout = []
        for index, (name, difference, floor, call) in enumerate(crossover):
            mark = ""
            if index and (crossover[index - 1][1] > 0) != (difference > 0):
                mark = "   <-- SIGN CHANGE from the previous ceiling"
            readout.append(f"  {name:>6s}: {difference:+6.1f}%  (floor {floor:.1f}% -> {call}){mark}")
        body.lines(readout)

    _verdict_tiles(report, crossover, arms)

    report.heading("summary", level=3)
    # Two charts, never two y-axes on one plot: throughput and fault counts are different scales,
    # and one plot would invent a relationship between them. They are stacked instead, sharing an
    # x axis, which is where the mechanism actually reads — the ceiling at which the lines cross is
    # the ceiling at which faults take over.
    _curves(report, cells, arms, ceilings)
    report.extend(body)

    notes = caveats(list(cells.values()))
    if notes:
        report.heading("caveats", level=3, folded=True).lines([f"  * {note}" for note in notes])
    if suite.notes:
        report.note(suite.notes)


def _curves(report: Report, cells: dict, arms: list[str], ceilings: list[float]) -> None:
    """Throughput and major faults against the ceiling, as two charts sharing an x axis.

    Ceilings run from loosest to tightest left to right, so the x axis is "increasing pressure" and
    the crossover reads as a crossing rather than as something to look up in a table. A void cell
    contributes no point: it measured nothing, so it must not draw a line through anywhere.
    """
    present = [ceiling for ceiling in ceilings if any((ceiling, arm) in cells for arm in arms)]
    if len(present) < 2:
        return
    labels = ["none" if ceiling == 0 else f"{ceiling:g}G" for ceiling in present]

    def values(arm: str, field: str) -> list[float | None]:
        out = []
        for ceiling in present:
            summary = cells.get((ceiling, arm))
            out.append(getattr(summary, field) if summary and summary.usable else None)
        return out

    for field, caption, unit, y_title in (
        ("qps", "Throughput as the memory ceiling falls", " qps", "throughput (qps)"),
        ("major_faults", "Major faults per rep as the memory ceiling falls", "", "major faults per rep"),
    ):
        report.chart(
            lines(
                labels,
                # Arms are one measurement across storage configurations, so they take the ordinal
                # arm ramp rather than a hue each — the report-wide rule, see `charts.arm_color`.
                [Series(arm, values(arm, field), arm=arm, tip={"arm": arm}) for arm in by_residency(arms)],
                caption,
                unit=unit,
                x_title="memory ceiling",
                y_title=y_title,
            ),
            caption,
        )


def _cells(loaded: list[Record]) -> dict[tuple[float, str], Summary]:
    cells = {}
    for key, group_records in group(loaded).items():
        dims = dict(key)
        ceiling = float(dims.get("ceiling_gb", 0) or 0)
        cells[(ceiling, dims.get("arm", "?"))] = summarise(group_records, ceiling or None)
    return cells


def _ceilings(suite: Suite, cells: dict) -> list[float]:
    """Ceilings in the order the suite declares them, plus any that only the results know about."""
    declared = [float(value) for value in suite.axes.get("ceiling_gb", [0])]
    seen = sorted({ceiling for ceiling, _ in cells}, reverse=True)
    return declared + [value for value in seen if value not in declared]


def _unfit_arms(out_dir: Path) -> dict[float, set[str]]:
    """Which (ceiling, arm) pairs were OOM-killed, from the dims in the runner's markers."""
    unfit: dict[float, set[str]] = {}
    for dims in unfit_cells(out_dir):
        ceiling = float(dims.get("ceiling_gb", 0) or 0)
        unfit.setdefault(ceiling, set()).add(dims.get("arm", "?"))
    return unfit


def _verdict_tiles(report: Report, crossover: list[tuple], arms: list[str]) -> None:
    """The ceiling at which the arms swap places — the one number this suite exists to find.

    A crossover is a SIGN CHANGE that also clears its floor on both sides. A sign change inside the
    floor is two unresolved cells in a row, not a crossing, and reporting it as one would hand an
    operator a memory budget the run never actually measured.
    """
    if not crossover or len(arms) < 2:
        return
    crossing = next(
        (
            crossover[index][0]
            for index in range(1, len(crossover))
            if (crossover[index - 1][1] > 0) != (crossover[index][1] > 0)
            and abs(crossover[index][1]) > crossover[index][2]
        ),
        None,
    )
    widest = max(crossover, key=lambda entry: abs(entry[1]))
    leader = arms[1] if widest[1] > 0 else arms[0]

    if crossing:
        status, reading = "good", (
            f"`{arms[1]}` and `{arms[0]}` swap places at the **{crossing}** ceiling, and the swap "
            f"clears its floor on both sides. That is the memory budget the deployment choice turns "
            f"on: above it the arms are interchangeable, below it they are not."
        )
    else:
        status, reading = "flat", (
            f"No resolved crossover: `{leader}` leads at every ceiling measured, or the sign changes "
            f"only inside the floor. The arms do not swap places anywhere in this range — extend the "
            f"ceilings downward before concluding they never do."
        )
    report.verdict(
        [
            ("crossover at", crossing or "none", "memory ceiling", status),
            ("widest gap", f"{widest[1]:+.1f}%", f"at {widest[0]}, floor ±{widest[2]:.1f}%", ""),
            ("leader there", leader, "of the two mapped arms", ""),
        ],
        reading,
    )
