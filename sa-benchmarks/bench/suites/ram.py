"""RAM scaling: one block per ceiling, and where the arms cross over.

The delta between arms is printed with the floor it has to clear — the wider of the two cells' own
slot spreads. Under a palindrome ordering each arm runs twice per ceiling, and the gap between an
arm's own two invocations is the honest limit on what that ceiling can resolve.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, lines
from ..config import Suite
from ..records import Record, Summary, delta_pct, group, noise_floor, summarise, unfit_cells, verdict
from ..report import Report, Table, band, caveats, count, gb, pct, qps


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = _cells(loaded)
    arms = [arm.name for arm in suite.arms]
    ceilings = _ceilings(suite, cells)
    unfit = _unfit_arms(out_dir)

    report.heading("summary", level=3)
    # Two charts, never two y-axes on one plot: throughput and fault counts are different scales,
    # and one plot would invent a relationship between them. They are stacked instead, sharing an
    # x axis, which is where the mechanism actually reads — the ceiling at which the lines cross is
    # the ceiling at which faults take over.
    _curves(report, cells, arms, ceilings)

    report.heading("per ceiling", level=3)
    table = Table(
        headers=["ceiling", "arm", "n", "qps", "band", "slots", "drift", "majflt/rep", "RSS GB"],
        aligns=["<", "<", ">", ">", ">", ">", ">", ">", ">"],
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

    report.table(table)

    if crossover:
        report.heading("crossover", level=3)
        lines = []
        for index, (name, difference, floor, call) in enumerate(crossover):
            mark = ""
            if index and (crossover[index - 1][1] > 0) != (difference > 0):
                mark = "   <-- SIGN CHANGE from the previous ceiling"
            lines.append(f"  {name:>6s}: {difference:+6.1f}%  (floor {floor:.1f}% -> {call}){mark}")
        report.lines(lines)

    notes = caveats(list(cells.values()))
    if notes:
        report.heading("caveats", level=3).lines([f"  * {note}" for note in notes])
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

    for field, caption, unit in (
        ("qps", "Throughput as the memory ceiling falls", " qps"),
        ("major_faults", "Major faults per rep as the memory ceiling falls", ""),
    ):
        report.chart(
            lines(
                labels,
                [Series(arm, values(arm, field), slot) for slot, arm in enumerate(arms)],
                caption,
                unit=unit,
                x_title="memory ceiling",
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
