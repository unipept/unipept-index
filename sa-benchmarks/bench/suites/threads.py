"""Thread sweep: one block per ceiling, each arm's curve against its own default-thread baseline.

Two questions are answered side by side, and they are not the same question. Down a column: does
oversubscription pay at this ceiling, and where does it peak. Across the arms: does the storage
choice still matter once faults are overlapped.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, by_residency, lines
from ..config import Suite
from ..records import Record, Summary, delta_pct, group, noise_floor, summarise, unfit_cells
from ..report import Report, Table, band, caveats, count, gb, pct, qps

DEFAULT_THREADS = "default"


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = _cells(loaded)
    arms = [arm.name for arm in suite.arms]
    ceilings = [str(value) for value in suite.axes.get("ceiling_gb", [0])]
    thread_counts = [str(value) for value in suite.axes.get("threads", [DEFAULT_THREADS])]
    unfit = {(dims.get("ceiling_gb", "0"), dims.get("threads", "?"), dims.get("arm", "?")) for dims in unfit_cells(out_dir)}

    for ceiling in ceilings:
        title = "unconstrained" if ceiling == "0" else f"{ceiling}G ceiling"
        report.heading(title, level=3)

        # Small multiples: one chart per ceiling, same series and same form each time, so the shape
        # of the curve can be compared across ceilings by eye. That comparison — does the optimum
        # move as the fault rate rises? — is the suite's whole question.
        caption = f"{title}: throughput against thread count"
        report.chart(
            lines(
                thread_counts,
                [
                    Series(
                        arm,
                        [
                            cells[(ceiling, threads, arm)].qps
                            if (ceiling, threads, arm) in cells and cells[(ceiling, threads, arm)].usable
                            else None
                            for threads in thread_counts
                        ],
                        arm=arm,
                        tip={"arm": arm},
                    )
                    for arm in by_residency(arms)
                ],
                caption,
                unit=" qps",
                y_title="throughput (qps)",
                x_title="RAYON_NUM_THREADS",
            ),
            caption,
        )

        headers = ["threads"]
        for arm in arms:
            headers += [f"{arm} qps", "vs dflt"]
        if len(arms) == 2:
            headers.append(f"{arms[1]} vs {arms[0]}")
        headers += [f"majflt {arm}" for arm in arms] + ["RSS GB"]
        table = Table(headers=headers, aligns=["<"] + [">"] * (len(headers) - 1))

        baselines = {arm: cells.get((ceiling, DEFAULT_THREADS, arm)) for arm in arms}
        best: dict[str, tuple[str, float]] = {}

        for threads in thread_counts:
            row = [threads]
            here: dict[str, Summary | None] = {}
            for arm in arms:
                summary = here[arm] = cells.get((ceiling, threads, arm))
                if summary is None:
                    row += ["did not fit" if (ceiling, threads, arm) in unfit else "-", "-"]
                    continue
                row.append("VOID" if summary.void_reason else qps(summary.qps))
                baseline = baselines.get(arm)
                if threads == DEFAULT_THREADS or not (baseline and baseline.usable and summary.usable):
                    row.append("base" if threads == DEFAULT_THREADS else "-")
                else:
                    row.append(pct(delta_pct(summary.qps, baseline.qps)))
                if summary.usable and (arm not in best or summary.qps > best[arm][1]):
                    best[arm] = (threads, summary.qps)

            if len(arms) == 2:
                left, right = here[arms[0]], here[arms[1]]
                if left and right and left.usable and right.usable:
                    difference = delta_pct(right.qps, left.qps)
                    floor = noise_floor(left, right)
                    row.append(f"{pct(difference)} (floor {floor:.1f}%)")
                else:
                    row.append("-")

            row += [count(here[arm].major_faults) if here.get(arm) else "-" for arm in arms]
            usable = [summary for summary in here.values() if summary]
            row.append(gb(max((summary.rss_gb for summary in usable), default=0)) if usable else "-")
            table.row(*row)

        report.table(table)
        if best:
            report.lines(
                [f"  best {arm:<8s} threads={threads} at {value:,.0f} qps" for arm, (threads, value) in best.items()]
            )
        _fault_flatness(report, cells, ceiling, arms, thread_counts)

    notes = caveats(list(cells.values()))
    if notes:
        report.heading("caveats", level=3, folded=True).lines([f"  * {note}" for note in notes])
    if suite.notes:
        report.note(suite.notes)


def _fault_flatness(report: Report, cells: dict, ceiling: str, arms: list[str], thread_counts: list[str]) -> None:
    """Major faults must stay flat across thread counts; if they move, something else moved too.

    Threads change how many faults are in flight, not how many there are. This makes that check an
    explicit line rather than something a reader has to eyeball across a column.
    """
    for arm in arms:
        values = [
            cells[(ceiling, threads, arm)].major_faults
            for threads in thread_counts
            if (ceiling, threads, arm) in cells and cells[(ceiling, threads, arm)].usable
        ]
        values = [value for value in values if value]
        if len(values) < 2:
            continue
        spread = (max(values) - min(values)) / min(values) * 100
        if spread > 5.0:
            report.warn(
                f"{arm}: major faults vary {spread:.1f}% across thread counts at this ceiling. "
                f"Threads change how many faults are in flight, not how many there are — something "
                f"other than the thread count differed between these cells."
            )


def _cells(loaded: list[Record]) -> dict[tuple[str, str, str], Summary]:
    cells = {}
    for key, group_records in group(loaded).items():
        dims = dict(key)
        ceiling = dims.get("ceiling_gb", "0")
        threads = dims.get("threads", DEFAULT_THREADS)
        cap = float(ceiling or 0)
        cells[(ceiling, threads, dims.get("arm", "?"))] = summarise(group_records, cap or None)
    return cells
