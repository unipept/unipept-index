"""Where the time goes inside a search, and how much cross-query batching still buys.

Everything here comes from an instrumented build (`metrics`), so the throughput column is perturbed
and is present only to show the batch curve's shape. The absolute number that matters is in
`defaults`.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, lines, stacked_rows
from ..config import Suite
from ..records import Record, Summary, delta_pct, group, median, summarise
from ..report import Report, Table, pct, qps


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = _cells(loaded)
    arms = [arm.name for arm in suite.arms]
    batches = [str(value) for value in suite.axes.get("tune.mlp_batch", [1])]

    report.warn(
        "instrumented build (metrics): these qps are perturbed by the counters and are NOT this "
        "version's throughput. Read the shape of the curve and the phase split; read `defaults` for "
        "the number."
    )

    report.heading("summary", level=3)
    # One chart for both arms rather than one each: the question is how the curve bends, and two
    # lines on one scale answer it where two charts would make the reader hold a shape in memory.
    report.chart(
        lines(
            batches,
            [
                Series(
                    name=arm,
                    values=[cells[(arm, batch)].qps if (arm, batch) in cells else None for batch in batches],
                    slot=slot,
                )
                for slot, arm in enumerate(arms)
            ],
            "Throughput against MLP batch size (instrumented build)",
            unit=" qps",
            x_title="peptides interleaved per task",
        ),
        "Throughput against MLP batch size (instrumented build)",
    )

    for arm in arms:
        report.heading(f"{arm}: MLP batch sweep", level=3)
        report.chart(*_phase_split(cells, arm, batches))
        table = Table(
            headers=["mlp_batch", "qps", f"vs B={batches[0]}", "search", "retrieval", "bounds", "iter", "accept%"],
            aligns=["<"] + [">"] * 7,
        )
        baseline = cells.get((arm, batches[0]))
        for batch in batches:
            summary = cells.get((arm, batch))
            if summary is None:
                table.row("scalar" if batch == "1" else batch, *["-"] * 7)
                continue
            phases = summary.phases
            table.row(
                "scalar" if batch == "1" else batch,
                qps(summary.qps),
                "base" if batch == batches[0] else pct(delta_pct(summary.qps, baseline.qps)) if baseline else "-",
                _share(phases["search_ns"], phases["total_ns"]),
                _share(phases["retrieval_ns"], phases["total_ns"]),
                _share(phases["bounds_ns"], phases["bounds_ns"] + phases["iter_ns"]),
                _share(phases["iter_ns"], phases["bounds_ns"] + phases["iter_ns"]),
                _accept(phases),
            )
        report.table(table)

        reference = cells.get((arm, batches[0]))
        if reference:
            report.lines(_regime(reference))

    if suite.notes:
        report.note(suite.notes)


def _phase_split(cells: dict, arm: str, batches: list[str]) -> tuple[str, str]:
    """Where the search phase's thread-time goes, as batching increases.

    Part-to-whole, so a stack — and the point is the shift: batching moves time out of the
    dependent binary-search chain, which nothing can prefetch, into the contiguous range scan,
    which readahead and prefetch distance both reach.
    """
    rows, bounds, iterate = [], [], []
    for batch in batches:
        summary = cells.get((arm, batch))
        if not summary or not summary.phases:
            continue
        total = summary.phases["bounds_ns"] + summary.phases["iter_ns"]
        if not total:
            continue
        rows.append("scalar" if batch == "1" else f"mlp {batch}")
        bounds.append(summary.phases["bounds_ns"] / 1e6)
        iterate.append(summary.phases["iter_ns"] / 1e6)

    caption = f"{arm}: search-phase thread-time, binary search against range scan"
    return (
        stacked_rows(
            rows,
            [Series("bounds (binary search)", bounds, 0), Series("iter (range scan)", iterate, 2)],
            caption,
            unit=" ms",
        ),
        caption,
    )


def _regime(summary: Summary) -> list[str]:
    """The one-line statement of which regime this dataset is in, and the parallelism factor."""
    phases = summary.phases
    total, search, retrieval = phases["total_ns"], phases["search_ns"], phases["retrieval_ns"]
    bounds, iterate = phases["bounds_ns"], phases["iter_ns"]
    lines = [
        f"  total {total / 1e6:.1f} ms   search {search / 1e6:.1f} ms ({_pctf(search, total)})"
        f"   retrieval {retrieval / 1e6:.1f} ms ({_pctf(retrieval, total)})"
    ]
    if bounds + iterate > 0:
        lines.append(
            f"  within search (summed thread-time, NOT comparable to wall clock): "
            f"bounds {bounds / 1e6:.0f} ms ({_pctf(bounds, bounds + iterate)}) "
            f"| iter {iterate / 1e6:.0f} ms ({_pctf(iterate, bounds + iterate)})"
        )
        if search:
            lines.append(f"  effective parallelism: {(bounds + iterate) / search:.0f}x")
    else:
        lines.append("  no phase split — this build has no `metrics`, so the counters read zero")
    if phases["matches_per_query"]:
        lines.append(f"  {phases['matches_per_query']:,.0f} matches/query")
    return lines


def _cells(loaded: list[Record]) -> dict[tuple[str, str], Summary]:
    cells = {}
    for key, group_records in group(loaded).items():
        dims = dict(key)
        summary = summarise(group_records)
        summary.phases = _phases(group_records)
        cells[(dims.get("arm", "?"), dims.get("tune.mlp_batch", "1"))] = summary
    return cells


def _phases(records: list[Record]) -> dict[str, float]:
    """Median of each timing/counter field across a cell's reps."""

    def field(name: str) -> float:
        return median(record.result.get(name, 0) for record in records)

    queries = median(record.result.get("amount_of_queries", 0) for record in records)
    suffixes = field("suffix_hit_count")
    return {
        "total_ns": field("total_duration_ns"),
        "search_ns": field("search_duration_ns"),
        "retrieval_ns": field("retrieval_duration_ns"),
        "bounds_ns": field("search_bounds_ns"),
        "iter_ns": field("match_iter_ns"),
        "examined": field("candidates_examined"),
        "accepted": field("candidates_accepted"),
        "matches_per_query": suffixes / queries if queries else 0.0,
    }


def _share(part: float, whole: float) -> str:
    return "-" if not whole else f"{part / whole * 100:.0f}%"


def _pctf(part: float, whole: float) -> str:
    return "-" if not whole else f"{part / whole * 100:.0f}%"


def _accept(phases: dict[str, float]) -> str:
    """Candidate acceptance rate; zero counters mean the build had no `metrics`."""
    examined = phases["examined"]
    if not examined:
        return "n/a"
    return f"{phases['accepted'] / examined * 100:.1f}%"
