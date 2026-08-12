"""Fallback analysis: one row per cell, no interpretation.

Used when a suite has no analysis module of its own — a new suite is runnable and readable before
anyone has decided what its headline number should be. It shows every dimension, so nothing is
silently dropped just because this module does not know what it means.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record, group, summarise
from ..report import Report, Table, band, count, gb, qps


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    grouped = group(loaded)
    dimension_names = sorted({key for cell in grouped for key, _ in cell})

    table = Table(
        headers=[*dimension_names, "n", "qps", "band", "majflt", "RSS GB"],
        aligns=["<"] * len(dimension_names) + [">"] * 5,
    )
    summaries = []
    for key in sorted(grouped):
        dims = dict(key)
        summary = summarise(grouped[key], _ceiling(dims))
        summaries.append(summary)
        table.row(
            *(dims.get(name, "-") for name in dimension_names),
            summary.n,
            qps(summary.qps),
            band(summary.band),
            count(summary.major_faults),
            gb(summary.rss_gb),
        )
    report.table(table)

    from ..report import caveats

    notes = caveats(summaries)
    if notes:
        report.note("Caveats:\n" + "\n".join(f"  * {note}" for note in notes))


def _ceiling(dims: dict) -> float | None:
    value = dims.get("ceiling_gb")
    return float(value) if value not in (None, "0") else None
