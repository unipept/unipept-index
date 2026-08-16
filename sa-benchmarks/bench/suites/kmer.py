"""What each k-mer table buys, read against attaching no table at all.

Every row is a delta from the `none` cell in the same length regime on the same backend, because
that is the only comparison the suite is for: absolute throughput here is `defaults`' job, and a
table's value only means something relative to not having one.

The resident cost travels with the win everywhere it is shown. A table is a space-for-probes trade —
on the full database the 6-mer is 3.06 GB against the 5-mer's 127 MB — and this suite runs
unconstrained, which is exactly the regime in which the cost is invisible in the throughput. A
report that showed the win without the bytes would recommend the 6-mer every time.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, by_residency, lines
from ..config import Suite
from ..records import Record, delta_pct
from ..report import Report, Table, band, gb, pct, qps
from .shared import (
    GRID_KEYS,
    by_cell,
    cell_band,
    floor_of,
    held_and_swept,
    kmer_label,
    column_label,
    context_text,
    knob_points,
    phase_switch,
    resolution_table,
    tips_for,
    varying,
)

#: Length regimes short to long; anything else follows alphabetically.
BUCKET_ORDER = ("summary", "mixed", "small", "medium", "large")

KMER_INDEX = GRID_KEYS.index("kmer_k")
NO_TABLE = 0

#: `KmerTable` is dense: `AMINO_ACID_COUNT ** k` entries of `(usize, usize)`. See
#: `sa-index/src/kmer_table.rs` — 24 amino acids including the ambiguity codes, 16 bytes a pair.
AMINO_ACIDS = 24
BYTES_PER_ENTRY = 16


def table_gb(k: int) -> float:
    """What attaching a k-mer table costs, in GB. 127 MB at k=5, 3.06 GB at k=6."""
    return 0.0 if k <= 0 else AMINO_ACIDS**k * BYTES_PER_ENTRY / 2**30


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    # Drift-corrected: this suite's process sweeps every table across every regime, which takes
    # long enough that the machine moves under it. The cadence is in `kmer.toml`.
    cells = by_cell(loaded, correct_drift=True)
    arms = [arm.name for arm in suite.arms]
    buckets = _ordered(sorted({key[0] for key in cells}))
    sizes = sorted({key[KMER_INDEX] for key in cells})
    # Everything else the suite varies — the two search options, and the query count that follows
    # tryptic. Each becomes a row of its own rather than three of the four being dropped on the
    # floor by a lookup that took the first match.
    extra = [name for name in varying(cells) if name not in ("kmer_k", "peptide_source")]
    points = knob_points(cells, GRID_KEYS, KMER_INDEX, extra)

    report.heading("summary", level=3)
    held_and_swept(report, loaded)
    resolution_table(report, loaded)
    if NO_TABLE not in sizes:
        report.warn(
            "no cell ran without a table, so there is nothing to read the tables against. Add 0 to "
            "`kmer` in the suite file — it is the reference row, not an extra one."
        )
        return

    _curve(report, points, buckets, sizes, arms, extra)

    report.heading("per cell", level=3, folded=True)
    panels = [
        (
            f"{source}{' · ' + context_text(extra, context) if context_text(extra, context) else ''}",
            [kmer_label(k) for k in sizes],
            (lambda i, arm, p=per_arm, z=sizes: p.get(arm, {}).get(z[i])),
        )
        for source in buckets
        for context, per_arm in sorted(((c, p) for (b, c), p in points.items() if b == source), key=str)
    ]
    phase_switch(report, panels, arms, "k-mer table attached", title="per cell")
    _cell_table(report, points, buckets, sizes, arms, extra)

    if suite.notes:
        report.note(suite.notes)


def _ordered(sources: list[str]) -> list[str]:
    known = [name for name in BUCKET_ORDER if name in sources]
    return known + [name for name in sources if name not in BUCKET_ORDER]


def _curve(report: Report, points: dict, buckets: list[str], sizes: list[int], arms: list[str], extra: list[str]) -> None:
    """One line per (context, arm): what the table does as k grows, as a percentage of no table.

    Normalised rather than absolute so the contexts share an axis — the regimes differ by two orders
    of magnitude in qps and tryptic differs again, and on a shared absolute axis every short-peptide
    line would sit flat against the bottom, which is exactly where the effect is. 100% is the
    no-table row, so a line that stays flat is a table that bought nothing.

    One panel per (regime, context) and one line per arm inside it — the arms take the ordinal ramp,
    so the panel title carries the regime rather than a hue, and no panel holds more than three
    lines. See `charts.facets` for why that is the shape rather than one axes carrying all of them.
    """
    from ..charts import facets, panel_min

    panels = []
    for source in buckets:
        for (bucket, context), per_arm in sorted(points.items(), key=str):
            if bucket != source:
                continue
            label_text = context_text(extra, context)
            series = []
            for arm in by_residency(arms):
                curve = per_arm.get(arm, {})
                reference = curve.get(NO_TABLE)
                if not reference or not reference["p50"]:
                    continue
                series.append(
                    Series(
                        arm,
                        [
                            curve[k]["p50"] / reference["p50"] * 100 if k in curve else None
                            for k in sizes
                        ],
                        arm=arm,
                        tip={
                            "peptides": source,
                            **{n: column_label(n, v) for n, v in zip(extra, context)},
                            "backend": arm,
                        },
                    )
                )
            if series:
                panels.append((f"{source}{' · ' + label_text if label_text else ''}", series))
    if not panels:
        return
    x_labels = [kmer_label(k) for k in sizes]
    caption = "Throughput as a percentage of no table, per k"
    report.figures(
        facets(
            panels,
            lambda name, series, frame, top, legend, bottom: lines(
                x_labels, series, name, unit="%", frame=frame,
                x_title="k-mer table attached", y_title="% of no table",
                y_max=top, y_min=bottom, legend=legend, baseline=100.0,
            ),
            axes="dashed rule: no table, at 100%",
            floor=panel_min,
        ),
        caption,
    )


def _cell_table(report: Report, points: dict, buckets: list[str], sizes: list[int], arms: list[str], extra: list[str]) -> None:
    # One column per swept coordinate, each its own chip group. A single joined `context` column
    # offers one chip per combination that occurred, which is a lookup rather than a filter. `file`
    # is one of them now: four sibling tables could not answer "every tryptic row, every regime".
    headers = ["file", "kmer", *extra, "arm", "qps", "band", "vs none", "floor", "table GB", "verdict"]
    table = Table(
        headers=headers,
        aligns=["<"] * (len(extra) + 3) + [">"] * 5 + ["<"],
        chips=["file", "kmer", *extra, "arm"],
        tips=tips_for(headers),
    )
    ordered = [
        (source, context, per_arm)
        for source in buckets
        for context, per_arm in sorted(((c, p) for (b, c), p in points.items() if b == source), key=str)
    ]
    for k in sizes:
        for source, context, per_arm in ordered:
            for arm in arms:
                curve = per_arm.get(arm, {})
                cell, reference = curve.get(k), curve.get(NO_TABLE)
                if not cell:
                    continue
                if k == NO_TABLE:
                    difference, floor, reading = float("nan"), float("nan"), "reference"
                else:
                    difference = delta_pct(cell["p50"], reference["p50"]) if reference else float("nan")
                    floor = floor_of(cell, reference)
                    reading = _reading(difference, floor, k)
                table.row(
                    source,
                    kmer_label(k),
                    *(column_label(name, value) for name, value in zip(extra, context)),
                    arm,
                    qps(cell["p50"]),
                    band(cell_band(cell)),
                    "-" if difference != difference else pct(difference),
                    "-" if floor != floor else band(floor),
                    "-" if k == NO_TABLE else gb(table_gb(k)),
                    reading,
                )
    report.table(table, raw=True)


def _reading(difference: float, floor: float, k: int) -> str:
    """The verdict, which is about the trade rather than about the sign of the delta.

    A table that wins by less than the floor has not been shown to win at all, and it is still
    holding its gigabytes — so it reads as a cost, not as a tie. That asymmetry is deliberate: the
    null action here is not attaching the table.
    """
    if difference != difference:
        return "no data"
    if difference <= -floor:
        return "SLOWER than no table"
    if abs(difference) <= floor:
        return f"unresolved — {table_gb(k):.2f} GB for nothing this run can show"
    return "the table pays"


