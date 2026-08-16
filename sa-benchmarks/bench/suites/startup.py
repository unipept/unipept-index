"""Startup cost per storage configuration.

One row per configuration, from the harness's `startup` section. `load` is work done before the
first query can be answered; `warmup` is the optional page-touch sweep. Preloading a structure moves
cost from the second into the first, which is why the totals are what get compared.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, stacked_rows
from ..config import Suite
from ..records import Record, group, summarise
from ..report import Report, Table, gb, seconds

#: Startup fields, in the order they are paid.
PHASES = (
    ("load_sa_ms", "sa"),
    ("load_proteins_ms", "proteins"),
    ("load_mapping_ms", "mapping"),
    ("kmer_table_ms", "kmer"),
)


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    summaries = {
        dict(key)["arm"]: summarise(cell) for key, cell in group(loaded).items() if "arm" in dict(key)
    }

    ordered = [arm.name for arm in suite.arms if arm.name in summaries]
    _verdict_tiles(report, summaries, ordered)

    # Part-to-whole per configuration, horizontal because the configuration names are long. What
    # this shows that the table cannot: preloading a structure does not remove its cost, it moves
    # the cost into a different segment of the same bar.
    report.heading("summary", level=3)
    report.chart(
        stacked_rows(
            ordered,
            [
                Series(
                    name=label,
                    values=[summaries[arm].startup.get(field, 0) / 1000 for arm in ordered],
                    slot=slot,
                    tip={"phase": label},
                )
                for slot, (field, label) in enumerate(PHASES + (("warmup_ms", "warmup"),))
            ],
            "Time before the first query can be answered",
            unit="s",
            x_title="time to first query (s)",
        ),
        "Time before the first query can be answered",
    )

    report.heading("per configuration", level=3, folded=True)
    table = Table(
        headers=["config", *(label for _, label in PHASES), "load", "warmup", "total", "RSS GB"],
        aligns=["<"] + [">"] * (len(PHASES) + 4),
    )
    # Suite order, not alphabetical: the arms are listed from fully preloaded to fully mapped, and
    # that ordering is what makes the trade visible down the column.
    for arm in suite.arms:
        summary = summaries.get(arm.name)
        if summary is None:
            table.row(arm.name, *["-"] * (len(PHASES) + 4))
            continue
        startup = summary.startup
        load = startup.get("load_total_ms")
        warmup = startup.get("warmup_ms")
        total = (load or 0) + (warmup or 0) if load is not None else None
        table.row(
            arm.name,
            *(seconds(startup.get(field)) for field, _ in PHASES),
            seconds(load),
            seconds(warmup),
            seconds(total),
            gb(summary.rss_gb),
        )
    report.table(table, raw=True)

    if suite.drop_caches:
        report.para("Page cache dropped before every configuration: these are cold-boot numbers.")
    else:
        report.para(
            "Warm numbers: after the first configuration the index files are in the page cache, so "
            "the later loads are memcpy rather than disk. Pass --cold for first-boot-after-deploy."
        )
    if suite.notes:
        report.note(suite.notes)


def _verdict_tiles(report: Report, summaries: dict, ordered: list[str]) -> None:
    """What each configuration costs before it can answer anything, as three numbers.

    The trade this suite exists to price is one sentence — preloading buys nothing at query time
    that it does not pay for at load time — and it is legible as the gap between the slowest and the
    fastest arm to first query. The per-phase breakdown underneath is where that gap comes from.
    """
    timed = [
        (arm, (summaries[arm].startup.get("load_total_ms") or 0) + (summaries[arm].startup.get("warmup_ms") or 0))
        for arm in ordered
        if summaries[arm].startup.get("load_total_ms") is not None
    ]
    if not timed:
        return
    fastest = min(timed, key=lambda entry: entry[1])
    slowest = max(timed, key=lambda entry: entry[1])
    rss = [(arm, summaries[arm].rss_gb) for arm in ordered if summaries[arm].rss_gb]
    spread = slowest[1] - fastest[1]

    report.verdict(
        [
            ("fastest to first query", fastest[0], f"{fastest[1] / 1000:.1f}s", "good"),
            ("slowest", slowest[0], f"{slowest[1] / 1000:.1f}s", ""),
            (
                "resident after load",
                f"{max(value for _, value in rss):.1f} GB" if rss else "",
                f"{min(value for _, value in rss):.1f} GB at the lightest" if rss else "",
                "",
            ),
        ],
        f"`{slowest[0]}` waits {spread / 1000:.1f}s longer than `{fastest[0]}` before its first "
        f"answer. What that buys at query time is what `defaults` and `ram` measure — preloading "
        f"does not remove a structure's cost, it moves it out of the query and into the load.",
    )
