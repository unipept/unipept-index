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

    # Part-to-whole per configuration, horizontal because the configuration names are long. What
    # this shows that the table cannot: preloading a structure does not remove its cost, it moves
    # the cost into a different segment of the same bar.
    report.heading("summary", level=3)
    ordered = [arm.name for arm in suite.arms if arm.name in summaries]
    report.chart(
        stacked_rows(
            ordered,
            [
                Series(
                    name=label,
                    values=[summaries[arm].startup.get(field, 0) / 1000 for arm in ordered],
                    slot=slot,
                )
                for slot, (field, label) in enumerate(PHASES + (("warmup_ms", "warmup"),))
            ],
            "Time before the first query can be answered",
            unit="s",
        ),
        "Time before the first query can be answered",
    )

    report.heading("per configuration", level=3)
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
    report.table(table)

    if suite.drop_caches:
        report.para("Page cache dropped before every configuration: these are cold-boot numbers.")
    else:
        report.para(
            "Warm numbers: after the first configuration the index files are in the page cache, so "
            "the later loads are memcpy rather than disk. Pass --cold for first-boot-after-deploy."
        )
    if suite.notes:
        report.note(suite.notes)
