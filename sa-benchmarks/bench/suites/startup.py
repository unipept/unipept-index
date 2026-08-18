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
from ..report import Report, Table, count, gb, seconds

#: Startup fields, in the order they are paid.
PHASES = (
    ("load_sa_ms", "sa"),
    ("load_proteins_ms", "proteins"),
    ("load_mapping_ms", "mapping"),
    ("kmer_table_ms", "kmer"),
)

#: Per-structure byte counters the page sweep reports, summed into what the whole sweep touched.
SWEPT_BYTES = ("warmup_sa_bytes", "warmup_proteins_bytes", "warmup_mapping_bytes")

#: Per-structure sweep durations. The three sweeps run concurrently, so the pass takes as long as
#: the slowest of them — summing these would divide the bytes by roughly three times the wall clock.
SWEPT_MS = ("warmup_sa_ms", "warmup_proteins_ms", "warmup_mapping_ms")


def swept_gb(startup: dict) -> float:
    return sum(startup.get(field, 0) for field in SWEPT_BYTES) / 2**30


def sweep_ms(startup: dict) -> float:
    """Wall time of the page sweep alone.

    Not `warmup_ms`: under `all:N` that also covers the peptides pushed through afterwards, which
    would drag the rate down by an arbitrary amount and make two suites' columns incomparable.
    """
    return max((startup.get(field, 0) for field in SWEPT_MS), default=0)


def sweep_gbps(startup: dict) -> float:
    """GB/s of the page sweep, or NaN when there was nothing to sweep or no time to measure it in."""
    gigabytes, milliseconds = swept_gb(startup), sweep_ms(startup)
    if not gigabytes or not milliseconds:
        return float("nan")
    return gigabytes / (milliseconds / 1000)


def sweep_rate(startup: dict) -> str:
    """Bytes per second of the page sweep — the column that separates an arm from its position.

    Elapsed time alone cannot say whether a sweep read from the device or from the page cache, and
    the two differ by an order of magnitude over the same bytes. Two arms that sweep the same
    structure must land at the same rate; if they do not, one of them was handed a warm cache and
    the pair is not comparable. An arm with nothing mapped sweeps no bytes at all, which is a
    different statement from sweeping them quickly and is printed as such.
    """
    if not swept_gb(startup):
        return "nothing mapped"
    rate = sweep_gbps(startup)
    return "-" if rate != rate else f"{rate:.2f}"


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
        headers=[
            "config",
            *(label for _, label in PHASES),
            "load",
            "warmup",
            "GB swept",
            "GB/s",
            "majflt",
            "total",
            "RSS GB",
        ],
        aligns=["<"] + [">"] * (len(PHASES) + 7),
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
            gb(swept_gb(startup)),
            sweep_rate(startup),
            count(startup.get("warmup_major_faults")),
            seconds(total),
            gb(summary.rss_gb),
        )
    report.table(table, raw=True)

    if suite.drop_caches:
        report.para(
            "Page cache dropped before every configuration: these are cold-boot numbers, and the "
            "`GB/s` column should read the same for two arms sweeping the same structure."
        )
    else:
        report.para(
            "WARM, and therefore ordered: the first configuration leaves the index in the page "
            "cache, so every later load and page sweep is memcpy rather than disk. Compare `GB/s` "
            "across the rows before reading `warmup` — an arm an order of magnitude faster per "
            "byte than another sweeping the same structure is measuring its position in the sweep. "
            "Pass --cold for numbers that compare arms."
        )
    _warn_if_sweeps_disagree(report, summaries, ordered)

    if suite.notes:
        report.note(suite.notes)


#: How far two arms' sweep bandwidths may diverge before the comparison is called void. Set well
#: above anything a real difference in what the arms sweep can produce — the arms cover overlapping
#: sections of the same three files — and well below the order of magnitude that separates a sweep
#: served by the device from one served by the page cache.
SWEEP_RATE_TOLERANCE = 2.0


def _warn_if_sweeps_disagree(report: Report, summaries: dict, ordered: list[str]) -> None:
    """Flags the failure this suite is most prone to: arms measuring the cache, not themselves.

    The arms sweep overlapping sections of the same three files, so their bytes-per-second should
    land within a small factor of each other whatever they are configured to hold. When they do not,
    the fast one was handed a page cache the slow one filled, and no `warmup` or `total` figure in
    the table means what it appears to. This is left to a check rather than to the reader because it
    is invisible in the column being read: a warm sweep looks exactly like a fast arm.
    """
    rates = [(arm, sweep_gbps(summaries[arm].startup)) for arm in ordered]
    rates = [(arm, rate) for arm, rate in rates if rate == rate]
    if len(rates) < 2:
        return

    slowest, fastest = min(rates, key=lambda e: e[1]), max(rates, key=lambda e: e[1])
    if fastest[1] <= slowest[1] * SWEEP_RATE_TOLERANCE:
        return

    report.warn(
        f"`{fastest[0]}` swept its pages at {fastest[1]:.2f} GB/s while `{slowest[0]}` managed "
        f"{slowest[1]:.2f} GB/s — {fastest[1] / slowest[1]:.1f}x apart over the same files. The two "
        f"were not handed the same page cache, so the `warmup` and `total` columns rank the arms by "
        f"their position in the sweep as much as by what they hold. Re-run with the cache dropped "
        f"before every configuration (`--cold`, which this suite sets by default) before quoting "
        f"any of it."
    )


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
