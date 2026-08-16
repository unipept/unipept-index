"""The production-defaults grid: preloaded against mmap, per search mode, per length regime.

Matrix-mode records already aggregate their reps, so each line here is one grid cell carrying its
own p10/p50/p90.

This suite is the regression gate, which decides everything about how it reports. It varies only the
two search options — the rest of the grid moved out to `kmer`, `mlp` and `validate` — so the tables
are four rows per length regime and the columns that would have read the same value in every row are
stated once, in prose, instead. The comparison against a baseline lives here too: with
`--baseline <session>`, every cell is diffed against the same cell in a previous run and only
movements that clear the wider of the two cells' bands are called changes.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, by_residency, grouped_columns
from ..config import Suite
from ..records import Record, delta_pct, load_dir, median
from ..report import Report, Table, band, count, pct, qps
from .shared import (
    GRID_KEYS,
    fmt_tune,
    by_cell,
    cell_band,
    floor_of,
    held_and_swept,
    phase_switch,
    kmer_label,
    ratio_of,
    tips_for,
    varying,
)

#: The order peptide files are reported in, whichever order they were measured in. It runs from the
#: whole-picture views to the length regimes, so a reader meets the summary before the detail:
#: `summary` is the cross-file overview, `mixed` the unbucketed 5..50 file, then the three buckets
#: short to long. A file not named here follows, alphabetically.
BUCKET_ORDER = ("summary", "mixed", "small", "medium", "large")

#: How a coordinate is printed in a table cell.
FORMATTERS = {"kmer_k": kmer_label, "amount_of_peptides": lambda n: f"{n:,}q"}


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = by_cell(loaded)
    # Instrumented arms are excluded from every throughput table, chart and verdict below. Their
    # counters perturb the very number those are about; they exist for `_inside_the_search`.
    arms = [arm.name for arm in suite.arms if not arm.metrics]
    instrumented = [arm.name for arm in suite.arms if arm.metrics]
    # Everything except the peptide file, which is not a column — it is what the tables are split by.
    columns = [key for key in varying(cells) if key != "peptide_source"]

    # The answer first, then the configuration it was measured at, then the grid behind it.
    _verdict_tiles(report, cells, arms)
    _summary(report, cells, arms, loaded)

    sources = _ordered(sorted({key[0] for key in cells}))
    report.heading("by length regime", level=3)
    # Figure first, grid folded underneath: the chart is what points at something, the grid is what
    # you open once it has. One grid for every regime rather than a section each — the regimes are
    # the comparison, so they belong on one scale in one figure.
    _regime_figure(report, cells, sources, arms, columns)
    _regime_table(report, cells, sources, arms, columns)

    _response_share(report, cells, arms)

    if instrumented:
        _inside_the_search(report, loaded, instrumented, columns)

    baseline = getattr(suite, "baseline", None)
    if baseline:
        _regressions(report, cells, baseline, arms, columns)

    if suite.notes:
        report.note(suite.notes)


def _ordered(sources: list[str]) -> list[str]:
    known = [name for name in BUCKET_ORDER if name in sources]
    return known + [name for name in sources if name not in BUCKET_ORDER]


def _label(key: str, value) -> str:
    return FORMATTERS.get(key, fmt_tune)(value)


def _coords(key: tuple, columns: list[str]) -> list[str]:
    """The cell's values for the columns being shown, in column order."""
    return [_label(name, key[GRID_KEYS.index(name)]) for name in columns]


def _summary(report: Report, cells: dict, arms: list[str], loaded: list[Record]) -> None:
    """What this version does at production defaults, before any of the grid detail."""
    report.heading("summary", level=3)
    held_and_swept(report, loaded)

    # The headline cell: production search options, whatever else this suite happens to hold. Found
    # by matching the two options rather than by naming a full coordinate tuple, so narrowing the
    # grid further does not silently empty this chart.
    buckets = _ordered(sorted({key[0] for key in cells}))
    production: dict[str, list[float | None]] = {arm: [] for arm in arms}
    for source in buckets:
        match = next(
            (key for key in sorted(cells) if key[0] == source and key[1] is True and key[2] is False),
            None,
        )
        for arm in arms:
            production[arm].append((cells.get(match) or {}).get(arm, {}).get("p50") if match else None)

    if not any(any(value for value in values) for values in production.values()):
        report.para("no cell at production defaults (equate_il on, tryptic off).")
        return

    caption = "Throughput at production defaults, per peptide length regime"
    report.chart(
        grouped_columns(
            buckets,
            # The arms are one ordinal ramp — how much is resident — not three categorical slots.
            [Series(arm, production[arm], arm=arm, tip={"backend": arm}) for arm in by_residency(arms)],
            caption,
            unit=" qps",
            x_title="peptide length regime",
            y_title="throughput (qps)",
        ),
        caption,
    )
    report.para(
        "Production defaults: equate_il on, tryptic off, at the tuning stated above. The rows below "
        "vary only the two search options; what the k-mer table and the MLP batch cost is measured "
        "by `kmer` and `mlp`."
    )


def _share_heatmap(report: Report, rows: list[tuple]) -> None:
    """The measured share as a grid: length regime down, search options across.

    The 48-row table under this is exhaustive and unreadable at a glance, and the thing it holds is
    a SHAPE — the share collapses as the answer grows, so a non-tryptic short-peptide request is
    almost entirely decode while a tryptic long one is almost entirely search. As a grid that is one
    look; as a column of percentages between 2% and 94% it is a lookup.

    Sequential, not diverging: this is a magnitude with no meaningful midpoint. Averaged over the
    arms, because the share is a property of the workload — the three storage arms agree on it to
    well within the noise, and the table below carries each of them.
    """
    from ..charts import sequential_heatmap

    buckets: dict[tuple[str, str], list[float]] = {}
    for key, _arm, cell in rows:
        phases = [cell.get(name) or 0.0 for name in ("search_ms", "retrieval_ms", "decode_ms", "serialise_ms")]
        if not sum(phases):
            continue
        # Two lines, drawn as two lines — `_text_lines` turns the newline into a real break. As one
        # run these headers were wider than their column and overlapped their neighbours, which is
        # most of why this figure could not be read.
        options = f"equate_il={_label('equate_il', key[1])}\ntryptic={_label('tryptic', key[2])}"
        buckets.setdefault((key[0], options), []).append((phases[0] + phases[1]) / sum(phases) * 100)
    if not buckets:
        return

    files = _ordered(sorted({source for source, _ in buckets}))
    columns = sorted({options for _, options in buckets})
    grid = {}
    for row, source in enumerate(files):
        for column, options in enumerate(columns):
            shares = buckets.get((source, options))
            if not shares:
                continue
            mean = sum(shares) / len(shares)
            grid[(row, column)] = (
                mean,
                f"{source} · {options}\n"
                f"measured share: {mean:.0f}%\n"
                f"the other {100 - mean:.0f}% is annotation decode and JSON serialisation,\n"
                f"which no suite in this report times",
            )
    caption = "Percent of a request spent in search + retrieval (the phases this report times)"
    report.chart(
        sequential_heatmap(
            columns,
            files,
            grid,
            caption,
            low=0.0,
            high=100.0,
            unit="%",
            x_title="search options",
            y_title="peptide length regime",
        ),
        caption,
    )


def _verdict_tiles(report: Report, cells: dict, arms: list[str]) -> None:
    """What this run found, before any of the grid it found it in.

    The lead number is the MEASURED SHARE, not a throughput. Every other figure in the report is
    search plus retrieval, which on a non-tryptic large request is a single-digit percentage of what
    a caller waits for — so it is the number that scales every verdict in every other suite, and it
    belongs where a reader meets it first rather than three folded sections down.
    """
    shares, leaders = [], []
    for per_arm in cells.values():
        for arm, cell in per_arm.items():
            if arm not in arms:
                continue
            phases = [cell.get(name) or 0.0 for name in ("search_ms", "retrieval_ms", "decode_ms", "serialise_ms")]
            if cell.get("decode_ms") and sum(phases):
                shares.append((phases[0] + phases[1]) / sum(phases) * 100)
        ranked = sorted(
            ((arm, cell) for arm, cell in per_arm.items() if arm in arms and cell.get("p50")),
            key=lambda entry: -entry[1]["p50"],
        )
        if len(ranked) >= 2:
            floor = floor_of(*(cell for _, cell in ranked))
            if delta_pct(ranked[0][1]["p50"], ranked[1][1]["p50"]) > floor:
                leaders.append(ranked[0][0])

    cell_count = len(cells)
    tiles = []
    if shares:
        tiles.append((
            "measured share",
            f"{min(shares):.0f}–{max(shares):.0f}%",
            "of a whole request",
            # Not "unresolved": this number is perfectly well resolved. It is amber because it
            # scales every verdict in every other suite, which is a different kind of warning.
            "warn:read this first" if min(shares) < 50 else "",
        ))
    if leaders:
        top = max(set(leaders), key=leaders.count)
        tiles.append(("fastest arm", top, f"ahead in {leaders.count(top)} of {cell_count} cells", "good"))
    else:
        tiles.append(("fastest arm", "none", f"0 of {cell_count} cells separate", "flat"))

    reading = (
        "No configuration separates from the others by more than its own noise floor — on this box "
        "that floor is wider than any difference between the storage arms."
        if not leaders
        else f"`{max(set(leaders), key=leaders.count)}` is ahead of the runner-up by more than the "
        f"floor in {leaders.count(max(set(leaders), key=leaders.count))} of {cell_count} cells."
    )
    if shares and min(shares) < 50:
        reading += (
            f" Read every other suite through the measured share: a knob that buys X% of search "
            f"buys as little as {min(shares):.0f}% of X to a caller."
        )
    report.verdict(tiles, reading)


def _regime_table(
    report: Report, cells: dict, sources: list[str], arms: list[str], columns: list[str]
) -> None:
    """Every cell, one row each, with the peptide file as its first column.

    One table rather than one per regime. The rows are identical either way; what changes is that a
    `file` chip filters to any one regime in a click, where four separate tables made "show me every
    tryptic row across the regimes" impossible in the page and a scroll in the markdown.

    The `verdict` column states the ratio's reading rather than leaving the reader to compare it
    against `noise` themselves. Two adjacent columns and a subtraction is exactly the step that gets
    skipped, and skipping it turns every 1.03x into a result.
    """
    # A ratio column per arm after the first, all against that first arm. With two arms this is the
    # one `ratio` column it always was; with three it stays readable, where a single ratio would
    # have to pick two of them and quietly drop the third.
    ratios = [f"{arm}/{arms[0]}" for arm in arms[1:]]
    headers = ["file", *columns, *arms, *ratios, "noise", "verdict"]
    table = Table(
        headers=headers,
        aligns=["<"] * (len(columns) + 1) + [">"] * (len(arms) + len(ratios) + 1) + ["<"],
        chips=["file", *columns],
        tips={
            **tips_for(headers),
            **{name: RATIO_TIP.format(name=name) for name in ratios},
        },
    )
    for source in sources:
        for key in sorted(key for key in cells if key[0] == source):
            values = [cells[key].get(arm) for arm in arms]
            table.row(
                source,
                *_coords(key, columns),
                *(qps(value["p50"]) if value else "-" for value in values),
                *(ratio_of([values[0], value]) for value in values[1:]),
                band(max((cell_band(value) for value in values if value), default=float("nan"))),
                _verdict(values, arms),
            )
    report.table(table, raw=True)


RATIO_TIP = (
    "`{name}` at the median. 1.00x means they measured the same. Only meaningful once its distance "
    "from 1.00 exceeds the noise column."
)


def _verdict(values: list[dict | None], arms: list[str]) -> str:
    """Which backend is ahead, once more than two of them can be.

    The comparison that matters with three arms is the leader against the RUNNER-UP, not against a
    fixed reference: a leader that beats the slowest arm by 30% while sitting inside the second's
    floor has not been shown to lead.
    """
    ranked = sorted(
        ((arm, cell) for arm, cell in zip(arms, values) if cell and cell["p50"]),
        key=lambda entry: -entry[1]["p50"],
    )
    if len(ranked) < 2:
        return "-"
    (best, top), (_, second) = ranked[0], ranked[1]
    difference = delta_pct(top["p50"], second["p50"])
    floor = floor_of(*(cell for _, cell in ranked))
    if difference <= floor:
        return f"cannot separate them (within the {floor:.1f}% floor)"
    return f"{best} ahead by {difference:.1f}%"


def _regime_figure(
    report: Report, cells: dict, sources: list[str], arms: list[str], columns: list[str]
) -> None:
    """One panel per length regime, all on one scale — see `shared.phase_switch`."""
    panels = []
    for source in sources:
        keys = sorted(key for key in cells if key[0] == source)
        if not keys:
            continue
        groups = [" · ".join(_coords(key, columns)) if columns else "default" for key in keys]
        # `keys` is bound per panel rather than closed over the loop variable, which would leave
        # every panel reading the last regime's cells.
        panels.append((source, groups, lambda index, arm, keys=keys: cells[keys[index]].get(arm)))
    phase_switch(
        report,
        panels,
        arms,
        " · ".join(columns) if columns else "configuration",
    )


def _regressions(
    report: Report, cells: dict, baseline_dir: Path, arms: list[str], columns: list[str]
) -> None:
    """Diffs every cell against the same cell in a previous session."""
    report.heading("regression check against the baseline", level=3)
    previous = by_cell(load_dir(baseline_dir))
    if not previous:
        report.warn(f"no records under {baseline_dir} — nothing to compare against")
        return

    # Against a baseline the columns are whichever coordinates differ in EITHER run: a cell that
    # only one of them has still has to be named in full, or two rows collapse into one label.
    shown = [key for key in varying({**cells, **previous}) if key != "peptide_source"] or columns
    table = Table(
        headers=["file", *shown, "arm", "base", "now", "delta", "verdict"],
        aligns=["<"] * (len(shown) + 1) + ["<", ">", ">", ">", "<"],
        tips=tips_for(["file", *shown, "arm", "base", "now", "delta", "verdict"]),
    )
    moved = 0
    for key in sorted(set(cells) & set(previous)):
        for arm in arms:
            now, base = cells[key].get(arm), previous[key].get(arm)
            if not (now and base):
                continue
            difference = delta_pct(now["p50"], base["p50"])
            floor = floor_of(now, base)
            changed = abs(difference) > floor
            moved += changed
            table.row(
                key[0],
                *_coords(key, shown),
                arm,
                qps(base["p50"]),
                qps(now["p50"]),
                pct(difference),
                "REGRESSION" if changed and difference < 0 else ("improvement" if changed else "unchanged"),
            )
    report.table(table)

    only_now = sorted(set(cells) - set(previous))
    only_before = sorted(set(previous) - set(cells))
    if only_now or only_before:
        report.para(
            f"{len(only_now)} cell(s) exist only in this run and {len(only_before)} only in the "
            f"baseline; those cannot be compared and are not counted above."
        )
    report.para(
        f"{moved} of {len(set(cells) & set(previous)) * len(arms)} comparable cells moved by more "
        f"than their own noise floor."
    )


# ---------------------------------------------------------------------------
# Inside the search — the instrumented arm
# ---------------------------------------------------------------------------

#: Column explanations for the counters. Only this section has them, so they live here.
PHASE_TIPS = {
    "search": "Share of the timed region spent finding matching suffixes.",
    "retrieval": "Share of the timed region spent resolving those suffixes to proteins.",
    "bounds": (
        "Share of search THREAD-time in the binary search — a dependent chain of probes that no "
        "prefetch can help and only the k-mer table shortens."
    ),
    "iter": (
        "Share of search thread-time scanning the matched suffix range. Contiguous, so readahead "
        "and the prefetch distance both reach it."
    ),
    "examined": "Candidate suffixes the range scan looked at.",
    "accept%": (
        "Candidates accepted as real matches, over candidates examined. A LOW rate with a low "
        "examined count means the work is already minimal and the hit rate is simply low; a low "
        "rate with a high examined count means exhaustive scanning. The two need opposite fixes, "
        "and this is the only column in the report that separates them."
    ),
    "parallelism": (
        "Search thread-time over search wall-time — how many cores the search phase actually kept "
        "busy. Well below the core count means the work did not spread."
    ),
}


def _inside_the_search(report: Report, loaded: list[Record], arms: list[str], columns: list[str]) -> None:
    """What the counters say, from the one instrumented arm.

    This was the `detail` suite. It is here because its two useful numbers are both about the search
    OPTIONS — the acceptance rate is the answer to "what is tryptic actually costing", and that
    question only exists in the tryptic rows this grid already runs — and because a separate suite
    meant a separate index load and a second, redundant MLP batch sweep that `mlp` measures better
    on an uninstrumented build.
    """
    report.heading("inside the search (instrumented)", level=3, folded=True)
    report.warn(
        f"the {', '.join(arms)} arm is built with `metrics`, whose counters perturb throughput by "
        "~2%. Its qps is shown for scale only and is excluded from every table above."
    )

    cells = _phase_cells(loaded, arms)
    _phase_split_chart(report, cells, columns)

    table = Table(
        headers=[*columns, "arm", "file", "qps", "search", "retrieval", "bounds", "iter", "examined", "accept%", "parallelism"],
        aligns=["<"] * len(columns) + ["<", "<"] + [">"] * 7,
        chips=[*columns, "arm", "file"],
        tips={**tips_for([*columns, "arm", "qps"]), **PHASE_TIPS},
    )
    for key, per_arm in sorted(cells.items(), key=str):
        source, coords = key[0], key[1:]
        for arm, phases in sorted(per_arm.items()):
            total = phases["total_ns"]
            thread_time = phases["bounds_ns"] + phases["iter_ns"]
            table.row(
                *(_label(name, value) for name, value in zip(columns, coords)),
                arm,
                source,
                qps(phases["qps"]),
                _share(phases["search_ns"], total),
                _share(phases["retrieval_ns"], total),
                _share(phases["bounds_ns"], thread_time),
                _share(phases["iter_ns"], thread_time),
                count(phases["examined"]),
                _accept(phases),
                f"{thread_time / phases['search_ns']:.1f}x" if phases["search_ns"] else "-",
            )
    report.table(table, raw=True)
    report.para(
        "`bounds` and `iter` are shares of search THREAD-time, which is summed across every rayon "
        "thread and so exceeds wall time by the parallelism factor in the last column — they are "
        "meaningful against each other, not against the clock. Which of the two dominates is what "
        "decides whether a prefetch distance ever had a phase to work on; see `prefetch`."
    )


def _phase_split_chart(report: Report, cells: dict, columns: list[str]) -> None:
    """Where search thread-time goes: the binary search, or the range scan.

    Drawn as shares summing to 100% rather than as two absolute bars, because the question is which
    of the two DOMINATES — that is what decides whether a prefetch distance ever had a phase to work
    on, and it is the split the README says tryptic inverts. Two bars of absolute nanoseconds put
    the four length regimes two orders of magnitude apart and hide the very ratio being asked about.
    """
    from ..charts import Series, facets, stacked_columns

    # One panel per length regime, not one axis carrying all sixteen cells. Every coordinate on one
    # x axis made a label like `large · false · true` sixteen times over in 618px: each one was
    # about three times its own slot, so all fifteen adjacent pairs overlapped and the axis was a
    # smear. The regime is the panel, which is how every other figure in this report is cut.
    per_regime: dict[str, tuple[list[str], list[float], list[float]]] = {}
    for key, per_arm in sorted(cells.items(), key=str):
        source, coords = key[0], key[1:]
        for _arm, phases in sorted(per_arm.items()):
            thread_time = phases["bounds_ns"] + phases["iter_ns"]
            if not thread_time:
                continue
            groups, bounds, iters = per_regime.setdefault(source, ([], [], []))
            groups.append(" · ".join(_label(name, value) for name, value in zip(columns, coords)))
            bounds.append(phases["bounds_ns"] / thread_time * 100)
            iters.append(phases["iter_ns"] / thread_time * 100)
    if not per_regime:
        return

    groups_of = {regime: groups for regime, (groups, _, _) in per_regime.items()}
    built = [
        (
            regime,
            [
                Series("bounds (binary search)", per_regime[regime][1], 0, tip={"phase": "bounds"}),
                Series("iter (range scan)", per_regime[regime][2], 1, tip={"phase": "iter"}),
            ],
        )
        for regime in _ordered(list(per_regime))
    ]

    caption = "Where search thread-time goes: binary search against range scan"
    report.figures(
        facets(
            built,
            lambda name, series, frame, top, legend: stacked_columns(
                groups_of[name], [""], series, name, unit="%",
                x_title=" · ".join(columns), y_title="share of search thread-time (%)",
                frame=frame, y_max=top, legend=legend,
            ),
            # Already shares of their own cell, so the axis is 100% by construction.
            extent=lambda series: 100.0,
        ),
        caption,
    )


def _phase_cells(loaded: list[Record], arms: list[str]) -> dict[tuple, dict[str, dict]]:
    """`(source, *coords) -> arm -> counters`, median across each cell's reps."""
    grouped: dict[tuple, dict[str, list[Record]]] = {}
    for record in loaded:
        arm = record.dims.get("arm", "?")
        if arm not in arms or record.config.get("sweep") == "drift":
            continue
        config = record.config
        key = (config.get("peptide_source"), config.get("equate_il"), config.get("tryptic"))
        grouped.setdefault(key, {}).setdefault(arm, []).append(record)

    out: dict[tuple, dict[str, dict]] = {}
    for key, per_arm in grouped.items():
        for arm, records in per_arm.items():
            def field(name: str) -> float:
                return median(record.result.get(name, 0) for record in records)

            out.setdefault(key, {})[arm] = {
                "qps": median(record.qps for record in records),
                "total_ns": field("total_duration_ns"),
                "search_ns": field("search_duration_ns"),
                "retrieval_ns": field("retrieval_duration_ns"),
                "bounds_ns": field("search_bounds_ns"),
                "iter_ns": field("match_iter_ns"),
                "examined": field("candidates_examined"),
                "accepted": field("candidates_accepted"),
            }
    return out


def _share(part: float, whole: float) -> str:
    return "-" if not whole else f"{part / whole * 100:.0f}%"


def _accept(phases: dict[str, float]) -> str:
    """Candidate acceptance rate; zero counters mean the build had no `metrics`."""
    examined = phases["examined"]
    if not examined:
        return "n/a"
    return f"{phases['accepted'] / examined * 100:.1f}%"


def _response_share(report: Report, cells: dict, arms: list[str]) -> None:
    """What fraction of a request the throughput in this report actually covers.

    `throughput_qps` is search plus retrieval, everywhere, deliberately — widening it would change
    what every suite means. But production does two more things before a client has bytes: it turns
    each `ProteinRef` into a `ProteinInfo`, decoding the functional annotations and allocating the
    accession, and it serialises the result to JSON. This suite times both.

    The number that matters is the last column. If decode and serialisation are most of a request,
    then a knob that buys 20% of search buys far less than 20% to a user, and every verdict in the
    rest of the report has to be read through that.
    """
    rows = [
        (key, arm, cell)
        for key, per_arm in sorted(cells.items(), key=str)
        for arm, cell in per_arm.items()
        if arm in arms and cell.get("decode_ms")
    ]
    if not rows:
        return

    report.heading("what a request actually costs", level=3)
    report.para(
        "**A server request has four phases; this report times the first two.** `search` finds the "
        "matching suffixes, `retrieval` turns them into protein references — and then `decode` "
        "unpacks each hit's annotations and `serialise` writes the JSON that goes back to the "
        "caller. Every throughput number in this report, in every suite, is phases 1-2 only. Only "
        "`defaults` pays to measure phases 3-4 at all, and this is where they are reported."
    )
    report.para(
        "The grid below is one number per cell: **what percentage of a request the timed phases "
        "were**. 90% means search and retrieval dominate and a knob that speeds them up reaches the "
        "caller almost undiluted. 10% means the request is nearly all decode and JSON, and the same "
        "knob buys a tenth of what it appears to. The share collapses as the answer grows — a short "
        "non-tryptic peptide matches half the database and spends its time serialising the result, "
        "while a long tryptic one matches almost nothing and spends it all in search."
    )
    _share_heatmap(report, rows)
    table = Table(
        headers=["file", "equate_il", "tryptic", "arm", "search", "retrieval", "decode", "serialise", "response KB", "measured share"],
        aligns=["<", "<", "<", "<", ">", ">", ">", ">", ">", ">"],
        chips=["file", "equate_il", "tryptic", "arm"],
        tips={
            **tips_for(["arm"]),
            "decode": (
                "Turning each ProteinRef into the ProteinInfo the server returns: an fa-compression "
                "decode of the functional annotations plus a String for the accession, per hit."
            ),
            "serialise": "Serialising the result to JSON, the last thing the server does.",
            "response KB": "How much JSON the request would have returned.",
            "measured share": (
                "Search plus retrieval over all four phases — the fraction of a request that every "
                "throughput figure in this report is measuring. The rest is real work no suite times."
            ),
        },
    )
    shares = []
    for key, arm, cell in rows:
        phases = [cell.get(name) or 0.0 for name in ("search_ms", "retrieval_ms", "decode_ms", "serialise_ms")]
        whole = sum(phases)
        if not whole:
            continue
        share = (phases[0] + phases[1]) / whole * 100
        shares.append(share)
        table.row(
            key[0],
            _label("equate_il", key[1]),
            _label("tryptic", key[2]),
            arm,
            *(f"{value:,.1f} ms" for value in phases),
            f"{(cell.get('response_bytes') or 0) / 1024:,.0f}",
            f"{share:.0f}%",
        )
    report.table(table, raw=True)
    if shares:
        low, high = min(shares), max(shares)
        report.warn(
            f"every throughput figure in this report — here and in every other suite — measures "
            f"search plus retrieval, which is {low:.0f}–{high:.0f}% of a request. The remainder is "
            f"annotation decoding and JSON serialisation, which the server does on every call. A "
            f"knob that buys X% of the measured part buys less than X% to a caller."
        )
