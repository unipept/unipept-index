"""The production-defaults grid: preloaded against mmap, per configuration.

Matrix-mode records already aggregate their reps, so each line here is one grid cell carrying its
own p10/p50/p90. Cells are keyed on `config` rather than `dims`, because a matrix invocation sweeps
the grid *inside* one process — every record from one arm shares that arm's dims, and it is the
config that distinguishes them.

The regression comparison lives here too: with `--baseline <session>`, every cell is diffed against
the same cell in a previous run and only movements that clear the wider of the two cells' bands are
called changes.
"""

from __future__ import annotations

from pathlib import Path

from ..charts import Series, grouped_columns, heatmap, sequential_heatmap
from ..config import Suite
from ..records import NOISE_FLOOR_PCT, Record, delta_pct, load_dir
from ..report import Report, Table, band, pct, qps

#: How a cell is identified across runs, and the order the columns appear in.
CELL_KEYS = ("peptide_source", "equate_il", "tryptic", "mlp_batch", "kmer_k")


#: The order peptide files are reported in, whichever order they were measured in. It runs from the
#: whole-picture views to the length regimes, so a reader meets the summary before the detail:
#: `summary` is the cross-file overview, `mixed` the unbucketed 5..50 file, then the three buckets
#: short to long. A file not named here follows, alphabetically.
BUCKET_ORDER = ("summary", "mixed", "small", "medium", "large")


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    cells = _by_cell(loaded)
    arms = [arm.name for arm in suite.arms]

    _summary(report, cells, arms, loaded)

    for source in _ordered(sorted({key[0] for key in cells})):
        report.heading(f"{source} peptides", level=3)
        _bucket_table(report, cells, source, arms)
        _bucket_figures(report, cells, source, arms)

    if len(arms) == 2:
        report.note(
            f"`ratio` is {arms[1]}/{arms[0]} at the median. `noise` is the wider of the two cells'\n"
            f"own p10..p90 half-spreads. A ratio whose distance from 1.00 is smaller than `noise` —\n"
            f"or smaller than the {NOISE_FLOOR_PCT}% measured full-database floor — says the run\n"
            f"cannot separate the two backends for that configuration, not that they are equal."
        )

    baseline = getattr(suite, "baseline", None)
    if baseline:
        _regressions(report, cells, baseline, arms)

    if suite.notes:
        report.note(suite.notes)


def _ordered(sources: list[str]) -> list[str]:
    known = [name for name in BUCKET_ORDER if name in sources]
    return known + [name for name in sources if name not in BUCKET_ORDER]


def _held(report: Report, loaded: list[Record]) -> None:
    """States the tuning these cells ran at, splitting it into what was swept and what was held.

    Everything here is read out of the records — both the values used and the shipped defaults they
    are compared against. Nothing about the knobs is written down on this side, so a field added to
    `SearchTuning` appears here the first time it is measured, with no change to this file. That is
    the point of carrying `tuning_defaults` in the record: a report that hardcoded them would
    quietly start lying the day a default is re-tuned, which is precisely when it matters.

    "The defaults sweep" is a claim about every knob, not only the ones on the axes. A run with one
    of the held knobs overridden looks identical in the tables and is measuring something else.
    """
    tunings = [record.config.get("tuning") for record in loaded if record.config.get("tuning")]
    defaults = next((record.config.get("tuning_defaults") for record in loaded), None)
    if not tunings:
        return

    swept, held, overridden, mixed = [], [], [], []
    for field in sorted({key for tuning in tunings for key in tuning}):
        values = {tuning.get(field) for tuning in tunings}
        default = (defaults or {}).get(field)
        if len(values) > 1:
            swept.append(f"{field} ∈ {sorted(map(_fmt_tune, values))}")
            if default is not None and default not in values:
                mixed.append(f"{field} (shipped {_fmt_tune(default)} not among the swept values)")
            continue
        value = next(iter(values))
        held.append(f"{field}={_fmt_tune(value)}")
        if default is not None and value != default:
            overridden.append(f"{field}={_fmt_tune(value)} (shipped {_fmt_tune(default)})")

    if swept:
        report.para("Swept across the cells below: " + ", ".join(swept) + ".")
    if held:
        report.para("Held at one value in every cell: " + ", ".join(held) + ".")
    if overridden:
        report.warn(
            "this run is NOT at the shipped tuning — " + ", ".join(overridden) + ". Its numbers "
            "describe that tuning, not the defaults."
        )
    if mixed:
        report.para(
            "note: " + ", ".join(mixed) + " — this sweep therefore does not include the "
            "configuration that ships."
        )


def _fmt_tune(value) -> str:
    """TOML/JSON booleans read better lowercase, matching how the knob is written in Rust."""
    return str(value).lower() if isinstance(value, bool) else str(value)


def _summary(report: Report, cells: dict, arms: list[str], loaded: list[Record]) -> None:
    """What this version does at production defaults, before any of the grid detail."""
    report.heading("summary", level=3)
    _held(report, loaded)
    buckets = _ordered(sorted({key[0] for key in cells}))
    production = {
        arm: [
            (cells.get((source, True, False, 16, 5)) or {}).get(arm, {}).get("p50")
            for source in buckets
        ]
        for arm in arms
    }
    if not any(any(value for value in values) for values in production.values()):
        report.para("no cell at production defaults (5-mer table, MLP batch 16, equate_il on).")
        return

    report.chart(
        grouped_columns(
            buckets,
            [Series(arm, production[arm], slot) for slot, arm in enumerate(arms)],
            "Throughput at production defaults, per peptide length regime",
            unit=" qps",
        ),
        "Throughput at production defaults, per peptide length regime",
    )
    report.para(
        "Production defaults: 5-mer table, MLP batch 16, equate_il on, tryptic off. Everything "
        "below varies one of those and shows what it costs."
    )


def _bucket_table(report: Report, cells: dict, source: str, arms: list[str]) -> None:
    """The full grid for one peptide file. First, because it is the exhaustive view."""
    table = Table(
        headers=["equate_il", "tryptic", "mlp_batch", "kmer", *arms, "ratio", "noise"],
        aligns=["<", "<", "<", "<"] + [">"] * len(arms) + [">", ">"],
    )
    for key in sorted(key for key in cells if key[0] == source):
        per_arm = cells[key]
        values = [per_arm.get(arm) for arm in arms]
        table.row(
            key[1],
            key[2],
            "scalar" if key[3] == 1 else key[3],
            _kmer(key[4]),
            *(qps(value["p50"]) if value else "-" for value in values),
            _ratio(values),
            band(max((_band(value) for value in values if value), default=float("nan"))),
        )
    report.table(table)


def _bucket_figures(report: Report, cells: dict, source: str, arms: list[str]) -> None:
    """Four grids per peptide file — one per (equate_il, tryptic) — three ways to colour them.

    Splitting on the two search options rather than folding them into one grid's rows is what makes
    each grid answer a single question: within one search mode, what do the two accelerators do?
    The switch then chooses whether the cells show each backend on its own (one hue, absolute
    throughput, a shared scale across the four grids) or the comparison between them (diverging,
    with everything inside its noise floor left neutral).
    """
    keys = [key for key in cells if key[0] == source]
    if not keys:
        return
    modes = sorted({(key[1], key[2]) for key in keys})
    columns = sorted({(key[3], key[4]) for key in keys})
    column_labels = [f"mlp {'scalar' if batch == 1 else batch} · {_kmer(kmer)}" for batch, kmer in columns]

    # One scale across all four grids of a file, so a cell in one is comparable with a cell in
    # another. Per-grid scaling would make every grid look the same.
    absolute = [
        value["p50"]
        for key in keys
        for value in cells[key].values()
        if value and value["p50"]
    ]
    low, high = (min(absolute), max(absolute)) if absolute else (0.0, 1.0)

    variants: list[tuple[str, list[str]]] = []
    for arm in arms:
        variants.append(
            (
                arm,
                [
                    sequential_heatmap(
                        column_labels,
                        [arm],
                        {
                            (0, c): (cells[(source, il, tr, batch, kmer)][arm]["p50"],
                                     f"{arm} · il={il} tryptic={tr} · mlp_batch={batch} {_kmer(kmer)}: "
                                     f"{cells[(source, il, tr, batch, kmer)][arm]['p50']:,.0f} qps")
                            for c, (batch, kmer) in enumerate(columns)
                            if (source, il, tr, batch, kmer) in cells and arm in cells[(source, il, tr, batch, kmer)]
                        },
                        f"equate_il={il}, tryptic={tr}",
                        low=low,
                        high=high,
                        unit=" qps",
                    )
                    for il, tr in modes
                ],
            )
        )

    if len(arms) == 2:
        variants.append(("ratio", [_mode_heatmap(cells, source, arms, il, tr, columns, column_labels) for il, tr in modes]))

    report.switch(f"{source}: colour cells by", variants, default="ratio" if len(arms) == 2 else arms[0])


def _mode_heatmap(
    cells: dict, source: str, arms: list[str], il: bool, tryptic: bool, columns: list, column_labels: list[str]
) -> str:
    """One (equate_il, tryptic) grid, coloured by which backend is ahead."""
    grid = {}
    for c, (batch, kmer) in enumerate(columns):
        pair = cells.get((source, il, tryptic, batch, kmer))
        if not pair or not all(arm in pair for arm in arms) or not pair[arms[0]]["p50"]:
            continue
        base, other = pair[arms[0]], pair[arms[1]]
        difference = delta_pct(other["p50"], base["p50"])
        floor = max(_band(base), _band(other), NOISE_FLOOR_PCT)
        verdict = (
            f"{arms[1] if difference > 0 else arms[0]} ahead"
            if abs(difference) > floor
            else "within the noise floor — this run cannot separate them"
        )
        grid[(0, c)] = (
            difference,
            floor,
            f"il={il} tryptic={tryptic} · mlp_batch={batch} {_kmer(kmer)}: "
            f"{arms[0]} {base['p50']:,.0f} vs {arms[1]} {other['p50']:,.0f} qps "
            f"({difference:+.1f}%, floor {floor:.1f}%) — {verdict}",
        )
    return heatmap(
        column_labels,
        [f"{arms[1]} vs {arms[0]}"],
        grid,
        f"equate_il={il}, tryptic={tryptic}",
        pos_label=arms[1],
        neg_label=arms[0],
    )


def _regressions(report: Report, cells: dict, baseline_dir: Path, arms: list[str]) -> None:
    """Diffs every cell against the same cell in a previous session."""
    report.heading("regression check against the baseline", level=3)
    previous = _by_cell(load_dir(baseline_dir))
    if not previous:
        report.warn(f"no records under {baseline_dir} — nothing to compare against")
        return

    table = Table(
        headers=["file", "equate_il", "tryptic", "mlp_batch", "kmer", "arm", "base", "now", "delta", "verdict"],
        aligns=["<", "<", "<", ">", ">", "<", ">", ">", ">", "<"],
    )
    moved = 0
    for key in sorted(set(cells) & set(previous)):
        for arm in arms:
            now, base = cells[key].get(arm), previous[key].get(arm)
            if not (now and base):
                continue
            difference = delta_pct(now["p50"], base["p50"])
            floor = max(_band(now), _band(base), NOISE_FLOOR_PCT)
            changed = abs(difference) > floor
            moved += changed
            table.row(
                key[0],
                key[1],
                key[2],
                "scalar" if key[3] == 1 else key[3],
                _kmer(key[4]),
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



def _by_cell(loaded: list[Record]) -> dict[tuple, dict[str, dict]]:
    """(config key) -> arm -> {p10, p50, p90}.

    A key field is looked up in `config` first and in `config.tuning` second, so a coordinate that
    is a searcher knob (`mlp_batch`) is found in the same way as one that is not (`equate_il`) —
    and moving another field into `SearchTuning` later would not break this.
    """
    cells: dict[tuple, dict[str, dict]] = {}
    for record in loaded:
        config = record.config
        tuning = config.get("tuning", {})
        key = tuple(config.get(name, tuning.get(name)) for name in CELL_KEYS)
        spread = record.spread()
        p10, p50, p90 = spread if spread else (record.qps, record.qps, record.qps)
        arm = record.dims.get("arm", "?")
        cells.setdefault(key, {})[arm] = {"p10": p10, "p50": p50, "p90": p90}
    return cells


def _band(value: dict | None) -> float:
    if not value or not value["p50"]:
        return float("nan")
    return (value["p90"] - value["p10"]) / 2 / value["p50"] * 100


def _ratio(values: list[dict | None]) -> str:
    if len(values) != 2 or not all(values) or not values[0]["p50"]:
        return "-"
    return f"{values[1]['p50'] / values[0]['p50']:.2f}x"


def _kmer(k) -> str:
    return {0: "none", 5: "5-mer", 6: "6-mer"}.get(k, str(k))
