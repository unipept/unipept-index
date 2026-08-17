"""Reporting every suite needs, that no single suite owns.

Two things live here. The first is stating what tuning a set of records ran at: both the suite that
holds the knobs fixed and the suite that sweeps them have to say so, and they have to say it the
same way — a report where "held at the shipped value" and "overridden" are worded differently
depending on which suite produced it is a report whose reader has to learn two vocabularies.

The second is reading a matrix-mode grid. `defaults`, `kmer` and `mlp` are the same experiment
pointed at different coordinates: one process per arm, a grid swept inside it, one row per cell with
both arms side by side. What differs between them is which coordinate is the subject and what the
extra column says, which is little enough that each suite builds its own table — so what is shared
is the primitives underneath, not a table builder with three modes.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..charts import FACET_PLANE, Series, by_residency, heatmap, lines, sequential_heatmap
from ..records import NOISE_FLOOR_PCT, Record, delta_pct
from ..report import Report, Table, band, pct, qps

#: Every coordinate a grid cell can differ by, in the order they read best as columns.
#:
#: `amount_of_peptides` is in here because a block that runs a shorter query stream is not
#: comparable with one that runs the full one — the per-rep setup amortises differently — so two
#: cells that differ only by it must stay two cells rather than silently becoming one.
GRID_KEYS = ("peptide_source", "equate_il", "tryptic", "kmer_k", "mlp_batch", "amount_of_peptides")

#: Coordinates the configuration table reports, in reading order. A superset of `GRID_KEYS`: the
#: table names everything a cell ran at, including coordinates no suite currently varies, because
#: "what was this measured at" is a different question from "what distinguishes these cells".
CONTEXT_KEYS = ("peptide_source", "amount_of_peptides", "kmer_k", "equate_il", "tryptic", "max_matches")


#: How a coordinate's value is printed in the configuration table.
_COORD_LABELS = {
    "kmer_k": lambda k: {0: "none", 5: "5-mer", 6: "6-mer"}.get(k, str(k)),
    "mlp_batch": lambda b: "scalar" if b == 1 else str(b),
    "amount_of_peptides": lambda n: f"{n:,} queries",
    "peptide_source": str,
}

#: Coordinates that are the workload rather than a setting, listed first and never called
#: "overridden" — there is no shipped default for which peptide file a run queried. Nor are they
#: "swept": every suite reports per length regime, so varying the file is how the results are
#: SPLIT, not what the suite is asking about. Saying otherwise would put two subjects in every
#: summary and hide the real one.
_WORKLOAD_KEYS = ("peptide_source", "amount_of_peptides")

#: Length regimes short to long, matching the order the sections appear in.
_BUCKET_ORDER = ("summary", "mixed", "small", "medium", "large")


def held_and_swept(report: Report, loaded: list[Record]) -> None:
    """One table saying what this run was configured as, and which single thing it varied.

    Every suite past `defaults` is "one variable moved against a fixed background", and the reader's
    first question is always which variable and what the background was. Prose could not answer it
    at a glance: four sentences of `field=value, field=value` is where the shipped defaults used to
    live, and nobody reads to the end of it to notice that one knob is not at its default.

    So: one row per setting, the swept one marked, and a `shipped` column that makes an override
    impossible to miss. Everything is read out of the records — both the values used and the
    defaults they are compared against — so a field added to `SearchTuning` appears here the first
    time it is measured with no change to this file. That is what `tuning_defaults` is carried for:
    a report that hardcoded the defaults would start lying the day one is re-tuned, which is exactly
    when it matters.

    Split in two, because the two halves are read at different times. What the suite VARIES is the
    suite's subject and stays open; what it HOLDS is the background, and in a full report it is the
    same seven rows of `held at the shipped default` under every one of ten suites. Folding it keeps
    it one click away rather than seventy lines of the page — and an OVERRIDDEN row is never folded,
    since a run that is not at the shipped tuning has to say so where it cannot be missed.
    """
    configs = [record.config for record in loaded if record.config]
    if not configs:
        return
    defaults = next((config.get("tuning_defaults") for config in configs if config.get("tuning_defaults")), {})

    rows: list[tuple[str, set, object]] = []
    for key in _WORKLOAD_KEYS + tuple(k for k in CONTEXT_KEYS if k not in _WORKLOAD_KEYS):
        values = {config.get(key) for config in configs if config.get(key) is not None}
        if values:
            rows.append((key, values, None))
    for field in sorted({key for config in configs for key in (config.get("tuning") or {})}):
        values = {(config.get("tuning") or {}).get(field) for config in configs}
        rows.append((field, {value for value in values if value is not None}, defaults.get(field)))

    def new_table() -> Table:
        return Table(
            headers=["setting", "this run", "shipped", "role"],
            aligns=["<", "<", ">", "<"],
        )

    varies, held = new_table(), new_table()
    overridden, mixed = [], []
    for name, values, default in rows:
        label = _COORD_LABELS.get(name, fmt_tune)
        order = _bucket_order(values) if name == "peptide_source" else _in_order(values)
        shown = ", ".join(label(value) for value in order)
        swept = False
        if len(values) > 1:
            workload = name in _WORKLOAD_KEYS
            swept = not workload
            role = "a panel of the figure each" if workload else "SWEPT — the subject of this suite"
            if default is not None and default not in values:
                mixed.append(name)
                role = "SWEPT, and the shipped value is not among them"
        elif default is None:
            role = "fixed"
        elif next(iter(values)) == default:
            role = "held at the shipped default"
        else:
            role = "OVERRIDDEN"
            overridden.append(name)
        # A row belongs above the fold if it is what this suite moved, or if it is at a value the
        # run was not supposed to be at. Everything else is background.
        subject = len(values) > 1 or role == "OVERRIDDEN"
        table = varies if subject else held
        # Bold on the page. This table's whole job is to answer "what does this suite vary", and a
        # role column that has to be read row by row answers it more slowly than a weight does.
        if swept or role == "OVERRIDDEN":
            table.strong.add(len(table.rows))
        table.row(name, shown, "-" if default is None else fmt_tune(default), role)
    if varies.rows:
        report.table(varies)
    if held.rows:
        report.table(held, raw=f"held at the shipped defaults ({len(held.rows)} settings)")

    if overridden:
        report.warn(
            "this run is NOT at the shipped tuning — " + ", ".join(overridden) + " was overridden. "
            "Its numbers describe that tuning, not the defaults."
        )
    if mixed:
        report.para(
            f"note: the swept values of {', '.join(mixed)} do not include the shipped one, so this "
            f"sweep does not measure the configuration that ships."
        )


def fmt_tune(value) -> str:
    """TOML/JSON booleans read better lowercase, matching how the knob is written in Rust."""
    return str(value).lower() if isinstance(value, bool) else str(value)


def _bucket_order(values: set) -> list:
    """Length regimes short to long; anything the order does not name follows alphabetically."""
    known = [name for name in _BUCKET_ORDER if name in values]
    return known + sorted(value for value in values if value not in _BUCKET_ORDER)


def _in_order(values: set):
    """Numeric knob values in numeric order, anything else alphabetically.

    Every knob so far is a batch size or a distance, and sorting those as strings puts 128 between
    1 and 16 — which reads as a shuffled list rather than as the ladder that was actually swept.
    """
    try:
        return sorted(values, key=float)
    except (TypeError, ValueError):
        return sorted(values, key=str)


# ---------------------------------------------------------------------------
# Reading a matrix grid
# ---------------------------------------------------------------------------


def by_cell(
    loaded: list[Record],
    keys: tuple[str, ...] = GRID_KEYS,
    *,
    correct_drift: bool = False,
) -> dict[tuple, dict[str, dict]]:
    """`(coordinate tuple) -> arm -> {p10, p50, p90, search_ms, retrieval_ms, floor}`.

    Cells are keyed on `config`, not on `dims`: a matrix invocation sweeps the grid *inside* one
    process, so every record from one arm shares that arm's dims and only the config tells them
    apart.

    Each key is looked up in `config` first and in `config.tuning` second, so a coordinate that is a
    searcher knob (`mlp_batch`) is found the same way as one that is not (`equate_il`) — and moving
    a field into `SearchTuning` later would not break this.

    The two phase timings come along because throughput alone cannot say WHERE a configuration
    spends itself. In matrix mode `result` is the median-qps rep, so these describe that rep rather
    than the pooled median — close enough to read the split from, and the only phase numbers an
    uninstrumented build produces at all.

    With `correct_drift`, throughput is rescaled against the reference cell the grid interleaved —
    see `drift_of`. Suites whose blocks set `base_every` want this; a suite without a cadence gets
    the numbers as measured either way, so passing it is harmless.
    """
    drift = drift_by_process(loaded) if correct_drift else {}
    cells: dict[tuple, dict[str, dict]] = {}
    for index, record in _positions(loaded):
        config = record.config
        if config.get("sweep") == DRIFT:
            # The reference repeats are the correction, not a measurement of anything.
            continue
        tuning = config.get("tuning", {})
        key = tuple(config.get(name, tuning.get(name)) for name in keys)
        spread = record.spread()
        p10, p50, p90 = spread if spread else (record.qps, record.qps, record.qps)
        result = record.result
        series = drift.get(_process_key(record))
        scale = series.scale_at(index) if series else 1.0
        cells.setdefault(key, {})[record.dims.get("arm", "?")] = {
            "p10": p10 * scale,
            "p50": p50 * scale,
            "p90": p90 * scale,
            "search_ms": _millis(result.get("search_duration_ns")),
            "retrieval_ms": _millis(result.get("retrieval_duration_ns")),
            # Only present where a cell opted into the response phase; None everywhere else, which
            # is what keeps them out of the stacked chart for suites that did not measure them.
            "decode_ms": _millis(result.get("decode_duration_ns")),
            "serialise_ms": _millis(result.get("serialise_duration_ns")),
            "response_bytes": result.get("response_bytes") or None,
            # What the process this cell came from could resolve at all, however tight this one
            # cell's own reps happened to be. `floor_of` folds it in.
            "floor": series.floor if series else float("nan"),
            # Kept so a caller can tell how far off the shipped tuning a cell is. A suite sweeping
            # one knob has to exclude the cells a PLANE block moved two knobs in, or its curve is
            # several curves overlaid and its "context" is whatever the other knob happened to be.
            "off_default": tuple(
                sorted(name for name, value in tuning.items() if (config.get("tuning_defaults") or {}).get(name) != value)
            ),
        }
    return cells


# ---------------------------------------------------------------------------
# Drift
# ---------------------------------------------------------------------------

#: The `sweep` name `grid._with_drift_cadence` gives the interleaved reference cells.
DRIFT = "drift"

#: What makes two records the same process. A rayon pool is built once per process and a cgroup
#: scope wraps one, so these are also exactly the coordinates that cannot change while it lives.
PROCESS_DIMS = ("arm", "threads", "ceiling_gb")


@dataclass
class DriftSeries:
    """One process's reference measurements, and what they say about what it could resolve.

    A process sweeping tens of cells takes long enough that the machine moves under it, and a single
    reference at the start cannot tell drift from effect. Measuring the reference repeatedly turns
    drift into a series: subtract its trend and what is left is the floor on what that process can
    say about anything.
    """

    #: (position in the process's record order, throughput) for each reference measurement.
    marks: list[tuple[int, float]]
    #: First mark against last, as a percentage. How far the machine moved over the whole process.
    change: float
    #: Worst deviation from the straight line through the marks, as a percentage.
    residual: float

    @property
    def floor(self) -> float:
        """The residual against the measured full-database floor, whichever is wider."""
        if self.residual != self.residual:
            return float("nan")
        return max(self.residual, NOISE_FLOOR_PCT)

    def scale_at(self, index: int) -> float:
        """What to multiply a cell at this position by to put it on the process's opening footing.

        Rescaled to the reference at the START of the process rather than divided by the local one:
        the correction has to follow the machine, whose speed is a function of when a cell ran, but
        the result still has to be a throughput. A dimensionless ratio would be unreadable beside
        every other table in this package.
        """
        reference = self.reference_at(index)
        anchor = self.marks[0][1] if self.marks else None
        return anchor / reference if reference and anchor else 1.0

    def reference_at(self, index: int) -> float | None:
        """The reference throughput at a position, interpolated between the marks around it.

        None when there is no cadence, in which case throughput is reported as measured and the
        report says the run has no drift correction rather than pretending to one.
        """
        marks = self.marks
        if not marks:
            return None
        if len(marks) == 1 or index <= marks[0][0]:
            return marks[0][1]
        if index >= marks[-1][0]:
            return marks[-1][1]
        for (left, low), (right, high) in zip(marks, marks[1:]):
            if left <= index <= right:
                span = right - left
                return low + (high - low) * ((index - left) / span) if span else low
        return marks[-1][1]


def _process_key(record: Record) -> tuple:
    return tuple(record.dims.get(name, "") for name in PROCESS_DIMS)


def _positions(loaded: list[Record]):
    """`(position within its process, record)`, in the order the process wrote them.

    Order is load-bearing here in a way it is nowhere else in this package: the cadence is only
    readable as a series, and a process's records are exactly one jsonl file written cell by cell as
    the sweep ran.
    """
    seen: dict[tuple, int] = {}
    for record in loaded:
        key = _process_key(record)
        index = seen.get(key, 0)
        seen[key] = index + 1
        yield index, record


def drift_by_process(loaded: list[Record]) -> dict[tuple, DriftSeries]:
    """Every process's reference series, keyed by `PROCESS_DIMS`."""
    marks: dict[tuple, list[tuple[int, float]]] = {}
    for index, record in _positions(loaded):
        if record.config.get("sweep") != DRIFT:
            continue
        spread = record.spread()
        marks.setdefault(_process_key(record), []).append((index, spread[1] if spread else record.qps))
    return {key: DriftSeries(found, *_drift_stats(found)) for key, found in marks.items()}


def _drift_stats(marks: list[tuple[int, float]]) -> tuple[float, float]:
    """(first-to-last change, worst deviation from the straight line through them), both percent.

    The residual is the WORST deviation, not the average one: the number is used to refuse to call
    something an effect, and an average would let a single badly-behaved stretch of the process hide
    inside a calm one.
    """
    if len(marks) < 2:
        return float("nan"), float("nan")
    first, last = marks[0][1], marks[-1][1]
    change = (last - first) / first * 100 if first else float("nan")

    span = marks[-1][0] - marks[0][0]
    worst = 0.0
    for index, value in marks:
        fitted = first + (last - first) * ((index - marks[0][0]) / span) if span else first
        if fitted:
            worst = max(worst, abs(value - fitted) / fitted * 100)
    return change, worst


def resolution_table(report: Report, loaded: list[Record]) -> None:
    """Per process: how far the machine drifted, and what is left once that is removed.

    Printed before any result, because it bounds every one of them. A process whose residual exceeds
    the effects below it has not measured them.
    """
    series = drift_by_process(loaded)
    if not series:
        return
    counts: dict[tuple, int] = {}
    for record in loaded:
        if record.config.get("sweep") != DRIFT:
            key = _process_key(record)
            counts[key] = counts.get(key, 0) + 1

    table = Table(
        headers=["arm", "threads", "ceiling", "cells", "marks", "drift", "residual", "floor"],
        aligns=["<", ">", ">", ">", ">", ">", ">", ">"],
        tips=tips_for(["arm", "drift", "residual", "floor"]),
    )
    for key in sorted(series, key=str):
        found = series[key]
        table.row(
            key[0] or "-",
            key[1] or "default",
            key[2] or "0",
            counts.get(key, 0),
            len(found.marks),
            pct(found.change) if found.marks else "no cadence",
            band(found.residual) if found.marks else "-",
            band(found.floor if found.floor == found.floor else NOISE_FLOOR_PCT),
        )
    report.table(table)
    # The four definitions this paragraph used to spell out now hang off the column headings that
    # need them, which is where a reader is when the question arises — and which stops the same 450
    # characters appearing under six suites in one report.
    report.para(
        f"Drift is measured and removed, not merely noted; what is left is each process's own "
        f"floor, against the {NOISE_FLOOR_PCT}% measured full-database floor. **Nothing below a "
        f"process's floor is a result from that process.** Hover any heading for its definition."
    )


def _millis(nanoseconds: float | None) -> float | None:
    """One rep's phase duration in milliseconds, or None when the counter is absent.

    Per rep rather than per query: at full-database speed one query's search is tens of
    microseconds and on swissprot under one, so a per-query figure is a column of leading zeros.
    Every cell in a suite queries the same stream length — the configuration table states it, and
    flags it as swept if it ever differs — so the reps are directly comparable as they stand.
    """
    return None if not nanoseconds else nanoseconds / 1e6


def varying(cells: dict[tuple, dict], keys: tuple[str, ...] = GRID_KEYS) -> list[str]:
    """The keys that actually take more than one value in these cells.

    A grid narrowed in its suite file should narrow its table too: a column that reads `16` in every
    row is not a coordinate, it is a caption, and it belongs in the prose above the table where it
    is stated once. `held_and_swept` already prints the tuning half of that.
    """
    return [
        name
        for index, name in enumerate(keys)
        if len({key[index] for key in cells}) > 1
    ]


def cell_band(value: dict | None) -> float:
    """Half a cell's p10..p90 spread, as a percentage of its median."""
    if not value or not value["p50"]:
        return float("nan")
    return (value["p90"] - value["p10"]) / 2 / value["p50"] * 100


def floor_of(*values: dict | None) -> float:
    """What a delta between these cells has to clear before it is an answer.

    The widest of: each cell's own band, the drift residual of the process each came from, and the
    measured full-database floor. Taking the widest is the same rule `records.noise_floor` applies,
    for the same reason — the point is to refuse to call something an effect, not to find the
    reading under which it qualifies. The residual belongs here and not only in the resolution
    table, because a process that could not hold still is one whose effects are unreadable however
    tight an individual cell's reps happened to be.
    """
    spreads = [
        value
        for cell in values
        if cell
        for value in (cell_band(cell), cell.get("floor", float("nan")))
    ]
    return max([NOISE_FLOOR_PCT, *(value for value in spreads if value == value)])


def ratio_of(values: list[dict | None]) -> str:
    """Second arm over first at the median, for the two-arm tables."""
    if len(values) != 2 or not all(values) or not values[0]["p50"]:
        return "-"
    return f"{values[1]['p50'] / values[0]['p50']:.2f}x"


def kmer_label(k) -> str:
    return {0: "none", 5: "5-mer", 6: "6-mer"}.get(k, str(k))


def batch_label(batch) -> str:
    """`mlp_batch = 1` is not a batch size, it is the scalar path, and reads better said so."""
    return "scalar" if batch == 1 else str(batch)


#: The arm a composition figure is drawn for, when only one can be. The deployed configuration, so
#: the split shown is the split production actually pays.
DEPLOYED_ARM = "pprot"

#: The k-mer size that ships, mirrored from `sa-builder`'s `--kmer-size` default.
#:
#: Every `SearchTuning` field is read out of the binary, so the driver cannot drift from what ships.
#: The table is the exception: it is a build-time artefact chosen by `sa-builder`, not a searcher
#: knob, so nothing in a record says which k is production. It is written here instead of at each
#: use, so changing the shipped table is one edit rather than a hunt — and the `kmer = [...]` lines
#: in `suites/*.toml` have to be moved with it, since those are what actually pin the background.
SHIPPED_KMER_K = 6


def phase_switch(
    report,
    panels: list[tuple[str, list[str], object]],
    arms: list[str],
    x_title: str,
    *,
    title: str = "by length regime",
    throughput: bool = True,
    default_reading: str | None = None,
) -> None:
    """The four readings of a suite, each as one small-multiple grid over the length regimes.

    `panels` is `(panel title, that panel's x groups, at)`, where `at(group_index, arm)` returns
    that cell's measurement dict — the only part each suite has to supply. `defaults`, `kmer` and
    `mlp` differ only in what sits on the x axis, so the figure is built once here rather than three
    times.

    This used to emit one switch PER regime, four figures deep, which is where 176 of the report's
    204 figures came from — the same chart drawn four times, each scaled to its own maximum, so the
    one comparison the split invites was the one it made impossible. Now the regimes are the panels
    of one grid on one scale, and the four readings stay what they were: a switch, because they are
    the same cells painted differently.

    The phase split is the reading throughput cannot give. A k-mer table shortens the binary
    search's probe chain and should move the SEARCH bar alone; `tryptic` changes how many candidates
    survive and should move RETRIEVAL. A change that lands in the wrong phase is measuring something
    other than what it claims to.

    `default_reading` names the tab that opens. Throughput otherwise, which is the reading every
    knob suite is about — `defaults` is the exception, because it is the only suite that measures
    decode and serialisation at all and the `time split` tab is where they are drawn.
    """
    from ..charts import Series, facets, grouped_columns, stacked_columns

    def per_arm(at, groups: list[str], field: str) -> list[Series]:
        return [
            Series(
                arm,
                [(at(index, arm) or {}).get(field) for index in range(len(groups))],
                arm=arm,
                tip={"backend": arm},
            )
            for arm in by_residency(arms)
        ]

    def column_grid(field: str, label: str, unit: str, y_title: str) -> list[str]:
        built = [(name, per_arm(at, groups, field)) for name, groups, at in panels]
        groups_of = {name: groups for name, groups, _ in panels}
        return facets(
            built,
            lambda name, series, frame, top, legend: grouped_columns(
                groups_of[name],
                series,
                name,
                unit=unit,
                x_title=x_title,
                y_title=y_title,
                frame=frame,
                y_max=top,
                legend=legend,
            ),
            # No axes note: every panel names its own axes now, and the note repeated them
            # verbatim above the grid.
        )

    # Throughput leads every time. It used to be dropped for suites whose curve grid already plotted
    # it, on the grounds that the curve says the same thing — but the curve is NORMALISED to the
    # shipped value, so it answers "did this knob help" and cannot answer "how fast was this cell",
    # and the reader who wants the second question had no picture of it at all. Two readings of one
    # measurement is what the switch is for.
    readings = [
        ("search_ms", "search time", " ms", "search time (ms)"),
        ("retrieval_ms", "retrieval time", " ms", "retrieval time (ms)"),
    ]
    if throughput:
        readings.insert(0, ("p50", "throughput", " qps", "throughput (qps)"))

    variants: list[tuple[str, list[str]]] = []
    for field, label, unit, y_title in readings:
        grid = column_grid(field, label, unit, y_title)
        if grid:
            variants.append((label, grid))

    # What a rep is made of. Search and retrieval are always there; decode and serialisation only
    # where a cell opted in, and those two are the half of a real request nothing else measures —
    # so where they ARE present the stack finally shows what share of a request the rest of the
    # report describes.
    #
    # Drawn for ONE arm. Composition is a property of the workload rather than of the storage
    # backend — the three arms' phase splits differ by less than the noise on any of them — so
    # putting all three in would spend a second categorical channel, on top of the four phases, to
    # redraw the same shape three times. Which backend is faster is what the other three readings
    # are for.
    split_arm = DEPLOYED_ARM if DEPLOYED_ARM in arms else (arms[0] if arms else "")
    phases = (
        ("search_ms", "search"),
        ("retrieval_ms", "retrieval"),
        ("decode_ms", "decode"),
        ("serialise_ms", "serialise"),
    )

    def split_panel(name, series, frame, top, legend):
        groups = next(groups for panel, groups, _ in panels if panel == name)
        return stacked_columns(
            groups, [split_arm], series, name,
            unit=" ms", frame=frame, y_max=top, legend=legend, share=True,
            x_title=x_title, y_title="share of time per rep (%)",
        )

    built_split = []
    for name, groups, at in panels:
        parts = [
            Series(label, [(at(index, split_arm) or {}).get(field) for index in range(len(groups))], slot)
            for slot, (field, label) in enumerate(phases)
        ]
        parts = [item for item in parts if any(value for value in item.values)]
        if parts:
            built_split.append((name, parts))
    if built_split:
        grid = facets(
            built_split,
            split_panel,
            axes=f"drawn for the {split_arm} arm — composition barely differs between them",
            # Every column is normalised to itself, so the scale is 100% by construction rather than
            # taken from the data — `stack_max` would hand back the tallest raw total in ms and
            # scale a chart whose columns are already percentages to it.
            extent=lambda series: 100.0,
        )
        if grid:
            variants.append(("time split", grid))

    if variants:
        offered = {label for label, _ in variants}
        # A named tab that this suite did not build is a caller error worth failing on, not a silent
        # fallback: it would open on throughput and read as if nothing had been asked for.
        if default_reading is not None and default_reading not in offered:
            raise ValueError(
                f"phase_switch: default_reading={default_reading!r} is not one of {sorted(offered)}"
            )
        report.switch(title, variants, default=default_reading or variants[0][0])


#: What each recurring column means, shown on hovering its header.
#:
#: One glossary rather than a paragraph per table. These four words appear in six suites and used to
#: be re-explained under each one, which meant six wordings of the same definition and a reader who
#: had to scroll away from the number to find the one that applied. Beside the column is where a
#: column's definition belongs.
COLUMN_TIPS = {
    "ratio": (
        "The second backend's throughput divided by the first, at the median. "
        "1.00x means they measured the same; below 1.00 the second one is slower. "
        "Only meaningful once its distance from 1.00 exceeds the noise column."
    ),
    "noise": (
        "Half the p10..p90 spread of the noisier of the two cells, as a percentage of its median. "
        "How steady the measurement was. A difference smaller than this is not a small effect — "
        "it is no answer from this run."
    ),
    "band": (
        "Half this cell's own p10..p90 spread, as a percentage of its median. "
        "How steady this one cell was across its reps."
    ),
    "floor": (
        "What a difference has to clear before it counts: the wider of the two cells' own bands, "
        f"and never below the {NOISE_FLOOR_PCT}% run-to-run noise floor measured on the full "
        "database. Taking the widest is deliberate — the point is to refuse to call something an "
        "effect, not to find the reading under which it qualifies."
    ),
    "gain": (
        "The best value's throughput against the shipped value's, as a percentage. "
        "Read it against the floor beside it: a gain inside the floor is not a gain."
    ),
    "qps": "Queries per second, at the median of this cell's timed reps.",
    "table GB": (
        "The k-mer table's own resident size, computed from k rather than measured: the table is "
        "dense, 24^k entries of 16 bytes. Process RSS cannot answer this, because one matrix "
        "process loads every k the grid names before it sweeps."
    ),
    "vs none": "This cell against the same cell with no k-mer table attached.",
    "vs shipped": "This cell against the same cell at the batch size the server actually passes.",
    "verdict": "What the delta and the floor beside it add up to, stated so it cannot be skimmed past.",
    "arm": "Which storage build produced this row.",
    "kmer": "The k-mer table attached for this row. `none` is the reference the others are read against.",
    "mlp_batch": (
        "How many independent peptide searches are interleaved per rayon task. "
        "1 is the scalar path; the shipped value is what the server passes on every request."
    ),
    "resolved": (
        "How many of this knob's contexts produced a difference clearing their own floor. "
        "The verdict is decided on these, not on the raw argmaxes in the winners column."
    ),
    "drift": (
        "The reference cell's first mark against its last — how far the machine moved over the "
        "whole process. Every throughput in this suite is already rescaled by that reference "
        "interpolated to its own position, so drift is removed here, not merely noted."
    ),
    "residual": (
        "What the reference still scatters by once its trend is removed — this process's own floor. "
        "An effect below it is not a small effect; it is no answer from this process."
    ),
}


def tips_for(headers) -> dict[str, str]:
    """The glossary entries that apply to one table's headers."""
    return {header: COLUMN_TIPS[header] for header in headers if header in COLUMN_TIPS}


# ---------------------------------------------------------------------------
# Knob suites
# ---------------------------------------------------------------------------

#: Length regimes short to long; anything else follows alphabetically.
BUCKET_ORDER = ("summary", "mixed", "small", "medium", "large")

#: Per-value formatting for the knobs that read better as words than as numbers.
KNOB_LABELS = {"mlp_batch": batch_label, "kmer_k": kmer_label}

#: How a context coordinate reads inside a compound label. `equate_il=false · tryptic=false ·
#: amount_of_peptides=10000` is three facts and sixty characters, and it prefixes every row of
#: every table in a suite that sweeps the search options — so each one gets a word instead.
CONTEXT_LABELS = {
    # Booleans stay booleans. `il` / `no-tryptic` read as jargon and made two columns that mean
    # opposite things look alike; `true` / `false` under a column already named `tryptic` says the
    # same thing in the reader's own words, and the page paints them as coloured pills.
    "amount_of_peptides": lambda value: f"{value:,}q",
    "kmer_k": kmer_label,
    "mlp_batch": lambda value: f"mlp {batch_label(value)}",
    "validate_batch": lambda value: f"validate {value}",
}


def context_label(name: str, value) -> str:
    formatter = CONTEXT_LABELS.get(name)
    return formatter(value) if formatter else f"{name}={fmt_tune(value)}"


def column_label(name: str, value) -> str:
    """The same value, for somewhere the coordinate is ALREADY named — a table column, a tip key.

    `context_label` writes `tryptic=true` because it joins coordinates into one string where nothing
    else would say which is which. Under a header that reads `tryptic` it says it twice, and it is
    not one of the spellings the page paints as a pill or that the search-mode control can drive.
    The formatters stay: `5-mer` is the value of `kmer_k`, not a restatement of its name.
    """
    formatter = CONTEXT_LABELS.get(name)
    return formatter(value) if formatter else fmt_tune(value)


def order_buckets(sources) -> list[str]:
    known = [name for name in BUCKET_ORDER if name in sources]
    return known + sorted(name for name in sources if name not in BUCKET_ORDER)


def shipped_value(loaded: list[Record], knob: str, fallback=None):
    """What the binary reports as this knob's default, so no suite restates it.

    Read out of `tuning_defaults` on the records rather than written down here: a report that
    hardcoded a default would start lying the day it is re-tuned, which is exactly when a re-tuning
    suite is being read.

    A coordinate OUTSIDE `SearchTuning` — `amount_of_peptides`, `max_matches` — has no entry there
    and falls through to `fallback`, which is the first swept value and is usually the wrong
    reference. Such a suite should pass `reference` to `knob_analysis` explicitly.
    """
    for record in loaded:
        value = (record.config.get("tuning_defaults") or {}).get(knob)
        if value is not None:
            return value
    return fallback


def knob_analysis(
    report: Report,
    suite,
    loaded: list[Record],
    *,
    knob: str,
    x_title: str = "",
    mechanism: str = "",
    reference=None,
) -> None:
    """The whole report for a suite whose subject is one `SearchTuning` knob.

    `mlp`, `validate` and any knob suite after them are the same analysis pointed at a different
    field, so it lives here once: the configuration table, the resolution table, a curve per context
    normalised to the shipped value, a verdict per context, and a folded section per length regime
    with the phase split.

    What each suite supplies is the knob's name and its mechanism. Nothing here knows the field list
    — a knob added to `SearchTuning` is analysable the day a suite names it.
    """
    label = KNOB_LABELS.get(knob, fmt_tune)
    keys = GRID_KEYS if knob in GRID_KEYS else GRID_KEYS + (knob,)
    index = keys.index(knob)
    cells = by_cell(loaded, keys, correct_drift=True)

    values = _in_order({key[index] for key in cells if key[index] is not None})
    if len(values) < 2:
        report.heading("summary", level=3)
        held_and_swept(report, loaded)
        report.warn(f"only one {knob} value ran — there is no curve to read.")
        return
    # `reference` is what every curve is read against. A `SearchTuning` field has one on the record
    # and needs nothing here; a coordinate that is not a knob — the query count, the match cutoff —
    # has no shipped value at all, and falling back to the first swept value picks the smallest.
    # That is how `stream` came to report "+2229.7%" against a ten-peptide call.
    shipped = reference if reference is not None else shipped_value(loaded, knob, values[0])
    if shipped not in values:
        shipped = values[-1] if reference is not None else values[0]

    # Everything that varies besides the knob and the peptide file, which is what sections the
    # report. A suite crossing its knob with `tryptic` gets two lines per regime rather than one
    # silently overwriting the other.
    # A knob the suite also crosses in a plane block varies in `cells` but is not a context — it is
    # the other axis of a figure this analysis does not draw. Excluded from both.
    swept_elsewhere = {
        name
        for per_arm in cells.values()
        for cell in per_arm.values()
        for name in cell.get("off_default", ())
        if name != knob
    }
    extra = [
        name
        for name in varying(cells, keys)
        if name not in (knob, "peptide_source") and name not in swept_elsewhere
    ]
    points = knob_points(cells, keys, index, extra, knob=knob)
    arms = [arm.name for arm in suite.arms]
    sources = order_buckets({source for source, _ in points})

    # The answer first, then how it was reached. The verdict table is built before anything is
    # emitted so its aggregate can lead the suite rather than close it.
    verdict_table, outcomes = _knob_verdicts(points, sources, arms, label, knob, shipped, extra)
    _knob_verdict_tiles(report, outcomes, label, knob, shipped)

    report.heading("summary", level=3)
    held_and_swept(report, loaded)
    resolution_table(report, loaded)
    _knob_curve(report, points, values, sources, arms, label, knob, x_title, shipped, extra)
    report.table(verdict_table)
    if mechanism:
        report.para(mechanism)

    # The working, in one folded section rather than one per regime. Every panel of the grid is a
    # (regime, context) pair on one scale, and the table below carries every regime with a `file`
    # chip — which is what turns "show me every tryptic row" from four scrolls into one click.
    report.heading("per cell", level=3, folded=True)
    x_labels = [label(value) for value in values]
    panels = [
        (
            f"{source}{' · ' + context_text(extra, context) if context_text(extra, context) else ''}",
            x_labels,
            (lambda i, arm, p=per_arm, v=values: p.get(arm, {}).get(v[i])),
        )
        for source in sources
        for context, per_arm in sorted(((c, p) for (s, c), p in points.items() if s == source), key=str)
    ]
    phase_switch(report, panels, arms, x_title or knob, title="per cell")
    _knob_table(report, points, sources, values, arms, label, knob, shipped, extra)

    if suite.notes:
        report.note(suite.notes)


def knob_points(cells: dict, keys: tuple, index: int, extra: list[str], knob: str | None = None) -> dict:
    """`(source, context label) -> arm -> {swept value: cell}`.

    With `knob`, cells that move any OTHER tuning field off its default are left out. A suite that
    sweeps one knob and also sweeps a plane has both kinds of cell in its records, and folding them
    together gives a curve that is several curves overlaid — the point at `validate_batch=64` would
    be whichever `mlp_batch` the plane happened to visit last.
    """
    points: dict[tuple[str, tuple], dict[str, dict]] = {}
    for key, per_arm in cells.items():
        # A TUPLE of the coordinate values, not a joined string: every one of them becomes a column
        # of its own, so the page's chips can filter on each independently. A single `context`
        # column reading `il · no-tryptic · 10,000q` is one chip group offering every combination
        # that occurred, which is the opposite of a filter.
        context = tuple(key[keys.index(name)] for name in extra)
        for arm, cell in per_arm.items():
            if knob and any(name != knob for name in cell.get("off_default", ())):
                continue
            points.setdefault((key[0], context), {}).setdefault(arm, {})[key[index]] = cell
    return points


def context_text(extra: list[str], values: tuple) -> str:
    """The same coordinates as one string, for a chart legend or a verdict row."""
    return " · ".join(context_label(name, value) for name, value in zip(extra, values))


def _knob_curve(report, points, values, sources, arms, label, knob, x_title, shipped, extra) -> None:
    """One panel per (regime, context); within a panel, one line per arm.

    Normalised to the shipped value because the length regimes differ by two orders of magnitude in
    absolute throughput and the question is the SHAPE of each curve, not which regime is fastest —
    `defaults` answers that. 100% is the value that ships and is drawn as a reference rule, so a
    knob that bought nothing is a flat line lying on it.

    Faceted because this was the report's worst chart: every (regime, context, arm) triple on one
    axes is twenty-four lines, three times the eight-slot ceiling the palette can distinguish, and
    the two paired lightnesses had to cover three arms — so a third of those lines were drawn in a
    colour already used. Eight panels of three is the tier where colour alone is comfortable for
    everyone, and the panels share one scale, so their shapes stay comparable.
    """
    from ..charts import facets, panel_min

    panels = []
    for source in sources:
        for (bucket, context), per_arm in sorted(points.items(), key=str):
            if bucket != source:
                continue
            label_text = context_text(extra, context)
            series = []
            for arm in by_residency(arms):
                curve = per_arm.get(arm, {})
                reference = curve.get(shipped)
                if not reference or not reference["p50"]:
                    continue
                series.append(
                    Series(
                        arm,
                        [
                            curve[value]["p50"] / reference["p50"] * 100 if value in curve else None
                            for value in values
                        ],
                        arm=arm,
                        tip={
                            "peptides": source,
                            **dict(zip(extra, (column_label(n, v) for n, v in zip(extra, context)))),
                            "backend": arm,
                        },
                    )
                )
            if series:
                panels.append((f"{source}{' · ' + label_text if label_text else ''}", series))
    if not panels:
        return

    x_labels = [label(value) for value in values]
    caption = f"Throughput as a percentage of {knob}={label(shipped)}"
    report.figures(
        facets(
            panels,
            lambda name, series, frame, top, legend, bottom: lines(
                x_labels, series, name, unit="%", frame=frame,
                x_title=x_title or knob, y_title=f"% of {knob}={label(shipped)}",
                y_max=top, y_min=bottom, legend=legend, baseline=100.0,
            ),
            axes=f"dashed rule: {knob}={label(shipped)}, the value that ships",
            floor=panel_min,
        ),
        caption,
    )


def _knob_verdicts(points, sources, arms, label, knob, shipped, extra) -> tuple[Table, list[tuple]]:
    """Per context: where the peak is, and whether the curve is a curve at all.

    Returns the table and the per-context outcomes behind it, rather than emitting — the aggregate
    of those outcomes is the suite's headline and has to be placed above everything, while the table
    itself belongs below the curve it summarises.

    One column per swept coordinate rather than one joined `context` column, so the page's chips
    filter on each independently — "show me every tryptic row" is a filter; "show me the
    `il · tryptic · 10,000q` rows" is a lookup.
    """
    headers = ["file", *extra, "arm", "shipped", "best", "gain", "floor", "reading"]
    table = Table(
        headers=headers,
        aligns=["<"] * (len(extra) + 2) + [">"] * 4 + ["<"],
        chips=["file", *extra, "arm"],
        tips=tips_for(headers),
    )
    #: One entry per context: (winning value, gain %, floor %, did it clear its floor).
    outcomes: list[tuple] = []
    for source in sources:
        for (bucket, context), per_arm in sorted(points.items(), key=str):
            if bucket != source:
                continue
            for arm in arms:
                curve = {value: cell for value, cell in per_arm.get(arm, {}).items() if cell and cell["p50"]}
                reference = curve.get(shipped)
                if not curve or not reference:
                    continue
                best = max(curve, key=lambda value: curve[value]["p50"])
                difference = delta_pct(curve[best]["p50"], reference["p50"])
                floor = floor_of(curve[best], reference)
                resolved = [
                    value
                    for value, cell in curve.items()
                    if abs(delta_pct(cell["p50"], reference["p50"])) > floor_of(cell, reference)
                ]
                outcomes.append((best, difference, floor, abs(difference) > floor and best != shipped))
                table.row(
                    source,
                    *(column_label(name, value) for name, value in zip(extra, context)),
                    arm,
                    label(shipped),
                    label(best),
                    pct(difference),
                    band(floor),
                    _knob_reading(best, shipped, difference, floor, resolved, label, knob),
                )
    return table, outcomes


def _knob_verdict_tiles(report, outcomes: list[tuple], label, knob: str, shipped) -> None:
    """The suite's answer, from the same per-context outcomes the table below prints.

    The rule is the README's, applied rather than left to the reader: a value is proposed as a new
    default only if it wins EVERY context that resolved at all. One that wins some and loses others
    is a deployment knob; one that wins none is a knob whose default stands.
    """
    if not outcomes:
        return
    cleared = [entry for entry in outcomes if entry[3]]
    winners = {entry[0] for entry in cleared}
    gains = sorted(abs(entry[1]) for entry in cleared)
    floors = sorted(entry[2] for entry in outcomes)
    median_gain = gains[len(gains) // 2] if gains else 0.0
    median_floor = floors[len(floors) // 2]

    if not cleared:
        status, reading = "flat", (
            f"**FLAT** — no value of `{knob}` clears its own noise floor in any context. "
            f"The shipped default stands, and on this evidence the knob is not worth having."
        )
        best_value = label(shipped)
    elif len(winners) == 1:
        winner = next(iter(winners))
        status = "good"
        best_value = f"{label(shipped)} → {label(winner)}"
        reading = (
            f"`{knob}={label(winner)}` wins every context that resolved "
            f"({len(cleared)} of {len(outcomes)}) — proposed as a change to `SearchTuning::default()`."
        )
    else:
        # Not "unresolved" — these contexts resolved perfectly well, they just disagree. The two
        # findings need different words or the amber dot means whichever the reader guesses.
        status = "warn:no single winner"
        best_value = " / ".join(label(value) for value in sorted(winners, key=str))
        reading = (
            f"No single value wins everywhere: {len(winners)} different values take "
            f"{len(cleared)} of {len(outcomes)} contexts. That is a **deployment knob** the server "
            f"should be able to set, not a new default."
        )

    report.verdict(
        [
            ("best value", best_value, f"shipped {label(shipped)}", ""),
            ("median gain", f"{median_gain:+.1f}%" if cleared else "none", f"floor ±{median_floor:.1f}%", ""),
            ("contexts", f"{len(cleared)} / {len(outcomes)}", "cleared the floor", status),
        ],
        reading,
    )


def _knob_reading(best, shipped, difference: float, floor: float, resolved: list, label, knob: str) -> str:
    if best == shipped:
        return "the shipped value is the peak"
    if abs(difference) <= floor:
        return "flat — nothing clears the floor against the shipped value"
    if len(resolved) == 1:
        return f"only {knob}={label(best)} clears the floor — one cell, not yet a curve"
    return f"{knob}={label(best)} wins, and {len(resolved)} values clear the floor"


def _knob_table(report, points, sources, values, arms, label, knob, shipped, extra) -> None:
    """Every cell of the sweep, one row each, with the peptide file as a column.

    One table for the whole suite rather than one per regime. `mlp` used to print this four times at
    42 rows apiece; the rows are the same either way, and as one chipped table a reader can hold the
    knob and vary the regime, which four sibling tables cannot do at all.
    """
    headers = ["file", knob, *extra, "arm", "qps", "band", "vs shipped", "floor"]
    table = Table(
        headers=headers,
        aligns=["<"] * (len(extra) + 3) + [">"] * 4,
        chips=["file", knob, *extra, "arm"],
        tips=tips_for(headers),
    )
    ordered = [
        (source, context, per_arm)
        for source in sources
        for context, per_arm in sorted(((c, p) for (s, c), p in points.items() if s == source), key=str)
    ]
    for value in values:
        for source, context, per_arm in ordered:
            for arm in arms:
                curve = per_arm.get(arm, {})
                cell, reference = curve.get(value), curve.get(shipped)
                if not cell:
                    continue
                difference = delta_pct(cell["p50"], reference["p50"]) if reference else float("nan")
                table.row(
                    source,
                    label(value),
                    *(column_label(name, coord) for name, coord in zip(extra, context)),
                    arm,
                    qps(cell["p50"]),
                    band(cell_band(cell)),
                    "reference" if value == shipped else ("-" if difference != difference else pct(difference)),
                    "-" if value == shipped else band(floor_of(cell, reference)),
                )
    report.table(table, raw=True)


# ---------------------------------------------------------------------------
# Knob planes
# ---------------------------------------------------------------------------


def knob_planes(report: Report, loaded: list[Record], keys: tuple[str, ...] | None = None) -> list[bool]:
    """Every two-knob plane in the records, as a heatmap against the shipped configuration.

    The question a plane answers and a pair of curves cannot: is the ridge diagonal? If the best
    value of one knob depends on the other, one-at-a-time sweeping reports the slice through the
    default point and never sees the rest.

    That question is now drawn as well as written: the ridge — the best cell in each row — is a
    polyline over the plane, so a peak that moves is a line that bends and a reader sees the
    interaction before reading a word. The text readout stays as the table-view twin.

    The planes are one facet grid on one diverging key rather than a column of full-width figures,
    for the same reason everything else here is: six of them, each with its own legend, is the same
    legend six times and six scroll-lengths between the first and the last.

    Cells that sit on an axis of the plane were measured as one-knob cells, so they live on that
    knob's curve rather than in the plane. Both are pulled in, which is what keeps the grid whole
    instead of leaving its edges blank.

    Returns one flag per plane: whether that plane's peak MOVED, which is the suite's finding.
    """
    grouped = _plane_cells(loaded, keys)
    if not grouped:
        return []
    planes, singles, bases = grouped
    figures: list[str] = []
    rates: list[str] = []
    readouts: list[str] = []
    moved: list[bool] = []

    for (knobs, context), points in sorted(planes.items(), key=str):
        base = bases.get(context)
        if base is None or not base["p50"]:
            continue
        rows = _in_order({key[1] for key in points})
        columns = _in_order({key[0] for key in points})
        grid = {}
        rate = {}
        for row, down in enumerate(rows):
            for column, across in enumerate(columns):
                cell = points.get((across, down)) or _axis_cell(singles, knobs, across, down, context, base)
                if cell is None or not cell["p50"]:
                    continue
                difference = delta_pct(cell["p50"], base["p50"])
                floor = floor_of(cell, base)
                hover = (
                    f"{knobs[0]}: {fmt_tune(across)}\n"
                    f"{knobs[1]}: {fmt_tune(down)}\n"
                    f"throughput: {cell['p50']:,.0f} qps\n"
                    f"at the shipped pair: {base['p50']:,.0f} qps\n"
                    f"difference: {difference:+.1f}% (floor {floor:.1f}%)"
                )
                grid[(row, column)] = (difference, floor, hover)
                rate[(row, column)] = (cell["p50"], hover)
        peaks = _peaks(grid, rows, columns)
        moved.append(len({column for _, column in peaks}) > 1)
        headers = [fmt_tune(value) for value in columns]
        side = [fmt_tune(value) for value in rows]
        figures.append(
            heatmap(
                headers,
                side,
                grid,
                context,
                pos_label="faster than shipped",
                neg_label="slower",
                frame=FACET_PLANE,
                x_title=knobs[0],
                y_title=knobs[1],
                # One key above the grid, drawn from the first plane only — every plane uses the
                # same diverging scale by construction.
                legend=not figures,
                ridge=peaks,
            )
        )
        # The same cells as the throughput they were measured at. The diverging view answers "is
        # this better than what ships", which is the suite's question, and answers it as a wall of
        # neutral cells whenever the knobs do nothing — correct, and unreadable as the only view.
        # This one just shows the numbers, on a ramp scaled to the plane so the shape is visible at
        # all; the key states the spread against the noise floor, so a shape made of scatter says so.
        values = [value for value, _ in rate.values()]
        if values:
            low, high = min(values), max(values)
            rates.append(
                sequential_heatmap(
                    headers,
                    side,
                    rate,
                    context,
                    low,
                    high,
                    " qps",
                    frame=FACET_PLANE,
                    ridge=peaks,
                    floor=base.get("floor"),
                    compact=True,
                    x_title=knobs[0],
                    y_title=knobs[1],
                )
            )
        readouts += _ridge(grid, rows, columns, knobs) + [""]

    if figures:
        knob_names = sorted({knobs for knobs, _ in planes})
        axes = (
            f"columns: {knob_names[0][0]} · rows: {knob_names[0][1]}"
            if len(knob_names) == 1
            else "one knob pair per plane — see the readout below"
        )
        # Throughput leads. It is the reading that always says something: the deltas collapse to one
        # neutral colour on a suite where nothing clears its floor, which is most of them, and a
        # reader who opens the section to a grey grid learns nothing until they have read the prose.
        variants = [("throughput", rates)] if rates else []
        variants.append(("vs shipped", figures))
        report.switch("cell", variants, default=variants[0][0])
        report.para(
            f"The same planes, twice ({axes}). **throughput** is what each pair measured, on a ramp "
            f"scaled to that plane — read the spread beside its key before its shape, since a plane "
            f"that varies by less than its floor has none. **vs shipped** is each cell against the "
            f"shipped pair, neutral wherever the difference is inside the noise. "
            f"The dark line is the ridge — the best cell in each row; a ridge that bends is the "
            f"interaction a one-at-a-time sweep cannot see."
        )
        report.lines(readouts[:-1])
    return moved


def _peaks(grid: dict, rows: list, columns: list) -> list[tuple[int, int]]:
    """The best cell in each row, as (row, column) — the ridge, for drawing and for counting."""
    out = []
    for row in range(len(rows)):
        row_cells = {column: grid[(row, column)][0] for column in range(len(columns)) if (row, column) in grid}
        if row_cells:
            out.append((row, max(row_cells, key=lambda column: row_cells[column])))
    return out


def _plane_cells(loaded: list[Record], keys: tuple[str, ...] | None):
    """`(planes, single-knob cells, reference cells)`, all keyed by context.

    Grouped straight off the records rather than through `by_cell`, because a plane's coordinates
    are the two knobs a cell moved — which is a property of the cell, not a fixed column list.
    """
    keys = keys or GRID_KEYS
    drift = drift_by_process(loaded)
    planes: dict[tuple, dict] = {}
    singles: dict[tuple, dict] = {}
    bases: dict[str, dict] = {}

    for index, record in _positions(loaded):
        config = record.config
        if config.get("sweep") == DRIFT:
            continue
        tuning = config.get("tuning") or {}
        defaults = config.get("tuning_defaults") or {}
        moved = tuple(sorted(name for name, value in tuning.items() if defaults.get(name) != value))
        spread = record.spread()
        p10, p50, p90 = spread if spread else (record.qps, record.qps, record.qps)
        series = drift.get(_process_key(record))
        scale = series.scale_at(index) if series else 1.0
        cell = {
            "p10": p10 * scale,
            "p50": p50 * scale,
            "p90": p90 * scale,
            "floor": series.floor if series else float("nan"),
            "tuning": tuning,
        }
        context = _plane_context(record, keys)
        if len(moved) == 2:
            planes.setdefault((moved, context), {})[(tuning[moved[0]], tuning[moved[1]])] = cell
        elif len(moved) == 1:
            singles[((moved[0], tuning[moved[0]]), context)] = cell
        else:
            bases[context] = cell
    return (planes, singles, bases) if planes else None


def _plane_context(record: Record, keys: tuple[str, ...]) -> str:
    """Everything about a cell except its tuning — what a plane's numbers must be read against.

    The two booleans are spelled out, `tryptic=false` as readily as `tryptic=true`. This used to
    append the coordinate's NAME when it was true and nothing at all when it was false, so the
    non-tryptic plane was identified by noticing an ABSENT word next to its twin — and a reader who
    saw only one of the two had no way to tell which they were looking at. It is also what lets the
    page tag the figure with the mode it was measured under; `charts._figure` reads this text, and a
    coordinate that goes unwritten cannot be filtered on.
    """
    config = record.config
    parts = [str(record.dims.get("arm", "?"))]
    for name in ("peptide_source", "kmer_k", "equate_il", "tryptic", "amount_of_peptides"):
        value = config.get(name)
        if value is None:
            continue
        if name == "peptide_source":
            parts.append(str(value))
        elif name == "kmer_k":
            parts.append(kmer_label(value))
        elif name == "amount_of_peptides":
            parts.append(f"{value:,}q")
        else:
            parts.append(context_label(name, bool(value)))
    return " · ".join(parts)


def _axis_cell(singles: dict, knobs: tuple, across, down, context: str, base: dict) -> dict | None:
    """The plane point that lies on an axis, recovered from the one-knob cells."""
    shipped = base["tuning"]
    if across == shipped.get(knobs[0]) and down == shipped.get(knobs[1]):
        return base
    if across == shipped.get(knobs[0]):
        return singles.get(((knobs[1], down), context))
    if down == shipped.get(knobs[1]):
        return singles.get(((knobs[0], across), context))
    return None


def _ridge(grid: dict, rows: list, columns: list, knobs: tuple) -> list[str]:
    """Where the best column sits in each row — a diagonal here is the interaction."""
    out = [f"  best {knobs[0]} per {knobs[1]}:"]
    peaks = []
    for row, peak in _peaks(grid, rows, columns):
        peaks.append(peak)
        out.append(
            f"    {knobs[1]}={fmt_tune(rows[row]):<6} -> {knobs[0]}={fmt_tune(columns[peak]):<6} "
            f"({grid[(row, peak)][0]:+.1f}% against the shipped pair)"
        )
    if len(set(peaks)) > 1:
        out.append(
            f"  the peak MOVES: the best {knobs[0]} is not the same at every {knobs[1]}, which is "
            f"the interaction a one-at-a-time sweep cannot see."
        )
    elif peaks:
        out.append(f"  the peak holds at one {knobs[0]} for every {knobs[1]} — no interaction here.")
    return out
