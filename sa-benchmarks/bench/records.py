"""Reading results back, and the statistics that decide whether a number means anything.

Each driver script used to carry its own copy of this arithmetic in a heredoc, keyed off the file
name. Since schema v8 the sweep coordinates travel inside the record, so grouping is a `dims`
lookup and there is exactly one implementation of each statistic.

Four of them are load-bearing, and none is optional decoration:

**band** — half the p10..p90 spread of a cell's reps, as a percent of its median. How steady one
cell was.

**slot spread** — the gap between a cell's own two invocations under a palindrome ordering. This is
the floor on what the experiment can resolve: a between-arm delta smaller than it is not a small
effect, it is *no answer*. The storage-axes run made the case — the real effects had 0.1-2.1%
spreads and the noise had 6.8-12.9%.

**rep drift** — the first quarter of a cell's reps against the last, in run order. A capped cell
starts with whatever the page sweep left in the cache and climbs; a large positive drift means it
never reached steady state, so its median understates that arm.

**void** — a capped cell whose RSS landed above its ceiling was never actually constrained. That is
how a bad probe once produced "zero major faults", and such a cell must be discarded rather than
read.
"""

from __future__ import annotations

import json
import signal
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

#: A capped cell may sit slightly above its ceiling (the harness's own allocations are charged too);
#: beyond this fraction the cap did not bind and the cell is void.
CEILING_TOLERANCE = 1.05

#: Measured run-to-run noise floor on the full database, from the reference matrix sweep. A delta
#: smaller than this is noise even when a cell's own band happens to look tighter.
NOISE_FLOOR_PCT = 3.9


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


@dataclass
class Record:
    """One JSONL line. `raw` is kept whole so a suite can reach a field this module never named."""

    raw: dict[str, Any]

    @property
    def dims(self) -> dict[str, str]:
        return self.raw.get("dims", {})

    @property
    def result(self) -> dict[str, Any]:
        return self.raw.get("result", {})

    @property
    def startup(self) -> dict[str, Any]:
        return self.raw.get("startup", {}) or {}

    @property
    def config(self) -> dict[str, Any]:
        return self.raw.get("config", {})

    @property
    def qps(self) -> float:
        return float(self.result.get("throughput_qps", 0.0))

    def spread(self) -> tuple[float, float, float] | None:
        """(p10, p50, p90) when this record already aggregates reps (matrix mode), else None."""
        stats = self.raw.get("stats")
        if not stats:
            return None
        return stats["qps_p10"], stats["qps_p50"], stats["qps_p90"]


def load_file(path: Path) -> list[Record]:
    records = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            records.append(Record(json.loads(line)))
        except json.JSONDecodeError as error:
            raise ValueError(f"{path}:{lineno}: invalid JSON — {error}") from None
    return records


def load_dir(path: Path) -> list[Record]:
    """Every record under a results directory, in file-name order."""
    records: list[Record] = []
    for jsonl in sorted(path.glob("*.jsonl")):
        records.extend(load_file(jsonl))
    return records


#: How a SIGKILL — under a cgroup ceiling, the OOM killer — is recorded in a marker.
#:
#: Mirrors `runner.OOM_STATUSES`; kept here too so a marker can be judged without importing the
#: runner. Two spellings because two things write one: `subprocess` reports a killed child as the
#: negated signal number, while the bash scripts this package replaces wrote the shell's
#: `128 + signal`. A marker from either is a result, not a crash.
OOM_EXIT = 137
OOM_STATUSES = (OOM_EXIT, -signal.SIGKILL)


def unfit_cells(path: Path) -> list[dict[str, str]]:
    """Dims of the cells that were OOM-killed under their ceiling.

    A cell killed under its ceiling produces no records, but "this arm cannot run at this ceiling"
    is an answer about the arm and has to reach the report. The runner therefore writes the cell's
    dims into the marker, so this needs no file-name parsing.

    The exit status is checked rather than trusted from the file extension. The runner now writes a
    marker only for an OOM, but a session started under an older driver has markers for every
    non-zero exit in it, and reading one of those would report a panic or a bad path as a fact about
    the arm's memory behaviour — a claim no measurement was made for. A marker that says otherwise
    is skipped; the cell then has neither records nor a marker, which is what "did not run" looks
    like everywhere else.
    """
    unfit = []
    for marker in sorted(path.glob("*.oom")):
        try:
            recorded = json.loads(marker.read_text())
        except (json.JSONDecodeError, OSError):
            # A truncated write. Better a cell with no dims than a silently dropped one.
            unfit.append({"label": marker.stem})
            continue
        if recorded.get("exit", OOM_EXIT) not in OOM_STATUSES:
            continue
        unfit.append(recorded.get("dims") or {"label": marker.stem})
    return unfit


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------

#: Dimensions that identify *when* a cell ran rather than *what* it measured, so grouping ignores
#: them: every slot of one configuration is the same cell, measured twice.
POSITIONAL_DIMS = ("slot",)


def group_key(dims: dict[str, str], ignore: Iterable[str] = POSITIONAL_DIMS) -> tuple:
    skip = set(ignore)
    return tuple(sorted((key, value) for key, value in dims.items() if key not in skip))


def group(records: Iterable[Record], ignore: Iterable[str] = POSITIONAL_DIMS) -> dict[tuple, list[Record]]:
    """Buckets records by their dims, ignoring the positional ones."""
    grouped: dict[tuple, list[Record]] = {}
    for record in records:
        grouped.setdefault(group_key(record.dims, ignore), []).append(record)
    return grouped


def dims_of(key: tuple) -> dict[str, str]:
    return dict(key)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def percentile(values: list[float], fraction: float) -> float:
    """Linear-interpolated percentile of an unsorted list. Mirrors the harness's own `percentile`."""
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = fraction * (len(ordered) - 1)
    low = int(rank)
    frac = rank - low
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * frac


def median(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.median(values) if values else 0.0


@dataclass
class Summary:
    """Everything a suite needs to know about one cell, pooled over its slots."""

    dims: dict[str, str]
    n: int
    qps: float
    p10: float
    p90: float
    #: Half the p10..p90 spread, as a percent of the median.
    band: float
    #: Gap between this cell's own slots, as a percent. NaN when it ran in only one slot.
    slot_spread: float
    #: First quarter of reps vs the last, as a percent. NaN when there are too few reps.
    drift: float
    major_faults: float
    minor_faults: float
    rss_gb: float
    startup: dict[str, float] = field(default_factory=dict)
    #: Set when a ceiling was requested but the cell's RSS shows it never bound.
    void_reason: str | None = None

    @property
    def usable(self) -> bool:
        return self.void_reason is None and self.n > 0


def summarise(records: list[Record], ceiling_gb: float | None = None) -> Summary:
    """Pools one cell's records into a `Summary`.

    Handles both record shapes: per-rep lines (single mode) and pre-aggregated lines carrying a
    `stats` spread (matrix mode). Mixing the two in one cell would be meaningless, so the
    pre-aggregated form wins when present.
    """
    dims = dims_of(group_key(records[0].dims)) if records else {}

    aggregated = [record for record in records if record.spread()]
    if aggregated:
        spreads = [record.spread() for record in aggregated]
        p10 = median(s[0] for s in spreads)
        qps = median(s[1] for s in spreads)
        p90 = median(s[2] for s in spreads)
        reps = sum(record.raw["stats"]["runs"] for record in aggregated)
        slot_spread, drift = _slot_spread(aggregated), float("nan")
    else:
        values = [record.qps for record in records]
        p10, qps, p90 = percentile(values, 0.10), median(values), percentile(values, 0.90)
        reps = len(values)
        slot_spread, drift = _slot_spread(records), _drift(records)

    rss_bytes = median(record.result.get("total_memory", 0) for record in records)
    summary = Summary(
        dims=dims,
        n=reps,
        qps=qps,
        p10=p10,
        p90=p90,
        band=(p90 - p10) / 2 / qps * 100 if qps else 0.0,
        slot_spread=slot_spread,
        drift=drift,
        major_faults=median(record.result.get("major_faults", 0) for record in records),
        minor_faults=median(record.result.get("minor_faults", 0) for record in records),
        rss_gb=rss_bytes / 2**30,
        startup=_startup(records),
    )

    if ceiling_gb and summary.rss_gb > ceiling_gb * CEILING_TOLERANCE:
        summary.void_reason = (
            f"RSS {summary.rss_gb:.0f} GB exceeds the {ceiling_gb:g} GB ceiling — the cap did not "
            f"bind (was the page cache dropped?), so this cell measures nothing"
        )
    return summary


def _slot_spread(records: list[Record]) -> float:
    """Gap between the medians of a cell's two slots, as a percent of their median."""
    by_slot: dict[str, list[float]] = {}
    for record in records:
        by_slot.setdefault(record.dims.get("slot", "a"), []).append(
            record.spread()[1] if record.spread() else record.qps
        )
    if len(by_slot) != 2:
        return float("nan")
    first, second = (median(values) for values in by_slot.values())
    centre = median([first, second])
    return abs(first - second) / centre * 100 if centre else float("nan")


def _drift(records: list[Record]) -> float:
    """First quarter of the reps against the last, in run order, within a single slot.

    Reps of different slots are not consecutive in time, so pooling them would compare a slot's
    early reps with another slot's late ones and report drift that is really the slot gap.
    """
    by_slot: dict[str, list[float]] = {}
    for record in records:
        by_slot.setdefault(record.dims.get("slot", "a"), []).append(record.qps)
    values = by_slot[sorted(by_slot)[0]]
    if len(values) < 4:
        return float("nan")
    quarter = max(1, len(values) // 4)
    first, last = median(values[:quarter]), median(values[-quarter:])
    return (last - first) / first * 100 if first else float("nan")


def _startup(records: list[Record]) -> dict[str, float]:
    """Startup is one value per invocation, repeated on every rep — take it, do not average it."""
    for record in records:
        if record.startup:
            return {key: float(value) for key, value in record.startup.items()}
    return {}


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def delta_pct(new: float, base: float) -> float:
    return (new - base) / base * 100 if base else float("nan")


def noise_floor(*summaries: Summary) -> float:
    """The floor below which a delta between these cells is not an answer.

    The widest of: each cell's own slot spread, each cell's own band, and the measured full-database
    noise floor. Taking the widest is deliberate — the point is to refuse to call something an
    effect, not to find the reading under which it qualifies.
    """
    candidates = [NOISE_FLOOR_PCT]
    for summary in summaries:
        if summary.slot_spread == summary.slot_spread:  # not NaN
            candidates.append(summary.slot_spread)
        candidates.append(summary.band)
    return max(candidates)


def verdict(delta: float, floor: float, better: str = "new", worse: str = "base") -> str:
    """Names the winner only when the delta clears the floor; otherwise says the run cannot tell."""
    if delta != delta:  # NaN
        return "no data"
    if abs(delta) <= floor:
        return f"unresolved (within the {floor:.1f}% floor)"
    return f"{better} wins" if delta > 0 else f"{worse} wins"
