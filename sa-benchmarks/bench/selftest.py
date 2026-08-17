"""Exercises every suite's analysis against fabricated records.

    python3 -m bench.selftest

The `ram` and `threads` suites need root, cgroup v2 and Linux, so on a development machine their
analysis code never runs — which is exactly the code most worth checking, since it is where the
crossover detection, the void-cell rule and the fault-flatness warning live. This builds records
with known shapes and asserts the analyses draw the right conclusions from them.

It is not a benchmark and touches no index. It checks that the reporting cannot silently lose a
cell: a cell that was killed under its ceiling, one whose cap did not bind, and one that never ran
must each appear in the output, distinguishable from each other.
"""

from __future__ import annotations

import json
import random
import re
import sys
import tempfile
from pathlib import Path

from . import config, records, rig
from .report import Report
from .suites import analysis_for


def write_cell(out_dir: Path, dims: dict, qps: float, majflt: int, rss_gb: float, reps: int = 40) -> None:
    """One cell's worth of per-rep records, with a little jitter so the spreads are non-degenerate."""
    out_dir.mkdir(parents=True, exist_ok=True)
    label = "__".join(f"{key}-{value}" for key, value in sorted(dims.items()))
    with (out_dir / f"{label}.jsonl").open("w") as handle:
        for _ in range(reps):
            handle.write(
                json.dumps(
                    {
                        "version": 8,
                        "label": label,
                        "commit": "selftest",
                        "suite": out_dir.name,
                        "dims": {key: str(value) for key, value in dims.items()},
                        "config": {},
                        "startup": {"load_total_ms": 1000, "warmup_ms": 500},
                        "result": {
                            "throughput_qps": qps * random.uniform(0.97, 1.03),
                            "major_faults": majflt,
                            "minor_faults": majflt * 10,
                            "total_memory": int(rss_gb * 2**30),
                            "amount_of_queries": 10_000,
                            "suffix_hit_count": 250_000,
                        },
                    }
                )
                + "\n"
            )


def build_ram(root: Path) -> Path:
    """A crossover: pprot ahead while resident, mmap ahead under pressure. Plus a void and an OOM."""
    out_dir = root / "ram"
    for ceiling, mmap_qps, pprot_qps, majflt in [
        (0, 100_000, 121_000, 0),
        (223, 98_000, 118_000, 200),
        (167, 60_000, 55_000, 24_000),
        (112, 31_000, 22_000, 51_000),
    ]:
        rss = 240 if ceiling == 0 else ceiling * 0.98
        for slot, arm, value in [
            ("a", "mmap", mmap_qps),
            ("b", "pprot", pprot_qps),
            ("c", "pprot", pprot_qps * 1.01),
            ("d", "mmap", mmap_qps * 0.99),
        ]:
            write_cell(out_dir, {"ceiling_gb": ceiling, "arm": arm, "slot": slot, "features": "mmap"}, value, majflt, rss)

    # A cell whose ceiling did not bind: RSS far above the cap. Must read as VOID, never as 90k qps.
    write_cell(out_dir, {"ceiling_gb": 78, "arm": "mmap", "slot": "a", "features": "mmap"}, 90_000, 100, 200)
    # And one the OOM killer took: no records at all, but still an answer about that arm.
    (out_dir / "ceiling_gb-78__arm-pprot__features-mmap__slot-b.oom").write_text(
        json.dumps({"exit": 137, "dims": {"ceiling_gb": "78", "arm": "pprot"}}) + "\n"
    )
    return out_dir


def build_threads(root: Path) -> Path:
    """Oversubscription costs ~10% unconstrained and pays ~60% under a ceiling, faults flat."""
    out_dir = root / "threads"
    for ceiling in (0, 167, 112):
        for threads, factor in [("default", 1.0), ("48", 1.3 if ceiling else 0.95), ("96", 1.6 if ceiling else 0.90)]:
            for arm, base in (("mmap", 60_000), ("pprot", 57_000)):
                write_cell(
                    out_dir,
                    {"ceiling_gb": ceiling, "threads": threads, "arm": arm, "features": "mmap"},
                    base * factor,
                    0 if ceiling == 0 else 24_000,
                    240 if ceiling == 0 else ceiling * 0.97,
                )
    return out_dir


# ---------------------------------------------------------------------------
# tuning: a process whose machine drifts, knobs that do and do not matter
# ---------------------------------------------------------------------------

#: The shipped tuning the fabricated records claim to have been built against.
TUNING_DEFAULTS = {
    "mlp_batch": 16,
    "validate_batch": 64,
    "prefetch_threshold": 32,
    "retrieval_prefetch_distance": 32,
}

#: How much slower the fabricated machine gets from the first cell of a process to the last.
#: Linear, so the drift cadence can cancel it exactly and any residual is a bug in the analysis
#: rather than in the fixture.
DRIFT_END = 0.70


def _grid_record(suite: str, dims: dict, config: dict, qps: float) -> dict:
    """One aggregated matrix-mode cell, in the shape `--grid-file` runs produce.

    The band is deliberately tight (±1%), so the floor stays at the measured full-database 3.9% and
    an effect either clears that or is honestly unresolved. A noisy fixture would let the analysis
    pass by being unable to decide anything.
    """
    queries = config.get("amount_of_peptides", 10_000)
    total_ns = queries / qps * 1e9
    return {
        "version": 10,
        "label": "selftest",
        "commit": "selftest",
        "suite": suite,
        "dims": {key: str(value) for key, value in dims.items()},
        "config": {"tuning_defaults": TUNING_DEFAULTS, "peptide_source": "mixed", **config},
        "startup": {"load_total_ms": 1000, "warmup_ms": 500},
        "result": {
            "throughput_qps": qps,
            "major_faults": 0,
            "minor_faults": 0,
            "total_memory": int(3.5 * 2**30),
            "amount_of_queries": queries,
            "suffix_hit_count": 250_000,
            # A plausible phase split summing to the wall time the throughput implies: retrieval is
            # a third of it. Present so the time charts have something to draw — a fixture without
            # them would let a broken chart pass by rendering nothing at all.
            "search_duration_ns": int(total_ns * 0.68),
            "retrieval_duration_ns": int(total_ns * 0.32),
            "total_duration_ns": int(total_ns),
        },
        "stats": {
            "runs": 8,
            "qps_min": qps * 0.98,
            "qps_p10": qps * 0.99,
            "qps_p50": qps,
            "qps_p90": qps * 1.01,
            "qps_max": qps * 1.02,
        },
    }


def _write_drifting(out_dir: Path, suite: str, plan: list[tuple[dict, float]], cadence: int) -> Path:
    """One process's jsonl, with the drift cadence interleaved and a linear slowdown applied.

    The slowdown is the point of the fixture. A real win placed late in a process that is getting
    slower only reads as a win if the cadence was actually applied — so a check that passes here
    could not pass against an analysis that merely *recorded* drift instead of removing it.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    dims = {"arm": "mmap", "threads": "default", "ceiling_gb": 0, "features": "mmap"}

    reference = {
        "tuning": dict(TUNING_DEFAULTS),
        "sweep": "drift",
        "kmer_k": 5,
        "equate_il": True,
        "tryptic": False,
        "amount_of_peptides": 10_000,
    }
    lines, marks = [], 0
    for index, cell in enumerate(plan):
        if index % cadence == 0:
            lines.append(({**reference, "grid_slot": f"z{marks}"}, 100_000.0))
            marks += 1
        lines.append(cell)
    lines.append(({**reference, "grid_slot": f"z{marks}"}, 100_000.0))

    with (out_dir / "ceiling_gb-0__threads-default__mmap__a.jsonl").open("w") as handle:
        for position, (config, true_qps) in enumerate(lines):
            slowdown = 1.0 - (1.0 - DRIFT_END) * (position / max(1, len(lines) - 1))
            handle.write(json.dumps(_grid_record(suite, dims, config, true_qps * slowdown)) + "\n")
    return out_dir


def build_validate(root: Path) -> Path:
    """The validate_batch curve, its tryptic twin, and a plane with a diagonal ridge.

    One claim per group:

      * `validate_batch = 128` is a real +15% at 10,000 queries, placed late in a process that is
        getting slower, so it reads as a win only under drift correction.
      * the 2,000-query tryptic block wins at a DIFFERENT value, and its throughput is deliberately
        implausible against the block above — a report that folded the two stream lengths into one
        curve would put the shipped value four times below every swept value and call it a
        landslide.
      * the `mlp_batch` x `validate_batch` plane has a diagonal ridge: the best batch rises with the
        validation batch, which is the interaction one-at-a-time sweeping structurally cannot see.
    """
    base = 100_000.0
    cells: list[tuple[dict, float]] = []

    def add(qps: float, *, amount: int = 10_000, tryptic: bool = False, sweep: str = "curve", **tune) -> None:
        cells.append((
            {
                "tuning": {**TUNING_DEFAULTS, **tune},
                "sweep": sweep,
                "grid_slot": "a",
                "kmer_k": 5,
                "equate_il": True,
                "tryptic": tryptic,
                "amount_of_peptides": amount,
                "peptide_source": "mixed",
            },
            qps,
        ))

    for tryptic in (False, True):
        add(base, tryptic=tryptic)
        for value in (16, 32, 256):
            add(base * 1.005, tryptic=tryptic, validate_batch=value)
    # Late, and therefore deepest into the fabricated slowdown: only a drift-corrected number can
    # still see this as a win.
    for tryptic in (False, True):
        add(base * 1.15, tryptic=tryptic, validate_batch=128)
    return _write_drifting(root / "validate", "validate", cells, cadence=12)


def build_mlp_validate(root: Path) -> Path:
    """A plane with a diagonal ridge: the best mlp_batch rises with the validation batch.

    That is the interaction one-at-a-time sweeping structurally cannot see, and the only reason this
    suite exists — so the fixture has it, and the check fails if the ridge readout stops saying so.
    """
    base = 100_000.0
    cells: list[tuple[dict, float]] = []

    def cell(tune: dict, qps: float) -> None:
        cells.append((
            {
                "tuning": {**TUNING_DEFAULTS, **tune},
                "sweep": "plane",
                "grid_slot": "a",
                "kmer_k": 5,
                "equate_il": True,
                "tryptic": False,
                "amount_of_peptides": 10_000,
                "peptide_source": "mixed",
            },
            qps,
        ))

    # The reference cell every plane point is read against. Without one the plane has no origin and
    # the whole grid is skipped — which is a real failure mode, so the fixture has to carry it.
    cell({}, base)
    for validate, peak in ((16, 4), (32, 4), (128, 16), (256, 64)):
        for mlp in (1, 4, 16, 64):
            cell({"mlp_batch": mlp, "validate_batch": validate}, base * (1.20 if mlp == peak else 0.95))
    return _write_drifting(root / "mlp_validate", "mlp_validate", cells, cadence=12)




def build_prefetch(root: Path) -> Path:
    """A 4x4 plane where nothing clears the floor — the null result this suite exists to state.

    Every pair within a percent of the shipped one. The suite has to report that as FLAT and as
    grounds for dropping both knobs, not as sixteen small effects.
    """
    base = 100_000.0
    cells = []
    for index, threshold in enumerate((8, 16, 32, 64)):
        for offset, distance in enumerate((8, 16, 32, 64)):
            cells.append((
                {
                    "tuning": {
                        **TUNING_DEFAULTS,
                        "prefetch_threshold": threshold,
                        "retrieval_prefetch_distance": distance,
                    },
                    "sweep": "plane",
                    "grid_slot": "a",
                    "kmer_k": 5,
                    "equate_il": True,
                    "tryptic": False,
                    "amount_of_peptides": 10_000,
                    "peptide_source": "mixed",
                },
                base * (1.0 + 0.002 * ((index + offset) % 3)),
            ))
    return _write_drifting(root / "prefetch", "prefetch", cells, cadence=8)


# ---------------------------------------------------------------------------
# defaults / kmer / mlp: the three narrow matrix suites
# ---------------------------------------------------------------------------


#: Counter fields belong in `result`, not `config`; the fixtures pass them alongside the config and
#: this is where they are moved across.
COUNTERS = (
    "search_bounds_ns",
    "match_iter_ns",
    "candidates_examined",
    "candidates_accepted",
    "response_duration_ns",
    "response_bytes",
)


def _write_grid(out_dir: Path, suite: str, arm: str, cells: list[tuple[dict, float]]) -> None:
    """One process's jsonl: an arm, and the cells it swept in the order it swept them."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dims = {"arm": arm, "features": "metrics" if arm.endswith("+m") else arm, "slot": "a"}
    with (out_dir / f"{arm.replace('+', '_')}__a.jsonl").open("w") as handle:
        for config, qps in cells:
            counters = {name: config.pop(name) for name in COUNTERS if name in config}
            record = _grid_record(suite, dims, {"tuning": dict(TUNING_DEFAULTS), **config}, qps)
            record["result"].update(counters)
            handle.write(json.dumps(record) + "\n")


#: Throughput per length regime, roughly the shape a real index gives: short peptides match far more
#: suffixes, so they are an order of magnitude slower.
BUCKETS = {"mixed": 1_200_000.0, "small": 400_000.0, "medium": 1_600_000.0, "large": 2_200_000.0}


#: What the instrumented arm's counters read. Tryptic examines far more candidates and accepts far
#: fewer of them, which is the shape the acceptance column exists to show — and the fixture has to
#: have it, or a report that dropped the column entirely would still pass.
def _counters(tryptic: bool) -> dict:
    examined = 900_000 if tryptic else 120_000
    return {
        "search_bounds_ns": 40_000_000,
        "match_iter_ns": 90_000_000 if tryptic else 20_000_000,
        "candidates_examined": examined,
        "candidates_accepted": int(examined * (0.02 if tryptic else 0.85)),
    }


#: The phase production runs after retrieval, shaped as the local index measures it: it dwarfs search
#: on a non-tryptic request and all but disappears under tryptic, which returns almost nothing. The
#: fixture has to have that asymmetry, or a report that lost the share column entirely would still
#: pass.
#:
#: Sized for the v13 measurement — production's shape, decode parallel across peptides — which is
#: several times smaller than the serial v12 number for the same work. A fixture carrying the old
#: serial ratio would let a silent revert to it pass.
def _response(tryptic: bool) -> dict:
    return {
        "response_duration_ns": 250_000 if tryptic else 11_000_000,
        "response_bytes": 92_000 if tryptic else 9_800_000,
    }


def build_defaults(root: Path, *, regressed: bool = False) -> Path:
    """The narrowed gate: equate_il x tryptic only, at one k and one batch.

    Two claims. A coordinate the suite holds fixed must NOT become a column — that is what makes
    the narrowing visible in the report rather than only in the run time. And mmap must read as
    behind preloaded on `small` by more than the floor, while the two are inside the floor
    everywhere else, so the suite is shown to separate a real gap from a tie rather than calling
    every difference a result.
    """
    out_dir = root / ("defaults-base" if regressed else "defaults")
    for arm in ("preloaded", "mmap", "pprot", "pprot+m"):
        cells = []
        for source, base in BUCKETS.items():
            for equate_il in (True, False):
                for tryptic in (True, False):
                    qps = base * (1.35 if tryptic else 1.0)
                    if arm == "mmap":
                        # A real gap in BOTH directions, and a tie in between: the sign has to
                        # survive, not just the magnitude.
                        qps *= {"small": 0.80, "large": 1.20}.get(source, 1.005)
                    elif arm == "pprot":
                        # The middle arm, and deliberately the middle number: with three arms the
                        # verdict has to compare the leader against the RUNNER-UP, so a third arm
                        # between them must not change who is reported as ahead.
                        qps *= {"small": 0.90, "large": 1.10}.get(source, 1.002)
                    # The thing a baseline diff has to catch: one cell, moved well past its floor.
                    if regressed and source == "medium" and equate_il and not tryptic:
                        qps *= 1.25
                    cells.append((
                        {
                            "peptide_source": source,
                            "kmer_k": 5,
                            "equate_il": equate_il,
                            "tryptic": tryptic,
                            "amount_of_peptides": 10_000,
                            "sweep": "defaults",
                            "grid_slot": "a",
                            **(_counters(tryptic) if arm.endswith("+m") else {}),
                            **_response(tryptic),
                        },
                        qps,
                    ))
        _write_grid(out_dir, "defaults", arm, cells)
    return out_dir


def build_kmer(root: Path) -> Path:
    """5-mer pays everywhere; the 6-mer only pays on `large`.

    That asymmetry is the point of the suite: on the short buckets the 6-mer must read as
    unresolved WITH its 2.85 GB attached to the verdict, because a table that cannot be shown to
    win is still holding its memory. A report that called it a tie would recommend it.
    """
    out_dir = root / "kmer"
    # Roughly the measured full-database shape: on the short bucket neither table can be shown to
    # help, and the 6-mer only pulls away from the 5-mer once peptides are long.
    gains = {
        "mixed": {0: 1.00, 5: 1.28, 6: 1.30},
        "small": {0: 1.00, 5: 1.02, 6: 1.03},
        "medium": {0: 1.00, 5: 1.30, 6: 1.32},
        "large": {0: 1.00, 5: 1.35, 6: 1.55},
    }
    for arm in ("preloaded", "mmap", "pprot"):
        cells = []
        for source, base in BUCKETS.items():
            for k, gain in gains[source].items():
                cells.append((
                    {
                        "peptide_source": source,
                        "kmer_k": k,
                        "equate_il": True,
                        "tryptic": False,
                        "amount_of_peptides": 10_000,
                        "sweep": "kmer",
                        "grid_slot": "a",
                    },
                    base * gain,
                ))
        _write_grid(out_dir, "kmer", arm, cells)
    return out_dir


#: The fabricated curve: rises to the shipped batch, then falls away inside the floor. `1` is the
#: scalar path. Nothing above 16 beats it, which is the reading that leaves the default alone.
MLP_CURVE = {1: 0.70, 4: 0.85, 8: 0.96, 16: 1.00, 32: 0.995, 64: 0.99, 128: 0.98}


def build_mlp(root: Path) -> Path:
    """A curve that peaks at the shipped batch for preloaded, and is flat for mmap.

    Both readings have to survive: "the shipped batch is the peak" is the answer that leaves the
    default alone, and a flat curve is the answer that says batching is not the constraint here.
    Neither is "no data", and a report that only knew how to name a winner would invent one.
    """
    out_dir = root / "mlp"
    for arm in ("preloaded", "mmap", "pprot"):
        cells = []
        for source, base in BUCKETS.items():
            for batch, gain in MLP_CURVE.items():
                cells.append((
                    {
                        "peptide_source": source,
                        "kmer_k": 5,
                        "equate_il": True,
                        "tryptic": False,
                        "amount_of_peptides": 10_000,
                        "sweep": "mlp",
                        "grid_slot": "a",
                        "tuning": {**TUNING_DEFAULTS, "mlp_batch": batch},
                    },
                    base * (1.0 if arm == "mmap" else gain),
                ))
        _write_grid(out_dir, "mlp", arm, cells)
    return out_dir


def _check_defaults_baseline(root: Path, repo) -> list[str]:
    """The gate's other half: a cell that moved past its floor must read as a REGRESSION.

    Untestable through `CHECKS`, which runs one directory — this needs two, and the comparison
    between them is the only reason the suite runs on every commit.
    """
    now = build_defaults(root)
    before = build_defaults(root, regressed=True)

    suite = config.load("defaults", repo)
    suite.baseline = before
    report = Report().heading("defaults")
    analysis_for("defaults")(report, suite, records.load_dir(now), now)
    text = report.to_text()

    failures = []
    for needle, why in (
        ("REGRESSION", "a cell that moved past its floor must read as a regression"),
        ("unchanged", "a cell that did not move must read as unchanged, not as a small delta"),
    ):
        if needle not in text:
            failures.append(f"defaults: {why} (expected {needle!r})")
        print(f"  {'ok  ' if needle in text else 'FAIL'} defaults {why}")
    return failures


def _check_grid() -> list[str]:
    """The expander: planes subsume the lines through them, and blocks dedup against each other."""
    from . import grid

    failures = []

    def check(label: str, ok: bool) -> None:
        if not ok:
            failures.append(f"grid: {label}")
        print(f"  {'ok  ' if ok else 'FAIL'} grid     {label}")

    knobs = {"tune": {"mlp_batch": [1, 16, 64], "validate_batch": [16, 64, 256]}}

    ofat = grid.tuning_points({"name": "t", "strategy": "ofat", **knobs}, TUNING_DEFAULTS)
    # 1 base + (3-1) per knob: the shipped value of each knob is already the base point.
    check("ofat spends 1 + sum(len-1) points, not their product", len(ofat) == 5)
    check("ofat emits no duplicate point", len({tuple(sorted(p.items())) for p in ofat}) == len(ofat))
    check(
        "an ofat point holds every other knob at its shipped value",
        all(all(p[k] == v for k, v in TUNING_DEFAULTS.items() if k not in ("mlp_batch", "validate_batch")) for p in ofat),
    )

    paired = grid.tuning_points(
        {"name": "t", "strategy": "pairs", "pairs": [{"axes": ["mlp_batch", "validate_batch"]}], **knobs},
        TUNING_DEFAULTS,
    )
    # The 3x3 plane already contains both lines through the default point, so the total is the
    # plane, not the plane plus the lines.
    check("a pair costs |A|x|B|, absorbing the lines it subsumes", len(paired) == 9)
    check("a pair emits no duplicate point", len({tuple(sorted(p.items())) for p in paired}) == len(paired))

    try:
        grid.tuning_points({"name": "t", "strategy": "pairs", "pairs": [{"axes": ["mlp_batch"]}], **knobs}, TUNING_DEFAULTS)
        check("a one-axis pair is rejected", False)
    except grid.GridError:
        check("a one-axis pair is rejected", True)
    try:
        grid.tuning_points({"name": "t", "strategy": "ofat", "tune": {"nope": [1]}, "pairs": [{"axes": ["nope", "x"]}]}, TUNING_DEFAULTS)
        check("a pair axis with no values is rejected", False)
    except grid.GridError:
        check("a pair axis with no values is rejected", True)
    try:
        grid.expand([{"name": "t", "strategy": "base", "arms": ["mmap"], "files": ["mixed"], "kmr": [5]}], TUNING_DEFAULTS)
        check("a misspelled context key is rejected", False)
    except grid.GridError:
        check("a misspelled context key is rejected", True)

    # Two blocks describing the same measurement must collapse: identity is what was measured, not
    # which block asked for it.
    common = {"arms": ["mmap"], "files": ["mixed"], "kmer": [5], "strategy": "base"}
    expanded = grid.expand(
        [{"name": "one", **common}, {"name": "two", **common}], TUNING_DEFAULTS, suite_defaults={"runs": 8, "amount": 10}
    )
    cells = expanded[("mmap", "default", 0)]
    check("two blocks describing the same cell collapse into one", len(cells) == 1)

    # ... but not when they measure it at different precision, which is a different measurement.
    expanded = grid.expand(
        [{"name": "one", **common}, {"name": "two", "amount": 2_000, **common}],
        TUNING_DEFAULTS,
        suite_defaults={"runs": 8, "amount": 10_000},
    )
    check("the same cell at two precisions stays two cells", len(expanded[("mmap", "default", 0)]) == 2)

    # Process grouping: a block naming one arm must not put cells in the other arm's process.
    expanded = grid.expand(
        [{"name": "m", "arms": ["mmap"], "threads": [48], "files": ["mixed"], "strategy": "base"}],
        TUNING_DEFAULTS,
    )
    check("cells land in the process their block named", list(expanded) == [("mmap", 48, 0)])

    # The cadence must interleave rather than replace, and must bracket the last stretch.
    expanded = grid.expand(
        [{"name": "k", "arms": ["mmap"], "files": ["mixed"], "strategy": "ofat", "base_every": 2, **knobs}],
        TUNING_DEFAULTS,
    )
    cells = expanded[("mmap", "default", 0)]
    drift = [cell for cell in cells if cell["sweep"] == "drift"]
    check("the drift cadence is interleaved, not substituted", len(cells) - len(drift) == 5)
    check("every drift mark has its own slot", len({cell["grid_slot"] for cell in drift}) == len(drift))
    check("the last stretch of cells is bracketed by a closing mark", cells[-1]["sweep"] == "drift")
    return failures


CHECKS = {
    "defaults": [
        ("equate_il", "the two search options must be the table's columns"),
        ("by length regime", "the regimes must be one figure on one scale, not a section each"),
        ("held at the shipped default", "held knobs must be stated against the value that ships"),
        ("a panel of the figure each", "the peptide file is how results are faceted, not what is swept"),
        ("preloaded ahead", "a backend gap past the floor must name the winner"),
        ("mmap ahead", "and must name the other one when the sign flips"),
        ("cannot separate them", "a gap inside the floor must say so, not read as a small effect"),
        ("search time", "the phase timings must each get a chart of their own"),
        ("time split", "the search/retrieval decomposition must be offered as a stacked chart"),
        ("inside the search", "the instrumented arm must get a section of its own"),
        ("accept%", "the candidate acceptance rate must reach the report"),
        ("perturb", "the instrumented arm's numbers must be labelled as perturbed"),
        ("what a request actually costs",
         "the two phases after retrieval must get a section of their own"),
        ("measured share",
         "the report must state what fraction of a request its throughput covers"),
        ("decode", "the annotation decode must reach the phase split"),
    ],
    "kmer": [
        ("per cell (throughput, search time, retrieval time, time split)",
         "the phase readings must survive as a switch over one faceted grid, throughput included"),
        ("the table pays", "a table that wins past its floor must be reported as paying"),
        ("for nothing this run can show", "a table inside the floor must be a cost, not a tie"),
        ("2.85 GB", "the price of an unresolved table must travel with the verdict"),
        ("reference", "the no-table row must be marked as what the others are read against"),
    ],
    "mlp": [
        ("the shipped value is the peak", "a curve peaking at the shipped value must leave it alone"),
        ("flat", "a curve with nothing above the floor must read as flat, not as a winner"),
        ("scalar", "mlp_batch=1 is the scalar path and must be named as one"),
    ],
    "validate": [
        ("validate_batch=128", "a real win late in a drifting process must survive the correction"),
        ("tryptic", "each swept search option must get a column of its own, not a joined label"),
        ("residual", "the resolution table must state what the process could resolve at all"),
    ],
    "mlp_validate": [
        ("the peak MOVES", "a diagonal ridge in a plane must be reported"),
        ("residual", "the resolution table must state what the process could resolve at all"),
    ],
    "prefetch": [
        ("FLAT", "a plane where nothing clears its floor must say so, not list small effects"),
        ("grounds for dropping both knobs", "a null result must be stated as the finding it is"),
        ("prefetch_threshold", "both knobs must be named as the suite's subject"),
    ],
    "ram": [
        ("VOID", "the cell whose cap did not bind must read as VOID, not as a throughput"),
        ("did not fit", "the OOM-killed cell must appear as did-not-fit"),
        ("(not run)", "a ceiling with no results must be shown, not omitted"),
        ("SIGN CHANGE", "the crossover between the arms must be flagged"),
        ("the cap did not bind", "the void cell must be explained in the caveats"),
    ],
    "threads": [
        ("best mmap", "each arm's best thread count must be named"),
        ("floor", "an arm-to-arm delta must never be shown without its noise floor"),
        ("unconstrained", "the no-ceiling block is the go/no-go and must be present"),
    ],
}


#: Cell classes the HTML page must apply, per suite. These are the ones that carry meaning: a VOID
#: cell rendered as an ordinary number is the whole failure mode the class exists to prevent, and it
#: is invisible in a passing markdown check.
#:
#: `ram` states its outcomes in words ("pprot wins", "did not fit"); `threads` states them as signed
#: deltas against each arm's own default-thread baseline. Both must survive into the page.
HTML_CLASSES = {"ram": ("void", "muted", "good"), "threads": ("pos", "neg")}

#: Page-only assertions, checked against the rendered HTML rather than the text report.
PAGE_CHECKS = {
    "defaults": (
        ("raw numbers", "the per-cell grid must fold away under the figure it belongs to"),
        ('class="pill t">true', "a boolean must render as a coloured boolean, not as jargon"),
        ('class="pill f">false', "and its false twin must be painted too"),
    ),
    "kmer": (("raw numbers", "the per-cell grid must fold away"),),
    "validate": (("raw numbers", "the per-cell grid must fold away"),),
}


def _check_html(name: str, report: Report) -> list[str]:
    """Renders the page and asserts it is well-formed, self-contained and marks the right cells."""
    from html.parser import HTMLParser

    from .html import render_page

    page = render_page(f"selftest — {name}", report, subtitle="selftest", statuses={name: "ok"})
    failures: list[str] = []

    class WellFormed(HTMLParser):
        VOID_TAGS = {"meta", "br", "hr", "input", "img", "link"}

        def __init__(self) -> None:
            super().__init__()
            self.stack: list[str] = []
            self.mismatched: list[str] = []

        def handle_starttag(self, tag: str, attrs: object) -> None:
            if tag not in self.VOID_TAGS:
                self.stack.append(tag)

        def handle_endtag(self, tag: str) -> None:
            if self.stack and self.stack[-1] == tag:
                self.stack.pop()
            else:
                self.mismatched.append(tag)

    parser = WellFormed()
    parser.feed(page)
    for label, problem in (("unclosed", parser.stack), ("mismatched", parser.mismatched)):
        if problem:
            failures.append(f"{name}: HTML is not well-formed ({label}: {problem[:4]})")
    print(f"  {'FAIL' if failures else 'ok  '} {name:<8} the page is well-formed")

    # These pages are scp'd off a server and opened from disk; a single external reference would
    # make the report render differently, or not at all, wherever it actually gets read.
    external = re.findall(r'(?:src|href)="(?:https?:)?//', page)
    if external:
        failures.append(f"{name}: page references {len(external)} external resource(s)")
    print(f"  {'FAIL' if external else 'ok  '} {name:<8} the page is self-contained")

    # No legend may show the same swatch twice. This is the check that would have caught the three
    # storage arms sharing two colours: `series_color` carries two lightnesses per hue, a third arm
    # was added, and every call site's `light=slot > 0` painted mmap and pprot identically in 138 of
    # the report's 183 legends — on the comparison the report exists to make. A fourth arm must fail
    # here rather than ship as a fourth invisible one.
    collisions = []
    for legend in re.findall(r'<div class="legend[^"]*">(.*?)</div>', page, re.S):
        swatches = re.findall(r"background:var\((--[a-z0-9-]+)\)", legend)
        if len(set(swatches)) < len(swatches):
            collisions.append(swatches)
    if collisions:
        failures.append(f"{name}: {len(collisions)} legend(s) show two series in one colour: {collisions[0]}")
    print(f"  {'FAIL' if collisions else 'ok  '} {name:<8} no two series in a legend share a colour")

    # A facet grid carries exactly one legend, above the panels. One per panel is the same row of
    # words repeated down the grid, which is what faceting was supposed to stop.
    for grid in re.findall(r'<div class="figgrid[^"]*">(.*?)(?=<div class="figgrid|</div></div>)', page, re.S):
        if grid.count('class="legend') > 1:
            failures.append(f"{name}: a facet grid draws {grid.count('class=\"legend')} legends")
            break
    print(f"  {'FAIL' if any('facet grid draws' in f for f in failures) else 'ok  '} "
          f"{name:<8} a facet grid carries one legend")

    # The column glossary has to survive into the page as a hoverable heading. It used to be a
    # paragraph under the table, and a paragraph that says what `floor` means is a paragraph the
    # reader has to scroll away from the number to find.
    if name in ("defaults", "kmer", "mlp"):
        hinted = 'class="r hint"' in page or 'class="l hint"' in page
        if not hinted:
            failures.append(f"{name}: no column heading carries its explanation on the page")
        print(f"  {'ok  ' if hinted else 'FAIL'} {name:<8} column headings explain themselves on hover")

        # And the swept setting has to be findable without reading the role column row by row.
        strong = "<tr class=strong>" in page
        if not strong:
            failures.append(f"{name}: the swept setting is not set apart in the configuration table")
        print(f"  {'ok  ' if strong else 'FAIL'} {name:<8} the swept setting is set in bold")

    # Page-only properties. The text and markdown forms can neither fold a table nor paint a
    # boolean, so these cannot be asserted on `to_text()` — and both are things a reader of the page
    # would notice immediately and a passing markdown check would never catch.
    for needle, why in PAGE_CHECKS.get(name, ()):
        present = needle in page
        if not present:
            failures.append(f"{name}: {why} (expected {needle!r} in the page)")
        print(f"  {'ok  ' if present else 'FAIL'} {name:<8} {why}")

    for css_class in HTML_CLASSES.get(name, ()):
        present = f' {css_class}"' in page
        if not present:
            failures.append(f"{name}: no cell was marked '{css_class}' in the HTML")
        print(f"  {'ok  ' if present else 'FAIL'} {name:<8} cells are marked '{css_class}'")
    return failures


#: Which columns become filter chips, and which are measurements that merely repeat. The chip
#: heuristic reads the data rather than being told by the suites, so it is the kind of thing that
#: silently starts offering six buttons for a `ratio` column the next time a suite adds one.
CATEGORY_CASES = [
    ("booleans", ["True", "False", "True", "False"], True),
    ("kmer", ["none", "5-mer", "none", "5-mer"], True),
    ("batch", ["scalar", "16", "scalar", "16"], True),
    ("threads", ["default", "48", "96", "default", "48", "96"], True),
    ("ceiling", ["none", "223G", "167G", "none", "223G", "167G"], True),
    ("status", ["ok", "ok", "ok", "skipped", "skipped"], True),
    ("arm", ["mmap", "pprot", "mmap", "pprot"], True),
    ("ratio", ["0.90x", "0.92x", "0.94x", "0.90x", "0.92x", "0.94x"], False),
    ("percent", ["77%", "79%", "80%", "77%", "79%", "80%"], False),
    ("seconds", ["0.7s", "0.8s", "1.4s", "0.7s", "0.8s", "1.4s"], False),
    ("plain numbers", ["0.8", "0.9", "0.8", "0.9"], False),
    ("wall clock", ["0.0 min", "0.1 min", "0.0 min", "0.3 min"], False),
    ("prose", ["cannot run suite ram: needs root (cgroup v2)"] * 2 + ["x", "y"], False),
    ("all distinct", ["a", "b", "c", "d"], False),
    # A "vs baseline" column: mostly deltas, with one `base` label that used to make it look
    # categorical because not every value was a number.
    ("vs baseline", ["base", "+9.1%", "+16.6%", "base", "+9.1%", "+16.6%"], False),
    ("rep count", ["80", "40", "80", "40"], False),
    ("qps", ["1,234,567", "987,654", "1,234,567", "987,654"], False),
]


def _check_categories() -> list[str]:
    from .html import _is_category

    failures = []
    for name, values, expected in CATEGORY_CASES:
        distinct = sorted({value for value in values if value and value != "-"})
        got = _is_category(values, distinct)
        if got != expected:
            failures.append(f"chips: '{name}' should {'' if expected else 'not '}be a filter group")
        print(f"  {'ok  ' if got == expected else 'FAIL'} chips    {name} -> {'chips' if got else 'no chips'}")
    return failures


def _check_charts() -> list[str]:
    """Renders every chart form and checks the SVG closes every element it opens.

    Worth its own check because the forms are not evenly exercised: a chart used only by the master
    run's headline appears in no single-suite report, so an unclosed tag in it would survive every
    other test here and only surface as a subtly broken page.
    """
    from html.parser import HTMLParser

    from .charts import Series, grouped_columns, heatmap, lines, stacked_columns, stacked_rows

    forms = {
        "grouped_columns": grouped_columns(
            ["small", "medium", "large"],
            [Series("preloaded", [320_000, 2_100_000, 1_080_000], 0), Series("mmap", [274_000, None, 1_228_000], 1)],
            "columns",
            unit=" qps",
        ),
        "lines": lines(
            ["1", "16", "128"],
            [Series("preloaded", [900_000, 1_200_000, 1_290_000], 0), Series("mmap", [880_000, None, 1_110_000], 1)],
            "lines",
            unit=" qps",
            x_title="batch",
        ),
        "stacked_columns": stacked_columns(
            ["il=true · tryptic=false", "il=true · tryptic=true"],
            ["preloaded", "mmap"],
            [Series("search", [0.7, 0.8, 1.4, None], 0), Series("retrieval", [0.3, 0.4, 0.2, 0.5], 1)],
            "stacked columns",
            unit=" µs",
        ),
        "stacked_rows": stacked_rows(
            ["preloaded", "mmap"],
            [Series("sa", [0.2, 0.0], 0), Series("proteins", [0.1, 0.0], 1), Series("mapping", [0.6, 0.0], 2)],
            "stack",
        ),
        "heatmap": heatmap(
            ["batch 16 · 5-mer", "scalar · none"],
            ["il=True tryptic=False", "il=False tryptic=True"],
            {
                (0, 0): (12.0, 3.9, "resolved, positive"),
                (0, 1): (-30.0, 3.9, "resolved, negative"),
                (1, 0): (1.0, 3.9, "inside the floor"),
                # (1, 1) deliberately absent: a cell with no data must draw as absent, not as zero.
            },
            "heat",
            pos_label="mmap",
            neg_label="preloaded",
        ),
    }

    class WellFormed(HTMLParser):
        VOID_TAGS = {"meta", "br", "hr", "input", "img", "link", "line", "use"}

        def __init__(self) -> None:
            super().__init__()
            self.stack: list[str] = []
            self.mismatched: list[str] = []

        def handle_startendtag(self, tag: str, attrs: object) -> None:
            pass  # self-closing: opens and closes in one go

        def handle_starttag(self, tag: str, attrs: object) -> None:
            if tag not in self.VOID_TAGS:
                self.stack.append(tag)

        def handle_endtag(self, tag: str) -> None:
            if self.stack and self.stack[-1] == tag:
                self.stack.pop()
            else:
                self.mismatched.append(tag)

    failures = []
    for name, svg in forms.items():
        parser = WellFormed()
        parser.feed(svg)
        broken = parser.stack or parser.mismatched
        if broken:
            failures.append(f"charts: {name} SVG is not well-formed (open {parser.stack}, bad {parser.mismatched})")
        if not svg:
            failures.append(f"charts: {name} rendered nothing")
        print(f"  {'FAIL' if broken or not svg else 'ok  '} charts   {name} renders and closes every element")

        # Every form must be hoverable. This is the check the previous round did not have: the
        # charts carried SVG <title> children, which the page's own root <title> shadowed, so a
        # chart that looked finished had no working hover at all and nothing said so.
        hoverable = 'data-tip="' in svg and 'class="mark' in svg
        if not hoverable:
            failures.append(f"charts: {name} has no hoverable marks — a chart with no hover is a picture")
        print(f"  {'ok  ' if hoverable else 'FAIL'} charts   {name} marks carry a hover payload")

    # A missing heatmap cell must not be painted as if it were a measured zero.
    if svg and "opacity=\".35\"" not in forms["heatmap"]:
        failures.append("charts: a heatmap cell with no data is not drawn as absent")
    return failures


def _check_unknown_knob() -> list[str]:
    """A `SearchTuning` field this driver has never heard of must still be swept and reported.

    That is the whole contract of the generic tuning path: adding a knob in Rust should need no
    change on this side. The check uses a made-up field name, so it passes only if nothing here is
    matching against a list of known knobs.
    """
    from .records import Record
    from .report import Report
    from .suites.shared import held_and_swept as _held

    made_up = "a_knob_this_driver_has_never_heard_of"
    loaded = [
        Record({"config": {"tuning": {made_up: value, "validate_batch": 64},
                           "tuning_defaults": {made_up: 7, "validate_batch": 64}}})
        for value in (7, 9)
    ]
    report = Report()
    _held(report, loaded)
    text = report.to_text()

    failures = []
    for needle, why in (
        (made_up, "a knob this driver has never heard of must appear in the configuration table"),
        ("SWEPT", "a knob varied across cells must be marked as the suite's subject"),
        ("held at the shipped default", "a knob left alone must say so against the shipped value"),
    ):
        if needle not in text:
            failures.append(f"tuning: {why}")
        print(f"  {'ok  ' if needle in text else 'FAIL'} tuning   {why}")

    # And an override of a knob it has never heard of must still be called out.
    overridden = [Record({"config": {"tuning": {made_up: 99}, "tuning_defaults": {made_up: 7}}})]
    report = Report()
    _held(report, overridden)
    flagged = "NOT at the shipped tuning" in report.to_text()
    if not flagged:
        failures.append("tuning: an overridden knob must be flagged even if unknown to the driver")
    print(f"  {'ok  ' if flagged else 'FAIL'} tuning   an overridden unknown knob is flagged")
    return failures


def main() -> int:
    random.seed(7)
    repo = rig.repo_root()
    failures: list[str] = _check_categories() + _check_charts() + _check_unknown_knob() + _check_grid()

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        failures += _check_defaults_baseline(root, repo)
        for name, build in (
            ("defaults", build_defaults),
            ("kmer", build_kmer),
            ("mlp", build_mlp),
            ("ram", build_ram),
            ("threads", build_threads),
            ("validate", build_validate),
            ("mlp_validate", build_mlp_validate),
            ("prefetch", build_prefetch),
        ):
            out_dir = build(root)
            suite = config.load(name, repo)
            report = Report().heading(name)
            analysis_for(name)(report, suite, records.load_dir(out_dir), out_dir)
            text = report.to_text()

            for needle, why in CHECKS[name]:
                status = "ok  " if needle in text else "FAIL"
                if needle not in text:
                    failures.append(f"{name}: {why} (expected {needle!r} in the report)")
                print(f"  {status} {name:<8} {why}")

            # Every rendering must survive the same report; a crash in one of them would otherwise
            # only surface at the end of a multi-hour master run.
            report.to_markdown()
            failures += _check_html(name, report)

    if failures:
        print("\n".join(["", "FAILED:"] + [f"  - {failure}" for failure in failures]))
        return 1
    print("\nall analyses render and reach the right conclusions on known inputs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
