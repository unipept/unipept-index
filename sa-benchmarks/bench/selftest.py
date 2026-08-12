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


CHECKS = {
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

    from .charts import Series, grouped_columns, heatmap, lines, stacked_rows

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
    from .suites.defaults import _held

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
        (f"{made_up} ∈", "a knob varied across cells must be reported as swept"),
        ("validate_batch=64", "a knob held at one value must be reported as held"),
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
    failures: list[str] = _check_categories() + _check_charts() + _check_unknown_knob()

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name, build in (("ram", build_ram), ("threads", build_threads)):
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
