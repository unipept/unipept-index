"""Rendering results, once, for both the terminal and `report.md`.

One renderer with two outputs. The suites build a `Report` out of headings, tables and notes; the
terminal gets it as plain text and the master run gets the same thing as markdown. A suite that
formatted its own strings would drift from its markdown twin exactly the way the ten driver scripts
drifted from each other.

Tables are fixed-width and render inside a fenced block in markdown, so numbers stay in columns in
both. Two formatting rules are enforced here rather than left to each suite:

* a delta is never shown without the floor it has to clear (`delta_cell`);
* a void or missing cell prints a reason, never a blank or a zero that reads as a measurement.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .records import Summary, verdict


@dataclass
class Table:
    """A fixed-width table. Column widths come from the content, alignment from the caller."""

    headers: list[str]
    #: One of "<" or ">" per column; defaults to right-aligned except the first column.
    aligns: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)
    #: Headers to offer as filter chips on the page. `None` lets the renderer decide from the data,
    #: which is right for a table whose columns it has never seen. Naming them is for the case where
    #: the guess is wrong in a specific way: a suite that varies one thing wants a filter for THAT
    #: thing, not a row of chips for every column that happens to repeat a few values.
    chips: list[str] | None = None
    #: Header -> what that column means, shown on hover. This is where a column glossary belongs:
    #: beside the column, at the moment of reading it, rather than in a paragraph under the table
    #: that has to be scrolled to and then mapped back onto the headings.
    tips: dict[str, str] = field(default_factory=dict)
    #: Row indices to set in bold on the page. For the one or two rows a table exists to point at —
    #: in the configuration table, the setting this suite actually sweeps. The fixed-width text and
    #: markdown forms cannot carry weight, so they say it in words instead; the page can do both.
    strong: set[int] = field(default_factory=set)

    def __post_init__(self) -> None:
        """Normalises `aligns` to exactly one entry per header.

        Short lists used to be silently tolerated, and the HTML renderer zips headers against them —
        so a table with eight headers and two aligns emitted two `<th>` and eight `<td>`, and every
        column below the second sat one or more columns to the left of its own heading. In the
        fixed-width text form the same table looked fine, which is why it survived: the bug was
        only visible on the page.
        """
        if not self.aligns:
            self.aligns = ["<"] + [">"] * (len(self.headers) - 1)
        elif len(self.aligns) != len(self.headers):
            fill = self.aligns[-1] if self.aligns else ">"
            self.aligns = (self.aligns + [fill] * len(self.headers))[: len(self.headers)]

    def row(self, *cells: Any) -> None:
        self.rows.append(["" if cell is None else str(cell) for cell in cells])

    def render(self) -> list[str]:
        widths = [len(header) for header in self.headers]
        for row in self.rows:
            for index, cell in enumerate(row):
                widths[index] = max(widths[index], len(cell))

        def line(cells: list[str]) -> str:
            return "  ".join(
                f"{cell:{align}{width}}" for cell, align, width in zip(cells, self.aligns, widths)
            ).rstrip()

        return [line(self.headers), "  ".join("-" * width for width in widths)] + [
            line(row) for row in self.rows
        ]


@dataclass
class Report:
    """An ordered list of blocks, renderable as text or as markdown."""

    blocks: list[tuple[str, Any]] = field(default_factory=list)

    def heading(self, text: str, level: int = 2, folded: bool = False) -> "Report":
        """A section heading. `folded` starts its section collapsed on the page.

        Declared by the suite rather than guessed from the text: which sections are worth showing
        first is a fact about what the suite is FOR — a summary is the answer and the per-regime
        detail is the working — and matching on heading strings would make that a naming
        convention nobody is told about.
        """
        self.blocks.append(("heading", (level, text, folded)))
        return self

    def para(self, text: str) -> "Report":
        self.blocks.append(("para", text.strip()))
        return self

    def table(self, table: Table, raw: bool | str = False) -> "Report":
        """A table. `raw` folds it away on the page, behind a "raw numbers" summary — or behind the
        string given, for a fold whose contents are not numbers at all.

        The page's charts are what a reader scans; the exhaustive per-cell grid underneath is what
        they open once a chart has pointed at something. Marking the grid `raw` puts the two in that
        order instead of making every figure scroll past forty rows to reach the next one.

        A table that states a CONCLUSION is never raw — the verdicts, the resolution floors, the
        configuration this run used. Those are the report; the cells behind them are the evidence.
        """
        self.blocks.append(("table", (table, raw)))
        return self

    def lines(self, lines: list[str]) -> "Report":
        """A pre-formatted block (a key/value list, a plan dump) kept monospace in both outputs."""
        self.blocks.append(("lines", lines))
        return self

    def chart(self, svg: str, caption: str) -> "Report":
        """An inline-SVG figure. Only the HTML page draws it.

        Charts are always an additional reading of numbers that are also printed in a table beside
        them, so the text and markdown outputs lose nothing by naming the figure and moving on —
        and the reader of a terminal is not shown a hole where something should be.
        """
        if svg:
            self.blocks.append(("chart", (svg, caption)))
        return self

    def verdict(self, tiles: list[tuple[str, str, str, str]], reading: str = "") -> "Report":
        """What this suite found, as a row of stat tiles above everything it found it from.

        A suite's answer is one or two numbers and a word — `6-mer, +4.1%, resolved in 20 of 24
        contexts`. It used to live in the `reading` column at the right edge of a 42-row
        table, which is the single most reliable way for a page to bury its own point: the reader
        who needs it most is the one who will not read that far.

        Each tile is `(label, value, sub, status)`, and `status` is one of `good` / `flat` / `warn` /
        `""`. The status is drawn as a dot AND spelled out — colour never carries it alone — and it
        uses the reserved status tokens rather than a series colour, because here it means good-or-
        bad rather than which-series.
        """
        tiles = [tile for tile in tiles if tile[1]]
        if tiles:
            self.blocks.append(("verdict", (tiles, reading.strip())))
        return self

    def figures(self, svgs: list[str], caption: str) -> "Report":
        """A small-multiple grid: several figures read as one, under one caption.

        Distinct from `chart` because the grid is the unit — the panels share a scale and a legend,
        and pulling one out of it says something the run did not measure.
        """
        svgs = [svg for svg in svgs if svg]
        if svgs:
            self.blocks.append(("figures", (svgs, caption)))
        return self

    def switch(self, title: str, variants: list[tuple[str, list[str]]], default: str = "") -> "Report":
        """One set of figures with several colourings, and a control to pick between them.

        The same grid answers different questions depending on what is painted into it — how fast
        each backend is on its own, or which of them is ahead — and they are the same cells either
        way. A switch keeps them one figure instead of three that have to be mentally aligned.
        """
        variants = [(label, [svg for svg in svgs if svg]) for label, svgs in variants]
        variants = [(label, svgs) for label, svgs in variants if svgs]
        if variants:
            self.blocks.append(("switch", (title, variants, default or variants[0][0])))
        return self

    def note(self, text: str) -> "Report":
        """The 'how to read this' prose. Carried through verbatim — it is the interpretation."""
        self.blocks.append(("note", text.strip("\n")))
        return self

    def warn(self, text: str) -> "Report":
        self.blocks.append(("warn", text.strip()))
        return self

    def extend(self, other: "Report") -> "Report":
        self.blocks.extend(other.blocks)
        return self

    # -- outputs

    def to_text(self) -> str:
        out: list[str] = []
        for kind, payload in self.blocks:
            if kind == "heading":
                level, text = payload[0], payload[1]
                marker = "==" if level <= 2 else "--"
                out += ["", f"{marker} {text} {marker}"]
            elif kind == "para":
                out += ["", payload]
            elif kind == "table":
                out += [""] + payload[0].render()
            elif kind == "lines":
                out += [""] + payload
            elif kind == "chart":
                out += ["", f"  [chart] {payload[1]} — see report.html"]
            elif kind == "verdict":
                tiles, reading = payload
                out += ["", "  " + "   ".join(f"{label}: {value}" for label, value, _, _ in tiles)]
                if reading:
                    out += [f"  {reading}"]
            elif kind == "figures":
                out += ["", f"  [figures] {payload[1]} — see report.html"]
            elif kind == "switch":
                title, variants, _ = payload
                shown = ", ".join(label for label, _ in variants)
                out += ["", f"  [figures] {title} ({shown}) — see report.html"]
            elif kind == "note":
                out += [""] + [f"  {line}" for line in payload.splitlines()]
            elif kind == "warn":
                out += ["", f"  !! {payload}"]
        return "\n".join(out).strip() + "\n"

    def to_markdown(self) -> str:
        out: list[str] = []
        for kind, payload in self.blocks:
            if kind == "heading":
                level, text = payload[0], payload[1]
                out += ["", f"{'#' * level} {text}"]
            elif kind == "para":
                out += ["", payload]
            elif kind in ("table", "lines"):
                body = payload[0].render() if kind == "table" else payload
                out += ["", "```text", *body, "```"]
            elif kind == "chart":
                out += ["", f"*[chart: {payload[1]} — see report.html]*"]
            elif kind == "verdict":
                # A tile row is text-shaped, so markdown keeps it rather than pointing at the page:
                # `report.md` is what gets pasted into a PR, and the verdict is the reason to paste.
                tiles, reading = payload
                out += ["", " · ".join(f"**{label}** {value}" for label, value, _, _ in tiles)]
                if reading:
                    out += ["", reading]
            elif kind == "figures":
                out += ["", f"*[figures: {payload[1]} — see report.html]*"]
            elif kind == "switch":
                title, variants, _ = payload
                shown = ", ".join(label for label, _ in variants)
                out += ["", f"*[figures: {title} ({shown}) — see report.html]*"]
            elif kind == "note":
                # Blockquote: the interpretation, visually separated from the numbers.
                out += [""] + [f"> {line}" if line else ">" for line in payload.splitlines()]
            elif kind == "warn":
                out += ["", f"**⚠ {payload}**"]
        return "\n".join(out).strip() + "\n"


# ---------------------------------------------------------------------------
# Cell formatting
# ---------------------------------------------------------------------------


def qps(value: float | None) -> str:
    return "-" if value is None else f"{value:,.0f}"


def pct(value: float, digits: int = 1) -> str:
    return "-" if value != value else f"{value:+.{digits}f}%"


def band(value: float) -> str:
    return "-" if value != value else f"±{value:.1f}%"


def gb(value: float | None) -> str:
    return "-" if value is None else f"{value:,.1f}"


def seconds(milliseconds: float | None) -> str:
    return "n/a" if milliseconds is None else f"{milliseconds / 1000:.1f}s"


def count(value: float | None) -> str:
    return "-" if value is None else f"{value:,.0f}"


def cell_qps(summary: Summary | None) -> str:
    """A cell's throughput, or why there is no number to show."""
    if summary is None:
        return "(not run)"
    if summary.void_reason:
        return "VOID"
    return qps(summary.qps)


def delta_cell(new: Summary | None, base: Summary | None) -> tuple[str, str]:
    """(delta, verdict) for one comparison, never showing a delta without its floor."""
    from .records import delta_pct, noise_floor

    if not (new and base and new.usable and base.usable):
        return "-", "no data"
    difference = delta_pct(new.qps, base.qps)
    floor = noise_floor(new, base)
    return pct(difference), verdict(difference, floor)


def caveats(summaries: list[Summary], drift_limit: float = 10.0) -> list[str]:
    """Everything about this run that weakens it: void cells, and cells still climbing.

    Collected rather than printed inline so the master report can gather them into one section —
    a caveat buried under a table is a caveat nobody reads.
    """
    notes: list[str] = []
    for summary in summaries:
        where = " ".join(f"{key}={value}" for key, value in sorted(summary.dims.items()))
        if summary.void_reason:
            notes.append(f"{where}: {summary.void_reason}")
        elif summary.drift == summary.drift and abs(summary.drift) > drift_limit:
            notes.append(
                f"{where}: reps drifted {summary.drift:+.1f}% from the first quarter to the last, "
                f"so this cell had not reached steady state and its median understates it"
            )
    return notes
