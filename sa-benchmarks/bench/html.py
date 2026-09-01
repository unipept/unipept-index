"""Rendering a `Report` as one self-contained, navigable HTML page.

A benchmark report is read by scrolling around it — comparing a cell to the same cell three suites
away, finding the one row that moved, checking whether a delta cleared its floor. Flat markdown
makes all of that a text search. This adds what actually helps:

* a sidebar built from the report's own headings, with a status dot per suite, so `ram` being
  skipped is visible before you scroll;
* real tables, with the cells that mean something — VOID, REGRESSION, did-not-fit — marked so they
  cannot be skimmed past;
* a filter box that hides non-matching rows across every table at once, which is how you get from a
  34-row grid to the four rows you care about;
* click-to-sort columns, parsing the numbers out of `1,234,567` and `+4.2%`.

Self-contained by construction: everything is inlined, because these files get scp'd off a
benchmark server and opened from disk with no network.

The renderer never changes what a suite said. It classifies cells by matching the strings the
suites already emit (see `CELL_MARKERS`), so there is still exactly one analysis behind the
terminal, the markdown and this page.
"""

from __future__ import annotations

import html
import re
from typing import Any

from .report import Report, Table

#: Whole-cell matches, checked first. Short words like "ok" are only safe as exact matches — as
#: substrings they would colour any cell that happened to contain them.
EXACT_MARKERS = {
    "ok": "good",
    "skipped": "warned",
    "planned": "muted",
    "failed": "bad",
    "base": "muted",
    "-": "muted",
}

#: Cell text -> CSS class. Matched case-sensitively as a substring, longest first, so "did not fit"
#: wins over a bare "not". These are the strings the suites emit for outcomes that must not be read
#: as ordinary numbers.
CELL_MARKERS = {
    "VOID": "void",
    "did not fit": "void",
    "(not run)": "muted",
    "no data": "muted",
    "REGRESSION": "bad",
    "improvement": "good",
    "unchanged": "neutral",
    "unresolved": "neutral",
    "wins": "good",
}

#: A signed percentage, which is how every delta and drift figure is formatted.
SIGNED = re.compile(r"^[+-]\d")


def render_page(title: str, report: Report, *, subtitle: str = "", statuses: dict[str, str] | None = None) -> str:
    """One complete HTML document for `report`."""
    sections, nav = _sections(report, statuses or {})
    return _DOCUMENT.format(
        title=html.escape(title),
        subtitle=html.escape(subtitle),
        style=_STYLE,
        script=_SCRIPT,
        nav="\n".join(nav),
        body="\n".join(sections),
    )


# ---------------------------------------------------------------------------
# Blocks -> HTML
# ---------------------------------------------------------------------------


def _sections(report: Report, statuses: dict[str, str]) -> tuple[list[str], list[str]]:
    """Renders every block, collecting sidebar entries from the headings as it goes.

    Each top-level heading opens a `<details>` that runs until the next one, so a master report is
    seven foldable suites rather than one column thousands of pixels tall. They start OPEN: a report
    that hides its findings by default is a worse report, and the collapse is for after you have
    seen what is there. The toolbar has expand/collapse-all, and the sidebar's script opens whatever
    section a link points into — a nav item that scrolls to a folded section would be a dead link.
    """
    body: list[str] = []
    nav: list[str] = []
    seen: set[str] = set()
    table_index = 0
    open_suite = False
    open_part = False

    def close_part() -> None:
        nonlocal open_part
        if open_part:
            body.append("</div></details>")
            open_part = False

    for kind, payload in report.blocks:
        if kind == "heading":
            level, text = payload[0], payload[1]
            folded = payload[2] if len(payload) > 2 else False
            anchor = _anchor(text, seen)
            status = statuses.get(text, "")
            if level <= 2:
                close_part()
                if open_suite:
                    body.append("</div></details>")
                body.append(
                    f'<details class="suite" open id="{anchor}">'
                    f'<summary class="h2">{html.escape(text)}</summary><div class="suitebody">'
                )
                open_suite = True
            elif level == 3:
                close_part()
                body.append(
                    f'<details class="part"{"" if folded else " open"} id="{anchor}">'
                    f'<summary class="h3">{html.escape(text)}</summary><div class="partbody">'
                )
                open_part = True
            else:
                body.append(
                    f'<h{min(level, 4)} id="{anchor}" class="h{level}">{html.escape(text)}</h{min(level, 4)}>'
                )
            nav.append(
                f'<a class="nav-l{min(level, 3)}{" nav-" + status if status else ""}" href="#{anchor}">'
                f'{"<span class=dot></span>" if status else ""}{html.escape(text)}</a>'
            )
        elif kind == "para":
            body.append(f"<p>{_inline(payload)}</p>")
        elif kind == "table":
            rendered = _table(payload[0], table_index)
            if payload[1]:
                label = payload[1] if isinstance(payload[1], str) else "raw numbers"
                rendered = (
                    f'<details class="rawtable"><summary>{html.escape(label)}</summary>'
                    f"<div>{rendered}</div></details>"
                )
            body.append(rendered)
            table_index += 1
        elif kind == "lines":
            body.append(f'<pre class="lines">{html.escape(chr(10).join(payload))}</pre>')
        elif kind == "chart":
            body.append(payload[0])
        elif kind == "verdict":
            body.append(_verdict(payload))
        elif kind == "figures":
            # The caption is rendered. `Report.figures` has always carried one and the markdown has
            # always printed it, but this branch dropped it, so every faceted grid on the page — the
            # knob curves, the phase splits, the response-size curves — was a wall of panels titled
            # only by their regime, with the sentence saying what was plotted reaching every
            # rendering except the one people read.
            grid = f'<div class="{_gridclass(payload[0])}">{"".join(payload[0])}</div>'
            caption = payload[1] if len(payload) > 1 else ""
            if caption:
                grid += f'<p class="figcap">{_inline(caption)}</p>'
            body.append(grid)
        elif kind == "switch":
            body.append(_switch(payload, len(body)))
        elif kind == "note":
            # Closed by default, and no longer titled "how to read this": the column glossaries that
            # used to live here are now tooltips on the headings they define, beside the number
            # rather than a scroll away. What is left is the suite's own prose — the mechanism
            # behind the measurement and what previous full-database runs found — which is worth
            # keeping and is not what anyone needs open while reading a table.
            body.append(
                '<details class="note"><summary>Suite notes</summary>'
                f"<div>{_prose(payload)}</div></details>"
            )
        elif kind == "warn":
            body.append(f'<div class="warn">{_inline(payload)}</div>')
    close_part()
    if open_suite:
        body.append("</div></details>")
    return body, nav


#: Values rendered as a coloured pill rather than plain text. Booleans are the case that matters:
#: `equate_il` and `tryptic` are two of the four coordinates of the defaults grid, and as bare
#: "True"/"False" they are the hardest column to scan for.
PILLS = {"True": "t", "False": "f", "true": "t", "false": "f"}

#: A column becomes a filter chip group when it has at least this many rows to filter and its
#: values are few enough to be categories rather than measurements.
_MIN_ROWS_FOR_CHIPS = 4
_MAX_CATEGORY_VALUES = 6
#: Longer than this and it is prose (a skip reason, a verdict sentence), not a category label.
_MAX_CATEGORY_LENGTH = 20

#: Telling a measurement from a label is the whole difficulty here, and a bare integer is genuinely
#: ambiguous: 16 is a batch size (a label) while 80 is a rep count (a measurement). Three patterns
#: resolve it between them.
#:
#: A number carrying a unit, or written with a decimal point or thousands separators, is a
#: measurement — that is what keeps a `ratio` column of 0.90x/0.92x/0.94x from becoming filter
#: buttons.
#:
#: The leading sign includes `±`, which is how `report.band()` writes a noise band. Without it a
#: `floor` or `noise` column reads as a set of labels and gets filter chips — turning the width of
#: the noise into something to click, in the one column that says how much of the table is not a
#: result.
_MEASURED = re.compile(
    r"^[+\-±]?[\d,]+(\.\d+)?\s*(%|x|s|ms|min|[KMGT]B?|G)$|^[+\-±]?[\d,]*\.\d+$|^\d{1,3}(,\d{3})+$",
    re.IGNORECASE
)
#: A signed number is a comparison against something else. Never a category, whatever else shares
#: the column — which is what a "vs baseline" column of `base, +9.1%, +16.6%` looks like.
_SIGNED_VALUE = re.compile(r"^[+-]\d")
#: A bare integer with no unit and no siblings that are labels: a count or a size, not a name.
_BARE_INTEGER = re.compile(r"^\d+$")


def _table(table: Table, index: int) -> str:
    head = "".join(
        f'<th class="{_align(align)}{" hint" if header in table.tips else ""}"'
        + (f' data-tip="{html.escape(table.tips[header], quote=True)}"' if header in table.tips else "")
        + f">{html.escape(header)}</th>"
        for header, align in zip(table.headers, table.aligns)
    )
    rows = []
    # NOT `index`: that is this table's number, and shadowing it here numbered every table by its
    # last ROW instead. Tables with the same row count then shared an id, and the chips are wired by
    # id — filtering a three-row table narrowed every other three-row table in the report, and the
    # "n of m rows" count reported one of them at random.
    for position, row in enumerate(table.rows):
        cells = "".join(
            f'<td class="{_align(align)}{_marker(cell)}" data-col="{column}">{_cell(cell)}</td>'
            for column, (cell, align) in enumerate(zip(row, table.aligns + ["<"] * len(row)))
        )
        rows.append(f'<tr{" class=strong" if position in table.strong else ""}>{cells}</tr>')

    table_id = f"t{index}"
    return (
        f"{_chips(table, table_id)}"
        f'<div class="tablewrap"><table class="grid" id="{table_id}"><thead><tr>'
        f"{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )


#: The default word beside a status dot. Spelled out because colour alone is never an encoding
#: here — and because "flat" and "unresolved" are different findings that a green and an amber dot
#: would not distinguish for a reader who cannot tell the two apart.
#:
#: A tile may override the word with `"<kind>:<word>"`, since the same amber means "no value won
#: everywhere" on a knob suite and "read this before the rest of the report" on `defaults`.
STATUS_WORDS = {
    "good": "resolved",
    "flat": "flat",
    "warn": "unresolved",
}


def _verdict(payload: tuple) -> str:
    """The suite's answer, as a KPI row above the evidence for it.

    Proportional figures on the values, not `tabular-nums`: these are standalone display numbers
    rather than a column, and equal-width digits make a short one look loose at this size.
    """
    tiles, reading = payload
    cells = []
    for label, value, sub, status in tiles:
        status, _, override = status.partition(":")
        word = override or STATUS_WORDS.get(status, "")
        badge = (
            f'<span class="tilestatus s-{status}"><span class="dot"></span>{word}</span>'
            if word
            else ""
        )
        cells.append(
            f'<div class="tile"><div class="tilelabel">{html.escape(label)}</div>'
            f'<div class="tileval">{html.escape(value)}</div>'
            f'<div class="tilesub">{html.escape(sub)}{badge}</div></div>'
        )
    body = f'<div class="kpi">{"".join(cells)}</div>'
    if reading:
        body += f'<p class="verdictread">{_inline(reading)}</p>'
    return body


def _gridclass(items: list[str]) -> str:
    """Two columns for a grid of four panels or fewer, auto-fit beyond that.

    Four is the count that goes wrong on its own: auto-fit at this width lays them out three and
    one, and the orphan on the second row reads as an afterthought rather than as the fourth of
    four. The leading legend row does not count as a panel.
    """
    panels = sum(1 for item in items if "facethead" not in item[:200])
    return "figgrid cols-2" if panels <= 4 else "figgrid"


def _switch(payload: tuple, index: int) -> str:
    """A set of figures with a control choosing which colouring is shown.

    All variants ship in the page and only their visibility changes, so switching is instant and
    works from a file with no network — the same reason everything else here is inlined.
    """
    title, variants, default = payload
    group = f"sw{index}"
    buttons = "".join(
        f'<button type="button" class="chip{" on" if label == default else ""}" '
        f'data-switch="{group}" data-variant="{html.escape(label, quote=True)}">{html.escape(label)}</button>'
        for label, _ in variants
    )
    panels = "".join(
        f'<div class="panel" data-switch="{group}" data-variant="{html.escape(label, quote=True)}"'
        f'{"" if label == default else " hidden"}>'
        f'<div class="{_gridclass(svgs)}">{"".join(svgs)}</div></div>'
        for label, svgs in variants
    )
    return (
        f'<div class="switch">'
        f'<div class="chips"><span class="chipgroup"><span class="label">{html.escape(title)}</span>'
        f"{buttons}</span></div>{panels}</div>"
    )


def _cell(cell: str) -> str:
    """Cell content: a pill for the values worth spotting at a glance, escaped text otherwise."""
    stripped = cell.strip()
    if stripped in PILLS:
        return f'<span class="pill {PILLS[stripped]}">{html.escape(stripped)}</span>'
    return html.escape(cell)


def _is_category(values: list[str], distinct: list[str]) -> bool:
    """Whether a column groups rows (filter it) or measures them (don't).

    Six ways a column fails to be a category, each one a chip group that would have been noise:
    too many values to be labels; every row different, so it identifies rather than groups; values
    long enough to be prose; all values measurements, which repeat only because they were rounded to
    the same few figures; any value signed, which makes the column a comparison; and all values bare
    integers, which makes it a count.

    A column survives when at least one of its values is a word — `scalar` beside 16, `default`
    beside 96, `none` beside 5-mer, True beside False. That is what a category looks like here.
    """
    if not 2 <= len(distinct) <= _MAX_CATEGORY_VALUES:
        return False
    if len(distinct) == len(values):
        return False
    if max(len(value) for value in distinct) > _MAX_CATEGORY_LENGTH:
        return False
    if any(_SIGNED_VALUE.match(value) for value in distinct):
        return False
    if all(_MEASURED.match(value) for value in distinct):
        return False
    return not all(_BARE_INTEGER.match(value) for value in distinct)


def _chips(table: Table, table_id: str) -> str:
    """Toggle buttons for this table's categorical columns.

    Derived from the data unless the table names them: a column whose values repeat a handful of
    times is a category worth filtering on, and one where every row differs is a measurement. That
    keeps most suites free of presentation concerns and means a new suite gets chips for free.

    `Table.chips` overrides the guess. The guess is right about what COULD be filtered and has no
    way to know what is worth filtering — a suite sweeping one coordinate ends up offering chips for
    every column that repeats, and a row of six chip groups is harder to use than none.

    Selections within a column are OR, across columns AND — so `tryptic: True` plus
    `kmer: 5-mer, 6-mer` reads the way it sounds.
    """
    if len(table.rows) < _MIN_ROWS_FOR_CHIPS:
        return ""

    groups = []
    for column, header in enumerate(table.headers):
        values = [row[column].strip() for row in table.rows if column < len(row)]
        distinct = sorted({value for value in values if value and value != "-"})
        wanted = header in table.chips if table.chips is not None else _is_category(values, distinct)
        if not wanted or len(distinct) < 2:
            continue
        buttons = "".join(
            f'<button type="button" class="chip" data-table="{table_id}" data-col="{column}" '
            f'data-value="{html.escape(value, quote=True)}">{html.escape(value)}</button>'
            for value in distinct
        )
        # `data-group` names the column this group filters, so the header's search-mode control can
        # find the tryptic chips and drive them through the ordinary filter path rather than
        # inventing a second one that the row count and the free-text box would not know about.
        groups.append(
            f'<span class="chipgroup" data-group="{html.escape(header, quote=True)}">'
            f'<span class="label">{html.escape(header)}</span>{buttons}</span>'
        )

    if not groups:
        return ""
    return (
        f'<div class="chips" data-for="{table_id}">{"".join(groups)}'
        f'<button type="button" class="chip clear" data-table="{table_id}">clear</button>'
        f'<span class="count" data-for="{table_id}"></span></div>'
    )


def _marker(cell: str) -> str:
    """The CSS class for one cell, from the strings the suites emit."""
    stripped = cell.strip()
    if stripped in EXACT_MARKERS:
        return " " + EXACT_MARKERS[stripped]
    for needle in sorted(CELL_MARKERS, key=len, reverse=True):
        if needle in cell:
            return " " + CELL_MARKERS[needle]
    if SIGNED.match(stripped):
        return " pos" if stripped.startswith("+") else " neg"
    return ""


def _align(align: str) -> str:
    return "r" if align == ">" else "l"


def _prose(text: str) -> str:
    """Note bodies are plain text with indentation and bullet lists; keep both readable."""
    parts = []
    for chunk in re.split(r"\n\s*\n", text):
        stripped = chunk.strip("\n")
        if re.match(r"^\s*\*", stripped):
            items = re.split(r"\n\s*\*\s", "\n" + stripped.strip())
            bullets = "".join(f"<li>{_inline(item)}</li>" for item in items if item.strip())
            parts.append(f"<ul>{bullets}</ul>")
        else:
            parts.append(f"<p>{_inline(stripped)}</p>")
    return "".join(parts)


def _inline(text: str) -> str:
    """Escapes, then restores the two bits of markup the suites use: `code` and **bold**."""
    escaped = html.escape(text.strip())
    escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    return re.sub(r"\s*\n\s*", " ", escaped)


def _anchor(text: str, seen: set[str]) -> str:
    base = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "section"
    anchor, index = base, 2
    while anchor in seen:
        anchor, index = f"{base}-{index}", index + 1
    seen.add(anchor)
    return anchor


# ---------------------------------------------------------------------------
# Shell
# ---------------------------------------------------------------------------

_STYLE = """
:root {
  --bg: #ffffff; --panel: #f6f7f9; --ink: #14171a; --muted: #6b7280; --line: #e3e6ea;
  --accent: #2563eb; --good: #047857; --bad: #b91c1c; --warnbg: #fef6e7; --warnline: #d98b12;
  --void: #7c3aed; --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  /* Chart palette: the validated default categorical order (blue, orange, aqua, yellow, magenta).
     Assigned to entities in fixed order and never cycled — the preloaded arm is --s1 everywhere. */
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100; --s5: #e87ba4;
  /* Paired light variants. A suite comparing the same thing across two backends uses one hue per
     thing and the two lightnesses for the backends, so the pairing is visible without reading the
     legend — five hues rather than ten, and small/preloaded sits next to small/mmap instead of
     next to whatever happened to be slot 2. */
  --s1-lt: #93bcee; --s2-lt: #f5b499; --s3-lt: #86dabd; --s4-lt: #f6d180; --s5-lt: #f4bdd2;
  /* The five storage arms, one hue each. They were a single-hue ordinal ramp — the arms ARE an
     ordinal axis, how much is resident — but three lightness steps of one blue are what the eye is
     worst at: telling mmap from pprot meant judging which of two blues was darker, on the
     comparison this report exists to make, in every legend on the page. Three hues are told apart
     at a glance and at any size. The residency order still governs which arm gets which slot and
     the order they are drawn in, so the assignment is stable; it is simply no longer the colour
     that has to carry it.
     Validated as a categorical palette (lightness band, chroma floor, CVD separation, normal-vision
     floor). Light-mode green is 2.74:1 on the surface, under the 3:1 gate — the documented relief
     applies: every chart on this page sits beside the table holding the same numbers, and every
     legend is labelled. Re-run the validator before substituting any hue. */
  /* Five arms now, and the two added for `ptext` and `pmap` come from the data-viz reference
     palette rather than being invented — slots 7 (violet) and 4 (yellow). Validated as a
     categorical palette IN DRAW ORDER (mmap, pprot, ptext, pmap, preloaded), which is the pairlist
     that matters for bars and lines: worst adjacent CVD ΔE 9.1 protan, normal-vision 22.9, both
     clear. `mmap`, `pprot` and `preloaded` keep the hues they have always had; only `ptext` moved,
     from a magenta that turned out to fail the dark-mode lightness band.
     Three of the five sit under 3:1 on the light surface — the documented relief applies, and it is
     already in force here: every chart sits beside the table holding the same numbers and every
     legend is labelled. Re-run `validate_palette.js` before substituting any hue.
     Violet is close to `--void`, which is table-cell text and never a chart mark. */
  --arm-1: #2a78d6; --arm-2: #eb6834; --arm-3: #1baf7a; --arm-4: #4a3aa7; --arm-5: #eda100;
  --grid: #e3e6ea; --axisink: #6b7280;
  /* Diverging: blue and red poles reading as opposite, three steps each, neutral gray midpoint. */
  /* Sequential: one hue, more-is-further-from-the-surface. */
  --seq-1: #cde2fb; --seq-2: #9ec5f4; --seq-3: #5598e7; --seq-4: #256abf; --seq-5: #104281;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0f1115; --panel: #171a20; --ink: #e6e8eb; --muted: #9aa3af; --line: #262b33;
    --accent: #60a5fa; --good: #34d399; --bad: #f87171; --warnbg: #2a2113; --warnline: #d98b12;
    --void: #c4b5fd;
    /* Dark is selected, not flipped: the same eight hues re-stepped for the dark surface, and
       validated as a set against it. The diverging arms run light-to-dark in the other direction,
       because "more" has to move away from the surface in both modes. */
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181;
    --s1-lt: #86b6ef; --s2-lt: #eaa084; --s3-lt: #74cbaf; --s4-lt: #e2bf6d; --s5-lt: #e9a1bb;
    /* The arm ramp re-stepped for the dark surface, and reversed with it: "more resident" has to
       move AWAY from the surface in both modes, so preloaded is the lightest step here. */
    --arm-1: #3987e5; --arm-2: #d95926; --arm-3: #199e70; --arm-4: #9085e9; --arm-5: #c98500;
    --grid: #262b33; --axisink: #9aa3af;
    --seq-1: #0d366b; --seq-2: #184f95; --seq-3: #2a78d6; --seq-4: #6da7ec; --seq-5: #b7d3f6;
  }
}
:root[data-theme="dark"] {
  --bg: #0f1115; --panel: #171a20; --ink: #e6e8eb; --muted: #9aa3af; --line: #262b33;
  --accent: #60a5fa; --good: #34d399; --bad: #f87171; --warnbg: #2a2113; --warnline: #d98b12;
  --void: #c4b5fd;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181;
  --s1-lt: #86b6ef; --s2-lt: #eaa084; --s3-lt: #74cbaf; --s4-lt: #e2bf6d; --s5-lt: #e9a1bb;
  --arm-1: #3987e5; --arm-2: #d95926; --arm-3: #199e70; --arm-4: #9085e9; --arm-5: #c98500;
  --grid: #262b33; --axisink: #9aa3af;
  --seq-1: #0d366b; --seq-2: #184f95; --seq-3: #2a78d6; --seq-4: #6da7ec; --seq-5: #b7d3f6;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--ink);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}
header {
  position: sticky; top: 0; z-index: 10; background: var(--bg);
  border-bottom: 1px solid var(--line); padding: 12px 20px;
  display: flex; gap: 16px; align-items: center; flex-wrap: wrap;
}
header h1 { font-size: 16px; margin: 0; font-weight: 650; }
header .sub { color: var(--muted); font-size: 13px; font-family: var(--mono); }
header .spacer { flex: 1; }
input[type=search], button {
  font: inherit; font-size: 13px; padding: 5px 10px; border-radius: 6px;
  border: 1px solid var(--line); background: var(--panel); color: var(--ink);
}
input[type=search] { min-width: 220px; }
button { cursor: pointer; }
.layout { display: grid; grid-template-columns: 250px minmax(0, 1fr); align-items: start; }
nav {
  position: sticky; top: 57px; max-height: calc(100vh - 57px); overflow-y: auto;
  padding: 16px 8px 40px; border-right: 1px solid var(--line);
}
nav a {
  display: flex; align-items: center; gap: 7px; text-decoration: none; color: var(--ink);
  padding: 4px 10px; border-radius: 6px; font-size: 13px;
}
nav a:hover { background: var(--panel); }
nav a.active { background: var(--panel); color: var(--accent); font-weight: 600; }
.nav-l1 { font-weight: 700; margin-top: 6px; }
.nav-l3 { padding-left: 26px !important; color: var(--muted) !important; font-size: 12.5px; }
nav .dot { width: 8px; height: 8px; border-radius: 50%; background: var(--good); flex: none; }
.nav-skipped .dot { background: var(--warnline); }
.nav-failed .dot { background: var(--bad); }
main { padding: 4px 28px 120px; min-width: 0; }
h1.h1 { font-size: 24px; margin: 24px 0 8px; }
h2.h2 { font-size: 19px; margin: 34px 0 6px; padding-top: 10px; border-top: 1px solid var(--line); }
h3.h3 { font-size: 15px; margin: 22px 0 4px; color: var(--muted); text-transform: uppercase;
        letter-spacing: .04em; }
p { margin: 10px 0; max-width: 78ch; }
code { font-family: var(--mono); font-size: .92em; background: var(--panel); padding: 1px 5px;
       border-radius: 4px; }
/* `overflow-x: auto` computes overflow-y to `auto` as well, which makes this element the containing
   scroller for anything sticky inside it. A `top` offset meant for the viewport would then measure
   from this box's own top edge and draw the header over the first rows. So: scroll in both axes
   here, and stick the header to this box at `top: 0`. Tall grids get a header that stays put while
   the rows scroll under it; short ones never scroll and it simply sits where it would anyway. */
.tablewrap {
  overflow: auto; max-height: 78vh; margin: 12px 0;
  border: 1px solid var(--line); border-radius: 8px;
}
table.grid { border-collapse: separate; border-spacing: 0; width: 100%; font-family: var(--mono); font-size: 13px; }
table.grid th, table.grid td { padding: 6px 12px; white-space: nowrap; border-bottom: 1px solid var(--line); }
table.grid th {
  position: sticky; top: 0; z-index: 1; background: var(--panel); text-align: left; cursor: pointer;
  user-select: none; font-weight: 600;
  /* The border travels with the sticky header; a `border-bottom` alone scrolls away from it. */
  box-shadow: inset 0 -1px 0 var(--line);
}
table.grid th:hover { color: var(--accent); }
/* The sort affordance sits on the OUTSIDE of the column: after the label in a left-aligned column,
   before it in a right-aligned one. Always appending it would push every numeric header left by
   the glyph's width while the numbers below stayed flush to the edge — the header and its column
   would then be out by ~10px in exactly the columns where lining up matters most. */
table.grid th::after { content: " ⇅"; opacity: .3; }
table.grid th.asc::after { content: " ↑"; opacity: 1; }
table.grid th.desc::after { content: " ↓"; opacity: 1; }
table.grid th.r::after { content: none; }
table.grid th.r::before { content: "⇅ "; opacity: .3; }
table.grid th.r.asc::before { content: "↑ "; opacity: 1; }
table.grid th.r.desc::before { content: "↓ "; opacity: 1; }
table.grid tbody tr:hover { background: var(--panel); }
/* The row a table exists to point at — the swept setting in a configuration table. Weight rather
   than colour: the colour classes already mean specific outcomes, and a fourth meaning on the same
   channel would make all of them harder to read. */
table.grid tr.strong td { font-weight: 700; }
table.grid tr:last-child td { border-bottom: 0; }
/* `table.grid th` above is (0,1,2) and would out-specify a bare `th.r` (0,1,1), which left every
   numeric HEADER aligned left while its values aligned right — the column and its own heading
   visibly out of step. Match the qualifier so the two rules compare on equal terms. */
table.grid td.r, table.grid th.r { text-align: right; }
td.good { color: var(--good); }
td.bad { color: var(--bad); font-weight: 650; }
td.void { color: var(--void); font-weight: 650; }
td.muted, td.neutral { color: var(--muted); }
td.warned { color: var(--warnline); font-weight: 650; }
.pill {
  /* No min-width and no centring: a padded, centred pill pushes its text off the column's
     alignment, so the header and the values below it stop lining up. The negative margin cancels
     the horizontal padding for the same reason — the tinted box may bleed into the cell's own
     padding, but the TEXT has to start where an unpilled cell's text would. */
  display: inline-block; padding: 1px 8px; margin: 0 -8px; border-radius: 999px; font-size: 12px;
}
/* A boolean is a state, not an outcome — `--good` here would collide with the green that means
   "improvement" three columns away. Accent for on, muted for off. */
.pill.t { background: color-mix(in srgb, var(--accent) 16%, transparent); color: var(--accent); font-weight: 600; }
.pill.f { background: color-mix(in srgb, var(--muted) 18%, transparent); color: var(--muted); }
.chips { display: flex; flex-wrap: wrap; gap: 6px 14px; align-items: center; margin: 16px 0 0; }
.chipgroup { display: flex; gap: 4px; align-items: center; }
.chipgroup .label {
  font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .05em;
  margin-right: 2px;
}
.chip {
  font-family: var(--mono); font-size: 12px; padding: 2px 10px; border-radius: 999px;
  border: 1px solid var(--line); background: var(--bg); color: var(--muted); cursor: pointer;
}
.chip:hover { border-color: var(--accent); color: var(--ink); }
.chip.on { background: var(--accent); border-color: var(--accent); color: var(--bg); font-weight: 600; }
.chip.clear { border-style: dashed; }
.count { font-size: 12px; color: var(--muted); font-family: var(--mono); }
/* Charts. Text wears text tokens, never the series colour; the mark beside it carries identity. */
.chartfig { margin: 16px 0 22px; max-width: 860px; }
.chartfig svg.chart { width: 100%; height: auto; display: block; }
/* Every hoverable mark. The cursor is the affordance — without it a chart looks like a picture,
   which is what "not interactive" means to someone using it. Labels and gridlines opt out of
   pointer events so they cannot steal the hover from the mark underneath. */
svg.chart .mark { cursor: pointer; }
svg.chart .mark:hover { filter: brightness(1.12) saturate(1.1); }
svg.chart .mark.dot:hover { filter: none; }
svg.chart text, svg.chart line { pointer-events: none; }
/* The tooltip: one element for the whole page, moved to whichever mark is under the cursor. */
#charttip {
  position: fixed; z-index: 40; pointer-events: none; opacity: 0; transition: opacity .08s;
  background: var(--ink); color: var(--bg); font: 12px/1.45 var(--mono);
  padding: 5px 9px; border-radius: 6px; max-width: 340px;
  /* One `setting: value` per line — the payload is newline-separated and set as textContent. */
  white-space: pre-line;
  box-shadow: 0 2px 10px rgba(0,0,0,.28);
}
#charttip.on { opacity: 1; }
.chartfig figcaption { font-size: 12.5px; color: var(--muted); margin-top: 2px; }
/* The caption for a whole faceted grid, which titles the panels above it collectively. */
p.figcap { font-size: 12.5px; color: var(--muted); margin: 2px 0 18px; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; margin-bottom: 6px; font-size: 12.5px; }
.legend .key { display: inline-flex; align-items: center; gap: 6px; color: var(--muted); }
.legend .swatch { width: 11px; height: 11px; border-radius: 3px; flex: none; }
.heatkey .swatch { width: 20px; height: 11px; border-radius: 2px; }
text.tick, text.axis { fill: var(--axisink); font-size: 11px;
  font-family: system-ui, -apple-system, sans-serif; }
text.axistitle { fill: var(--axisink); font-size: 11px; font-style: italic;
  font-family: system-ui, -apple-system, sans-serif; }
/* A heatmap names both axes on one line, and a knob name can be long — `columns: prefetch_threshold ·
   rows: retrieval_prefetch_distance`, the widest pair this ever carried, is the case to size for.
   A step down keeps a pair that long inside a 360-wide plane. */
text.axistitle.gridaxis { font-size: 9.5px; }
text.val { fill: var(--ink); font-size: 11px; font-family: var(--mono); }
text.cellval { fill: var(--ink); font-size: 11px; font-family: var(--mono); }
text.cellval.muted { fill: var(--axisink); }
/* Ink for a value sitting on a dark fill — chosen against the fill, not against the page. */
text.cellval.on-fill { fill: #ffffff; }
:root[data-theme="dark"] text.cellval.on-fill { fill: #0f1115; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) text.cellval.on-fill { fill: #0f1115; }
}
.legend .swatch.seq { width: 22px; border-radius: 0; }
.legend .key.ramp { gap: 0; }
.legend .key.ramp .swatch:first-child { border-radius: 3px 0 0 3px; }
.legend .key.ramp .swatch:last-child { border-radius: 0 3px 3px 0; }
/* The verdict row: what the suite found, above what it found it from. */
.kpi {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
  gap: 10px; margin: 14px 0 6px;
}
.tile {
  background: var(--panel); border: 1px solid var(--line); border-radius: 8px;
  padding: 10px 12px 9px;
}
.tilelabel {
  font-size: 11px; letter-spacing: .04em; text-transform: uppercase; color: var(--muted);
}
/* Proportional figures, not tabular: a standalone display number, not a column. */
.tileval {
  font-size: 26px; line-height: 1.15; margin: 2px 0 1px; color: var(--ink);
  font-variant-numeric: proportional-nums; font-weight: 600;
}
.tilesub { font-size: 12px; color: var(--muted); display: flex; flex-wrap: wrap; gap: 4px 10px; }
/* Status wears the RESERVED status tokens, never a series colour, and always with its word beside
   the dot — the dot is the glance, the word is what makes it readable without colour. */
.tilestatus { display: inline-flex; align-items: center; gap: 5px; }
.tilestatus .dot { width: 7px; height: 7px; border-radius: 50%; background: currentColor; flex: none; }
.tilestatus.s-good { color: var(--good); }
.tilestatus.s-warn { color: var(--warnline); }
.tilestatus.s-flat { color: var(--muted); }
.verdictread { margin: 2px 0 16px; color: var(--ink); }
/* A switch: one control, several colourings of the same figures. */
.switch { margin: 14px 0 24px; }
.switch .chips { margin-top: 0; }
.figgrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(330px, 1fr)); gap: 4px 22px; }
/* Four panels want 2x2, not a row of three with one orphan underneath — a grid whose last row is
   one panel wide reads as "and also this", which is the opposite of what a small multiple says. */
.figgrid.cols-2 { grid-template-columns: repeat(auto-fit, minmax(min(100%, 430px), 1fr)); }
/* One legend for the whole grid, spanning it, above every panel. The series are identical in each
   panel by construction, so a legend per panel is the same row of words repeated down the grid. */
.facethead {
  grid-column: 1 / -1; display: flex; flex-wrap: wrap; align-items: baseline;
  gap: 4px 18px; margin: 2px 0 2px;
}
.facethead .legend { margin-bottom: 0; }
.facethead .axesnote { font-size: 12px; color: var(--muted); font-style: italic; }
/* A facet panel: a titled small multiple, not a shrunken figure. Its caption names only the panel
   — the grid's own caption carries the measurement and the units once. */
.chartfig.facet { margin: 2px 0 8px; }
.chartfig.facet figcaption {
  order: -1; margin: 0 0 1px; font-size: 12px; font-weight: 600; color: var(--ink);
}
.chartfig.facet { display: flex; flex-direction: column; }
/* Same 860px cap as a standalone figure. `max-width: none` used to be here so several figures could
   share a row, but a switch holding ONE figure then stretched it to the full column — the same
   chart wider than its neighbours purely because of how it was placed. The cap is harmless in the
   multi-figure case, where the grid track is already the narrower constraint. */
.figgrid .chartfig { margin: 8px 0 4px; max-width: 860px; }
.nofigs {
  grid-column: 1 / -1; margin: 4px 0 12px; padding: 10px 12px;
  border: 1px dashed var(--line); border-radius: 6px;
  color: var(--muted); font-size: 13px; font-family: var(--mono);
}
[hidden] { display: none !important; }
/* Section separation. A suite is a real card — its own surface, border and radius — not a rule
   across the page. Ten suites separated only by a hairline read as one continuous column, and the
   question "am I still in `kmer`?" was answered by scrolling back to the last heading. A card has
   edges, so a suite ends somewhere visible. */
details.suite {
  background: var(--bg); border: 1px solid var(--line); border-radius: 12px;
  margin: 26px 0; padding: 4px 22px 18px;
}
details.suite > summary.h2 {
  border-top: 0; margin-top: 0; padding-top: 14px;
}
h2.h2 {
  border-top: 2px solid var(--line); margin-top: 44px; padding-top: 18px;
}
/* The suite's own title bar, so the card has a head rather than starting at its first paragraph. */
details.suite > summary.h2 {
  border-bottom: 1px solid var(--line); padding-bottom: 12px;
}
details.suite:not([open]) { padding-bottom: 4px; }
details.suite:not([open]) > summary.h2 { border-bottom: 0; padding-bottom: 2px; }
/* A suite folds. The marker is drawn rather than native, so it can sit ahead of the title at a
   size that matches the heading instead of the browser's 10px triangle. */
details.suite > summary.h2 {
  font-size: 19px; font-weight: 650; margin-bottom: 6px; cursor: pointer;
  list-style: none; display: flex; align-items: baseline; gap: 10px;
}
details.suite > summary.h2::-webkit-details-marker { display: none; }
details.suite > summary.h2::before {
  content: "▾"; color: var(--muted); font-size: 13px; width: 12px; flex: none;
}
details.suite:not([open]) > summary.h2::before { content: "▸"; }
details.suite:not([open]) > summary.h2 { color: var(--muted); }
details.suite > summary.h2:hover { color: var(--accent); }
details.suite:first-of-type > summary.h2 { margin-top: 8px; border-top: 0; }
/* A subsection folds the same way, one step quieter. Folded by default for the per-regime detail,
   open for the summaries — declared by the suite, see `Report.heading`. */
/* A part is set apart within the card by a tinted band across its title, so the eye finds the
   boundaries between "summary", "the plane" and "per cell" without reading them. A dashed hairline
   was doing that job and disappeared next to the tables' own borders. */
details.part > summary.h3 {
  margin: 26px -22px 0; padding: 11px 22px; border-top: 1px solid var(--line);
  border-bottom: 1px solid var(--line); background: var(--panel);
  color: var(--ink); font-size: 15.5px; font-weight: 650; cursor: pointer;
  list-style: none; display: flex; align-items: baseline; gap: 9px;
}
details.part > summary.h3::-webkit-details-marker { display: none; }
details.part > summary.h3::before {
  content: "▾"; color: var(--muted); font-size: 12px; width: 11px; flex: none;
}
details.part:not([open]) > summary.h3::before { content: "▸"; }
details.part:not([open]) > summary.h3::after {
  content: "show"; font-size: 11px; font-weight: 500; color: var(--muted);
  border: 1px solid var(--line); border-radius: 999px; padding: 0 7px;
}
details.part:not([open]) > summary.h3 { color: var(--muted); }
details.part > summary.h3:hover { color: var(--accent); }
details.part > .partbody { padding-top: 4px; }
details.suite > .suitebody > details.part:first-child > summary.h3 { margin-top: 16px; }
/* A column heading that explains itself on hover. The dotted underline is the affordance — without
   it the explanation exists but nothing says to go looking for it. */
table.grid th.hint { text-decoration: underline dotted var(--muted); text-underline-offset: 3px; }
table.grid th.hint:hover { color: var(--accent); }
.legend.shades { margin-top: 2px; opacity: .85; }
/* The first heading inside a suite needs no rule of its own — the summary above it is the divider. */
details.suite > .suitebody > h3.h3:first-child { border-top: 0; margin-top: 14px; padding-top: 0; }
details.part { margin: 0; }
h3.h3 {
  margin-top: 30px; padding-top: 14px; border-top: 1px dashed var(--line);
  color: var(--ink); text-transform: none; letter-spacing: 0; font-size: 15.5px; font-weight: 650;
}
h4 { font-size: 13px; margin: 20px 0 2px; color: var(--muted); text-transform: uppercase;
     letter-spacing: .05em; }
td.pos { color: var(--good); }
td.neg { color: var(--bad); }
pre.lines {
  font-family: var(--mono); font-size: 13px; background: var(--panel); padding: 12px 14px;
  border-radius: 8px; overflow-x: auto; border: 1px solid var(--line); margin: 12px 0;
}
/* The exhaustive per-cell grid, folded away under the figure it belongs to. Charts are what a
   reader scans; this is what they open once a chart has pointed somewhere. */
details.rawtable { margin: 6px 0 18px; }
details.rawtable > summary {
  cursor: pointer; font-size: 12px; color: var(--muted); list-style: none;
  display: inline-flex; align-items: baseline; gap: 6px;
  border: 1px solid var(--line); border-radius: 999px; padding: 2px 12px;
}
details.rawtable > summary::-webkit-details-marker { display: none; }
details.rawtable > summary::before { content: "▸"; font-size: 10px; }
details.rawtable[open] > summary::before { content: "▾"; }
details.rawtable > summary:hover { color: var(--accent); border-color: var(--accent); }
details.rawtable .tablewrap { margin-top: 8px; }

details.note {
  margin: 14px 0; border-left: 3px solid var(--accent); background: var(--panel);
  border-radius: 0 8px 8px 0; padding: 8px 16px; max-width: 88ch;
}
details.note summary {
  cursor: pointer; font-weight: 600; font-size: 13px; color: var(--muted);
  text-transform: uppercase; letter-spacing: .04em;
}
details.note p { font-size: 14px; }
details.note ul { padding-left: 20px; margin: 8px 0; }
details.note li { margin: 5px 0; max-width: 78ch; }
.warn {
  margin: 14px 0; padding: 10px 14px; background: var(--warnbg);
  border-left: 3px solid var(--warnline); border-radius: 0 8px 8px 0; max-width: 88ch;
}
tr.hidden { display: none; }
"""

_SCRIPT = """
// Every `localStorage` access goes through these two. The report is written to be opened from
// disk — it gets copied off a benchmark server and double-clicked — and a `file://` page has no
// storage origin in every browser: Safari throws `SecurityError` on the first access. Unguarded,
// and being the first statement in this script, that throw took the whole script with it: no
// theme, no folding, no filtering, no sorting, no tooltips, no sidebar, and nothing on screen to
// say why. A remembered preference is worth having and worth nothing if it costs the page.
const recall = (key) => { try { return localStorage.getItem(key); } catch (e) { return null; } };
const remember = (key, value) => { try { localStorage.setItem(key, value); } catch (e) {} };

// Theme: follow the system until someone picks, then remember the pick.
const saved = recall('bench-theme');
if (saved) document.documentElement.dataset.theme = saved;
document.getElementById('theme').onclick = () => {
  const dark = document.documentElement.dataset.theme === 'dark'
    || (!document.documentElement.dataset.theme
        && matchMedia('(prefers-color-scheme: dark)').matches);
  const next = dark ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  remember('bench-theme', next);
};

// Chart hover. One tooltip element for the page, following the cursor over any `[data-tip]` mark.
// Native SVG title elements were the previous mechanism and read as no interaction at all: about a
// second of delay, OS styling, and shadowed by the root title naming the figure.
const tip = document.getElementById('charttip');
function moveTip(event) {
  const mark = event.target.closest && event.target.closest('[data-tip]');
  if (!mark) { tip.classList.remove('on'); return; }
  tip.textContent = mark.getAttribute('data-tip');
  tip.classList.add('on');
  // Flip to the other side of the cursor near the right or bottom edge, so the tooltip is never
  // clipped by the viewport and never sits under the pointer.
  const box = tip.getBoundingClientRect();
  const x = event.clientX + 14 + box.width > innerWidth ? event.clientX - 14 - box.width : event.clientX + 14;
  const y = event.clientY + 18 + box.height > innerHeight ? event.clientY - 12 - box.height : event.clientY + 18;
  tip.style.left = Math.max(4, x) + 'px';
  tip.style.top = Math.max(4, y) + 'px';
}
addEventListener('mousemove', moveTip, {passive: true});
addEventListener('mouseleave', () => tip.classList.remove('on'), {passive: true});
addEventListener('scroll', () => tip.classList.remove('on'), {passive: true});

// Folding. Sections start open; the button flips whichever state most of them are in, so one
// click always does something visible. A nav link into a folded section would scroll to a closed
// element, so anything targeted by the hash — or containing it — is opened first.
const foldables = () => [...document.querySelectorAll('details.suite, details.part')];
const fold = document.getElementById('fold');
fold.onclick = () => {
  const anyOpen = foldables().some(d => d.open);
  foldables().forEach(d => { d.open = !anyOpen; });
  fold.textContent = anyOpen ? 'expand all' : 'collapse all';
};
function revealHash() {
  const id = decodeURIComponent(location.hash.slice(1));
  if (!id) return;
  const target = document.getElementById(id);
  if (!target) return;
  for (let node = target; node; node = node.parentElement) {
    if (node.tagName === 'DETAILS') node.open = true;
  }
  target.scrollIntoView({block: 'start'});
}
addEventListener('hashchange', revealHash);
revealHash();

// Filtering, from two sources at once: the header's free-text box (applies to every table) and each
// table's own category chips. A row survives when it passes both. Within a column the selected
// chips are OR'd; across columns they are AND'd, which is how "tryptic: True" plus
// "kmer: 5-mer, 6-mer" reads out loud.
const filter = document.getElementById('filter');

function applyFilters() {
  const query = filter.value.trim().toLowerCase();
  document.querySelectorAll('table.grid').forEach(table => {
    // Selected chips for this table, grouped by the column they filter.
    const wanted = new Map();
    document.querySelectorAll(`.chip.on[data-table="${table.id}"]`).forEach(chip => {
      const column = chip.dataset.col;
      if (!wanted.has(column)) wanted.set(column, new Set());
      wanted.get(column).add(chip.dataset.value);
    });

    let shown = 0;
    const rows = [...table.tBodies[0].rows];
    rows.forEach(row => {
      let visible = query === '' || row.textContent.toLowerCase().includes(query);
      if (visible) {
        for (const [column, values] of wanted) {
          const cell = row.querySelector(`td[data-col="${column}"]`);
          if (!cell || !values.has(cell.textContent.trim())) { visible = false; break; }
        }
      }
      row.classList.toggle('hidden', !visible);
      if (visible) shown++;
    });

    const count = document.querySelector(`.count[data-for="${table.id}"]`);
    if (count) count.textContent = shown === rows.length ? '' : `${shown} of ${rows.length} rows`;
  });
}

filter.oninput = applyFilters;

// Search mode. Every matrix suite sweeps `tryptic`, and it splits cells rather than pooling them,
// so most of this page is a tryptic figure beside its non-tryptic twin. This picks one workload.
//
// Figures are hidden by the tag `charts._figure` put on them; tables are narrowed by turning on the
// chip they already had, so the row counts and the free-text box stay correct and a reader can
// still override one table by hand. A figure with no tag is not a duplicate of anything and is
// never hidden — it carries tryptic as an axis of its own, or does not vary in it at all.
const trypticBar = document.getElementById('tryptic');

function applyTryptic(mode) {
  trypticBar.querySelectorAll('.chip').forEach(chip => {
    chip.classList.toggle('on', chip.dataset.tryptic === mode);
  });
  document.querySelectorAll('figure[data-tryptic]').forEach(figure => {
    figure.hidden = mode !== '' && figure.dataset.tryptic !== mode;
  });
  // A grid can empty out completely: a suite only draws the contexts that resolved, and those can
  // all be one search mode. Left alone that is a heading, prose describing planes, and no planes —
  // which reads as a broken page rather than as a filtered one. Say which it is, and take the
  // shared legend down with the panels it describes.
  document.querySelectorAll('.figgrid').forEach(grid => {
    const figures = [...grid.querySelectorAll('figure[data-tryptic]')];
    const empty = figures.length > 0 && figures.every(figure => figure.hidden);
    grid.querySelectorAll('.facethead').forEach(head => { head.hidden = empty; });
    let note = grid.querySelector('.nofigs');
    if (!note) {
      note = document.createElement('p');
      note.className = 'nofigs';
      grid.append(note);
    }
    note.hidden = !empty;
    note.textContent = figures.length + ' figure(s) here, none of them '
      + (mode === 'true' ? 'tryptic' : 'non-tryptic') + '. Switch the search mode to see them.';
  });
  // The column carries bare booleans in every suite that has it, which is also what the page paints
  // as a pill. Nothing is selected for "both", which is the chip state that means "do not narrow".
  document.querySelectorAll('.chipgroup[data-group="tryptic"]').forEach(group => {
    group.querySelectorAll('.chip').forEach(chip => {
      chip.classList.toggle('on', mode !== '' && chip.dataset.value === mode);
    });
  });
  remember('bench-tryptic', mode);
  applyFilters();
}

// Non-tryptic on a first visit: it halves the page, and it is the workload the other suites hold
// their defaults against. The pick is remembered, empty string meaning "show both".
const savedTryptic = recall('bench-tryptic');
trypticBar.querySelectorAll('.chip').forEach(chip => {
  chip.onclick = () => applyTryptic(chip.dataset.tryptic);
});
applyTryptic(savedTryptic === null ? 'false' : savedTryptic);

// Figure switches: one control, several colourings of the same cells. Every variant is already in
// the page; only visibility changes.
document.querySelectorAll('.chip[data-switch]').forEach(chip => {
  chip.onclick = () => {
    const group = chip.dataset.switch, want = chip.dataset.variant;
    document.querySelectorAll(`.chip[data-switch="${group}"]`)
      .forEach(other => other.classList.toggle('on', other === chip));
    document.querySelectorAll(`.panel[data-switch="${group}"]`)
      .forEach(panel => { panel.hidden = panel.dataset.variant !== want; });
  };
});

// The header's search-mode chips are excluded: they carry `data-tryptic` and already have a handler
// that drives figures and tables together, and this loop runs later, so it would replace it.
document.querySelectorAll('.chip:not([data-switch]):not([data-tryptic])').forEach(chip => {
  chip.onclick = () => {
    if (chip.classList.contains('clear')) {
      document.querySelectorAll(`.chip.on[data-table="${chip.dataset.table}"]`)
        .forEach(other => other.classList.remove('on'));
    } else {
      chip.classList.toggle('on');
    }
    applyFilters();
  };
});

// Sort: numeric when the column parses as numbers, text otherwise. Strips the separators the
// reports use (1,234,567 / +4.2% / ±3.9% / 1.07x) so the numbers compare as numbers.
const asNumber = s => {
  const cleaned = s.replace(/[,%±x]/g, '').replace(/[()]/g, '').trim();
  const value = parseFloat(cleaned);
  return (cleaned !== '' && !isNaN(value)) ? value : null;
};
document.querySelectorAll('table.grid th').forEach(th => {
  th.onclick = () => {
    const table = th.closest('table');
    const index = [...th.parentNode.children].indexOf(th);
    const desc = !th.classList.contains('desc');
    table.querySelectorAll('th').forEach(o => o.classList.remove('asc', 'desc'));
    th.classList.add(desc ? 'desc' : 'asc');
    const body = table.querySelector('tbody');
    [...body.rows]
      .sort((a, b) => {
        const x = a.cells[index]?.textContent ?? '', y = b.cells[index]?.textContent ?? '';
        const nx = asNumber(x), ny = asNumber(y);
        // Non-numeric cells (VOID, did not fit) sort to the bottom rather than interleaving.
        if (nx === null && ny === null) return x.localeCompare(y) * (desc ? -1 : 1);
        if (nx === null) return 1;
        if (ny === null) return -1;
        return (nx - ny) * (desc ? -1 : 1);
      })
      .forEach(row => body.appendChild(row));
  };
});

// Sidebar follows the scroll position: highlight the topmost heading currently on screen.
//
// Only the elements the sidebar actually links to. `main [id]` also matched each table, whose id
// is a filter target with no nav entry, so a table scrolling into view could win `find` and blank
// the highlight — `links.get('t0')` is undefined and the `?.` swallows it.
//
// The section elements are still included, and have to be: a `<details>` carries the id its nav
// link points at. That an open one encloses its subsections is why `find` reads in document order
// and takes the outermost — the suite highlights until its first part is reached.
const headings = [...document.querySelectorAll('main details[id], main h4[id]')];
const links = new Map([...document.querySelectorAll('nav a')].map(a => [a.hash.slice(1), a]));
const onScreen = new Set();
const spy = new IntersectionObserver(entries => {
  for (const entry of entries) {
    if (entry.isIntersecting) onScreen.add(entry.target.id);
    else onScreen.delete(entry.target.id);
  }
  const current = headings.find(h => onScreen.has(h.id));
  links.forEach(a => a.classList.remove('active'));
  if (current) links.get(current.id)?.classList.add('active');
}, { rootMargin: '-60px 0px -70% 0px' });
headings.forEach(h => spy.observe(h));
"""

_DOCUMENT = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{style}</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <span class="sub">{subtitle}</span>
  <span class="spacer"></span>
  <input type="search" id="filter" placeholder="filter rows…" aria-label="Filter table rows">
  <span class="chipgroup" id="tryptic" role="group" aria-label="Search mode">
    <span class="label">tryptic</span>
    <button type="button" class="chip" data-tryptic="false">non-tryptic</button>
    <button type="button" class="chip" data-tryptic="true">tryptic</button>
    <button type="button" class="chip" data-tryptic="">both</button>
  </span>
  <button id="fold" type="button">collapse all</button>
  <button id="theme" type="button">theme</button>
</header>
<div class="layout">
  <nav>{nav}</nav>
  <main>{body}</main>
</div>
<div id="charttip" role="status" aria-live="polite"></div>
<script>{script}</script>
</body>
</html>
"""
