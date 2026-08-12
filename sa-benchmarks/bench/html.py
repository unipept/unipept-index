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
    """Renders every block, collecting sidebar entries from the headings as it goes."""
    body: list[str] = []
    nav: list[str] = []
    seen: set[str] = set()
    table_index = 0

    for kind, payload in report.blocks:
        if kind == "heading":
            level, text = payload
            anchor = _anchor(text, seen)
            status = statuses.get(text, "")
            body.append(f'<h{min(level, 4)} id="{anchor}" class="h{level}">{html.escape(text)}</h{min(level, 4)}>')
            nav.append(
                f'<a class="nav-l{min(level, 3)}{" nav-" + status if status else ""}" href="#{anchor}">'
                f'{"<span class=dot></span>" if status else ""}{html.escape(text)}</a>'
            )
        elif kind == "para":
            body.append(f"<p>{_inline(payload)}</p>")
        elif kind == "table":
            body.append(_table(payload, table_index))
            table_index += 1
        elif kind == "lines":
            body.append(f'<pre class="lines">{html.escape(chr(10).join(payload))}</pre>')
        elif kind == "chart":
            body.append(payload[0])
        elif kind == "switch":
            body.append(_switch(payload, len(body)))
        elif kind == "note":
            # Open by default: the interpretation is the point, not an appendix. Collapsible so it
            # can be folded away once read.
            body.append(
                '<details class="note" open><summary>How to read this</summary>'
                f"<div>{_prose(payload)}</div></details>"
            )
        elif kind == "warn":
            body.append(f'<div class="warn">{_inline(payload)}</div>')
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
_MEASURED = re.compile(r"^[\d,]+(\.\d+)?\s*(%|x|s|ms|min|[KMGT]B?|G)$|^[\d,]*\.\d+$|^\d{1,3}(,\d{3})+$", re.IGNORECASE)
#: A signed number is a comparison against something else. Never a category, whatever else shares
#: the column — which is what a "vs baseline" column of `base, +9.1%, +16.6%` looks like.
_SIGNED_VALUE = re.compile(r"^[+-]\d")
#: A bare integer with no unit and no siblings that are labels: a count or a size, not a name.
_BARE_INTEGER = re.compile(r"^\d+$")


def _table(table: Table, index: int) -> str:
    head = "".join(
        f'<th class="{_align(align)}">{html.escape(header)}</th>'
        for header, align in zip(table.headers, table.aligns)
    )
    rows = []
    for row in table.rows:
        cells = "".join(
            f'<td class="{_align(align)}{_marker(cell)}" data-col="{column}">{_cell(cell)}</td>'
            for column, (cell, align) in enumerate(zip(row, table.aligns + ["<"] * len(row)))
        )
        rows.append(f"<tr>{cells}</tr>")

    table_id = f"t{index}"
    return (
        f"{_chips(table, table_id)}"
        f'<div class="tablewrap"><table class="grid" id="{table_id}"><thead><tr>'
        f"{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )


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
        f'{"" if label == default else " hidden"}><div class="figgrid">{"".join(svgs)}</div></div>'
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

    Derived from the data rather than declared by the suites: a column whose values repeat a handful
    of times is a category worth filtering on, and one where every row differs is a measurement.
    That keeps the suites free of presentation concerns and means a new suite gets chips for free.

    Selections within a column are OR, across columns AND — so `tryptic: True` plus
    `kmer: 5-mer, 6-mer` reads the way it sounds.
    """
    if len(table.rows) < _MIN_ROWS_FOR_CHIPS:
        return ""

    groups = []
    for column, header in enumerate(table.headers):
        values = [row[column].strip() for row in table.rows if column < len(row)]
        distinct = sorted({value for value in values if value and value != "-"})
        if not _is_category(values, distinct):
            continue
        buttons = "".join(
            f'<button type="button" class="chip" data-table="{table_id}" data-col="{column}" '
            f'data-value="{html.escape(value, quote=True)}">{html.escape(value)}</button>'
            for value in distinct
        )
        groups.append(f'<span class="chipgroup"><span class="label">{html.escape(header)}</span>{buttons}</span>')

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
  --grid: #e3e6ea; --axisink: #6b7280;
  /* Diverging: blue and red poles reading as opposite, three steps each, neutral gray midpoint. */
  --div-mid: #f0efec;
  --div-pos-1: #b7d3f6; --div-pos-2: #5598e7; --div-pos-3: #1c5cab;
  --div-neg-1: #f7c9c8; --div-neg-2: #e88b8a; --div-neg-3: #a82b2a;
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
    --grid: #262b33; --axisink: #9aa3af;
    --div-mid: #383835;
    --div-pos-1: #184f95; --div-pos-2: #2a78d6; --div-pos-3: #9ec5f4;
    --div-neg-1: #7a1f1e; --div-neg-2: #cf4a49; --div-neg-3: #f0a9a8;
  --seq-1: #0d366b; --seq-2: #184f95; --seq-3: #2a78d6; --seq-4: #6da7ec; --seq-5: #b7d3f6;
    --seq-1: #0d366b; --seq-2: #184f95; --seq-3: #2a78d6; --seq-4: #6da7ec; --seq-5: #b7d3f6;
  }
}
:root[data-theme="dark"] {
  --bg: #0f1115; --panel: #171a20; --ink: #e6e8eb; --muted: #9aa3af; --line: #262b33;
  --accent: #60a5fa; --good: #34d399; --bad: #f87171; --warnbg: #2a2113; --warnline: #d98b12;
  --void: #c4b5fd;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181;
  --grid: #262b33; --axisink: #9aa3af;
  --div-mid: #383835;
  --div-pos-1: #184f95; --div-pos-2: #2a78d6; --div-pos-3: #9ec5f4;
  --div-neg-1: #7a1f1e; --div-neg-2: #cf4a49; --div-neg-3: #f0a9a8;
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
table.grid th::after { content: " ⇅"; opacity: .3; }
table.grid th.asc::after { content: " ↑"; opacity: 1; }
table.grid th.desc::after { content: " ↓"; opacity: 1; }
table.grid tbody tr:hover { background: var(--panel); }
table.grid tr:last-child td { border-bottom: 0; }
td.r, th.r { text-align: right; }
td.good { color: var(--good); }
td.bad { color: var(--bad); font-weight: 650; }
td.void { color: var(--void); font-weight: 650; }
td.muted, td.neutral { color: var(--muted); }
td.warned { color: var(--warnline); font-weight: 650; }
.pill {
  /* No min-width and no centring: a padded, centred pill pushes its text off the column's
     alignment, so the header and the values below it stop lining up. */
  display: inline-block; padding: 1px 8px; border-radius: 999px; font-size: 12px;
}
.pill.t { background: color-mix(in srgb, var(--good) 18%, transparent); color: var(--good); }
.pill.f { background: color-mix(in srgb, var(--muted) 20%, transparent); color: var(--muted); }
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
.chartfig figcaption { font-size: 12.5px; color: var(--muted); margin-top: 2px; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; margin-bottom: 6px; font-size: 12.5px; }
.legend .key { display: inline-flex; align-items: center; gap: 6px; color: var(--muted); }
.legend .swatch { width: 11px; height: 11px; border-radius: 3px; flex: none; }
.heatkey .swatch { width: 20px; height: 11px; border-radius: 2px; }
text.tick, text.axis { fill: var(--axisink); font-size: 11px;
  font-family: system-ui, -apple-system, sans-serif; }
text.axistitle { fill: var(--axisink); font-size: 11px; font-style: italic;
  font-family: system-ui, -apple-system, sans-serif; }
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
/* A switch: one control, several colourings of the same figures. */
.switch { margin: 14px 0 24px; }
.switch .chips { margin-top: 0; }
.figgrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(330px, 1fr)); gap: 4px 22px; }
.figgrid .chartfig { margin: 8px 0 4px; max-width: none; }
[hidden] { display: none !important; }
/* Section separation: a suite is a card, and each figure/table group inside it is set apart. */
h2.h2 {
  border-top: 2px solid var(--line); margin-top: 44px; padding-top: 18px;
}
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
// Theme: follow the system until someone picks, then remember the pick.
const saved = localStorage.getItem('bench-theme');
if (saved) document.documentElement.dataset.theme = saved;
document.getElementById('theme').onclick = () => {
  const dark = document.documentElement.dataset.theme === 'dark'
    || (!document.documentElement.dataset.theme
        && matchMedia('(prefers-color-scheme: dark)').matches);
  const next = dark ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  localStorage.setItem('bench-theme', next);
};

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

document.querySelectorAll('.chip:not([data-switch])').forEach(chip => {
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
const headings = [...document.querySelectorAll('main [id]')];
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
  <button id="theme" type="button">theme</button>
</header>
<div class="layout">
  <nav>{nav}</nav>
  <main>{body}</main>
</div>
<script>{script}</script>
</body>
</html>
"""
