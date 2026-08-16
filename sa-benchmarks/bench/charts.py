"""Inline SVG charts for the HTML report.

Generated as SVG text rather than drawn by a library, for the same reason the page inlines
everything else: these reports get scp'd off a benchmark server and opened from disk, where a CDN
script is a blank rectangle. No runtime, no dependencies, no layout shift.

Colours come from CSS custom properties defined once in the page (`--s1`…`--s5`, the diverging
poles, the chrome tokens), so a chart re-themes with the rest of the page instead of baking in hex.
The palette is the validated default from the data-viz reference: categorical slots blue/orange/
aqua/yellow/magenta, which pass the lightness, chroma, CVD-separation and normal-vision gates in
both modes; diverging blue↔red with a neutral gray midpoint.

Two rules from that reference matter more than the rest here:

* **Never a dual axis.** Throughput and fault counts are different scales, so they are different
  charts, stacked one above the other, never two y-axes on one plot.
* **Every chart sits beside the table holding the same numbers.** Three light-mode slots fall below
  3:1 against the surface, and the table view is the documented relief — so a chart is always an
  additional reading of data that is also printed, never the only place a number appears.

One rule of this project's own: a difference below the measured noise floor is not a result. The
heatmap paints those cells the neutral midpoint rather than a faint tint, so the picture cannot
claim more than the statistics do.

Marks carry their value in a `data-tip` attribute rather than an SVG `<title>` child. Native
`<title>` tooltips wait about a second, are styled by the OS, and are shadowed by the root `<title>`
that names the figure — hovering a bar looked like nothing happening. The page's script reads
`data-tip` and draws its own tooltip at the cursor, instantly, and highlights the mark under it.
Every number in a chart is also printed in the table beside it, which is the path that needs no
pointer at all.
"""

from __future__ import annotations

import html
import math
import re
from dataclasses import dataclass, field, replace

#: Categorical slots, in fixed order. Assigned to entities and never cycled or reassigned by rank —
#: the preloaded arm is slot 1 in every chart in the report.
SERIES_VARS = ("--s1", "--s2", "--s3", "--s4", "--s5")

#: Mark specs from the reference: bars capped so the band keeps some air, 2px lines, markers big
#: enough to hover, and a 2px surface gap doing the separating between touching marks.
GAP = 2
LINE_W = 2


@dataclass(frozen=True)
class Frame:
    """The geometry one chart is drawn into.

    A facet is not a full-size chart made smaller. Scaling a 720-wide viewBox into a 330px grid cell
    halves every label with it, and an 84px left margin that carries a tick column and a rotated axis
    title at full size is a quarter of the plot at facet size. So a facet gets its own viewBox —
    narrower, with margins and marks re-cut for it — and the text renders at very nearly the size it
    does in a full-width figure.
    """

    width: int
    height: int
    #: Left padding carries the tick labels AND a rotated y-axis title; bottom carries the category
    #: labels AND the x-axis title, so both are two lines deep.
    pad_l: int
    pad_r: int
    pad_t: int
    pad_b: int
    bar_max: int
    dot_r: int
    #: Invisible hover target around each dot. A mark you must land on dead-centre is not a hover
    #: layer; the hit area has to be comfortably bigger than what is drawn.
    hit_r: int
    #: Extra CSS class on the <svg>, so the stylesheet can step the type down for facets.
    css: str = ""

    @property
    def plot_w(self) -> float:
        return self.width - self.pad_l - self.pad_r

    @property
    def plot_h(self) -> float:
        return self.height - self.pad_t - self.pad_b


#: A figure that owns its own width.
FULL = Frame(width=720, height=255, pad_l=84, pad_r=18, pad_t=14, pad_b=50, bar_max=24, dot_r=4, hit_r=12)
#: One panel of a small-multiple grid. Same form, same rules, re-cut for a ~330px cell.
FACET = Frame(
    width=360, height=200, pad_l=48, pad_r=10, pad_t=12, pad_b=40,
    bar_max=14, dot_r=3, hit_r=10, css=" facet",
)
#: A plane in a grid of planes. Wider left margin than `FACET` — a heatmap's row labels are knob
#: values, not tick numbers — and no bottom band, because it has no x axis to title.
FACET_PLANE = Frame(
    width=360, height=200, pad_l=58, pad_r=8, pad_t=10, pad_b=8,
    bar_max=14, dot_r=3, hit_r=10, css=" facet",
)

WIDTH, PAD_L, PAD_R, PAD_T, PAD_B = FULL.width, FULL.pad_l, FULL.pad_r, FULL.pad_t, FULL.pad_b
BAR_MAX, DOT_R, HIT_R = FULL.bar_max, FULL.dot_r, FULL.hit_r


def series_color(slot: int, light: bool = False) -> str:
    """The one place a categorical slot becomes a colour.

    The HUE encodes what is being compared — a search phase, a knob, a length regime — and two
    marks of the same hue are the same measurement seen two ways. Slots are assigned to entities in
    fixed order and never cycled or reassigned by rank, so filtering a series out never repaints the
    survivors.

    `light` is the paired tint, for a genuine PAIR of things: two lightnesses of one hue read as
    related where two hues read as unrelated. It is not a general second axis — it has exactly two
    steps. The storage arms outgrew it (there are three) and now use `arm_color` instead.
    """
    return f"var({SERIES_VARS[slot % len(SERIES_VARS)]}{'-lt' if light else ''})"


#: The storage arms, least resident first: `mmap` maps everything, `pprot` maps everything except
#: the protein store, `preloaded` owns it all in RAM.
#:
#: The order is what fixes each arm's hue and the order they are drawn in, so an arm keeps its
#: colour whatever a filter leaves on screen. It no longer sets a lightness ramp — see `arm_color`.
ARM_ORDER = ("mmap", "pprot", "preloaded")


def by_residency(arms: list[str]) -> list[str]:
    """The arms in residency order, least resident first.

    One fixed order everywhere, so `mmap` is the first series in every legend on the page and the
    reader learns the sequence once. Declaration order varied between suites, which put the arms in
    a different order in different figures of the same report. An arm this does not name keeps its
    place at the end.
    """
    known = [arm for arm in ARM_ORDER if arm in arms]
    return known + [arm for arm in arms if arm not in ARM_ORDER]


def arm_color(arm: str) -> str:
    """A storage arm's own hue.

    One hue per arm, not three lightness steps of one. The arms are an ordinal axis — how much is
    resident — and were painted as a ramp to show it, but reading a ramp means judging which of two
    blues is darker, and telling `mmap` from `pprot` is the single comparison this report exists to
    make. Three hues separate at a glance, at facet size, and in a screenshot. The ordering still
    lives in `ARM_ORDER`, which is where it can be read exactly rather than estimated from a tint.

    An arm the scale does not name falls back to the neutral axis ink rather than to a hue that
    already means another arm.
    """
    if arm not in ARM_ORDER:
        return "var(--axisink)"
    return f"var(--arm-{ARM_ORDER.index(arm) + 1})"


@dataclass
class Series:
    """One named line or bar group. `slot` fixes its hue and `light` its lightness — see
    `series_color` for what each of the two encodes."""

    name: str
    values: list[float | None]
    slot: int = 0
    light: bool = False
    #: Set instead of `slot`/`light` when this series IS a storage arm. The arms are an ordinal
    #: scale of their own (see `arm_color`), so they do not spend a categorical hue — which leaves
    #: the hue channel free for whatever else the chart is comparing.
    arm: str = ""
    #: Label -> value pairs naming what this series is, for the hover. `Series("mixed · preloaded")`
    #: tells the legend enough and the tooltip nothing: a reader hovering a mark wants
    #: `peptides: mixed` and `backend: preloaded` on their own lines, not one compound string they
    #: have to parse back into two facts. Falls back to the name when a caller supplies none.
    tip: dict[str, str] = field(default_factory=dict)

    @property
    def color(self) -> str:
        return arm_color(self.arm) if self.arm else series_color(self.slot, self.light)


# ---------------------------------------------------------------------------
# Frame
# ---------------------------------------------------------------------------


def _tip(*pairs) -> str:
    """The hover payload for one mark: `setting: value`, one per line, in a fixed order.

    Ordered and labelled rather than run together, because a tooltip reading
    `preloaded · mixed: 1,555,098 qps` asks the reader to work out which fragment is the backend,
    which is the workload and which is the measurement — every time, for every mark.
    """
    lines_out = []
    for pair in pairs:
        # A caller with prose rather than fields — the heatmap's cell explanation — passes one
        # string, which is already a sentence and is not improved by being split into columns.
        if isinstance(pair, str):
            lines_out.append(pair)
        elif pair[1] is not None:
            lines_out.append(f"{pair[0]}: {pair[1]}")
    text = "\n".join(lines_out)
    # `&#10;` rather than a raw newline: attribute values keep literal newlines in most parsers but
    # not reliably, and an entity survives every one of them.
    return 'data-tip="' + html.escape(text, quote=True).replace("\n", "&#10;") + '"'


def _series_pairs(item: "Series", fallback: str = "series") -> list[tuple[str, str]]:
    return list(item.tip.items()) if item.tip else [(fallback, item.name)]


def _axis_pairs(x_title: str, group: str) -> list[tuple[str, str]]:
    """The x coordinate as one line per field.

    A grouped chart whose x axis is a compound — `equate_il · tryptic` against `False · True` —
    should say so as two lines, not as one pair the reader has to re-split by counting separators.
    Zipped only when both sides have the same number of parts; anything else stays a single pair.
    """
    names = [part.strip() for part in x_title.split("·")] if x_title else []
    values = [part.strip() for part in group.split("·")]
    if len(names) > 1 and len(names) == len(values):
        return list(zip(names, values))
    return [(x_title or "group", group)]


def _metric(y_title: str) -> str:
    """The measurement's name, taken from the y-axis title with its unit dropped.

    The axis reads `throughput (qps)` and the tooltip line should read `throughput: 1.2M qps` — the
    unit belongs to the value, and repeating it in the label gives `throughput (qps): 1.2M qps`.
    """
    return re.sub(r"\s*\([^)]*\)\s*$", "", y_title).strip() or "value"


def _open(frame: Frame, title: str) -> list[str]:
    return [
        f'<svg class="chart{frame.css}" viewBox="0 0 {frame.width} {frame.height}" role="img" '
        f'aria-label="{html.escape(title)}" preserveAspectRatio="xMidYMid meet">',
    ]


def _nice_ticks(top: float, count: int = 4) -> list[float]:
    """Round tick values covering 0..top, so the axis reads 0 / 500k / 1M rather than 0 / 437k.

    Trimmed to the last tick the data actually reaches. Always emitting `count` intervals of the
    chosen step overshoots whenever the step had to round up — 2.6M of data drew an axis to 4M and
    spent a third of every panel on empty air above the tallest bar.
    """
    if top <= 0:
        return [0.0]
    raw = top / count
    magnitude = 10 ** (len(str(int(raw))) - 1) if raw >= 1 else 0.1
    for multiple in (1, 2, 2.5, 5, 10):
        step = magnitude * multiple
        if step * count >= top:
            break
    needed = max(1, math.ceil(top / step - 1e-9))
    return [step * index for index in range(needed + 1)]


def _nice_span(low: float, high: float, count: int = 4) -> list[float]:
    """Round ticks covering low..high — for a scale that does NOT start at zero.

    Legal here and nowhere else in this module. A bar encodes its value as a LENGTH, so a bar chart
    that starts anywhere but zero lies about ratios. A line indexed to a reference encodes its value
    as a POSITION against that reference, and forcing it to zero spends most of the plot on empty
    space below the data: the knob curves live between 70% and 130% of the shipped value, and drawn
    from zero every one of them is a flat line across the top third of its panel.
    """
    if high <= low:
        return [low, high or 1.0]
    raw = (high - low) / count
    magnitude = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    step = magnitude
    for multiple in (1, 2, 2.5, 5, 10):
        step = magnitude * multiple
        if step * count >= high - low:
            break
    # Rounded outwards at BOTH ends, the way `_nice_ticks` rounds up at its top. Trimming to the
    # last tick below `high` instead — which is what this did — ends the axis under the data, and
    # `lines` takes its domain from `ticks[-1]`: every value above that last tick was then drawn
    # above the plot, some of them at a negative y, outside the viewBox and clipped by it.
    start = math.floor(low / step) * step
    needed = max(1, math.ceil((high - start) / step - 1e-9))
    return [start + step * index for index in range(needed + 1)]


def _fmt(value: float) -> str:
    """Axis and label numbers: compact where the magnitude allows, never more precision than read."""
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M".replace(".0M", "M")
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    if value >= 10:
        return f"{value:.0f}"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _grid(frame: Frame, ticks: list[float], scale) -> list[str]:
    """Hairline solid gridlines and left-hand tick labels, both recessive.

    Bare numbers: the unit is on the axis title, and repeating it down the ticks gives a column of
    `0 qps / 1M qps / 2M qps` where `0 / 1M / 2M` says the same thing and leaves the plot wider.
    """
    out = []
    for tick in ticks:
        y = scale(tick)
        out.append(
            f'<line x1="{frame.pad_l}" y1="{y:.1f}" x2="{frame.width - frame.pad_r}" y2="{y:.1f}" '
            f'stroke="var(--grid)" stroke-width="1" />'
        )
        out.append(
            f'<text x="{frame.pad_l - 6}" y="{y + 4:.1f}" text-anchor="end" class="tick">'
            f"{_fmt(tick)}</text>"
        )
    return out


def _text_lines(x: float, y: float, text: str, cls: str, line_h: float = 11.5, anchor: str = "middle") -> str:
    """One `<text>` carrying an embedded newline as real lines.

    SVG has no wrapping and collapses a newline inside `<text>` to a space, so a two-part label
    written with `\\n` rendered as one long run that overlapped its neighbours. Each line becomes a
    `<tspan>` re-anchored at the same x, which is the only way a multi-line label stacks.
    """
    parts = text.split("\n")
    spans = "".join(
        f'<tspan x="{x:.1f}"{"" if index == 0 else f" dy={line_h}"}>{html.escape(part)}</tspan>'
        for index, part in enumerate(parts)
    )
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" class="{cls}">{spans}</text>'


def _axis_titles(frame: Frame, x_title: str, y_title: str) -> list[str]:
    """Names both axes. A chart whose y axis is bare leaves "1.2M of what?" to the caption.

    The y title is rotated into the left margin rather than floated above the plot, so it cannot be
    mistaken for a series label, and it reads bottom-to-top as every other rotated axis does.

    Facets name them too. `FACET` is cut with the room — its left and bottom bands are two lines
    deep for exactly this — and the alternative, one note above a grid of eight panels, is a label
    the reader has to scroll back to on every panel after the first. A panel that cannot say what
    its own y axis is is not a chart, it is a shape.
    """
    out = []
    if x_title:
        out.append(
            f'<text x="{(frame.pad_l + frame.width - frame.pad_r) / 2:.0f}" y="{frame.height - 6}" '
            f'text-anchor="middle" class="axistitle">{html.escape(x_title)}</text>'
        )
    if y_title:
        mid = (frame.pad_t + frame.height - frame.pad_b) / 2
        out.append(
            f'<text x="12" y="{mid:.0f}" text-anchor="middle" class="axistitle" '
            f'transform="rotate(-90 12 {mid:.0f})">{html.escape(y_title)}</text>'
        )
    return out


def _legend(series: list[Series]) -> str:
    """Always present for two or more series; identity never rests on colour alone."""
    if len(series) < 2:
        return ""
    items = "".join(
        f'<span class="key"><span class="swatch" style="background:{item.color}"></span>'
        f"{html.escape(item.name)}</span>"
        for item in series
    )
    return f'<div class="legend">{items}</div>'


#: A caption that names the search mode it was measured under, which every chart split by `tryptic`
#: does — the coordinate is part of the panel title (`mixed · tryptic=false`) and of a plane's
#: context line. Matching the caption rather than threading the value through eight suites keeps the
#: tag in ONE place; what makes it reliable is that both spellings are explicit, which is why
#: `shared._plane_context` no longer writes a bare `tryptic` for true and nothing for false.
_TRYPTIC = re.compile(r"\btryptic=(true|false)\b")


def _figure(svg: str, caption: str, legend: str = "", css: str = "") -> str:
    """One chart, captioned, and tagged with the search mode if its caption names one.

    `data-tryptic` is what the page's header control filters on. A caption that does NOT name the
    coordinate gets no attribute and is never hidden — correct for the charts that carry tryptic as
    an axis of their own rather than as a split across figures, which would otherwise vanish under a
    filter they are not a duplicate of.
    """
    found = _TRYPTIC.search(caption)
    tag = f' data-tryptic="{found.group(1)}"' if found else ""
    return (
        f'<figure class="chartfig{css}"{tag}>{legend}{svg}'
        f"<figcaption>{html.escape(caption)}</figcaption></figure>"
    )


def _ticks_for(values: list[float], y_max: float | None) -> list[float]:
    """The y scale, either this chart's own or one imposed from outside.

    `y_max` is what makes a small-multiple grid honest. Four panels that each scale to their own
    maximum invite exactly the comparison they cannot support — the peptide length regimes differ by
    two orders of magnitude, so four independently scaled panels draw four identical-looking charts
    of wildly different numbers.
    """
    return _nice_ticks(y_max if y_max is not None else max(values))


def panel_max(series: list[Series]) -> float:
    """The tallest mark in a panel whose series sit side by side."""
    values = [value for item in series for value in item.values if value is not None]
    return max(values) if values else 0.0


def stack_max(series: list[Series]) -> float:
    """The tallest mark in a panel whose series stack — the column total, not the largest segment.

    Passing `panel_max` for a stacked form would scale the axis to the biggest single phase and let
    every column overflow the top of its own plot.
    """
    width = max((len(item.values) for item in series), default=0)
    return max(
        (sum((item.values[i] or 0) for item in series) for i in range(width)),
        default=0.0,
    )


def panel_min(series: list[Series]) -> float:
    """The lowest mark in a panel. Only meaningful for a form whose scale may leave zero."""
    values = [value for item in series for value in item.values if value is not None]
    return min(values) if values else 0.0


def facets(
    panels: list[tuple[str, list[Series]]],
    draw,
    *,
    axes: str = "",
    extent=panel_max,
    floor=None,
) -> list[str]:
    """A small-multiple grid: one legend, then one panel per group, all on ONE scale.

    `panels` is `(panel title, that panel's series)`; `draw(title, series, frame, y_max, legend)`
    renders one panel with whatever form the caller wants. The return value is the list of figure
    strings `Report.switch` expects, so a grid drops straight into an existing switch variant.

    Two properties, and both are the point:

    * **One scale.** The maximum is taken across every panel and imposed on all of them. Panels that
      each pick their own maximum draw four charts of identical shape from numbers two orders of
      magnitude apart, which is worse than not drawing them — the grid's whole promise is that the
      panels are comparable.
    * **One legend.** The series are the same in every panel by construction, so a legend per panel
      is the same row of words repeated down the grid, in the space the panels needed.

    `floor` opts a grid out of a zero baseline — pass `panel_min` for an indexed line grid, where
    the data sits in a narrow band around a reference and a zero-based axis would draw every panel
    as the same flat line near the top. Never pass it for bars.

    This is also the answer to a series count that has outgrown the palette. Eight is the token
    ceiling and this report had a chart carrying twenty-four; faceting turns it into eight panels of
    three, and three is the tier where colour alone is comfortable for everyone.
    """
    live = [(title, series) for title, series in panels if any(
        value is not None for item in series for value in item.values
    )]
    if not live:
        return []
    top = max((extent(series) for _, series in live), default=0.0)
    if top <= 0:
        return []
    bottom = min((floor(series) for _, series in live), default=0.0) if floor else None
    legend = _legend(live[0][1])
    figures = [
        draw(title, series, FACET, top, False) if bottom is None
        else draw(title, series, FACET, top, False, bottom)
        for title, series in live
    ]
    figures = [figure for figure in figures if figure]
    if not figures:
        return []
    head = f'<div class="facethead">{legend}{f"<span class=axesnote>{html.escape(axes)}</span>" if axes else ""}</div>'
    return [head, *figures]


# ---------------------------------------------------------------------------
# Grouped columns — compare a handful of entities across a few categories
# ---------------------------------------------------------------------------


def grouped_columns(
    groups: list[str],
    series: list[Series],
    caption: str,
    unit: str = "",
    x_title: str = "",
    y_title: str = "",
    *,
    frame: Frame = FULL,
    y_max: float | None = None,
    legend: bool = True,
) -> str:
    """Columns grouped by category, coloured by `series_color`.

    The form for "tell distinct series apart across a few categories" — production-default
    throughput per backend per peptide bucket, which is the report's headline question.

    No value printed on the cap. Every bar carrying its own number turns the chart into a worse
    copy of the table beside it, and the numbers crowd first at exactly the widths where there are
    most bars to compare. The number lives on hover; the table underneath is the exhaustive view.
    """
    values = [value for item in series for value in item.values if value is not None]
    if not values:
        return ""
    ticks = _ticks_for(values, y_max)
    top = ticks[-1] or 1

    def scale(value: float) -> float:
        return frame.pad_t + frame.plot_h - (value / top) * frame.plot_h

    out = _open(frame, caption)
    out += _grid(frame, ticks, scale)

    band = frame.plot_w / max(len(groups), 1)
    bar_w = max(min(frame.bar_max, (band - band / 4) / max(len(series), 1) - GAP), 1.5)
    for index, group in enumerate(groups):
        centre = frame.pad_l + band * (index + 0.5)
        span = (bar_w + GAP) * len(series)
        for position, item in enumerate(series):
            value = item.values[index] if index < len(item.values) else None
            if value is None:
                continue
            x = centre - span / 2 + position * (bar_w + GAP)
            y = scale(value)
            out.append(
                _column(
                    x,
                    y,
                    bar_w,
                    scale(0) - y,
                    item.color,
                    _tip(
                        *_series_pairs(item),
                        *_axis_pairs(x_title, group),
                        (_metric(y_title), f"{value:,.0f}{unit}"),
                    ),
                )
            )
        out.append(
            f'<text x="{centre:.1f}" y="{frame.height - frame.pad_b + 18}" text-anchor="middle" '
            f'class="axis">{html.escape(group)}</text>'
        )
    out += _axis_titles(frame, x_title, y_title)
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series) if legend else "", frame.css)


def stacked_columns(
    groups: list[str],
    bars: list[str],
    series: list[Series],
    caption: str,
    unit: str = "",
    x_title: str = "",
    y_title: str = "",
    *,
    frame: Frame = FULL,
    y_max: float | None = None,
    legend: bool = True,
    share: bool = False,
) -> str:
    """Columns grouped by category, each split into the parts that sum to it.

    `series` are the PARTS (search, retrieval), each carrying one value per (group, bar) pair,
    flattened as `groups x bars` in that order. `bars` name the sub-columns within a group — the two
    backends — so one figure shows both the total and its composition per configuration.

    Vertical, unlike `stacked_rows`, because here the categories are short labels on an x axis
    rather than long row names.

    Each segment carries its own hover, and hovering the stack's parts is how the total is read —
    nothing is printed on the bars.

    `bars` is for a genuine PAIR sharing a group — two lightnesses of one hue read as related where
    two hues read as unrelated. It has exactly two steps, so a third sub-bar has nowhere to go and
    would silently repeat the second; passing more than two is refused rather than drawn wrong. The
    way to show three is to put them on the x axis and facet whatever was there before.

    `share` normalises every column to its own total, so all of them reach 100% and the chart answers
    the part-to-whole question it was drawn for. Composition is a RATIO, and columns whose totals
    differ by two orders of magnitude — which is what the length regimes do — leave the short ones a
    few pixels tall, exactly where the split is least visible. `unit` keeps describing the raw
    values, which stay in the hover beside their percentage: the axis becomes a share, the numbers
    behind it are not thrown away.
    """
    if len(bars) > 2:
        raise ValueError(
            f"stacked_columns: {len(bars)} sub-bars, but the light variant carries only two. "
            "Put them on the x axis and facet the groups instead."
        )
    columns = [(group, bar) for group in groups for bar in bars]
    totals = [sum((item.values[i] or 0) for item in series) for i in range(len(columns))]
    if not totals or max(totals) <= 0:
        return ""

    def plotted(value: float, column: int) -> float:
        """What the segment is drawn as: its share of its own column, or the value itself."""
        total = totals[column]
        return (value / total * 100.0) if share and total else value

    ticks = _nice_ticks(100.0) if share else _ticks_for(totals, y_max)
    top = ticks[-1] or 1

    def scale(value: float) -> float:
        return frame.pad_t + frame.plot_h - (value / top) * frame.plot_h

    out = _open(frame, caption)
    out += _grid(frame, ticks, scale)

    band = frame.plot_w / max(len(groups), 1)
    bar_w = max(min(frame.bar_max, (band - band / 4) / max(len(bars), 1) - GAP), 1.5)
    for index, group in enumerate(groups):
        centre = frame.pad_l + band * (index + 0.5)
        span = (bar_w + GAP) * len(bars)
        for position, bar in enumerate(bars):
            column = index * len(bars) + position
            x = centre - span / 2 + position * (bar_w + GAP)
            cursor = scale(0)
            for item in series:
                value = item.values[column] if column < len(item.values) else None
                if not value or value <= 0:
                    continue
                segment = (scale(0) - scale(plotted(value, column)))
                cursor -= segment
                # Under `share` the hover is the only place the absolute number survives, so it
                # carries both: what this phase cost, and what fraction of the rep that was.
                reading = (
                    f"{value:,.2f}{unit} ({plotted(value, column):.1f}%)"
                    if share
                    else f"{value:,.2f}{unit}"
                )
                out.append(
                    f'<rect class="mark" x="{x:.1f}" y="{cursor:.1f}" width="{bar_w:.1f}" '
                    f'height="{max(segment - GAP, 0.5):.1f}" '
                    f'fill="{series_color(item.slot, item.light or position > 0)}" rx="2" '
                    + _tip(
                        *_axis_pairs(x_title, group),
                        ("backend", bar) if bar else ("", None),
                        ("phase", item.name),
                        (_metric(y_title), reading),
                        ("total", f"{totals[column]:,.2f}{unit}"),
                    )
                    + " />"
                )
        out.append(
            f'<text x="{centre:.1f}" y="{frame.height - frame.pad_b + 18}" text-anchor="middle" '
            f'class="axis">{html.escape(group)}</text>'
        )
    out += _axis_titles(frame, x_title, y_title)
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series) if legend else "", frame.css)


def _column(x: float, y: float, width: float, height: float, color: str, tip: str) -> str:
    """A complete column: 4px rounded data-end, square foot on the baseline, its own hover payload."""
    radius = min(4.0, width / 2, max(height, 0.1))
    if height <= radius:
        return (
            f'<rect class="mark" {tip} x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" '
            f'height="{max(height, 0.5):.1f}" fill="{color}" />'
        )
    return (
        f'<path class="mark" {tip} '
        f'd="M{x:.1f},{y + height:.1f} V{y + radius:.1f} A{radius},{radius} 0 0 1 {x + radius:.1f},{y:.1f} '
        f'H{x + width - radius:.1f} A{radius},{radius} 0 0 1 {x + width:.1f},{y + radius:.1f} '
        f'V{y + height:.1f} Z" fill="{color}" />'
    )


# ---------------------------------------------------------------------------
# Lines — an ordered scale on the x axis
# ---------------------------------------------------------------------------


def lines(
    x_labels: list[str],
    series: list[Series],
    caption: str,
    unit: str = "",
    x_title: str = "",
    y_title: str = "",
    *,
    frame: Frame = FULL,
    y_max: float | None = None,
    y_min: float | None = None,
    legend: bool = True,
    baseline: float | None = None,
) -> str:
    """Multi-line chart over an ordered x scale (batch size, memory ceiling, thread count).

    Points carry their own hover title, which is the per-mark tooltip layer, and nothing is printed
    on the plot: with six lines an endpoint label per line lands in a column of overlapping text
    right where the lines converge. The table underneath is the exhaustive view.

    `baseline` draws one reference rule — the value that ships, at 100%. A knob curve is read as
    distance from that line, and without it drawn the reader has to find the shipped value on the x
    axis and track back to its height, on every panel of the grid.
    """
    values = [value for item in series for value in item.values if value is not None]
    if not values:
        return ""
    if y_min is None:
        ticks = _ticks_for(values, y_max)
    else:
        ticks = _nice_span(y_min, y_max if y_max is not None else max(values))
    top, base = ticks[-1], ticks[0]
    # `y_span`, not `span`: the x axis below claims that name, and a `y_of` closing over the x
    # extent silently drew every panel's data squashed against its own baseline.
    y_span = (top - base) or 1

    def y_of(value: float) -> float:
        return frame.pad_t + frame.plot_h - ((value - base) / y_span) * frame.plot_h

    steps = max(len(x_labels) - 1, 1)
    # An inset at each end, so the first and last dots and their hover rings sit inside the frame
    # rather than half on the axis.
    inset = frame.hit_r + 4
    span = frame.plot_w - 2 * inset

    def x_of(index: int) -> float:
        return frame.pad_l + inset + span * index / steps

    out = _open(frame, caption)
    out += _grid(frame, ticks, y_of)

    if baseline is not None and base <= baseline <= top:
        # Dashed, unlike the grid, because this one IS a threshold — the single case the report's
        # no-dashed-lines rule exists to leave room for.
        out.append(
            f'<line x1="{frame.pad_l}" y1="{y_of(baseline):.1f}" '
            f'x2="{frame.width - frame.pad_r}" y2="{y_of(baseline):.1f}" '
            f'stroke="var(--muted)" stroke-width="1" stroke-dasharray="4 3" />'
        )

    for index, label in enumerate(x_labels):
        out.append(
            f'<text x="{x_of(index):.1f}" y="{frame.height - frame.pad_b + 18}" '
            f'text-anchor="middle" class="axis">{html.escape(label)}</text>'
        )
    out += _axis_titles(frame, x_title, y_title)

    for item in series:
        points = [(x_of(i), y_of(v)) for i, v in enumerate(item.values) if v is not None]
        if len(points) > 1:
            path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points))
            out.append(
                f'<path d="{path}" fill="none" stroke="{item.color}" stroke-width="{LINE_W}" '
                f'stroke-linejoin="round" stroke-linecap="round" />'
            )
        for i, value in enumerate(item.values):
            if value is None:
                continue
            # 2px surface ring so dots stay legible where the lines cross.
            out.append(
                f'<circle cx="{x_of(i):.1f}" cy="{y_of(value):.1f}" r="{frame.dot_r}" '
                f'fill="{item.color}" stroke="var(--bg)" stroke-width="2" />'
            )
            # An invisible hit target around the dot: an 8px mark you must land on dead-centre is
            # not a hover layer. This one meets the ~24px minimum without changing what is drawn.
            out.append(
                f'<circle class="mark dot" cx="{x_of(i):.1f}" cy="{y_of(value):.1f}" '
                f'r="{frame.hit_r}" fill="transparent" '
                + _tip(
                    *_series_pairs(item),
                    *_axis_pairs(x_title, x_labels[i]),
                    (_metric(y_title), f"{value:,.0f}{unit}"),
                )
                + " />"
            )
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series) if legend else "", frame.css)


# ---------------------------------------------------------------------------
# Horizontal stacked bars — part-to-whole with long category names
# ---------------------------------------------------------------------------


def stacked_rows(
    rows: list[str],
    series: list[Series],
    caption: str,
    unit: str = "s",
    x_title: str = "",
    row_title: str = "configuration",
    *,
    frame: Frame = FULL,
) -> str:
    """One horizontal stacked bar per row. Horizontal because the row names are long.

    Segments are separated by a 2px gap in the surface colour rather than a stroke, and nothing is
    printed on the plot — each segment's hover carries its own value and the row's total, the same
    rule the other forms follow.
    """
    totals = [sum((item.values[i] or 0) for item in series) for i in range(len(rows))]
    if not totals or max(totals) <= 0:
        return ""
    top = max(totals)
    row_h = 26
    # Height follows the row count rather than the frame's: a three-row chart in a 255px box is
    # mostly empty box.
    frame = replace(frame, height=frame.pad_t + len(rows) * row_h + (40 if x_title else 28))
    height = frame.height
    PAD_L, WIDTH, PAD_R, PAD_T = frame.pad_l, frame.width, frame.pad_r, frame.pad_t
    plot_w = frame.plot_w - 56

    out = _open(frame, caption)
    for index, row in enumerate(rows):
        y = PAD_T + index * row_h
        x = PAD_L
        out.append(
            f'<text x="{PAD_L - 8}" y="{y + row_h / 2 + 4:.1f}" text-anchor="end" class="axis">'
            f"{html.escape(row)}</text>"
        )
        for item in series:
            value = item.values[index] or 0
            if value <= 0:
                continue
            width = value / top * plot_w
            out.append(
                f'<rect class="mark" x="{x:.1f}" y="{y + 4:.1f}" width="{max(width - GAP, 0.5):.1f}" '
                f'height="{row_h - 10}" fill="{item.color}" rx="2" '
                + _tip(
                    (row_title, row),
                    ("phase", item.name),
                    (_metric(x_title), f"{value:,.2f}{unit}"),
                    ("total", f"{totals[index]:,.2f}{unit}"),
                )
                + " />"
            )
            x += width
    if x_title:
        out.append(
            f'<text x="{(PAD_L + WIDTH - PAD_R) / 2:.0f}" y="{height - 4}" text-anchor="middle" '
            f'class="axistitle">{html.escape(x_title)}</text>'
        )
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series))


# ---------------------------------------------------------------------------
# Diverging heatmap — which of two arms is ahead, across a grid of configurations
# ---------------------------------------------------------------------------

#: Bin edges for |delta| as a percentage, above the per-cell noise floor. Discrete steps rather than
#: a continuous ramp: the eye reads three levels reliably and cannot be fooled into ranking two
#: cells that differ by a percent.
HEAT_BINS = (10.0, 25.0)


#: Steps of the one-hue sequential ramp, for absolute magnitude. Light-to-dark in light mode and
#: the other way in dark: "more" always moves away from the surface.
SEQ_STEPS = 5


def sequential_heatmap(
    columns: list[str],
    rows: list[str],
    cells: dict[tuple[int, int], tuple[float, str]],
    caption: str,
    low: float,
    high: float,
    unit: str = "",
    *,
    frame: Frame = FULL,
    legend: bool = True,
    ridge: list[tuple[int, int]] | None = None,
    floor: float | None = None,
    compact: bool = False,
    x_title: str = "",
    y_title: str = "",
) -> str:
    """A grid of absolute magnitudes on one hue, light to dark.

    Used for "how fast is this backend here", where there is no polarity to encode — the diverging
    form is for the comparison between two backends, and using it here would invent a midpoint that
    means nothing.

    `low`/`high` are the caller's, which lets one scale serve several grids that must be compared —
    the four grids of one peptide file — or one grid be scaled to itself, which is what a knob plane
    needs: its cells are one backend at one context, and a scale wide enough to hold every backend
    would paint the whole plane a single step and hide the shape it was drawn to show.

    That second use is only honest if the reader can see how wide the scale is, because a ramp
    stretched over a plane that varies by less than its own noise paints five confident steps of
    nothing. `floor` puts that comparison in the key — the plane's spread against the noise floor —
    and is why this form may be scaled per grid at all. `compact` labels the values in the shortened
    form the key already uses, for cells too narrow to hold a full number.
    """
    if not cells:
        return ""
    span = (high - low) or 1.0

    def paint(value: float) -> tuple[str, bool]:
        step = 1 + min(SEQ_STEPS - 1, int((value - low) / span * SEQ_STEPS))
        return f"var(--seq-{step})", step >= 4

    def label(value: float) -> str:
        # Compact cells carry a bare number, for the reason the tick labels do: the unit is already
        # on the key at both ends of the ramp, and `2.3M qps` in all sixteen cells spends the width
        # that told them apart on a word that never changes.
        return _fmt(value) if compact else f"{value:,.0f}{unit}"

    return _grid_svg(
        columns,
        rows,
        {key: (paint(value), label(value), hover) for key, (value, hover) in cells.items()},
        caption,
        _seq_key(low, high, unit, floor) if legend else "",
        frame,
        ridge,
        x_title,
        y_title,
    )


def heatmap(
    columns: list[str],
    rows: list[str],
    cells: dict[tuple[int, int], tuple[float, float, str]],
    caption: str,
    pos_label: str,
    neg_label: str,
    *,
    frame: Frame = FULL,
    legend: bool = True,
    ridge: list[tuple[int, int]] | None = None,
    x_title: str = "",
    y_title: str = "",
) -> str:
    """A grid of signed deltas, diverging from a neutral midpoint.

    `cells` maps (row, column) to (delta percent, noise floor percent, hover text). A positive delta
    takes the blue pole and is labelled `pos_label`; a negative one takes red and `neg_label`. A
    delta that does not clear its own floor is painted the neutral midpoint — the same rule the
    tables apply, so the picture and the numbers cannot disagree.
    """
    if not cells:
        return ""
    painted = {}
    for key, (delta, floor, hover) in cells.items():
        resolved = abs(delta) > floor
        painted[key] = ((_heat_color(delta, floor), _heat_is_dark(delta, floor)), f"{delta:+.0f}%", hover)
        if not resolved:
            # A cell inside its floor keeps its number but wears the muted ink, so it reads as
            # "measured, not resolved" rather than as a small effect.
            painted[key] = ((("var(--div-mid)"), False), f"{delta:+.0f}%", hover)

    # The key reads as a scale, ends inward, with the neutral step named rather than left to be
    # guessed at — it is the step that says "this run cannot tell these two apart".
    key = "".join(
        f'<span class="key"><span class="swatch" style="background:{fill}"></span>{label}</span>'
        for fill, label in (
            ("var(--div-neg-3)", f"{html.escape(neg_label)} by &gt;25%"),
            ("var(--div-neg-1)", "past the floor"),
            ("var(--div-mid)", "within the noise floor"),
            ("var(--div-pos-1)", "past the floor"),
            ("var(--div-pos-3)", f"{html.escape(pos_label)} by &gt;25%"),
        )
    )
    return _grid_svg(
        columns,
        rows,
        painted,
        caption,
        f'<div class="legend heatkey">{key}</div>' if legend else "",
        frame,
        ridge,
        x_title,
        y_title,
    )


def _grid_svg(
    columns: list[str],
    rows: list[str],
    cells: dict[tuple[int, int], tuple[tuple[str, bool], str, str]],
    caption: str,
    key: str,
    frame: Frame = FULL,
    ridge: list[tuple[int, int]] | None = None,
    x_title: str = "",
    y_title: str = "",
) -> str:
    """The grid itself, shared by both heatmap flavours.

    `cells` maps (row, column) to ((fill, needs-light-ink), label, hover text). A missing key draws
    as absent — a hatched blank — because a configuration that was never measured must not look
    like one that measured zero.

    `ridge` is the best cell per column, as (row, column) pairs. Drawn as one polyline across the
    plane, it answers the question the plane exists for at a glance: a straight ridge means the two
    knobs are separable and a one-at-a-time sweep would have found the same answer; a ridge that
    bends is the interaction such a sweep cannot see.
    """
    if not cells:
        return ""
    # Both axes are named on the figure, so the plane makes sense without its caption. That costs a
    # line above the column labels and a rotated title in the left margin, and the frame grows by
    # exactly that much rather than stealing it from the cells.
    head_lines = max((column.count("\n") + 1 for column in columns), default=1)
    PAD_R = frame.pad_r
    PAD_L = frame.pad_l
    # Both axes are named in ONE line across the top, rather than as a centred x title and a rotated
    # y title. A plane is short — three rows is ~130px — and a knob's name is long:
    # `retrieval_prefetch_distance` rotated into that height runs off both ends of the figure. On
    # one line it fits the full width, and it can say which way each axis runs, which a rotated
    # word beside the row labels still leaves the reader to infer.
    axis_line = " · ".join(
        part for part in (
            f"columns: {x_title}" if x_title else "",
            f"rows: {y_title}" if y_title else "",
        ) if part
    )
    PAD_T = frame.pad_t + (14 if axis_line else 0)
    cell_w = (frame.width - PAD_L - PAD_R) / max(len(columns), 1)
    cell_h = 34 if frame is FULL else 26
    label_band = 12 + (head_lines - 1) * 11.5
    frame = replace(frame, height=int(PAD_T + label_band + 10 + len(rows) * cell_h + 10), pad_t=PAD_T)
    height = frame.height

    out = _open(frame, caption)
    if axis_line:
        out.append(
            f'<text x="2" y="{frame.pad_t - 15:.1f}" text-anchor="start" '
            f'class="axistitle gridaxis">{html.escape(axis_line)}</text>'
        )
    for index, column in enumerate(columns):
        out.append(
            _text_lines(PAD_L + cell_w * (index + 0.5), PAD_T + 10, column, "axis")
        )
    # First row sits below the column labels, however many lines those took.
    top = PAD_T + label_band + 10
    for r, row in enumerate(rows):
        y = top + r * cell_h
        out.append(
            f'<text x="{PAD_L - 8}" y="{y + cell_h / 2 + 4:.1f}" text-anchor="end" class="axis">'
            f"{html.escape(row)}</text>"
        )
        for c in range(len(columns)):
            entry = cells.get((r, c))
            x = PAD_L + cell_w * c
            if entry is None:
                out.append(
                    f'<rect x="{x + 1:.1f}" y="{y + 1:.1f}" width="{cell_w - GAP:.1f}" '
                    f'height="{cell_h - GAP}" fill="var(--grid)" opacity=".35" rx="3" />'
                )
                out.append(
                    f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 + 4:.1f}" text-anchor="middle" '
                    f'class="cellval muted">not run</text>'
                )
                continue
            (fill, dark_fill), label, hover = entry
            out.append(
                f'<rect class="mark" x="{x + 1:.1f}" y="{y + 1:.1f}" width="{cell_w - GAP:.1f}" '
                f'height="{cell_h - GAP}" fill="{fill}" rx="3" ' + _tip(hover) + " />"
            )
            # Ink chosen against the fill it sits on, not against the page. `pointer-events: none`
            # in the stylesheet keeps the label from stealing the hover off the cell under it.
            out.append(
                f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 + 4:.1f}" text-anchor="middle" '
                f'class="cellval{" on-fill" if dark_fill else ""}">{html.escape(label)}</text>'
            )
    if ridge and len(ridge) > 1:
        points = " ".join(
            f"{PAD_L + cell_w * (c + 0.5):.1f},{top + r * cell_h + cell_h / 2:.1f}"
            for r, c in ridge
        )
        # A surface-coloured casing under the line, so it stays visible over both poles of the
        # diverging fill without a colour of its own that could read as another series.
        out.append(
            f'<polyline points="{points}" fill="none" stroke="var(--bg)" stroke-width="5" '
            f'stroke-linejoin="round" stroke-linecap="round" opacity=".75" />'
        )
        out.append(
            f'<polyline points="{points}" fill="none" stroke="var(--ink)" stroke-width="2" '
            f'stroke-linejoin="round" stroke-linecap="round" />'
        )
    out.append("</svg>")
    return _figure("\n".join(out), caption, key, frame.css)


def _seq_key(low: float, high: float, unit: str, floor: float | None = None) -> str:
    """The ramp's ends, and — where the caller scaled it to its own data — how wide that scale is.

    A ramp says "more" and "less"; it cannot say "by how much", and stretched across a grid whose
    values differ by a percent it draws five convincing steps of noise. Spelling the spread out
    beside the noise floor is what lets a reader tell the two apart at a glance: a spread inside the
    floor means the ramp is measurement scatter, whatever shape it appears to have.
    """
    steps = "".join(
        f'<span class="swatch seq" style="background:var(--seq-{step})"></span>' for step in range(1, SEQ_STEPS + 1)
    )
    # A process with too few reference cells to resolve anything reports its floor as NaN, and every
    # comparison against NaN is false — which would quietly print "past the floor" for a plane that
    # could not measure one. No floor, no verdict.
    note = ""
    if floor is not None and math.isfinite(floor) and low > 0:
        spread = (high - low) / low * 100.0
        verdict = "inside the floor" if spread <= floor else "past the floor"
        note = f'<span class="key">spread {spread:.1f}% · floor ±{floor:.1f}% — {verdict}</span>'
    return (
        f'<div class="legend heatkey"><span class="key">{_fmt(low)}{unit}</span>'
        f'<span class="key ramp">{steps}</span>'
        f'<span class="key">{_fmt(high)}{unit}</span>{note}</div>'
    )


def _heat_is_dark(delta: float, floor: float) -> bool:
    """Whether this cell's fill is dark enough that the value needs light ink."""
    return abs(delta) > floor and 1 + sum(abs(delta) > edge for edge in HEAT_BINS) >= 3


def _heat_color(delta: float, floor: float) -> str:
    """Neutral inside the noise floor, then three steps out along whichever pole applies."""
    if abs(delta) <= floor:
        return "var(--div-mid)"
    step = 1 + sum(abs(delta) > edge for edge in HEAT_BINS)
    return f"var(--div-{'pos' if delta > 0 else 'neg'}-{step})"
