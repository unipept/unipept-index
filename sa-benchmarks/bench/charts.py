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
"""

from __future__ import annotations

import html
from dataclasses import dataclass

#: Categorical slots, in fixed order. Assigned to entities and never cycled or reassigned by rank —
#: the preloaded arm is slot 1 in every chart in the report.
SERIES_VARS = ("--s1", "--s2", "--s3", "--s4", "--s5")

WIDTH = 720
PAD_L, PAD_R, PAD_T, PAD_B = 66, 18, 14, 34

#: Mark specs from the reference: bars capped so the band keeps some air, 2px lines, markers big
#: enough to hover, and a 2px surface gap doing the separating between touching marks.
BAR_MAX = 24
GAP = 2
LINE_W = 2
DOT_R = 4
#: Invisible hover target around each dot. A mark you must land on dead-centre is not a hover layer;
#: the hit area has to be comfortably bigger than what is drawn.
HIT_R = 12


@dataclass
class Series:
    """One named line or bar group. `slot` fixes its colour for the whole report."""

    name: str
    values: list[float | None]
    slot: int = 0

    @property
    def color(self) -> str:
        return f"var({SERIES_VARS[self.slot % len(SERIES_VARS)]})"


# ---------------------------------------------------------------------------
# Frame
# ---------------------------------------------------------------------------


def _open(height: int, title: str) -> list[str]:
    return [
        f'<svg class="chart" viewBox="0 0 {WIDTH} {height}" role="img" '
        f'aria-label="{html.escape(title)}" preserveAspectRatio="xMidYMid meet">',
        f"<title>{html.escape(title)}</title>",
    ]


def _nice_ticks(top: float, count: int = 4) -> list[float]:
    """Round tick values covering 0..top, so the axis reads 0 / 500k / 1M rather than 0 / 437k."""
    if top <= 0:
        return [0.0]
    raw = top / count
    magnitude = 10 ** (len(str(int(raw))) - 1) if raw >= 1 else 0.1
    for multiple in (1, 2, 2.5, 5, 10):
        step = magnitude * multiple
        if step * count >= top:
            break
    return [step * index for index in range(count + 1)]


def _fmt(value: float) -> str:
    """Axis and label numbers: compact where the magnitude allows, never more precision than read."""
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M".replace(".0M", "M")
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    if value >= 10:
        return f"{value:.0f}"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _grid(height: int, ticks: list[float], scale, unit: str = "") -> list[str]:
    """Hairline solid gridlines and left-hand tick labels, both recessive."""
    out = []
    for tick in ticks:
        y = scale(tick)
        out.append(
            f'<line x1="{PAD_L}" y1="{y:.1f}" x2="{WIDTH - PAD_R}" y2="{y:.1f}" '
            f'stroke="var(--grid)" stroke-width="1" />'
        )
        out.append(
            f'<text x="{PAD_L - 8}" y="{y + 4:.1f}" text-anchor="end" class="tick">{_fmt(tick)}{unit}</text>'
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


def _figure(svg: str, caption: str, legend: str = "") -> str:
    return (
        f'<figure class="chartfig">{legend}{svg}'
        f"<figcaption>{html.escape(caption)}</figcaption></figure>"
    )


# ---------------------------------------------------------------------------
# Grouped columns — compare a handful of entities across a few categories
# ---------------------------------------------------------------------------


def grouped_columns(groups: list[str], series: list[Series], caption: str, unit: str = "") -> str:
    """Columns grouped by category, one colour per series.

    The form for "tell distinct series apart across a few categories" — production-default
    throughput per backend per peptide bucket, which is the report's headline question.
    """
    values = [value for item in series for value in item.values if value is not None]
    if not values:
        return ""
    ticks = _nice_ticks(max(values))
    height = 250
    plot_h = height - PAD_T - PAD_B
    top = ticks[-1] or 1

    def scale(value: float) -> float:
        return PAD_T + plot_h - (value / top) * plot_h

    out = _open(height, caption)
    out += _grid(height, ticks, scale, unit)

    band = (WIDTH - PAD_L - PAD_R) / max(len(groups), 1)
    bar_w = min(BAR_MAX, (band - 24) / max(len(series), 1) - GAP)
    for index, group in enumerate(groups):
        centre = PAD_L + band * (index + 0.5)
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
                    f"{item.name} · {group}: {value:,.0f}{unit}",
                )
            )
            # Value on the cap: few enough columns that every one can carry its number.
            out.append(f'<text x="{x + bar_w / 2:.1f}" y="{y - 6:.1f}" text-anchor="middle" class="val">{_fmt(value)}</text>')
        out.append(
            f'<text x="{centre:.1f}" y="{height - PAD_B + 20}" text-anchor="middle" class="axis">'
            f"{html.escape(group)}</text>"
        )
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series))


def _column(x: float, y: float, width: float, height: float, color: str, hover: str) -> str:
    """A complete column: 4px rounded data-end, square foot on the baseline, its own hover title."""
    marker = f"<title>{html.escape(hover)}</title>"
    radius = min(4.0, width / 2, max(height, 0.1))
    if height <= radius:
        return (
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{max(height, 0.5):.1f}" '
            f'fill="{color}">{marker}</rect>'
        )
    return (
        f'<path d="M{x:.1f},{y + height:.1f} V{y + radius:.1f} A{radius},{radius} 0 0 1 {x + radius:.1f},{y:.1f} '
        f'H{x + width - radius:.1f} A{radius},{radius} 0 0 1 {x + width:.1f},{y + radius:.1f} '
        f'V{y + height:.1f} Z" fill="{color}">{marker}</path>'
    )


# ---------------------------------------------------------------------------
# Lines — an ordered scale on the x axis
# ---------------------------------------------------------------------------


def lines(x_labels: list[str], series: list[Series], caption: str, unit: str = "", x_title: str = "") -> str:
    """Multi-line chart over an ordered x scale (batch size, memory ceiling, thread count).

    Points carry their own hover title, which is the per-mark tooltip layer; the table underneath
    is the exhaustive view.
    """
    values = [value for item in series for value in item.values if value is not None]
    if not values:
        return ""
    ticks = _nice_ticks(max(values))
    height = 260
    plot_h = height - PAD_T - PAD_B
    top = ticks[-1] or 1

    def y_of(value: float) -> float:
        return PAD_T + plot_h - (value / top) * plot_h

    steps = max(len(x_labels) - 1, 1)
    # Leave room on the right for the endpoint label, which rides outside the last point rather
    # than being clipped by the frame.
    span = WIDTH - PAD_L - PAD_R - 72

    def x_of(index: int) -> float:
        return PAD_L + 20 + span * index / steps

    out = _open(height, caption)
    out += _grid(height, ticks, y_of, unit)

    for index, label in enumerate(x_labels):
        out.append(
            f'<text x="{x_of(index):.1f}" y="{height - PAD_B + 20}" text-anchor="middle" class="axis">'
            f"{html.escape(label)}</text>"
        )
    if x_title:
        out.append(
            f'<text x="{WIDTH - PAD_R}" y="{height - 4}" text-anchor="end" class="axistitle">'
            f"{html.escape(x_title)}</text>"
        )

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
                f'<circle cx="{x_of(i):.1f}" cy="{y_of(value):.1f}" r="{DOT_R}" fill="{item.color}" '
                f'stroke="var(--bg)" stroke-width="2" />'
            )
            # An invisible hit target around the dot: an 8px mark you must land on dead-centre is
            # not a hover layer. This one meets the ~24px minimum without changing what is drawn.
            out.append(
                f'<circle cx="{x_of(i):.1f}" cy="{y_of(value):.1f}" r="{HIT_R}" fill="transparent">'
                f"<title>{html.escape(item.name)} · {html.escape(x_labels[i])}: {value:,.0f}{unit}</title></circle>"
            )
        # Label the endpoint only — a number on every point is chaos and goes unread.
        last = [(i, v) for i, v in enumerate(item.values) if v is not None]
        if last:
            i, value = last[-1]
            out.append(
                f'<text x="{x_of(i) + 9:.1f}" y="{y_of(value) + 4:.1f}" class="val">{_fmt(value)}</text>'
            )
    out.append("</svg>")
    return _figure("\n".join(out), caption, _legend(series))


# ---------------------------------------------------------------------------
# Horizontal stacked bars — part-to-whole with long category names
# ---------------------------------------------------------------------------


def stacked_rows(rows: list[str], series: list[Series], caption: str, unit: str = "s") -> str:
    """One horizontal stacked bar per row. Horizontal because the row names are long.

    Segments are separated by a 2px gap in the surface colour rather than a stroke, and interior
    segments are never labelled inline — the legend and the hover title carry them.
    """
    totals = [sum((item.values[i] or 0) for item in series) for i in range(len(rows))]
    if not totals or max(totals) <= 0:
        return ""
    top = max(totals)
    row_h = 26
    height = PAD_T + len(rows) * row_h + 28
    plot_w = WIDTH - PAD_L - PAD_R - 56

    out = _open(height, caption)
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
                f'<rect x="{x:.1f}" y="{y + 4:.1f}" width="{max(width - GAP, 0.5):.1f}" height="{row_h - 10}" '
                f'fill="{item.color}" rx="2">'
                f"<title>{html.escape(row)} · {html.escape(item.name)}: {value:,.2f}{unit}</title></rect>"
            )
            x += width
        out.append(
            f'<text x="{x + 8:.1f}" y="{y + row_h / 2 + 4:.1f}" class="val">{totals[index]:,.1f}{unit}</text>'
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
) -> str:
    """A grid of absolute magnitudes on one hue, light to dark.

    Used for "how fast is this backend here", where there is no polarity to encode — the diverging
    form is for the comparison between two backends, and using it here would invent a midpoint that
    means nothing.

    `low`/`high` are supplied by the caller rather than taken per grid, so the four grids of one
    peptide file share a scale and can be compared with each other.
    """
    if not cells:
        return ""
    span = (high - low) or 1.0

    def paint(value: float) -> tuple[str, bool]:
        step = 1 + min(SEQ_STEPS - 1, int((value - low) / span * SEQ_STEPS))
        return f"var(--seq-{step})", step >= 4

    return _grid_svg(
        columns,
        rows,
        {key: (paint(value), f"{value:,.0f}{unit}", hover) for key, (value, hover) in cells.items()},
        caption,
        _seq_key(low, high, unit),
    )


def heatmap(
    columns: list[str],
    rows: list[str],
    cells: dict[tuple[int, int], tuple[float, float, str]],
    caption: str,
    pos_label: str,
    neg_label: str,
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
    return _grid_svg(columns, rows, painted, caption, f'<div class="legend heatkey">{key}</div>')


def _grid_svg(
    columns: list[str],
    rows: list[str],
    cells: dict[tuple[int, int], tuple[tuple[str, bool], str, str]],
    caption: str,
    key: str,
) -> str:
    """The grid itself, shared by both heatmap flavours.

    `cells` maps (row, column) to ((fill, needs-light-ink), label, hover text). A missing key draws
    as absent — a hatched blank — because a configuration that was never measured must not look
    like one that measured zero.
    """
    if not cells:
        return ""
    cell_w = (WIDTH - PAD_L - PAD_R) / max(len(columns), 1)
    cell_h = 34
    height = PAD_T + 22 + len(rows) * cell_h + 10

    out = _open(height, caption)
    for index, column in enumerate(columns):
        out.append(
            f'<text x="{PAD_L + cell_w * (index + 0.5):.1f}" y="{PAD_T + 12}" text-anchor="middle" '
            f'class="axis">{html.escape(column)}</text>'
        )
    for r, row in enumerate(rows):
        y = PAD_T + 22 + r * cell_h
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
                f'<rect x="{x + 1:.1f}" y="{y + 1:.1f}" width="{cell_w - GAP:.1f}" height="{cell_h - GAP}" '
                f'fill="{fill}" rx="3"><title>{html.escape(hover)}</title></rect>'
            )
            # Ink chosen against the fill it sits on, not against the page.
            out.append(
                f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 + 4:.1f}" text-anchor="middle" '
                f'class="cellval{" on-fill" if dark_fill else ""}">{html.escape(label)}</text>'
            )
    out.append("</svg>")
    return _figure("\n".join(out), caption, key)


def _seq_key(low: float, high: float, unit: str) -> str:
    steps = "".join(
        f'<span class="swatch seq" style="background:var(--seq-{step})"></span>' for step in range(1, SEQ_STEPS + 1)
    )
    return (
        f'<div class="legend heatkey"><span class="key">{_fmt(low)}{unit}</span>'
        f'<span class="key ramp">{steps}</span>'
        f"<span class=\"key\">{_fmt(high)}{unit}</span></div>"
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
