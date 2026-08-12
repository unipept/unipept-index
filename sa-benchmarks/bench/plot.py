#!/usr/bin/env python3
"""
Plot sa-benchmarks .jsonl result files.

Each file produces one line on the chart.  X = run number (1-based),
Y = the chosen result field.  The legend label comes from the first
record's `label` field.

Usage (file mode):
    python plot.py [OPTIONS] FILE [FILE ...]

Usage (directory mode):
    python plot.py [OPTIONS] --dirs DIR1 DIR2

    Plots only the .jsonl files that appear in both directories.
    Files from the same basename share a color; the second directory
    gets a lighter tint of that color.

Supported --stat values:
    mean        arithmetic mean
    median      median (P50)
    P<N>        Nth percentile, e.g. P95, P99

Examples:
    python plot.py results/run-a.jsonl results/run-b.jsonl
    python plot.py --field total_duration_ns --output compare.png a.jsonl b.jsonl
    python plot.py --stat mean --stat P95 --hide-data a.jsonl b.jsonl
    python plot.py --stat P95 --dirs before/ after/
"""

import argparse
import glob
import itertools
import json
import os
import statistics
import sys

import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

_STAT_LINESTYLES = ["--", "-.", ":", (0, (3, 1, 1, 1, 1, 1))]


def load_records(path: str) -> list[dict]:
    records = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: {path}:{lineno}: skipping invalid JSON — {e}", file=sys.stderr)
    return records


def extract_values(records: list[dict], field: str, path: str) -> list:
    values = []
    for i, rec in enumerate(records, 1):
        val = rec.get("result", {}).get(field)
        if val is None:
            print(
                f"Warning: {path} run {i}: field '{field}' not found in result — plotting None",
                file=sys.stderr,
            )
        values.append(val)
    return values


def compute_stat(values: list, stat: str, path: str):
    """Return the scalar stat value for *values*, or None on error."""
    clean = [v for v in values if v is not None]
    if not clean:
        return None
    s = stat.lower()
    if s == "mean":
        return statistics.mean(clean)
    if s == "median":
        return statistics.median(clean)
    if s.startswith("p") and s[1:].isdigit():
        n = int(s[1:])
        if not 1 <= n <= 99:
            print(f"Warning: {path}: percentile must be between 1 and 99, got {n} — skipping", file=sys.stderr)
            return None
        return statistics.quantiles(clean, n=100)[n - 1]
    print(f"Warning: {path}: unknown stat '{stat}' — skipping", file=sys.stderr)
    return None


def lighten_color(color, factor: float = 0.5):
    """Return a lighter variant of *color* by blending toward white."""
    r, g, b = mcolors.to_rgb(color)
    return 1 - factor * (1 - r), 1 - factor * (1 - g), 1 - factor * (1 - b)


def find_common_files(dirs: list[str]) -> list[tuple[str, list[str]]]:
    """Return [(basename, [path_in_dir1, path_in_dir2, ...]), ...] for files present in all dirs."""
    dir_maps = []
    for d in dirs:
        mapping = {os.path.basename(p): p for p in glob.glob(os.path.join(d, "*.jsonl"))}
        dir_maps.append(mapping)

    common = set(dir_maps[0].keys())
    for m in dir_maps[1:]:
        common &= m.keys()

    skipped = set(dir_maps[0].keys()) - common
    for basename in sorted(skipped):
        print(f"Warning: '{basename}' not present in all directories — skipping", file=sys.stderr)

    return [(basename, [m[basename] for m in dir_maps]) for basename in sorted(common)]


def plot_series(ax, x, values, color, label, hide_data, stats, stat_linestyles):
    """Plot one data series (and its stat lines) onto *ax*."""
    if not hide_data:
        ax.plot(x, values, marker="o", label=label, color=color)

    for stat, linestyle in zip(stats, stat_linestyles):
        val = compute_stat(values, stat, label)
        if val is not None:
            ax.axhline(val, linestyle=linestyle, color=color)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot sa-benchmarks .jsonl results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("files", nargs="*", help=".jsonl benchmark result files")
    parser.add_argument(
        "--dirs",
        nargs=2,
        metavar="DIR",
        help="Compare matching .jsonl files across two directories instead of listing files",
    )
    parser.add_argument(
        "--field",
        default="throughput_qps",
        help="Result field to plot on the Y axis (default: throughput_qps)",
    )
    parser.add_argument(
        "--output",
        help="Save the plot to this file instead of showing it interactively",
    )
    parser.add_argument("--title", help="Plot title (defaults to the field name)")
    parser.add_argument(
        "--stat",
        metavar="STAT",
        action="append",
        help="Draw a horizontal line per file for this stat (mean, median, P95, …). "
             "Can be repeated for multiple stats.",
    )
    parser.add_argument(
        "--hide-data",
        action="store_true",
        help="Hide the individual data point lines; only stat lines are shown (if --stat is given)",
    )
    args = parser.parse_args()

    if args.dirs and args.files:
        parser.error("Provide either positional files or --dirs, not both")
    if not args.dirs and not args.files:
        parser.error("Provide either positional files or --dirs DIR1 DIR2")

    stats = args.stat or []
    stat_linestyles = [_STAT_LINESTYLES[i % len(_STAT_LINESTYLES)] for i in range(len(stats))]

    fig, ax = plt.subplots(figsize=(10, 6))
    plotted = 0
    color_cycle = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    if args.dirs:
        dir_names = [os.path.basename(os.path.normpath(d)) for d in args.dirs]
        common_files = find_common_files(args.dirs)

        for basename, paths in common_files:
            base_color = next(color_cycle)
            colors = [base_color, lighten_color(base_color)]

            for path, color, dir_name in zip(paths, colors, dir_names):
                records = load_records(path)
                if not records:
                    print(f"Warning: {path} is empty — skipping", file=sys.stderr)
                    continue
                label = f"{records[0].get('label', basename)} ({dir_name})"
                values = extract_values(records, args.field, path)
                x = list(range(1, len(values) + 1))
                plot_series(ax, x, values, color, label, args.hide_data, stats, stat_linestyles)

            plotted += 1
    else:
        for path in args.files:
            records = load_records(path)
            if not records:
                print(f"Warning: {path} is empty — skipping", file=sys.stderr)
                continue

            label = records[0].get("label", path)
            values = extract_values(records, args.field, path)
            x = list(range(1, len(values) + 1))
            color = next(color_cycle)
            plot_series(ax, x, values, color, label, args.hide_data, stats, stat_linestyles)
            plotted += 1

    if plotted == 0:
        print("Error: no data to plot", file=sys.stderr)
        sys.exit(1)

    # Add one legend entry per stat type (gray, showing the linestyle)
    handles, labels = ax.get_legend_handles_labels()
    for stat, linestyle in zip(stats, stat_linestyles):
        handles.append(mlines.Line2D([], [], color="gray", linestyle=linestyle, label=stat))
        labels.append(stat)
    ax.legend(handles=handles, labels=labels, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.set_xlabel("Run")
    ax.set_ylabel(args.field)
    ax.set_title(args.title or args.field)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.subplots_adjust(right=0.75)

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
