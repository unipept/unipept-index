#!/usr/bin/env python3
"""
Plot sa-benchmarks .jsonl result files.

Each file produces one line on the chart.  X = run number (1-based),
Y = the chosen result field.  The legend label comes from the first
record's `label` field.

Usage:
    python plot.py [--field FIELD] [--output PATH] [--title TITLE] FILE [FILE ...]

Examples:
    python plot.py results/run-a.jsonl results/run-b.jsonl
    python plot.py --field total_duration_ns --output compare.png a.jsonl b.jsonl
    python plot.py --field search_duration_ns --title "Search latency" a.jsonl
"""

import argparse
import json
import sys

import matplotlib.pyplot as plt


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot sa-benchmarks .jsonl results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("files", nargs="+", help=".jsonl benchmark result files")
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
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(10, 6))
    plotted = 0

    for path in args.files:
        records = load_records(path)
        if not records:
            print(f"Warning: {path} is empty — skipping", file=sys.stderr)
            continue

        label = records[0].get("label", path)
        values = extract_values(records, args.field, path)

        x = list(range(1, len(values) + 1))
        ax.plot(x, values, marker="o", label=label)
        plotted += 1

    if plotted == 0:
        print("Error: no data to plot", file=sys.stderr)
        sys.exit(1)

    ax.set_xlabel("Run")
    ax.set_ylabel(args.field)
    ax.set_title(args.title or args.field)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
