"""Per-suite analysis: how each suite's numbers become a table, and how to read that table.

Each module here exposes one function:

    analyse(report: Report, suite: Suite, records: list[Record], out_dir: Path) -> None

It appends to `report`; it never prints. That is what lets the terminal and `report.md` show the
same analysis, and what stops the master run from growing a second, divergent copy of every table.

The prose each module attaches with `report.note(...)` is carried over from the driver scripts these
suites replace. It is not commentary — it is the part that says which column to read first and which
differences are not answers, and it was the most easily lost thing in this consolidation.
"""

from __future__ import annotations

from importlib import import_module
from typing import Callable


def analysis_for(suite_name: str) -> Callable:
    """The `analyse` function for a suite, or a fallback that at least shows the raw cells."""
    try:
        return import_module(f"{__name__}.{suite_name}").analyse
    except ModuleNotFoundError:
        from .generic import analyse

        return analyse
