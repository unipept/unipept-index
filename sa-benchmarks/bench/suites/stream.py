"""How throughput depends on the number of peptides in one call.

`shared.knob_analysis` pointed at `amount_of_peptides` — the one coordinate a caller chooses per
request, and the only thing left in the harness that varies within a process.

The knee is what this suite is for: the smallest call size that reaches within the floor of the
saturated figure. Every other suite measures at 10,000 peptides per call, where the machine is
fully occupied; nothing else says how far below that a realistic request lands.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record
from ..report import Report
from .shared import knob_analysis

MECHANISM = (
    "`search_all_matching_suffixes` splits its input with `par_chunks(MLP_BATCH)`, so a call of N "
    "peptides yields N/16 rayon tasks — 625 at the 10,000 every other suite measures, and "
    "three at a 50-peptide request. Below the knee the shortfall is idle cores, not slower code, "
    "which is why the fix is on the caller's side: batch the requests."
)


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    knob_analysis(
        report,
        suite,
        loaded,
        knob="amount_of_peptides",
        x_title="peptides per call",
        mechanism=MECHANISM,
        # The saturated call — the size every other suite measures at — so the curve reads "a
        # 100-peptide call reaches 34% of saturated throughput" rather than "+2229.7% against a
        # ten-peptide call", which is what comparing against the smallest swept value gave.
        reference=suite.defaults.get("amount"),
    )
