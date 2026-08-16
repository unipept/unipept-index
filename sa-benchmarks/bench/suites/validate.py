"""The candidate-validation batch, one knob against a fixed background.

`shared.knob_analysis`, the same analysis `mlp` uses, pointed at a different field. The plane this
knob sits in — `mlp_batch` x `validate_batch` — is its own suite, `mlp_validate`, because a plane is
a different question from a curve and mixing them put the plane's cells on the curve.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record
from ..report import Report
from .shared import knob_analysis

MECHANISM = (
    "`validate_batch` holds candidate suffixes in flight while their text is compared, so it exists "
    "to cover one memory latency and stops paying once it does — a cliff, not a peak. It is the "
    "second of two batches: `mlp_batch` interleaves whole peptide searches, this one interleaves "
    "candidates within a single range scan. Both consume the same line-fill buffers, which is why "
    "`mlp_validate` crosses them."
)


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    knob_analysis(
        report,
        suite,
        loaded,
        knob="validate_batch",
        x_title="validate_batch (candidates per two-pass batch)",
        mechanism=MECHANISM,
    )

