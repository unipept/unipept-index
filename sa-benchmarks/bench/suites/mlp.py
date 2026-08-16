"""The cross-query MLP batch curve, per length regime and per backend.

The curve, not the winner, is the output. A peak is a number anyone can read off a table; the shape
around it is what says whether memory-level parallelism was the constraint at all — a curve that
rises and plateaus, one that is flat, and one that rises then falls call for three different next
moves, and only the last of them is "change the default".

Everything below `analyse` is `shared.knob_analysis`, which `validate` and any later knob suite use
the same way: the three differ in which field is the subject, not in how a knob is read.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record
from ..report import Report
from .shared import knob_analysis

MECHANISM = (
    "A suffix-array probe is a random read whose address depends on the previous probe's result, so "
    "one search is a dependent chain of cache misses with nothing to overlap. Batching B searches "
    "per rayon task hands the memory system B independent chains at once, and the win is however "
    "much of that latency the hardware can then hide. A flat curve therefore says something OTHER "
    "than memory-level parallelism is the wall — on the full database that has been the measured "
    "answer, because a 242 GB working set misses the TLB on every access and there are only a "
    "handful of hardware page walkers per core."
)


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    knob_analysis(
        report,
        suite,
        loaded,
        knob="mlp_batch",
        x_title="mlp_batch (peptides interleaved per task)",
        mechanism=MECHANISM,
    )
