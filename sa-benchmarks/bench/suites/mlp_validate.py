"""The two batch knobs as a plane, and the one question a pair of curves cannot answer.

`mlp` and `validate` each measure one of these knobs against a fixed background. Their curves are
slices through this surface, taken at the other knob's shipped value — which is the whole picture
only if the two do not interact. They plausibly do: both hold in-flight misses, and the misses
outstanding at once go as their product against a fixed line-fill budget per core.

So the output is a ridge readout, not a throughput table: does the best `mlp_batch` stay the same at
every `validate_batch`? A yes makes the two curves sufficient. A no makes them each one slice, and
the pair has to be chosen together.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Suite
from ..records import Record
from ..report import Report
from .shared import held_and_swept, knob_planes, resolution_table


def analyse(report: Report, suite: Suite, loaded: list[Record], out_dir: Path) -> None:
    # The planes are drawn into a scratch report first: this suite's verdict IS a property of the
    # ridges, so it cannot be stated until they are computed, and it has to be printed before them.
    planes = Report()
    moved = knob_planes(planes, loaded)
    _verdict_tiles(report, moved)

    report.heading("summary", level=3)
    held_and_swept(report, loaded)
    resolution_table(report, loaded)
    report.para(
        "One plane per context. Each cell is that pair of batch sizes against the shipped pair, and "
        "a cell inside its floor is painted the neutral midpoint rather than a faint tint — so a "
        "plane that looks flat IS the finding, not a rendering failure."
    )

    report.heading("the planes", level=3)
    report.extend(planes)
    report.para(
        "The ridge drawn on each plane is the result. If the best `mlp_batch` is not the same at "
        "every `validate_batch`, the two knobs interact and neither `mlp` nor `validate` can be "
        "read on its own; if it holds, both curves are sufficient and this suite need not be run "
        "again until the search path changes."
    )

    if suite.notes:
        report.note(suite.notes)


def _verdict_tiles(report: Report, moved: list[bool]) -> None:
    """Whether the two knobs are separable — the only thing this suite exists to decide.

    The answer is binary and it licenses (or withdraws) the reading of two other suites, so it is a
    sentence rather than six grids to compare by eye.
    """
    if not moved:
        return
    bent = sum(1 for flag in moved if flag)
    if bent:
        status, reading = "warn", (
            f"The ridge bends in {bent} of {len(moved)} planes: the best `mlp_batch` is not the same "
            f"at every `validate_batch`. `mlp` and `validate` are each ONE SLICE through this "
            f"surface, taken at the other knob's shipped value — their curves cannot be read "
            f"independently, and the pair has to be chosen together."
        )
    else:
        status, reading = "good", (
            f"The ridge holds straight in all {len(moved)} planes: the best `mlp_batch` is the same "
            f"at every `validate_batch`. The two knobs are separable, `mlp` and `validate` are each "
            f"sufficient on their own, and this suite need not run again until the search path "
            f"changes."
        )
    report.verdict(
        [
            ("interaction", "yes" if bent else "no", "between the two batches", status),
            ("ridge bends in", f"{bent} / {len(moved)}", "planes measured", ""),
        ],
        reading,
    )
