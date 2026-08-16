"""How far along a session is, on one line that stays at the bottom of the terminal.

A full-database `run.sh all` is an overnight job. Everything it prints is per-cell — a timestamped
label as each starts — which answers "is it alive" and nothing else: not how much is left, not
whether it will be done before morning. This adds the missing half.

Three decisions worth stating:

* **Progress is measured in timed queries, not cells.** Cells differ by orders of magnitude — a
  matrix suite's process sweeps a whole grid, and a tryptic cell at a fifth the query count is a
  fifth of the cost. Counting cells would make the bar lurch; weighting by the queries each cell
  will actually run makes it roughly linear in time, which is what makes the ETA worth printing.
* **The ETA is measured, never assumed.** It extrapolates from the queries this session has already
  completed and the seconds they took, so it is empty until the first cell finishes and it adapts to
  the machine instead of to a throughput constant baked in here. Cells skipped as already-complete
  are deliberately excluded from that basis: they cost no time, and letting them into the average
  would predict the rest of the run at infinite speed.
* **Nothing is lost when the output is not a terminal.** Under `nohup`, in CI, or piped to a file,
  the bar degrades to one ordinary line per cell — the same information, without the escape codes
  that would otherwise fill the log.
"""

from __future__ import annotations

import shutil
import sys
import time
from dataclasses import dataclass, field


@dataclass
class Progress:
    """Session-wide progress across every suite that will run.

    Constructed with the preflight's per-suite estimates, so the bar is meaningful from the first
    cell rather than filling in as suites start. `begin_suite` replaces an estimate with the exact
    weight once the arms are built and the grid can be expanded for real.
    """

    #: Suite name -> timed queries expected from it, from the preflight plan.
    estimates: dict[str, int] = field(default_factory=dict)
    stream = sys.stderr

    def __post_init__(self) -> None:
        self.total = sum(self.estimates.values())
        self.done = 0
        self.suite = ""
        self.suite_remaining = 0
        self.cells_done = 0
        self.cells_total = 0
        #: Only the cells this process actually executed, which is what the ETA extrapolates from.
        self.ran_weight = 0
        self.ran_seconds = 0.0
        self.started = time.monotonic()
        self.suite_started = time.monotonic()
        self.begun: set[str] = set()
        #: Suite name -> wall clock seconds it took, filled in as each one ends. The session's own
        #: record of where its time went: the report has the same numbers, but only once the run is
        #: over, and an interrupted session never gets that far.
        self.durations: dict[str, float] = {}
        self._drawn = False
        self.enabled = self.stream.isatty()

    # -- lifecycle

    def begin_suite(self, name: str, weight: int, cells: int) -> None:
        """Starts a suite, correcting the session total with its now-exact weight.

        The preflight's estimate for a matrix suite is an upper bound until its arms are built (the
        shipped tuning defaults are read out of a binary), so the total is restated here rather than
        left to describe a grid nobody ran.
        """
        self.begun.add(name)
        self.total += weight - self.estimates.get(name, 0)
        self.estimates[name] = weight
        self.suite = name
        self.suite_remaining = weight
        self.suite_started = time.monotonic()
        self.cells_done = 0
        self.cells_total = cells
        self.draw()

    def end_suite(self) -> None:
        """Closes a suite's bar.

        Whatever it did not run leaves the session total: a suite skipped for want of root, or one
        aborted part-way, must not leave a shared bar permanently short of 100%.

        The last bar is left on screen as a line of its own rather than erased, so the scrollback
        keeps one line per suite saying what it cost — including how long that suite took, which is
        the number a session is planned by next time.
        """
        self.total -= max(self.suite_remaining, 0)
        self.suite_remaining = 0
        if self.suite:
            # `+=`: a suite reached twice in one session (an `--only` rerun, a resumed suite that a
            # later pass extends) has spent the sum of both, not the last one's time.
            self.durations[self.suite] = self.durations.get(self.suite, 0.0) + (
                time.monotonic() - self.suite_started
            )
        self.close()
        self.suite = ""

    def drop_suite(self, name: str) -> None:
        """Removes a suite that never started — blocked on privileges, or failed while building.

        A no-op once it has begun, since `end_suite` has already reconciled what it did not run.
        """
        if name not in self.begun:
            self.total -= self.estimates.pop(name, 0)

    def cell_done(self, weight: int, *, ran: bool, seconds: float = 0.0) -> None:
        self.done += weight
        self.suite_remaining = max(self.suite_remaining - weight, 0)
        self.cells_done += 1
        if ran:
            self.ran_weight += weight
            self.ran_seconds += seconds
        if self.enabled:
            self.draw()
        elif self.total > 0:
            # No terminal to redraw in, so the same line is printed once per cell instead. A log
            # that says only "cell 41 started" cannot be read for how much is left.
            print(f"  progress: {self.line()}")
            sys.stdout.flush()

    # -- output

    def echo(self, message: str = "") -> None:
        """Prints above the bar. Hand this to the runner as its `echo`."""
        self.clear()
        print(message)
        sys.stdout.flush()
        self.draw()

    def draw(self) -> None:
        # Only while a suite is running. Before the first one there is nothing to say (the builds
        # are not measured), and after the last one the bar has already been closed onto a line of
        # its own — redrawing it there would put it back in front of the report.
        if not self.suite or not self.enabled or self.total <= 0:
            return
        self.stream.write("\r\033[K" + self._bar())
        self.stream.flush()
        self._drawn = True

    def clear(self) -> None:
        if self._drawn:
            self.stream.write("\r\033[K")
            self.stream.flush()
            self._drawn = False

    def close(self) -> None:
        """Leaves the final bar on screen as a line of its own, rather than being overwritten."""
        if self._drawn:
            self.stream.write("\r\033[K" + self._bar() + "\n")
            self.stream.flush()
            self._drawn = False

    def line(self) -> str:
        """The same information as one plain line, for a non-terminal stream."""
        return self._bar(width=0)

    # -- rendering

    def _bar(self, width: int | None = None) -> str:
        fraction = min(self.done / self.total, 1.0) if self.total else 0.0
        parts = [f"{fraction * 100:5.1f}%"]
        if self.suite:
            # The suite's own clock beside the session's: under `all` the interesting question at
            # any moment is how long THIS suite has been going, and the closing line then states
            # what it cost outright.
            parts.append(
                f"{self.suite} {self.cells_done}/{self.cells_total} in "
                f"{_duration(time.monotonic() - self.suite_started)}"
            )
        parts.append(f"elapsed {_duration(time.monotonic() - self.started)}")
        parts.append(f"eta {self._eta()}")
        text = "  ·  ".join(parts)

        if width is None:
            # Whatever is left of the terminal after the text, clamped to something still readable
            # as a bar; a narrow window drops it entirely rather than wrapping the line.
            width = max(0, min(30, shutil.get_terminal_size((100, 24)).columns - len(text) - 6))
        if width < 8:
            return text
        filled = int(round(fraction * width))
        return f"[{'#' * filled}{'.' * (width - filled)}] {text}"

    def _eta(self) -> str:
        if self.ran_weight <= 0 or self.ran_seconds <= 0:
            return "?"
        remaining = max(self.total - self.done, 0)
        return "~" + _duration(remaining * self.ran_seconds / self.ran_weight)


def _duration(seconds: float) -> str:
    if seconds < 90:
        return f"{seconds:.0f}s"
    if seconds < 5400:
        return f"{seconds / 60:.1f}m"
    return f"{seconds / 3600:.1f}h"
