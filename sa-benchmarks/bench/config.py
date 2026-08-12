"""Suite definitions: what to build, what to sweep, and in what order to run it.

A suite is a TOML file under `sa-benchmarks/suites/`. It names the *arms* (feature sets to build)
and the *axes* (values to sweep), and this module expands them into the flat list of cells the
runner executes. Nothing here touches the machine — `--dry-run` stops after expansion, which is how
a multi-hour sweep gets eyeballed before it is committed to.

Axes are not free-form. Every axis name below has a defined effect on the harness invocation, and an
unrecognised name is an error rather than a silently inert extra dimension:

    ceiling_gb   cgroup v2 MemoryMax for this cell, in GB. 0 = unconstrained.
    threads      RAYON_NUM_THREADS. "default" leaves it unset.
    mlp_batch    --mlp-batch: peptides interleaved per rayon task. 1 = scalar.
    peptides     profile peptide-file key -> --peptide-file.
    kmer_table   profile k-mer-table key -> --kmer-table-file. "none" attaches no table.
    equate_il    --equate-il true|false.
    tryptic      --tryptic true|false.

Ordering decides how cells are interleaved in time, which is not cosmetic: it is the only defence
against machine drift being read as an arm effect.

    sequential   Each cell runs once, arms adjacent within an axis combination. Use when the arm
                 comparison is the delicate one and should span the smallest possible time window
                 (thread_arm_matrix.sh's reasoning: six thread counts each give an independent
                 estimate of the same arm gap, and agreement across them is the evidence).
    palindrome   Each arm runs twice per axis combination, forward then reversed (a b b a for two
                 arms, a b c c b a for three). Every arm gets one early and one late slot, so drift
                 that is linear in slot position cannot land entirely on one arm. The gap between an
                 arm's own two slots is then the floor on what the experiment can resolve — see
                 `records.slot_spread`. "abba" is accepted as an alias.

Caveat inherited from ab_storage_axes.sh (2026-08-09): a palindrome only cancels drift that is
LINEAR in slot position. If the per-invocation table shows a monotone trend, run a second pass or
discard the first invocation as warmup rather than reading the pooled medians.
"""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Any

#: Axis names with a defined effect, and a one-line description used in error messages.
#:
#: `tune.<field>` is handled separately: it sweeps any field of the searcher's `SearchTuning`, and
#: the set of those is defined in Rust, not here. That is deliberate — a knob added to that struct
#: becomes sweepable and reportable without this file changing, and a misspelled one is rejected by
#: the harness (which knows the real field list) rather than being silently accepted here.
KNOWN_AXES = {
    "ceiling_gb": "cgroup MemoryMax in GB (0 = unconstrained)",
    "threads": "RAYON_NUM_THREADS ('default' = unset)",
    "peptides": "profile peptide-file key",
    "kmer_table": "profile k-mer-table key ('none' = no table)",
    "equate_il": "true | false",
    "tryptic": "true | false",
}

#: Prefix for a SearchTuning field used as a sweep axis, e.g. `tune.mlp_batch = [1, 16]`.
TUNE_PREFIX = "tune."

ORDERINGS = ("sequential", "palindrome", "abba")


class ConfigError(Exception):
    """A suite file is malformed, or asks for something that cannot be expanded."""


@dataclass(frozen=True)
class Arm:
    """One build of the harness: a name and the cargo features that produce it."""

    name: str
    features: tuple[str, ...]

    @property
    def feature_string(self) -> str:
        """Comma-separated features, i.e. what goes after `--features`. Empty = default build."""
        return ",".join(self.features)


@dataclass(frozen=True)
class Cell:
    """One harness invocation: an arm, a point in the axis space, and its slot in the ordering."""

    suite: str
    arm: Arm
    #: Axis name -> value, exactly as written in the suite file.
    axes: dict[str, Any]
    #: Position within a palindrome ("a", "b", ...); always "a" under `sequential`.
    slot: str
    #: Ordinal position in the full run order, for reading drift down a column.
    order: int

    @property
    def label(self) -> str:
        """Stable, filesystem-safe name. Also the harness's `--label`, so `<label>.jsonl`."""
        parts = [f"{key}-{self.axes[key]}" for key in sorted(self.axes)]
        parts += [self.arm.name, self.slot]
        return re.sub(r"[^A-Za-z0-9_.-]", "_", "__".join(parts))

    @property
    def dims(self) -> dict[str, str]:
        """What gets written into every record's `dims`, as strings.

        The arm and its feature list are in here deliberately: a record must say which binary
        produced it without anyone having to remember how the sweep was invoked.
        """
        dims = {key: str(value) for key, value in self.axes.items()}
        dims["arm"] = self.arm.name
        dims["features"] = self.arm.feature_string
        dims["slot"] = self.slot
        return dims

    def describe(self) -> str:
        axes = " ".join(f"{key}={self.axes[key]}" for key in sorted(self.axes))
        return f"{self.arm.name:<12} {axes}  [slot {self.slot}]"


@dataclass
class Suite:
    """A parsed suite file."""

    name: str
    description: str
    #: "single" (one invocation per cell) or "matrix" (one invocation per arm, sweeping a grid
    #: inside the process). Matrix mode exists because index loads dominate everything else at
    #: full-database scale, so a process per grid cell would spend the whole run loading.
    mode: str
    arms: list[Arm]
    axes: dict[str, list[Any]]
    ordering: str
    needs_root: bool
    drop_caches: bool
    #: Harness settings shared by every cell (runs, amount, warmup, peptides, ...).
    defaults: dict[str, Any] = field(default_factory=dict)
    #: `mode = "matrix"` only: the grid the harness sweeps in-process.
    matrix: dict[str, Any] = field(default_factory=dict)
    #: Build the arms with sa-index's `metrics` feature. Perturbs throughput; see `detail`.
    metrics: bool = False
    #: Prose printed under this suite's tables, explaining how to read them.
    notes: str = ""
    #: Results directory of a previous run of this suite, for the regression comparison. Set from
    #: `--baseline` rather than from the suite file: it names a past run, not a property of the suite.
    baseline: Path | None = None

    def expand(self) -> list[Cell]:
        """Expands arms x axes into the ordered list of cells this suite runs."""
        axis_names = sorted(self.axes)
        combinations = [
            dict(zip(axis_names, values))
            for values in product(*(self.axes[name] for name in axis_names))
        ] or [{}]

        cells: list[Cell] = []
        for combination in combinations:
            for slot, arm in self._slots():
                cells.append(
                    Cell(
                        suite=self.name,
                        arm=arm,
                        axes=combination,
                        slot=slot,
                        order=len(cells),
                    )
                )
        return cells

    def _slots(self) -> list[tuple[str, Arm]]:
        """The arm order within one axis combination, as (slot, arm) pairs."""
        if self.ordering == "sequential":
            return [("a", arm) for arm in self.arms]
        # palindrome / abba: forward then reversed, so every arm holds one early and one late slot.
        forward = list(enumerate(self.arms))
        ordered = forward + list(reversed(forward))
        return [(chr(ord("a") + position), arm) for position, (_, arm) in enumerate(ordered)]


def suites_dir(repo: Path) -> Path:
    return repo / "sa-benchmarks" / "suites"


def available(repo: Path) -> list[str]:
    return sorted(path.stem for path in suites_dir(repo).glob("*.toml"))


def load(name: str, repo: Path) -> Suite:
    """Loads and validates `suites/<name>.toml`."""
    path = suites_dir(repo) / f"{name}.toml"
    if not path.exists():
        raise ConfigError(f"no suite '{name}' at {path}\n  available: {', '.join(available(repo))}")

    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    mode = raw.get("mode", "single")
    if mode not in ("single", "matrix"):
        raise ConfigError(f"{path}: mode must be 'single' or 'matrix', got '{mode}'")

    arms = [_arm(entry, path) for entry in raw.get("arms", [])]
    if not arms:
        raise ConfigError(f"{path}: no [[arms]] — a suite must build at least one configuration")
    if len({arm.name for arm in arms}) != len(arms):
        raise ConfigError(f"{path}: two arms share a name; arm names become file names")

    axes = _flatten_tune(raw.get("axes", {}))
    for axis, values in axes.items():
        if axis.startswith(TUNE_PREFIX):
            if len(axis) == len(TUNE_PREFIX):
                raise ConfigError(f"{path}: axis '{axis}' names no tuning field")
        elif axis not in KNOWN_AXES:
            known = "\n    ".join(f"{key:<12} {why}" for key, why in KNOWN_AXES.items())
            raise ConfigError(
                f"{path}: unknown axis '{axis}'. An axis with no defined effect would sweep "
                f"nothing while still multiplying the run time.\n  known axes:\n    {known}\n"
                f"    {TUNE_PREFIX}<field>  any SearchTuning field "
                f"(see `sa-benchmarks --help-tuning`)"
            )
        if not isinstance(values, list) or not values:
            raise ConfigError(f"{path}: axis '{axis}' must be a non-empty list")

    ordering = raw.get("ordering", "sequential")
    if ordering not in ORDERINGS:
        raise ConfigError(f"{path}: ordering must be one of {ORDERINGS}, got '{ordering}'")
    if ordering == "abba":
        ordering = "palindrome"

    suite = Suite(
        name=name,
        description=raw.get("description", "").strip(),
        mode=mode,
        arms=arms,
        axes=axes,
        ordering=ordering,
        needs_root=bool(raw.get("needs_root", False)),
        drop_caches=bool(raw.get("drop_caches", False)),
        defaults=raw.get("defaults", {}),
        matrix=raw.get("matrix", {}),
        metrics=bool(raw.get("metrics", False)),
        notes=raw.get("notes", "").strip(),
    )

    if suite.drop_caches and not suite.needs_root:
        raise ConfigError(f"{path}: drop_caches needs root, so needs_root must also be true")
    if any(value != 0 for value in axes.get("ceiling_gb", [])) and not suite.needs_root:
        raise ConfigError(f"{path}: a non-zero ceiling_gb needs root (cgroup v2), so set needs_root")
    if mode == "matrix" and axes:
        raise ConfigError(
            f"{path}: mode='matrix' sweeps its grid inside one process, so it cannot also have "
            f"[axes] — put the grid under [matrix]"
        )
    return suite


def _flatten_tune(axes: dict) -> dict:
    """Turns the `[axes.tune]` table into flat `tune.<field>` axis names.

    TOML reads `tune.mlp_batch = [...]` as a nested table, not as a key with a dot in it, so a
    suite writes the natural thing:

        [axes.tune]
        mlp_batch = [1, 16]

    and everything downstream sees one flat axis named `tune.mlp_batch`. Flattening here rather
    than teaching the expander about nesting keeps cells, dims and labels one level deep.
    """
    flat = {key: value for key, value in axes.items() if key != "tune"}
    for field, values in (axes.get("tune") or {}).items():
        flat[f"{TUNE_PREFIX}{field}"] = values
    return flat


def _arm(entry: dict, path: Path) -> Arm:
    try:
        name = entry["name"]
    except KeyError:
        raise ConfigError(f"{path}: an [[arms]] entry has no name") from None
    features = entry.get("features", [])
    if not isinstance(features, list):
        raise ConfigError(f"{path}: arm '{name}': features must be a list")
    return Arm(name=name, features=tuple(features))
