"""Suite definitions: what to build, what to sweep, and in what order to run it.

A suite is a TOML file under `sa-benchmarks/suites/`. It names the *arms* (feature sets to build)
and the *axes* (values to sweep), and this module expands them into the flat list of cells the
runner executes. Nothing here touches the machine — `--dry-run` stops after expansion, which is how
a multi-hour sweep gets eyeballed before it is committed to.

Axes are not free-form. Every axis name below has a defined effect on the harness invocation, and an
unrecognised name is an error rather than a silently inert extra dimension:

    ceiling_gb   cgroup v2 MemoryMax for this cell, in GB. 0 = unconstrained.
    threads      RAYON_NUM_THREADS. "default" leaves it unset.
    peptides     profile peptide-file key -> --peptide-file.
    kmer         k of the k-mer table to attach. 0 attaches none.
    equate_il    --equate-il true|false.
    tryptic      --tryptic true|false.

`kmer` is spelled the same here as in a `[[sweep]]` block, and means the same thing: the integer k.
It used to be `kmer_table = "k6"` on this side and `kmer = [6]` on the other — the same fact in two
vocabularies, which meant anything asking about both (`preflight._table_notes`) carried two parsers
for one question. The profile key is derived from k where the file is actually needed, which is what
the `k<N>` naming convention in `[kmer_tables]` is for. `kmer_table` is still accepted and
normalised to `kmer`, with a deprecation note in the preflight; see `_normalise_kmer`.

A matrix-mode suite is different: it loads the index once and sweeps in-process, so an axis there
would cost a whole index load per value. Only `threads` and `ceiling_gb` may be axes (see
`PROCESS_AXES` — neither can change while a process lives), and everything else it varies goes in
`[[sweep]]` blocks, which `bench.grid` expands into a grid file the harness reads.

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
#: There used to be a `tune.<field>` form that swept any field of the searcher's `SearchTuning`.
#: That struct is gone — its fields are compile-time constants now, because no sweep could separate
#: their values from noise — so the searcher exposes nothing to vary at runtime and every axis left
#: here is either the workload, the build, or the machine.
KNOWN_AXES = {
    "ceiling_gb": "cgroup MemoryMax in GB (0 = unconstrained)",
    "threads": "RAYON_NUM_THREADS ('default' = unset)",
    "peptides": "profile peptide-file key",
    "kmer": "k of the k-mer table to attach (0 = none)",
    "equate_il": "true | false",
    "tryptic": "true | false",
}

#: Axes that decide which *process* a cell runs in, rather than what it measures inside one.
#:
#: These are the only axes a matrix-mode suite may have. Rayon's global pool is built once per
#: process before any searcher exists and a cgroup scope wraps a process, so neither can change
#: while an index stays loaded. Everything else a matrix suite varies belongs in a `[[sweep]]`
#: block, where it costs a cell rather than a whole index load.
PROCESS_AXES = ("threads", "ceiling_gb")

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
        """Comma-separated features, i.e. what goes after `--features`. Empty = default build.

        This string is what every record reports as the binary that produced it, so two arms that
        differ only in their features stay distinguishable in the records.
        """
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
    #: `mode = "matrix"` only: `[[sweep]]` blocks, which `bench.grid` expands into the grid file the
    #: harness sweeps in-process. This is the only way a matrix suite says what it measures — the
    #: harness has no built-in grid of its own, so a sweep cannot be described in two places.
    sweeps: list[dict[str, Any]] = field(default_factory=list)
    #: Prose printed under this suite's tables, explaining how to read them.
    notes: str = ""
    #: Results directory of a previous run of this suite, for the regression comparison. Set from
    #: `--baseline` rather than from the suite file: it names a past run, not a property of the suite.
    baseline: Path | None = None
    #: Spellings this file used that still work but should be changed. Carried rather than printed,
    #: because `load` has no output of its own and a warning printed from here would arrive before
    #: the preflight block it belongs in — see `preflight._plan_suite`.
    deprecations: list[str] = field(default_factory=list)

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

    axes = raw.get("axes", {})
    defaults = raw.get("defaults", {})
    deprecations = _normalise_kmer(axes, path, axis=True) + _normalise_kmer(defaults, path)

    for axis, values in axes.items():
        if axis not in KNOWN_AXES:
            known = "\n    ".join(f"{key:<12} {why}" for key, why in KNOWN_AXES.items())
            hint = (
                "\n  `tune.<field>` was removed with `SearchTuning`: the searcher's performance "
                "parameters are compile-time constants now."
                if axis.startswith("tune")
                else ""
            )
            raise ConfigError(
                f"{path}: unknown axis '{axis}'. An axis with no defined effect would sweep "
                f"nothing while still multiplying the run time.\n  known axes:\n    {known}{hint}"
            )
        if not isinstance(values, list) or not values:
            raise ConfigError(f"{path}: axis '{axis}' must be a non-empty list")

    ordering = raw.get("ordering", "sequential")
    if ordering not in ORDERINGS:
        raise ConfigError(f"{path}: ordering must be one of {ORDERINGS}, got '{ordering}'")
    if ordering == "abba":
        ordering = "palindrome"

    drop_caches = bool(raw.get("drop_caches", False))
    capped = any(value != 0 for value in axes.get("ceiling_gb", []))

    suite = Suite(
        name=name,
        description=raw.get("description", "").strip(),
        mode=mode,
        arms=arms,
        axes=axes,
        ordering=ordering,
        # Derived, not declared. A cgroup ceiling and `drop_caches` both need root by construction,
        # and a suite whose ceiling list is the thing being edited should not also have to keep a
        # boolean in sync with it — a mismatch there fails the run rather than the edit.
        needs_root=bool(raw.get("needs_root", False)) or drop_caches or capped,
        drop_caches=drop_caches,
        defaults=defaults,
        sweeps=raw.get("sweep", []),
        notes=raw.get("notes", "").strip(),
        deprecations=deprecations,
    )

    if mode == "matrix":
        _check_matrix(suite, path)
    elif suite.sweeps:
        raise ConfigError(
            f"{path}: [[sweep]] blocks are expanded into a grid the harness sweeps in one process, "
            f"which only mode='matrix' does. In single mode use [axes] instead."
        )
    return suite


def _normalise_kmer(table: dict, path: Path, axis: bool = False) -> list[str]:
    """Rewrites a deprecated `kmer_table = "k6"` into `kmer = 6`, in place.

    `axis` says whether this table may hold several values. `[axes]` may; `[defaults]` names one
    value per setting, and a list there reaches consumers that treat it as a scalar.

    One vocabulary downstream. `[[sweep]]` blocks have always named the integer k, because a matrix
    process keeps a pool of tables keyed by k and swaps them per cell; `[axes]` and `[defaults]`
    named the profile key instead, so the same fact was spelled two ways and every reader of both
    needed two parsers. The profile key is now derived from k at the one place the file is opened.

    Accepted rather than rejected, because a suite file is also a record of past runs and breaking
    every archived one to rename a key is a poor trade. The `k<N>` convention is what makes the
    rewrite possible at all — the same convention that lets a missing table be rebuilt from its
    name — so a table named anything else is an error here rather than a silent misreading.
    """
    if "kmer_table" not in table:
        return []
    value = table.pop("kmer_table")

    def k_of(name) -> int:
        if name in (None, "", "none"):
            return 0
        name = str(name)
        if not (name.startswith("k") and name[1:].isdigit()):
            raise ConfigError(
                f"{path}: kmer_table = '{name}' cannot be read as a k. Names are 'k<N>' (k5, k6) or "
                f"'none'. This spelling is deprecated in favour of `kmer = <N>` — write that instead."
            )
        return int(name[1:])

    if isinstance(value, list) and not axis:
        # A list is a sweep, and only an axis can sweep. In `[defaults]` it would be accepted here
        # and then fail far away: `runner._kmer_args` does `int(settings.get("kmer"))` and
        # `preflight._table_notes` indexes it, and the `TypeError` from either escapes
        # `_plan_suite`, which catches only `ConfigError` — so one mis-typed default takes down a
        # session that was meant to report it and carry on.
        raise ConfigError(
            f"{path}: kmer_table = {value!r} lists several tables in [defaults], which names one "
            f"value per setting. List them in [axes] to sweep them, or pick one here."
        )

    table["kmer"] = [k_of(item) for item in value] if isinstance(value, list) else k_of(value)
    shown = table["kmer"]
    # `kmer = [5, 6]` is only valid advice for an axis; in [defaults] it is the shape just rejected.
    return [
        f"`kmer_table` is deprecated: write `kmer = {shown}` instead (same meaning, and the same "
        f"spelling a [[sweep]] block uses)"
    ]


def _check_matrix(suite: Suite, path: Path) -> None:
    """Matrix mode may vary processes through `[axes]` and cells through `[[sweep]]`, nothing else."""
    if not suite.sweeps:
        raise ConfigError(
            f"{path}: mode='matrix' needs at least one [[sweep]] block saying what to sweep. "
            f"(The harness used to carry a built-in grid selected by a [matrix] table; it no longer "
            f"does, because a grid described half in TOML and half in Rust could only be widened by "
            f"editing both.)"
        )

    stray = sorted(set(suite.axes) - set(PROCESS_AXES))
    if stray:
        raise ConfigError(
            f"{path}: axis '{stray[0]}' cannot be a matrix-mode axis. Matrix mode loads the index "
            f"once and sweeps in-process, so an axis here would cost a whole index load per value.\n"
            f"  process axes (one process each): {', '.join(PROCESS_AXES)}\n"
            f"  everything else goes in a [[sweep]] block, where it costs one cell"
        )


def _arm(entry: dict, path: Path) -> Arm:
    try:
        name = entry["name"]
    except KeyError:
        raise ConfigError(f"{path}: an [[arms]] entry has no name") from None
    features = entry.get("features", [])
    if not isinstance(features, list):
        raise ConfigError(f"{path}: arm '{name}': features must be a list")
    return Arm(name=name, features=tuple(features))
