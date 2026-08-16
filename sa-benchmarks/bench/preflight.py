"""What this machine is about to be asked to do, and whether it can honestly do it.

Every check in `rig` exists because its absence has invalidated a benchmark before. This module is
where they are all asked *at once, before anything runs* — and where the run says out loud how big
it is. Two failures motivate it:

* **A sweep that dies four hours in.** The peptide-supply check, the cgroup probe and the profile
  validation each used to fire when the cell that needed them was reached. Under `all`, the suites
  that need root run last, so a session could spend an evening on `defaults` and then discover it
  could never have run `ram`. Everything askable up front is asked up front.
* **A sweep whose size was a surprise.** "16 processes" and "38 million queries" are the difference
  between a coffee break and an overnight run, and the moment to learn which one it is is before
  starting, not from the progress bar an hour later.

The output is a block of checks and a plan table. A `FAIL` aborts before the first build; a `warn`
is printed and the run proceeds, because most of what invalidates a comparison (a dirty tree, a
co-tenant job) is the operator's call and not something to refuse over.

The counts here are the same expansion the runner executes — `select` is what both call — so the
plan cannot describe a different run from the one that happens. The one caveat is a matrix suite's
grid before its arms are built: the shipped tuning defaults are read out of a binary, so the cell
count is an upper bound by at most one cell per knob until then (see `runner.tuning_defaults`). The
progress bar re-totals from the real grid once the arms exist.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from . import config, rig
from .config import Cell, Suite
from .profile import INDEX_FILES, Profile
from .runner import Runner

OK, WARN, FAIL = "ok", "warn", "FAIL"

#: Free scratch space below which a long session is likely to end in ENOSPC rather than a report.
#: Records are small; logs under `--matrix` are not, and a full-database session writes both for
#: every cell of every suite.
SCRATCH_FLOOR_GB = 5.0


@dataclass
class Check:
    """One answerable question about the machine, with the answer and how bad it is."""

    name: str
    status: str
    detail: str

    def line(self) -> str:
        return f"  {self.status:<4}  {self.name:<12} {self.detail}"


@dataclass
class SuitePlan:
    """How big one suite is, and whether it can run here."""

    name: str
    #: "run"  — will execute
    #: "skip" — optional and blocked on this machine; the report will say so in its section
    #: "warn" — this suite will fail, but only this one: a malformed suite file, or a `--baseline`
    #:          that holds nothing for it (normal for a suite added since that session). The rest of
    #:          the session is unaffected, so it is not worth refusing to start over.
    #: "FAIL" — would waste the session: a short peptide file, a bucket the profile lacks, a
    #:          required suite this machine cannot run.
    status: str
    processes: int = 0
    grid_cells: int = 0
    queries: int = 0
    #: Cells that already have results (or a did-not-fit marker) in the session directory. A resumed
    #: session skips these, so they are neither work nor missing data.
    done: int = 0
    reason: str = ""
    notes: list[str] = field(default_factory=list)
    #: What this suite needs from the machine, so the privilege checks can be asked only when
    #: something in the session actually depends on them.
    needs_root: bool = False
    needs_cgroup: bool = False
    #: Timed queries still to run — what the progress bar apportions its bar over.
    weight: int = 0


@dataclass
class Preflight:
    machine: list[Check]
    suites: list[SuitePlan]

    @property
    def ok(self) -> bool:
        return not any(check.status == FAIL for check in self.machine) and not any(
            plan.status == FAIL for plan in self.suites
        )

    @property
    def failures(self) -> list[str]:
        return [f"{check.name}: {check.detail}" for check in self.machine if check.status == FAIL] + [
            f"{plan.name}: {plan.reason}" for plan in self.suites if plan.status == FAIL
        ]

    @property
    def weights(self) -> dict[str, int]:
        """Suite name -> timed queries left to run, for the progress bar's initial total."""
        return {plan.name: plan.weight for plan in self.suites if plan.status == "run"}


# ---------------------------------------------------------------------------
# The cell list, shared with the runner
# ---------------------------------------------------------------------------


def results_dir(session: Path, name: str) -> Path:
    """Where a suite's records go.

    A standalone run gets `<session>/results`; under `all` the session is shared and each suite
    writes into its own subdirectory. Defined once, because preflight has to count the cells that
    are already there and would otherwise be counting a different directory from the one the runner
    writes to.
    """
    return session / "results" if session.name == name else session / name


def select(suite: Suite, runner: Runner, only: str | None) -> list[Cell]:
    """The cells this invocation will run: the suite's expansion, filtered by `--only`.

    Both the plan and the run go through here, so the number in the preflight table is the number
    of cells that execute rather than a second, drifting estimate of it.
    """
    cells = suite.expand()
    if only:
        wanted = [token.strip() for token in only.split(",") if token.strip()]
        cells = [cell for cell in cells if any(token in cell.label for token in wanted)]
        if not cells:
            raise config.ConfigError(f"--only '{only}' matched no cell in suite '{suite.name}'")
    return runner.applicable(cells)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def check(
    names: list[str],
    args,
    repo: Path,
    profile: Profile,
    session: Path,
    optional: set[str] | None = None,
) -> Preflight:
    """Everything askable before the first build, for one suite or for a whole session."""
    optional = optional or set()
    suites: list[SuitePlan] = []
    for name in names:
        suites.append(_plan_suite(name, args, repo, profile, session, name in optional))
    return Preflight(machine=_machine(repo, profile, suites), suites=suites)


def _machine(repo: Path, profile: Profile, plans: list[SuitePlan]) -> list[Check]:
    """The checks that are about the box rather than about a suite — asked once per session."""
    facts = rig.host_facts()
    checks: list[Check] = []

    state = rig.git_state(repo)
    checks.append(Check("host", OK, f"{facts['cpu']} · {facts['cores']} cores · {facts['kernel']} {facts['arch']}"))
    checks.append(Check("profile", OK, f"{profile.name} (python {facts['python']})"))
    checks.append(Check("commit", OK, f"{state.short} ({state.branch})"))
    checks.append(
        Check("tree", WARN, "DIRTY — these numbers are not attributable to that commit")
        if state.dirty
        else Check("tree", OK, "clean")
    )

    cores = int(facts["cores"] or 1)
    load = rig.load_average()[0]
    checks.append(
        Check(
            "load",
            WARN if load > cores * 0.25 else OK,
            f"{facts['load']} (1 / 5 / 15 min) on {cores} cores"
            + (" — a co-tenant job invalidates every comparison in this run" if load > cores * 0.25 else ""),
        )
    )

    checks.append(_toolchain_check(facts))

    checks += _privilege_checks(plans, facts)
    checks += _storage_checks(profile, facts)
    return checks


def _toolchain_check(facts: dict[str, str]) -> Check:
    """Can the arms be built — asked along the route `bench.build` actually takes.

    Under `sudo`, PATH is replaced by `secure_path` and a rustup toolchain lives in the invoking
    user's `~/.cargo/bin`, so a probe of root's PATH answers "no toolchain" for a box that builds
    fine. Builds are handed back to that user under a login shell, and so is this.

    Never a `FAIL`, deliberately, however it comes out. The two mistakes are not symmetric: a false
    negative here would refuse to start a session that would have worked, while a false positive
    costs the seconds it takes cargo to say so itself — and the builds run before the first cell, so
    a real missing toolchain still stops the run before anything is measured.
    """
    missing = [tool for tool in ("cargo", "rustc") if facts[tool] == "?"]
    if not missing:
        return Check("toolchain", OK, f"{facts['rustc']} · {facts['cargo']}")

    named = " and ".join(missing)
    # Runnable here but not along the build's route: the toolchain exists, the build just will not
    # find it. Worth separating, because the fix is the invoking user's shell profile, not rustup.
    if rig.dropping_privileges() and all(rig.tool_version_here(tool) != "?" for tool in missing):
        user = os.environ.get("SUDO_USER", "the invoking user")
        return Check(
            "toolchain",
            WARN,
            f"{named} runs here but not under `bash -lc` as {user}, which is how builds run — "
            f"add ~/.cargo/bin to that account's login profile if the build fails",
        )
    if rig.is_root() and not rig.dropping_privileges():
        return Check(
            "toolchain",
            WARN,
            f"{named} not runnable as root (sudo's secure_path drops ~/.cargo/bin) — prefer "
            f"'sudo ./sa-benchmarks/run.sh ...' from your own account, which hands builds back to you",
        )
    return Check("toolchain", WARN, f"{named} not runnable — the first build will fail if this is real")


def _privilege_checks(plans: list[SuitePlan], facts: dict[str, str]) -> list[Check]:
    """Root and cgroups, asked only when something in this session actually needs them.

    Never fatal on its own: what a missing privilege costs is already the affected suite's status —
    `skip` when it is optional, `FAIL` when it is not — so this reports the machine's capability and
    lets the plan table say who needed it.
    """
    checks: list[Check] = []
    wanted = [plan for plan in plans if plan.needs_root]
    capped = any(plan.needs_cgroup for plan in plans)

    if rig.is_root():
        detail = "running as root" + ("" if rig.dropping_privileges() else " with no SUDO_USER")
        status = OK if rig.dropping_privileges() else WARN
        if status == WARN:
            detail += " — builds and git run as root and leave root-owned files in target/"
        checks.append(Check("privileges", status, detail))
    else:
        checks.append(
            Check(
                "privileges",
                WARN if wanted else OK,
                "not root" + (f" — {len(wanted)} suite(s) need it: {', '.join(p.name for p in wanted)}" if wanted else ""),
            )
        )

    # The cgroup probe is only meaningful under root (systemd-run --scope needs it to bind), and
    # `setenv_reaches_child` actually starts a scope, so it is not asked of a session that has no
    # capped suite in it.
    if capped:
        detail = f"{facts['cgroup']} · systemd-run {facts['systemd_run']}"
        if rig.is_root() and rig.has_systemd_run() and not rig.setenv_reaches_child():
            checks.append(Check("cgroups", WARN, detail + " · --setenv does NOT reach the child"))
        else:
            checks.append(Check("cgroups", OK, detail))
    return checks


def _storage_checks(profile: Profile, facts: dict[str, str]) -> list[Check]:
    """The dataset and the room to write about it."""
    checks: list[Check] = []

    index_bytes = sum((profile.index_dir / name).stat().st_size for name in INDEX_FILES)
    ram_gb = float(facts["ram_gb"]) if facts["ram_gb"].replace(".", "").isdigit() else 0.0
    index_gb = index_bytes / 2**30
    detail = f"{index_gb:,.2f} GB in {profile.index_dir}"
    if ram_gb and index_gb > ram_gb * 0.9:
        # Not fatal — that comparison is half of what `ram` measures — but a preloaded arm that does
        # not fit is measuring the swap device or is about to be OOM-killed, so it is said up front.
        checks.append(
            Check("index", WARN, detail + f" — larger than this box's {ram_gb:,.0f} GB of RAM; a preloaded arm will not fit")
        )
    else:
        checks.append(Check("index", OK, detail))

    swap = facts["swap_gb"]
    checks.append(Check("memory", OK, f"{facts['ram_gb']} GB RAM, swap {swap if not swap.isdigit() else swap + ' GB'}"))

    counts = ", ".join(f"{name} {_lines(path):,}" for name, path in sorted(profile.peptides.items()))
    checks.append(Check("peptides", OK, f"{len(profile.peptides)} file(s): {counts} lines"))
    checks.append(
        Check("kmer tables", OK, ", ".join(sorted(profile.kmer_tables)) or "none in this profile — a suite that asks for one builds it in-process")
    )

    free_gb = shutil.disk_usage(profile.scratch).free / 2**30
    checks.append(
        Check(
            "scratch",
            WARN if free_gb < SCRATCH_FLOOR_GB else OK,
            f"{free_gb:,.1f} GB free at {profile.scratch}"
            + (f" — below {SCRATCH_FLOOR_GB:g} GB; a long session writes a log per cell" if free_gb < SCRATCH_FLOOR_GB else ""),
        )
    )
    return checks


def _plan_suite(name: str, args, repo: Path, profile: Profile, session: Path, optional: bool) -> SuitePlan:
    """One suite's size, and the reasons it might not run.

    Everything raised here is turned into a status rather than propagated: a session must be able to
    report that its third suite is misconfigured while still saying how big the other four are.
    """
    from .__main__ import apply_overrides, needs_cgroup

    try:
        suite = config.load(name, repo)
        apply_overrides(suite, args)
    except (config.ConfigError, rig.RigError) as error:
        # `warn`, not `FAIL`: a suite file that will not load, or a `--baseline` session that holds
        # nothing for this suite, costs this suite and no other. Under `all` that is exactly what
        # happened before there was a preflight — the suite is reported as failed and the session
        # carries on — and refusing to start eleven suites over one of them would be a worse trade.
        return SuitePlan(name, WARN, reason=str(error))

    out_dir = results_dir(session, name)
    runner = Runner(suite, profile, {}, out_dir, echo=lambda _: None)
    try:
        cells = select(suite, runner, args.only)
    except config.ConfigError as error:  # --only matched nothing here; it may still match elsewhere
        return SuitePlan(name, WARN, reason=str(error))

    plan = SuitePlan(name, "run", needs_root=suite.needs_root, needs_cgroup=needs_cgroup(suite))
    plan.processes = len(cells)
    for cell in cells:
        plan.grid_cells += len(runner.grid_for(cell)) if suite.sweeps else 1
        plan.queries += runner.weight(cell)
        if runner.completed(cell):
            plan.done += 1
        else:
            plan.weight += runner.weight(cell)

    if not cells:
        plan.status = WARN
        plan.reason = "expands to no cells — nothing in it would be measured"
        return plan

    # A suite the machine cannot run needs no notes about what it would have consumed; the reason it
    # is not running is the only thing to say about it.
    blocking = rig.blockers(suite.needs_root, needs_cgroup(suite))
    if blocking:
        plan.status = "skip" if optional else FAIL
        plan.reason = "; ".join(blocking)
        plan.weight = 0
        return plan

    plan.notes += _supply_notes(runner, cells, profile, plan)
    plan.notes += _table_notes(suite, profile)

    if plan.status != "run":
        return plan
    if plan.done == plan.processes:
        plan.reason = "already complete in this session — every cell would be skipped"
    elif plan.done:
        plan.reason = f"resuming: {plan.done} of {plan.processes} cells already have results"
    return plan


def _supply_notes(runner: Runner, cells: list[Cell], profile: Profile, plan: SuitePlan) -> list[str]:
    """The peptide-supply check, hours before the runner would reach it.

    A file shorter than what the plan consumes does not fail the harness — it shortens the later
    reps, so the run completes and measures something else. That is the one check here that turns a
    plan from `run` into `FAIL`.
    """
    supply: list[str] = []
    for name, needed in sorted(runner.requirements(cells).items()):
        try:
            path = profile.peptide_file(name)
        except Exception as error:  # ProfileError: the suite names a bucket this profile lacks
            plan.status = FAIL
            plan.reason = str(error)
            continue
        have = _lines(path)
        if have < needed:
            plan.status = FAIL
            plan.reason = (
                f"peptide file '{name}' has {have:,} lines, needs >= {needed:,} "
                f"(runs x amount) — generate more with sa-benchmarks/src/generate_peptides.rs, "
                f"or lower runs/amount"
            )
        else:
            supply.append(f"{name} {needed:,}/{have:,}")
    return [f"queries (needs/has): {' · '.join(supply)}"] if supply else []


def _table_notes(suite: Suite, profile: Profile) -> list[str]:
    """k-mer tables a suite asks for that this profile has to build instead of load.

    Not a failure — the same table is built in-process — but it is minutes of startup per process,
    which is worth knowing before rather than after. A profile with no tables at all says so once,
    in the machine checks; repeating it under every suite would bury the case that matters, which is
    the profile that has the 5-mer and not the 6-mer.
    """
    if not profile.kmer_tables:
        return []
    wanted: set[int] = set()
    for block in suite.sweeps:
        wanted |= {int(k) for k in block.get("kmer", []) if int(k)}
    table = suite.defaults.get("kmer_table")
    if table and str(table) != "none" and str(table).startswith("k") and str(table)[1:].isdigit():
        wanted.add(int(str(table)[1:]))

    built = sorted(k for k in wanted if not profile.kmer_table(f"k{k}"))
    if not built:
        return []
    return [
        "kmer: no pre-built table for "
        + ", ".join(f"k={k}" for k in built)
        + " in this profile — built in-process, paid at every process startup"
    ]


def _lines(path: Path) -> int:
    """Lines in a peptide file. `rig.count_lines` caches, so asking per suite costs one `wc -l`."""
    return rig.count_lines(path)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render(preflight: Preflight, *, header: str = "preflight") -> list[str]:
    lines = [f"== {header} ==", ""]
    lines += [check.line() for check in preflight.machine]

    lines += ["", "  suite         processes   grid cells        queries   status"]
    totals = [0, 0, 0]
    for plan in preflight.suites:
        note = "" if plan.status == "run" and not plan.reason else f"   {plan.reason}"
        lines.append(
            f"  {plan.name:<12}  {plan.processes:>9}   {plan.grid_cells:>10}   {plan.queries:>12,}   "
            f"{plan.status}{note}"
        )
        for note in plan.notes:
            lines.append(f"       {note}")
        if plan.status == "run":
            totals = [a + b for a, b in zip(totals, (plan.processes, plan.grid_cells, plan.queries))]

    if len(preflight.suites) > 1:
        lines.append(
            f"  {'total':<12}  {totals[0]:>9}   {totals[1]:>10}   {totals[2]:>12,}   "
            f"{sum(1 for plan in preflight.suites if plan.status == 'run')} suite(s) to run"
        )

    lines.append("")
    if not preflight.ok:
        lines.append("  NOT OK — nothing was run:")
        lines += [f"    - {reason}" for reason in preflight.failures]
        return lines

    warned = [plan for plan in preflight.suites if plan.status == WARN]
    if warned:
        # Said twice deliberately: the table row is easy to miss, and a suite that will produce
        # nothing is the kind of thing that is only noticed once the report is missing a section.
        lines.append(
            "  this suite cannot produce a report:"
            if len(preflight.suites) == 1
            else "  ok to run, but these suites will fail — the rest of the session is unaffected:"
        )
        lines += [f"    - {plan.name}: {plan.reason}" for plan in warned]
    elif any(check.status == WARN for check in preflight.machine):
        lines.append("  ok to run, with the warnings above")
    else:
        lines.append("  ok to run")
    return lines
