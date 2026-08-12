"""`run.sh all`: every suite in one session, one report characterising this commit.

This is an orchestrator, not a sixth suite. It reuses each suite's own `analyse()` — the master
report must never re-implement a table or reword a narrative, or the two would drift apart exactly
the way the ten driver scripts drifted from each other.

Three things it owns that no single suite can:

* **One session.** All suites share `<scratch>/<commit>-<timestamp>/`, one built binary per arm, and
  one preflight. Resumable at cell granularity like anything else, so an interrupted overnight run
  restarts where it stopped.
* **Explicit degradation.** Without root there is no cgroup ceiling and no `drop_caches`, so `ram`
  and `threads` cannot run. They are skipped, and their sections say why. A suite missing from a
  report reads as "nothing to say"; a suite that says "not run — needs root" reads as what it is.
* **The verdict.** With `--baseline`, what moved since the last session, judged against the floor
  each comparison has to clear rather than against zero.
"""

from __future__ import annotations

import json
import time
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

from . import config, records, rig
from .html import render_page
from .profile import Profile, ProfileError, load as load_profile
from .report import Report, Table


@dataclass
class SuiteOutcome:
    name: str
    status: str  # "ok" | "skipped" | "failed"
    report: Report | None = None
    reason: str = ""
    seconds: float = 0.0
    #: One line naming what this suite was run with, for provenance.
    settings: str = ""


def run_all(args, repo: Path) -> int:
    """Runs every suite listed in `suites/all.toml` and writes report.md + report.json."""
    from .__main__ import run_suite

    try:
        plan = _load_plan(repo)
        profile = load_profile(args.profile, repo)
    except (ProfileError, config.ConfigError, rig.RigError) as error:
        print(f"error: {error}")
        return 1

    state = rig.git_state(repo)
    session = args.out or (profile.scratch / f"{state.short}-{time.strftime('%Y%m%d-%H%M%S')}")
    session.mkdir(parents=True, exist_ok=True)

    skip = {name.strip() for name in (args.skip or "").split(",") if name.strip()}
    names = [name for name in plan["order"] if name not in skip]
    optional = set(plan.get("optional", []))

    print(f"== session {session} ==")
    print(f"== suites: {', '.join(names)} ==\n")

    outcomes: list[SuiteOutcome] = []
    for name in names:
        outcome = _run_one(run_suite, name, args, repo, session, optional)
        outcomes.append(outcome)
        print(f"-- {name}: {outcome.status}{(' — ' + outcome.reason) if outcome.reason else ''}\n")

    if args.dry_run:
        return 0

    report = _assemble(plan, profile, repo, state, session, outcomes)
    # Three renderings of one report object: the page to read, the markdown to paste into a PR, and
    # the JSON a later run consumes as its baseline.
    (session / "report.html").write_text(
        render_page(
            f"Suffix-array index — {state.short}",
            report,
            subtitle=f"{state.branch} · {profile.name} · {time.strftime('%Y-%m-%d %H:%M')}",
            statuses={outcome.name: outcome.status for outcome in outcomes},
        )
    )
    (session / "report.md").write_text(report.to_markdown())
    (session / "report.json").write_text(json.dumps(_machine_readable(state, outcomes), indent=2) + "\n")

    print(report.to_text())
    print(f"\nreport: {session / 'report.html'}")
    print(f"        {session / 'report.md'}")
    return 0


def _run_one(run_suite, name: str, args, repo: Path, session: Path, optional: set[str]) -> SuiteOutcome:
    """Runs one suite, converting a can't-run into a recorded skip when the suite is optional."""
    print(f"===== {name} =====")
    started = time.monotonic()
    try:
        _, report = run_suite(name, args, repo, session=session)
    except rig.RigError as error:
        if name in optional:
            return SuiteOutcome(name, "skipped", reason=str(error), seconds=time.monotonic() - started)
        return SuiteOutcome(name, "failed", reason=str(error), seconds=time.monotonic() - started)
    except (ProfileError, config.ConfigError) as error:
        # A misconfigured suite must not take the whole session down: the other four still answer
        # their questions, and the report will say this one did not.
        return SuiteOutcome(name, "failed", reason=str(error), seconds=time.monotonic() - started)

    elapsed = time.monotonic() - started
    if report is None:
        # `--dry-run` plans but produces nothing; calling that "ok" would read as a completed suite.
        return SuiteOutcome(name, "planned", seconds=elapsed)
    outcome = SuiteOutcome(name, "ok", report=report, seconds=elapsed)
    outcome.settings = _settings_line(name, args, repo)
    return outcome


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def _assemble(
    plan: dict, profile: Profile, repo: Path, state: rig.GitState, session: Path, outcomes: list[SuiteOutcome]
) -> Report:
    report = Report().heading(f"Suffix-array index — {state.short}", level=1)

    _provenance(report, profile, state, session, outcomes)

    report.heading("Suites")
    table = Table(headers=["suite", "status", "wall clock", "note"], aligns=["<", "<", ">", "<"])
    for outcome in outcomes:
        table.row(outcome.name, outcome.status, f"{outcome.seconds / 60:.1f} min", outcome.reason)
    report.table(table)

    for outcome in outcomes:
        report.heading(outcome.name, level=2)
        if outcome.status == "ok" and outcome.report:
            # The suite's own analysis, verbatim. Its top-level heading is dropped so it does not
            # nest a duplicate title under this one.
            report.blocks.extend(outcome.report.blocks[1:])
        else:
            report.para(f"**not run** — {outcome.reason or 'no results'}")

    report.heading("Caveats")
    caveats = [
        f"{outcome.name}: {outcome.reason}" for outcome in outcomes if outcome.status != "ok"
    ]
    if state.dirty:
        caveats.append("the working tree was dirty")
    report.lines([f"  * {line}" for line in caveats] or ["  (none)"])
    return report


def _provenance(
    report: Report, profile: Profile, state: rig.GitState, session: Path, outcomes: list[SuiteOutcome]
) -> None:
    """Everything needed to reproduce this run, or to decide it cannot be compared with another.

    Four groups, in the order a sceptical reader wants them: what code, on what data, on what
    machine, with what settings.
    """
    report.heading("Provenance")
    facts = rig.host_facts()
    table = Table(headers=["what", "value"], aligns=["<", "<"])

    table.row("commit", f"{state.commit}  ({state.branch})")
    table.row("tree", "DIRTY — these numbers are not attributable to the commit above" if state.dirty else "clean")
    table.row("toolchain", f"{facts['rustc']} · {facts['cargo']}")
    table.row("driver", f"python {facts['python']}")
    table.row("", "")

    for label, value in profile.describe():
        table.row(label, value)
    table.row("", "")

    table.row("cpu", f"{facts['cpu']} · {facts['cores']} logical cores · {facts['arch']}")
    swap = facts["swap_gb"]
    table.row("memory", f"{facts['ram_gb']} GB RAM, swap {swap if not swap.isdigit() else swap + ' GB'}")
    table.row("kernel", facts["kernel"])
    table.row("cgroups", f"{facts['cgroup']} · systemd-run {facts['systemd_run']}")
    # Spelled out because "13.64" beside a core count is meaningless without knowing which it is.
    table.row("load average", f"{facts['load']} (1 / 5 / 15 min) at the start of the session")
    table.row("", "")

    table.row("session", str(session))
    table.row("started", time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    for outcome in outcomes:
        if outcome.settings:
            table.row(f"  {outcome.name}", outcome.settings)
    report.table(table)

    if state.dirty:
        report.warn("working tree was dirty — these numbers are not attributable to that commit")
    cores = int(facts["cores"] or 1)
    if rig.load_average()[0] > cores * 0.25:
        report.warn(
            "the box was busy when this started — a co-tenant job invalidates every comparison in "
            "this report"
        )









def _settings_line(name: str, args, repo: Path) -> str:
    """What this suite was actually run with — the other half of reproducing a number."""
    try:
        suite = config.load(name, repo)
    except config.ConfigError:
        return ""
    defaults = dict(suite.defaults)
    if args.runs is not None:
        defaults["runs"] = args.runs
    if args.amount is not None:
        defaults["amount"] = args.amount

    parts = [
        f"{defaults.get('runs', '?')} reps x {defaults.get('amount', '?'):,} peptides",
        f"arms {'/'.join(arm.name for arm in suite.arms)}",
    ]
    for key, label in (("peptides", "queries"), ("kmer_table", "kmer"), ("warmup", "warmup")):
        if defaults.get(key):
            parts.append(f"{label} {defaults[key]}")
    if suite.metrics:
        parts.append("instrumented (metrics)")
    if suite.axes:
        parts.append("axes " + ", ".join(f"{axis}={values}" for axis, values in sorted(suite.axes.items())))
    if suite.mode == "matrix":
        # A matrix suite's grid lives under [matrix], not [axes] — without this its line would say
        # nothing about what was actually swept.
        grid = suite.matrix
        parts.append(f"files {'/'.join(grid.get('files', []))}")
        parts.append(f"mlp batches {grid.get('batches', [])}")
        parts.append(f"6-mer {'included' if grid.get('kmer6') else 'excluded'}")
    return " · ".join(parts)


def _machine_readable(state: rig.GitState, outcomes: list[SuiteOutcome]) -> dict:
    """`report.json`: the same facts without the prose, so a later run can use this as a baseline."""
    return {
        "commit": state.commit,
        "branch": state.branch,
        "dirty": state.dirty,
        "suites": {
            outcome.name: {
                "status": outcome.status,
                "reason": outcome.reason,
                "seconds": round(outcome.seconds, 1),
                "settings": outcome.settings,
            }
            for outcome in outcomes
        },
    }


def _load_plan(repo: Path) -> dict:
    path = config.suites_dir(repo) / "all.toml"
    if not path.exists():
        raise config.ConfigError(f"no master plan at {path}")
    with path.open("rb") as handle:
        plan = tomllib.load(handle)
    if "order" not in plan:
        raise config.ConfigError(f"{path}: missing 'order'")
    missing = [name for name in plan["order"] if not (config.suites_dir(repo) / f"{name}.toml").exists()]
    if missing:
        raise config.ConfigError(f"{path}: order names suites that do not exist: {', '.join(missing)}")
    return plan
