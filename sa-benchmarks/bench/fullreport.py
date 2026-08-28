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

from . import config, preflight, records, rig
from .html import render_page
from .profile import Profile, ProfileError, load as load_profile
from .progress import Progress
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
        plan = load_plan(repo)
        profile = load_profile(args.profile, repo)
    except (ProfileError, config.ConfigError, rig.RigError) as error:
        print(f"error: {error}")
        return 1

    state = rig.git_state(repo)
    if args.report_only:
        # Re-rendering names the session it re-renders. Falling back to a fresh timestamped
        # directory would silently produce an empty report instead of the one that was asked for.
        if not args.out:
            print("error: --report-only needs --out <session dir> naming the session to re-render")
            return 1
        if not args.out.is_dir():
            print(f"error: no session at {args.out}")
            return 1
    session = args.out or (profile.scratch / f"{state.short}-{time.strftime('%Y%m%d-%H%M%S')}")
    session.mkdir(parents=True, exist_ok=True)
    facts = _open_session(session, report_only=args.report_only, dry_run=args.dry_run)

    skip = {name.strip() for name in (args.skip or "").split(",") if name.strip()}
    names = [name for name in plan["order"] if name not in skip]
    optional = set(plan.get("optional", []))

    print(f"== session {session} ==")
    print(f"== suites: {', '.join(names)} ==\n")

    # One preflight for the whole session, before the first build. Suites that need root run last,
    # so without this a session can spend an evening on `defaults` and only then discover it was
    # never going to be able to run `ram` — and nobody learns how big the run is until it is over.
    progress = None
    if not args.report_only:
        flight = preflight.check(names, args, repo, profile, session, optional=optional)
        print("\n".join(preflight.render(flight, header=f"preflight — {len(names)} suites")))
        print()
        # Under --dry-run the checks are the deliverable, not a gate: a plan is usually being read on
        # a machine other than the one that will run it.
        if not flight.ok and not args.dry_run:
            print("error: preflight failed; nothing was run")
            return 1
        if not args.dry_run:
            progress = Progress(flight.weights)

    echo = progress.echo if progress else print
    outcomes: list[SuiteOutcome] = []
    for name in names:
        outcome = _run_one(run_suite, name, args, repo, session, optional, progress)
        outcomes.append(outcome)
        # Under a second means the suite never started (blocked, misconfigured); a clock there
        # would read as "it ran, and took no time".
        clock = f" in {outcome.seconds / 60:.1f} min" if outcome.seconds >= 1 else ""
        echo(f"-- {name}: {outcome.status}{clock}{(' — ' + outcome.reason) if outcome.reason else ''}\n")
        if not args.dry_run:
            _write_timings(session, outcomes)
    if progress:
        progress.close()

    if args.dry_run:
        return 0

    report = _assemble(plan, profile, repo, state, session, outcomes, facts)
    # Three renderings of one report object: the page to read, the markdown to paste into a PR, and
    # the JSON that records what this session was (see `_machine_readable` — not the baseline).
    (session / "report.html").write_text(
        render_page(
            f"Suffix-array index — {state.short}",
            report,
            # The session's date, not the render's: a `--report-only` redraw a week later must not
            # retitle the page with the day somebody happened to look at it.
            subtitle=f"{state.branch} · {profile.name} · {facts['started_short']}",
            statuses={outcome.name: outcome.status for outcome in outcomes},
        )
    )
    (session / "report.md").write_text(report.to_markdown())
    (session / "report.json").write_text(json.dumps(_machine_readable(state, outcomes), indent=2) + "\n")

    print(report.to_text())
    print(f"\nreport: {session / 'report.html'}")
    print(f"        {session / 'report.md'}")
    return 0


def _run_one(
    run_suite,
    name: str,
    args,
    repo: Path,
    session: Path,
    optional: set[str],
    progress: Progress | None = None,
) -> SuiteOutcome:
    """Runs one suite, converting a can't-run into a recorded skip when the suite is optional."""
    echo = progress.echo if progress else print
    echo(f"===== {name} =====")
    started = time.monotonic()
    try:
        # `checked`: the session-wide preflight above already cleared this suite, and re-running it
        # per suite would re-print the same block a dozen times.
        _, report = run_suite(name, args, repo, session=session, echo=echo, progress=progress, checked=True)
    except rig.RigError as error:
        # A suite that never started still holds its share of the bar; without this the session
        # could only ever reach the fraction the suites that did run were worth.
        if progress:
            progress.drop_suite(name)
        if name in optional:
            return SuiteOutcome(name, "skipped", reason=str(error), seconds=time.monotonic() - started)
        return SuiteOutcome(name, "failed", reason=str(error), seconds=time.monotonic() - started)
    except (ProfileError, config.ConfigError) as error:
        # A misconfigured suite must not take the whole session down: the other four still answer
        # their questions, and the report will say this one did not.
        if progress:
            progress.drop_suite(name)
        return SuiteOutcome(name, "failed", reason=str(error), seconds=time.monotonic() - started)

    elapsed = time.monotonic() - started
    if report is None:
        # `--dry-run` plans but produces nothing; calling that "ok" would read as a completed suite.
        return SuiteOutcome(name, "planned", seconds=elapsed)
    # A re-render's own duration is not what the suite cost, and printing it under "wall clock"
    # would quietly replace a six-minute sweep with the second it took to redraw its page.
    outcome = SuiteOutcome(name, "ok", report=report, seconds=0.0 if args.report_only else elapsed)
    outcome.settings = _settings_line(name, args, repo)
    return outcome


def _open_session(session: Path, *, report_only: bool, dry_run: bool) -> dict:
    """The host facts and wall-clock start for this session, sampled ONCE and persisted.

    These used to be read while the report was being assembled, which is the wrong end of the run:
    assembly happens after the last suite, so `started` printed the FINISH, and the load average
    was sampled seconds after `threads` had finished oversubscribing the box to 96 rayon threads.
    On the 2dfa6517b7 session that read 61.33 on twelve cores and tripped the co-tenant warning
    across a run that was clean — the harness accusing a neighbour of its own last cell.

    So the sample point is session start, and the answer is written to `session.json` next to
    `timings.json`. A `--report-only` redraw reads that file rather than the machine doing the
    redraw, which is usually a laptop and never the box the numbers came from.

    A resumed session keeps the facts from the invocation that created it: `started` is when the
    session began, not when its last suite was picked up, and the load average belongs to the same
    moment. `sampled` records which of the two happened, because a session predating this file has
    nothing to read back and its provenance has to say so rather than quietly report the redraw.
    """
    path = session / "session.json"
    if path.exists():
        try:
            facts = json.loads(path.read_text())
            facts.setdefault("sampled", "session start")
            return facts
        except (OSError, json.JSONDecodeError):
            pass  # Unreadable is the same as absent: sample now and say where it came from.

    now = time.localtime()
    facts = {
        "started": time.strftime("%Y-%m-%d %H:%M:%S %Z", now),
        "started_short": time.strftime("%Y-%m-%d %H:%M", now),
        "load": list(rig.load_average()),
        "host": rig.host_facts(),
        "sampled": "this re-render — the session predates session.json" if report_only else "session start",
    }
    # A dry run plans a session it does not start; stamping a start time on the directory would
    # date the real run that follows to whenever somebody last planned it.
    if not (report_only or dry_run):
        path.write_text(json.dumps(facts, indent=2) + "\n")
    return facts


def _write_timings(session: Path, outcomes: list[SuiteOutcome]) -> None:
    """Per-suite wall clock, rewritten after every suite finishes.

    `report.json` carries the same numbers, but only once the whole session is over — and the
    sessions whose timings are most worth having are the ones that were interrupted. Written after
    each suite so the file is always current for the suites that have run.
    """
    payload = {
        "total_minutes": round(sum(outcome.seconds for outcome in outcomes) / 60, 1),
        "suites": {
            outcome.name: {"status": outcome.status, "minutes": round(outcome.seconds / 60, 1)}
            for outcome in outcomes
        },
    }
    (session / "timings.json").write_text(json.dumps(payload, indent=2) + "\n")


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def _assemble(
    plan: dict,
    profile: Profile,
    repo: Path,
    state: rig.GitState,
    session: Path,
    outcomes: list[SuiteOutcome],
    facts: dict,
) -> Report:
    report = Report().heading(f"Suffix-array index — {state.short}", level=1)

    _provenance(report, profile, state, session, outcomes, facts)

    report.heading("Suites")
    table = Table(
        headers=["suite", "status", "wall clock", "share", "note"], aligns=["<", "<", ">", ">", "<"]
    )
    total = sum(outcome.seconds for outcome in outcomes)
    for outcome in outcomes:
        clock = "-" if outcome.seconds < 1 else f"{outcome.seconds / 60:.1f} min"
        share = f"{outcome.seconds / total * 100:.0f}%" if total and outcome.seconds >= 1 else "-"
        table.row(outcome.name, outcome.status, clock, share, outcome.reason)
    table.row("total", "", f"{total / 60:.1f} min", "", "")
    report.table(table)
    # Said because it decides what to cut when a session has to fit in a night: a suite's clock is
    # its own cells plus whichever arms it was the first to need.
    report.note(
        "Wall clock is per suite and includes the arm builds that suite paid for; every later suite "
        "reuses those binaries, so the first suite in the order carries the build cost for all of "
        "them. The same numbers, per suite, are in `timings.json` and `report.json`."
    )

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
    report: Report,
    profile: Profile,
    state: rig.GitState,
    session: Path,
    outcomes: list[SuiteOutcome],
    session_facts: dict,
) -> None:
    """Everything needed to reproduce this run, or to decide it cannot be compared with another.

    Four groups, in the order a sceptical reader wants them: what code, on what data, on what
    machine, with what settings.

    Every machine fact here was sampled by `_open_session` when the session began, and is read back
    from `session.json` rather than measured now — see that function for why the difference matters.
    """
    report.heading("Provenance")
    facts = session_facts["host"]
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
    # Spelled out because "13.64" beside a core count is meaningless without knowing which it is,
    # and `sampled` because the whole point of the number is WHEN it was taken.
    table.row("load average", f"{facts['load']} (1 / 5 / 15 min), sampled at {session_facts['sampled']}")
    table.row("", "")

    table.row("session", str(session))
    table.row("started", session_facts["started"])
    for outcome in outcomes:
        if outcome.settings:
            table.row(f"  {outcome.name}", outcome.settings)
    report.table(table)

    if state.dirty:
        report.warn("working tree was dirty — these numbers are not attributable to that commit")
    cores = int(facts["cores"] or 1)
    # The session-start load, never a fresh reading: at assembly time the box is still carrying
    # whatever the last suite left behind, and `threads` deliberately ends by oversubscribing it.
    if session_facts["load"][0] > cores * 0.25:
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
    for key, label in (("peptides", "queries"), ("kmer", "kmer"), ("warmup", "warmup")):
        if defaults.get(key):
            parts.append(f"{label} {defaults[key]}")
    if suite.measure:
        parts.append("instrumented (measure)")
    if suite.axes:
        parts.append("axes " + ", ".join(f"{axis}={values}" for axis, values in sorted(suite.axes.items())))
    if suite.mode == "matrix":
        # A matrix suite's grid lives in its [[sweep]] blocks, not in [axes] — without this its line
        # would say nothing about what was actually swept.
        parts.append("files " + "/".join(defaults.get("files", [])))
        parts.append("sweeps " + ", ".join(block.get("name", "?") for block in suite.sweeps))
    return " · ".join(parts)


def _machine_readable(state: rig.GitState, outcomes: list[SuiteOutcome]) -> dict:
    """`report.json`: what this session WAS, without the prose — commit, and each suite's outcome.

    Not the baseline, despite the name being the obvious candidate for one. `--baseline` names a
    session DIRECTORY and every suite reads its own jsonl out of it (see `__main__.apply_overrides`
    and `defaults._regressions`), because a comparison is per cell and this file holds no cells. It
    is provenance and status: which commit, which suites ran, what each cost.
    """
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


def load_plan(repo: Path) -> dict:
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
