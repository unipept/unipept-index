"""Command line: `./sa-benchmarks/run.sh <suite> [options]`.

    run.sh defaults              the production-defaults sweep (the regression gate)
    run.sh all --check           can this box run it, and how big is it? nothing is built or written
    run.sh ram --dry-run         ... plus the per-cell plan, without touching the index
    run.sh all                   every suite in one session, into one report.md
    run.sh all --baseline <dir>  ... and diff it against a previous session
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from . import build, config, preflight, records, rig
from .html import render_page
from .profile import ProfileError, load as load_profile
from .progress import Progress
from .report import Report
from .runner import Runner
from .suites import analysis_for


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)

    try:
        repo = rig.repo_root()
    except rig.RigError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    if args.check:
        return run_checks(args, repo)

    if args.suite == "all":
        from .fullreport import run_all

        return run_all(args, repo)

    try:
        status, report = run_suite(args.suite, args, repo)
    except (ProfileError, config.ConfigError, rig.RigError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    if report is not None:
        print(report.to_text())
    return status


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="run.sh", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("suite", help="suite name, or 'all' for every suite plus a combined report")
    parser.add_argument("--profile", default="local", help="machine profile (default: local)")
    parser.add_argument("--out", type=Path, help="session directory (default: <scratch>/<suite>)")
    parser.add_argument(
        "--check",
        action="store_true",
        help="run the preflight checks and print the plan table, then exit — nothing is built, "
        "nothing is written, exit status says whether the run would be able to start",
    )
    parser.add_argument("--dry-run", action="store_true", help="the checks plus the per-cell plan, then exit")
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="re-render an existing session's report from its jsonl, without running anything",
    )
    parser.add_argument("--runs", type=int, help="override the suite's timed reps per cell")
    parser.add_argument("--amount", type=int, help="override the suite's peptides per rep")
    parser.add_argument("--only", help="comma-separated substrings; run only matching cells")
    parser.add_argument("--skip", help="comma-separated suite names to skip (with 'all')")
    parser.add_argument("--baseline", type=Path, help="a previous session directory to compare against")
    parser.add_argument(
        "--cold",
        action="store_true",
        help="drop the page cache before every cell (needs root) — first-boot rather than warm-restart numbers",
    )
    return parser


def run_checks(args, repo: Path) -> int:
    """`--check`: is this machine in a state to run what was asked, and how big is it?

    Everything the preflight prints at the start of a real run, and nothing else — no build, no
    cells, no session directory. Cheap enough to run before walking away from the terminal, and its
    exit status is the answer, so it can gate an overnight run from a shell script:

        ./sa-benchmarks/run.sh all --check && nohup ./sa-benchmarks/run.sh all &

    Under `all` the session is the one a run started now would create, so nothing is reported as
    already complete. Pass `--out <dir>` to ask the question of a session being resumed instead.
    """
    from .fullreport import load_plan

    try:
        profile = load_profile(args.profile, repo)
        if args.suite == "all":
            plan = load_plan(repo)
            skip = {name.strip() for name in (args.skip or "").split(",") if name.strip()}
            names = [name for name in plan["order"] if name not in skip]
            optional = set(plan.get("optional", []))
        else:
            names, optional = [args.suite], set()
    except (ProfileError, config.ConfigError, rig.RigError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    state = rig.git_state(repo)
    session = args.out or (
        profile.scratch / f"{state.short}-{time.strftime('%Y%m%d-%H%M%S')}"
        if args.suite == "all"
        else profile.scratch / args.suite
    )
    flight = preflight.check(names, args, repo, profile, session, optional=optional)
    print("\n".join(preflight.render(flight, header=f"preflight — {', '.join(names)}")))
    return 0 if flight.ok else 1


def run_suite(
    name: str,
    args,
    repo: Path,
    *,
    session: Path | None = None,
    echo=print,
    progress: Progress | None = None,
    checked: bool = False,
) -> tuple[int, Report | None]:
    """Runs one suite end to end.

    Returns `(exit status, report)`. The report is handed back rather than printed so the master
    run can splice the same object into `report.md` — one analysis, two destinations.

    `checked` says the caller has already run the preflight for this suite (`run.sh all` runs one
    for the whole session, before the first build, which is the only way a session can refuse
    something it will not be able to do four hours from now).
    """
    profile = load_profile(args.profile, repo)
    suite = config.load(name, repo)
    apply_overrides(suite, args)

    # A standalone run gets its own session dir; under `all` the session is shared, and each suite
    # writes into its own subdirectory beside the one `bin/` every suite builds into.
    session = session or args.out or (profile.scratch / name)
    out_dir = preflight.results_dir(session, name)
    bin_dir = session / "bin"

    runner = Runner(suite, profile, {}, out_dir, echo=echo, progress=progress)

    if not (checked or args.report_only):
        flight = preflight.check([name], args, repo, profile, session)
        echo("\n".join(preflight.render(flight)))
        # A dry run reports what would block it and still prints the plan — being told the plan is
        # most of what it is for, and a machine that cannot run a suite is often not the machine the
        # sweep is being planned for.
        if not flight.ok and not args.dry_run:
            raise rig.RigError(f"preflight failed for suite '{name}' — nothing was run")
        if progress is None and not args.dry_run:
            progress = runner.progress = Progress(flight.weights)
            # Everything the run prints has to go above the bar rather than through it.
            echo = runner.echo = progress.echo

    cells = preflight.select(suite, runner, args.only)

    if args.dry_run:
        echo("")
        echo("\n".join(runner.plan(cells)))
        return 0, None

    if args.report_only:
        # Everything a suite concludes is already in its jsonl — `build_report` reads nothing else.
        # So changing how the report LOOKS never needs the sweep run again, which on the full
        # database is the difference between iterating on the page and iterating overnight.
        if not records.load_dir(out_dir):
            raise rig.RigError(f"nothing to re-render: no records under {out_dir}")
        report = build_report(suite, out_dir, [], None)
    else:
        blocking = rig.blockers(suite.needs_root, needs_cgroup(suite))
        if blocking:
            raise rig.RigError(f"cannot run suite '{name}': " + "; ".join(blocking))

        echo("")
        runner.check_supply(cells)
        runner.binaries = build.build_arms(suite, repo, bin_dir, echo=echo)

        started = time.monotonic()
        try:
            results = runner.run(cells)
        finally:
            # Whatever this suite did not get through must leave the session total, or a bar shared
            # with the suites after it would never reach 100%.
            if progress:
                progress.end_suite()
        elapsed = time.monotonic() - started

        report = build_report(suite, out_dir, results, elapsed)

    # A single-suite run gets the same page the master run produces, so there is no reason to run
    # `all` just to get something readable.
    state = rig.git_state(repo)
    page = out_dir / "report.html"
    page.write_text(
        render_page(
            f"{suite.name} — {state.short}",
            report,
            subtitle=f"{state.branch} · {args.profile} · {time.strftime('%Y-%m-%d %H:%M')}",
            statuses={f"{suite.name} — {suite.description or 'results'}": "ok"},
        )
    )
    (out_dir / "report.md").write_text(report.to_markdown())
    echo(f"\nreport: {page}")
    return 0, report


def build_report(suite: config.Suite, out_dir: Path, results, elapsed: float | None) -> Report:
    """Loads what the run produced and hands it to the suite's own analysis."""
    loaded = records.load_dir(out_dir)
    report = Report().heading(f"{suite.name} — {suite.description or 'results'}")
    if not loaded:
        return report.warn(f"no records in {out_dir}; every cell failed or was skipped")

    analysis_for(suite.name)(report, suite, loaded, out_dir)

    failures = [result for result in results if result.status == "failed"]
    if failures:
        report.warn(
            f"{len(failures)} cell(s) failed outright: "
            + ", ".join(f"{result.cell.label} ({result.detail})" for result in failures)
        )
    unfit = [result for result in results if result.status == "did-not-fit"]
    if unfit:
        report.para(
            "did not fit under its ceiling: " + ", ".join(result.cell.label for result in unfit)
        )
    # `elapsed is None` is a re-render of an earlier session: this process did not time anything,
    # and printing its own zero would overwrite what that run actually cost.
    report.para(
        f"raw jsonl in {out_dir}"
        if elapsed is None
        else f"wall clock {elapsed / 60:.1f} min; raw jsonl in {out_dir}"
    )
    return report


def apply_overrides(suite: config.Suite, args, session: Path | None = None) -> None:
    """`--runs` / `--amount` exist so a smoke run can be small without editing the suite file.

    They are caps, not assignments, wherever a `[[sweep]]` block set its own: a block that asks for
    fewer queries because its cells are slow should keep asking for fewer, but one that asks for
    more must not quietly ignore `--amount 300` and turn a smoke run into a real one.
    """
    if args.runs is not None:
        suite.defaults["runs"] = args.runs
    if args.amount is not None:
        suite.defaults["amount"] = args.amount
    for block in suite.sweeps:
        for key, cap in (("runs", args.runs), ("amount", args.amount)):
            if cap is not None and block.get(key) is not None:
                block[key] = min(block[key], cap)
    if getattr(args, "cold", False):
        # Dropping the cache is only meaningful, and only possible, with root.
        suite.drop_caches = True
        suite.needs_root = True
    if args.baseline:
        # `--baseline` names a session directory; a suite compares against its own results inside it.
        candidates = [args.baseline / suite.name, args.baseline / "results", args.baseline]
        suite.baseline = next((path for path in candidates if path.is_dir()), None)
        if suite.baseline is None:
            raise config.ConfigError(
                f"--baseline {args.baseline} holds no results for suite '{suite.name}' "
                f"(looked for {suite.name}/, results/, and the directory itself)"
            )


def needs_cgroup(suite: config.Suite) -> bool:
    """Only a memory ceiling needs cgroup v2 and `systemd-run`.

    A pinned thread count used to count too, which made every thread sweep need root on a box where
    setting `RAYON_NUM_THREADS` would have done. The runner now sets it on the child's environment
    directly and reaches for the scope only when there is a ceiling to impose.
    """
    return any(value != 0 for value in suite.axes.get("ceiling_gb", []))


if __name__ == "__main__":
    sys.exit(main())
