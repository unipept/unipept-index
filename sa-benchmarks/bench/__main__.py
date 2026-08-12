"""Command line: `./sa-benchmarks/run.sh <suite> [options]`.

    run.sh defaults              the production-defaults sweep (the regression gate)
    run.sh ram --dry-run         plan a sweep without touching the index
    run.sh all                   every suite in one session, into one report.md
    run.sh all --baseline <dir>  ... and diff it against a previous session
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from . import build, config, records, rig
from .html import render_page
from .profile import ProfileError, load as load_profile
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
    parser.add_argument("--dry-run", action="store_true", help="print the plan and exit")
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


def run_suite(
    name: str,
    args,
    repo: Path,
    *,
    session: Path | None = None,
    echo=print,
) -> tuple[int, Report | None]:
    """Runs one suite end to end.

    Returns `(exit status, report)`. The report is handed back rather than printed so the master
    run can splice the same object into `report.md` — one analysis, two destinations.
    """
    profile = load_profile(args.profile, repo)
    suite = config.load(name, repo)
    _apply_overrides(suite, args)

    # A standalone run gets its own session dir; under `all` the session is shared, and each suite
    # writes into its own subdirectory beside the one `bin/` every suite builds into.
    session = session or args.out or (profile.scratch / name)
    out_dir = session / "results" if session.name == name else session / name
    bin_dir = session / "bin"

    cells = suite.expand()
    if args.only:
        wanted = [token.strip() for token in args.only.split(",") if token.strip()]
        cells = [cell for cell in cells if any(token in cell.label for token in wanted)]
        if not cells:
            raise config.ConfigError(f"--only '{args.only}' matched no cell in suite '{name}'")

    runner = Runner(suite, profile, {}, out_dir, echo=echo)

    if args.dry_run:
        echo("\n".join(runner.plan(cells)))
        blocking = rig.blockers(suite.needs_root, _needs_cgroup(suite))
        if blocking:
            echo("\nWOULD BE SKIPPED on this machine:")
            for reason in blocking:
                echo(f"  - {reason}")
        return 0, None

    blocking = rig.blockers(suite.needs_root, _needs_cgroup(suite))
    if blocking:
        raise rig.RigError(f"cannot run suite '{name}': " + "; ".join(blocking))

    echo("\n".join(_provenance(profile, repo)))
    runner.preflight(cells)
    runner.binaries = build.build_arms(suite, repo, bin_dir, echo=echo)

    started = time.monotonic()
    results = runner.run(cells)
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


def build_report(suite: config.Suite, out_dir: Path, results, elapsed: float) -> Report:
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
    report.para(f"wall clock {elapsed / 60:.1f} min; raw jsonl in {out_dir}")
    return report


def _apply_overrides(suite: config.Suite, args, session: Path | None = None) -> None:
    """`--runs` / `--amount` exist so a smoke run can be small without editing the suite file."""
    if args.runs is not None:
        suite.defaults["runs"] = args.runs
    if args.amount is not None:
        suite.defaults["amount"] = args.amount
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


def _needs_cgroup(suite: config.Suite) -> bool:
    return any(value != 0 for value in suite.axes.get("ceiling_gb", [])) or any(
        str(value) != "default" for value in suite.axes.get("threads", [])
    )


def _provenance(profile, repo: Path) -> list[str]:
    """The short form printed before a single-suite run starts; the report carries the full table."""
    state = rig.git_state(repo)
    facts = rig.host_facts()
    lines = [
        f"commit    : {state.describe()}",
        *(f"{label:<10}: {value}" for label, value in profile.describe()),
        f"host      : {facts['cpu']} · {facts['cores']} cores · {facts['ram_gb']} GB RAM",
        f"load      : {facts['load']} (1 / 5 / 15 min)",
    ]
    load_1min = rig.load_average()[0]
    cores = int(facts["cores"] or 1)
    if load_1min > cores * 0.25:
        lines.append(
            "  !! the box is busy — a co-tenant job invalidates every comparison in this run; abort"
        )
    if state.dirty:
        lines.append("  !! working tree is dirty — results are not attributable to this commit")
    warning = rig.warn_if_root_without_sudo_user()
    if warning:
        lines.append(f"  !! {warning}")
    return lines


if __name__ == "__main__":
    sys.exit(main())
