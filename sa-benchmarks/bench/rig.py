"""The machine, and whether it can honestly run a given suite.

Every check here exists because its absence has silently invalidated a benchmark before. They fail
loudly and up front rather than producing a clean-looking, meaningless run:

* a cgroup ceiling that does not bind, because the page cache was never dropped;
* `--setenv` not reaching the child, so every cell in a thread sweep ran at the default count;
* a peptide file shorter than `runs * amount`, which shortens the later reps instead of failing;
* a co-tenant job on the box, which invalidates any comparison made while it runs.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


class RigError(Exception):
    """The machine cannot run what was asked for."""


# ---------------------------------------------------------------------------
# Identity: repo, commit, privileges
# ---------------------------------------------------------------------------


def repo_root() -> Path:
    out = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True, check=False
    )
    if out.returncode != 0:
        raise RigError("not inside a git repository")
    return Path(out.stdout.strip())


@dataclass
class GitState:
    commit: str
    short: str
    branch: str
    dirty: bool

    def describe(self) -> str:
        suffix = "  ** DIRTY — results are not attributable to this commit **" if self.dirty else ""
        return f"{self.short} ({self.branch}){suffix}"


def git_state(repo: Path) -> GitState:
    """Commit, branch and dirtiness of `repo`.

    Runs as the invoking user: git refuses a user-owned repository when called as root ("dubious
    ownership"), and a sweep that needs root would otherwise lose its provenance entirely.
    """

    def git(*args: str) -> str:
        result = subprocess.run(
            as_user(["git", "-C", str(repo), *args]), capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"

    # `status --porcelain`, not `diff --quiet`: the latter only sees unstaged edits, so a tree with
    # everything staged but nothing committed would report clean and the numbers would be attributed
    # to a commit that does not contain them.
    status = subprocess.run(
        as_user(["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"]),
        capture_output=True,
        text=True,
        check=False,
    )
    dirty = bool(status.stdout.strip())
    commit = git("rev-parse", "HEAD")
    return GitState(
        commit=commit, short=commit[:10], branch=git("rev-parse", "--abbrev-ref", "HEAD"), dirty=dirty
    )


def is_root() -> bool:
    return os.geteuid() == 0


def dropping_privileges() -> bool:
    """True when `as_user` will actually hand work back to the invoking user."""
    return is_root() and bool(os.environ.get("SUDO_USER"))


def as_user(command: list[str]) -> list[str]:
    """Wraps `command` so it runs as the invoking user even when the sweep runs as root.

    Builds and git must not run as root: sudo resets PATH so cargo is not found, root-owned
    artefacts in `target/` break the next ordinary build there, and git rejects a user-owned repo.
    """
    if dropping_privileges():
        return ["sudo", "-u", os.environ["SUDO_USER"], "-H", *command]
    return command


def warn_if_root_without_sudo_user() -> str | None:
    if is_root() and not os.environ.get("SUDO_USER"):
        return (
            "running as root with no SUDO_USER: builds and git checks will run as root, leaving "
            "root-owned files in target/. Prefer 'sudo ./sa-benchmarks/run.sh ...' from your account."
        )
    return None


# ---------------------------------------------------------------------------
# Capabilities
# ---------------------------------------------------------------------------


def has_cgroup_memory() -> bool:
    try:
        return "memory" in Path("/sys/fs/cgroup/cgroup.controllers").read_text().split()
    except OSError:
        return False


def has_systemd_run() -> bool:
    return shutil.which("systemd-run") is not None


def setenv_reaches_child() -> bool:
    """Proves `systemd-run --setenv` actually reaches the child before hours are spent on it.

    `--setenv` is the portable spelling; `-E` only exists in newer systemd. Without this probe a
    thread sweep runs every cell at the default thread count and looks entirely plausible.
    """
    if not has_systemd_run():
        return False
    result = subprocess.run(
        ["systemd-run", "--scope", "--quiet", "--setenv=RAYON_NUM_THREADS=42", "--", "env"],
        capture_output=True,
        text=True,
        check=False,
    )
    return "RAYON_NUM_THREADS=42" in result.stdout.splitlines()


def blockers(needs_root: bool, needs_cgroup: bool) -> list[str]:
    """Why this machine cannot run a suite — empty means it can.

    Returned rather than raised so the master run can skip a suite, say so in the report, and carry
    on with the rest instead of aborting the session.
    """
    reasons: list[str] = []
    if needs_root and not is_root():
        reasons.append("needs root (cgroup ceilings + drop_caches)")
    if needs_cgroup:
        if not has_cgroup_memory():
            reasons.append("cgroup v2 memory controller not available at /sys/fs/cgroup")
        if not has_systemd_run():
            reasons.append("systemd-run not found")
        elif is_root() and not setenv_reaches_child():
            reasons.append("systemd-run --setenv does not reach the child")
    return reasons


# ---------------------------------------------------------------------------
# State the run depends on
# ---------------------------------------------------------------------------


def drop_caches() -> None:
    """Drops the page cache. Root only.

    Not hygiene: cgroup v2 charges a page-cache page to the cgroup that FIRST faults it in. A page
    left resident by the previous cell is reused without being charged again, so without dropping,
    a constrained cell silently runs with more memory than its ceiling.
    """
    if not is_root():
        raise RigError("drop_caches needs root")
    subprocess.run(["sync"], check=True)
    Path("/proc/sys/vm/drop_caches").write_text("3\n")


def load_average() -> tuple[float, float, float]:
    return os.getloadavg()


def host_facts() -> dict[str, str]:
    """Everything about this machine that a reader would need to reproduce or discount a run.

    A benchmark number means nothing without the box it came from, and "the box" is more than a core
    count: the CPU model sets the memory hierarchy the whole index is tuned around, the kernel
    decides how `MemoryMax` and `drop_caches` behave, and the toolchain decides what was compiled.
    Anything that cannot be read on this platform reports "?" rather than being left out — a missing
    row is a fact about the run too.
    """
    facts = {
        "cpu": _cpu_model(),
        "cores": str(os.cpu_count() or 0),
        "kernel": f"{platform.system()} {platform.release()}",
        "arch": platform.machine(),
        "rustc": _first_line(["rustc", "--version"]),
        "cargo": _first_line(["cargo", "--version"]),
        "python": platform.python_version(),
    }
    facts.update(_memory())
    one, five, fifteen = load_average()
    facts["load"] = f"{one:.2f} / {five:.2f} / {fifteen:.2f}"
    facts["cgroup"] = (
        "v2 memory controller available" if has_cgroup_memory() else "no cgroup v2 memory controller"
    )
    facts["systemd_run"] = "present" if has_systemd_run() else "absent"
    return facts


def _cpu_model() -> str:
    if platform.system() == "Darwin":
        return _first_line(["sysctl", "-n", "machdep.cpu.brand_string"])
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "?"


def _memory() -> dict[str, str]:
    """Total RAM and swap in GB. The ceilings a `ram` sweep imposes are only readable against these."""
    try:
        meminfo = dict(
            (line.split(":", 1)[0], line.split(":", 1)[1].strip())
            for line in Path("/proc/meminfo").read_text().splitlines()
        )
        return {
            "ram_gb": f"{int(meminfo['MemTotal'].split()[0]) / 2**20:.0f}",
            "swap_gb": f"{int(meminfo['SwapTotal'].split()[0]) / 2**20:.0f}",
        }
    except (OSError, KeyError, ValueError):
        pass
    if platform.system() == "Darwin":
        total = _first_line(["sysctl", "-n", "hw.memsize"])
        if total.isdigit():
            # macOS has no fixed swap partition to report; the sweeps that care are Linux-only.
            return {"ram_gb": f"{int(total) / 2**30:.0f}", "swap_gb": "dynamic"}
    return {"ram_gb": "?", "swap_gb": "?"}


def _first_line(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return result.stdout.strip().splitlines()[0] if result.returncode == 0 and result.stdout.strip() else "?"
    except (OSError, IndexError):
        return "?"


def count_lines(path: Path) -> int:
    result = subprocess.run(["wc", "-l", str(path)], capture_output=True, text=True, check=True)
    return int(result.stdout.split()[0])


def check_peptide_supply(path: Path, needed: int) -> None:
    """The harness consumes `amount` lines per rep, sequentially.

    A short file does not fail — it shortens the later reps, so a run looks complete and measures
    something else. This is the check that turns that into an error.
    """
    have = count_lines(path)
    if have < needed:
        raise RigError(
            f"{path} has {have:,} lines, need >= {needed:,} (runs x amount). "
            f"Generate more with sa-benchmarks/src/generate_peptides.rs, or lower runs/amount."
        )
