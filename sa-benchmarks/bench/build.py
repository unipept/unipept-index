"""Building one harness binary per arm, up front.

Two rules, both learned the hard way:

**Never alternate builds and timed runs.** Every arm is built before the first cell runs, so a
compile cannot land in the middle of a sweep and charge its CPU time to whichever cell was unlucky.

**Prove the binaries differ.** Every build writes the same `target/release/sa-benchmarks`, so each
is copied out immediately. More importantly, a feature that is declared but not forwarded through
every crate compiles fine and produces a byte-identical binary — a silently meaningless arm that
reports plausible numbers for a configuration that was never built. Comparing the copies pairwise is
what catches it.
"""

from __future__ import annotations

import subprocess
from filecmp import cmp as file_cmp
from pathlib import Path

from .config import Arm, Suite
from .rig import RigError, as_user, dropping_privileges


def build_arms(suite: Suite, repo: Path, bin_dir: Path, echo=print) -> dict[str, Path]:
    """Builds every arm of `suite` into `bin_dir`, returning arm name -> binary path.

    Resumable: an arm whose binary is already present is not rebuilt, so an interrupted session
    restarts without paying for the builds again.
    """
    bin_dir.mkdir(parents=True, exist_ok=True)
    binaries: dict[str, Path] = {}

    for arm in suite.arms:
        target = bin_dir / arm.name
        features = _features(arm, suite)
        # The manifest is what makes sharing one bin/ across suites safe: two suites may both call
        # an arm "mmap", and reusing a binary built from a different feature set would measure the
        # wrong configuration under the right name.
        manifest = bin_dir / f"{arm.name}.features"
        if target.exists() and manifest.exists() and manifest.read_text().strip() == features:
            echo(f"  skip build {arm.name} (already in {bin_dir})")
            binaries[arm.name] = target
            continue
        echo(f"== build {arm.name} (features: {features or 'none'}) ==")
        _cargo_build(repo, features)
        # Copy immediately: the next arm's build overwrites this exact path.
        target.write_bytes((repo / "target" / "release" / "sa-benchmarks").read_bytes())
        target.chmod(0o755)
        manifest.write_text(f"{features}\n")
        binaries[arm.name] = target

    _assert_distinct(suite, binaries)
    return binaries


def _features(arm: Arm, suite: Suite) -> str:
    """The feature string for this arm, with `measure` folded in when either level asks for it."""
    features = list(arm.features)
    if (suite.measure or arm.measure) and "measure" not in features:
        features.append("measure")
    return ",".join(features)


def _cargo_build(repo: Path, features: str) -> None:
    command = ["cargo", "build", "--release", "-q", "-p", "sa-benchmarks", "--no-default-features"]
    if features:
        command += ["--features", features]

    if dropping_privileges():
        # A login shell, because cargo lives on the user's PATH and sudo's secure_path drops it.
        # `cd` inside the shell rather than passing cwd=, which sudo would evaluate as root.
        wrapped = as_user(
            ["bash", "-lc", f'cd "$1" && {subprocess.list2cmdline(command)}', "_", str(repo)]
        )
        result = subprocess.run(wrapped, check=False)
    else:
        result = subprocess.run(command, cwd=repo, check=False)

    if result.returncode != 0:
        raise RigError(f"cargo build failed for features '{features or 'none'}'")


def _assert_distinct(suite: Suite, binaries: dict[str, Path]) -> None:
    """Every pair of arms must differ in bytes; identical arms mean a feature never took effect."""
    names = sorted(binaries)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            if file_cmp(binaries[left], binaries[right], shallow=False):
                left_features = _features(next(a for a in suite.arms if a.name == left), suite)
                right_features = _features(next(a for a in suite.arms if a.name == right), suite)
                raise RigError(
                    f"arms '{left}' and '{right}' produced byte-identical binaries\n"
                    f"  '{left}'  = features {left_features or 'none'}\n"
                    f"  '{right}' = features {right_features or 'none'}\n"
                    f"  A feature is probably not forwarded through every crate, which would make "
                    f"one of these arms a meaningless duplicate of the other."
                )
