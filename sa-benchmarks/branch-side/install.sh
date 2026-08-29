#!/usr/bin/env bash
# Copies the baseline suite definitions into a branch checkout. See README.md next to this script.
#
#   bash sa-benchmarks/branch-side/install.sh <branch-tree>
#   bash sa-benchmarks/branch-side/install.sh <branch-tree> --uninstall
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BRANCH="${1:?usage: install.sh <branch-tree> [--uninstall]}"
MODE="${2:-install}"

[ -d "$BRANCH/sa-benchmarks/suites" ] || { echo "ERROR: $BRANCH does not look like the branch checkout"; exit 1; }

copy () {  # src dst
  if [ "$MODE" = "--uninstall" ]; then
    rm -f "$2" && echo "  removed $2"
    return
  fi
  # Refuse rather than clobber: these are untracked in the branch tree, so an overwrite is
  # unrecoverable, and someone editing a suite there is the expected case rather than a mistake.
  if [ -e "$2" ] && ! cmp -s "$1" "$2"; then
    echo "  SKIP $2 — exists and differs; delete it first if you want this version"
    return
  fi
  cp "$1" "$2" && echo "  $2"
}

for f in baseline.toml baseline_startup.toml; do
  copy "$HERE/suites/$f" "$BRANCH/sa-benchmarks/suites/$f"
done
for f in baseline.py baseline_startup.py; do
  copy "$HERE/bench-suites/$f" "$BRANCH/sa-benchmarks/bench/suites/$f"
done

if [ "$MODE" != "--uninstall" ]; then
  echo
  echo "installed. These four files are UNTRACKED in the branch tree — 'git status' will show them,"
  echo "and '--uninstall' removes them again."
fi
