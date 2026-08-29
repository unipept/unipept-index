#!/usr/bin/env bash
# Stages this tree's TWO c4f0f30 arms into a benchmark session, so the branch's driver runs them as
# ordinary arms.
#
# Two, not one: at this commit the storage choice is a runtime bool covering all three structures,
# so `mmap` and `preloaded` are the only configurations that exist — and both are built here,
# because the driver cannot build either from the branch tree.
#
# The driver skips an arm whose binary is already in `<session>/bin` with a matching `.features`
# manifest; that is its resume path, and putting the binaries there is what turns "runs in another
# tree whose records get merged afterwards" into "one session, one ordering, one report". That
# matters because `baseline` uses palindrome ordering, which interleaves the arms so
# between-invocation drift is measured rather than absorbed into an arm difference. Separate runs
# cannot be interleaved, and merging them afterwards silently discards it.
#
# Re-run whenever this tree's harness changes: the driver skips a staged arm without checking how
# old it is, and a stale binary writes perfectly well-formed records.
#
#   bash sa-benchmarks/stage.sh <session-dir>
set -euo pipefail

WORKTREE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SESSION="${1:?usage: stage.sh <session-dir>}"
# c4f0f30 has no rust-toolchain.toml, so it takes the box default. CARGO=... to override.
CARGO="${CARGO:-cargo}"

BIN="$SESSION/bin"
mkdir -p "$BIN"

stage () {  # arm-name  cargo-features
  local arm="$1" features="$2"
  echo "building $arm (features: ${features:-none}) ..."
  if [ -n "$features" ]; then
    ( cd "$WORKTREE" && $CARGO build --release -q -p sa-benchmarks --features "$features" )
  else
    ( cd "$WORKTREE" && $CARGO build --release -q -p sa-benchmarks )
  fi
  cp "$WORKTREE/target/release/sa-benchmarks" "$BIN/$arm"
  chmod 755 "$BIN/$arm"
  # The manifest the driver compares against the arm's feature string in the suite file.
  printf '%s\n' "$features" > "$BIN/$arm.features"
  echo "  staged -> $BIN/$arm"
}

stage c4f0f30-preloaded ""
stage c4f0f30-mmap      "mmap"

# The one failure this cannot detect later: two arms that compiled to the same bytes would mean the
# feature never reached the code, and the records would be plausible and wrong. The driver checks
# this too, but only after the index has been loaded.
if cmp -s "$BIN/c4f0f30-preloaded" "$BIN/c4f0f30-mmap"; then
  echo "ERROR: the two arms are byte-identical — the 'mmap' feature did not take effect" >&2
  exit 1
fi

echo
echo "both arms staged, and they differ."
echo "  from commit : $(cd "$WORKTREE" && git rev-parse HEAD)"
echo "  the driver will now report 'skip build' for both and run them as ordinary arms."
