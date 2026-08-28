#!/usr/bin/env bash
# Stages this tree's fbc9328 harness into a benchmark session, so the branch's driver runs it as an
# ordinary arm.
#
# The driver builds every declared arm from the tree it is invoked in, and it cannot build this one —
# fbc9328 is a different commit with a different sa-index API. What it CAN do is skip an arm whose
# binary is already in `<session>/bin` with a matching `.features` manifest; that is its resume path.
# Putting the ported binary there is what turns "a run in another tree whose records get merged
# afterwards" into "one session, one ordering, one report".
#
# That distinction is not cosmetic. `baseline` uses palindrome ordering, which interleaves the arms
# so between-invocation drift is measured rather than absorbed into an arm difference. Two separate
# runs cannot be interleaved, and merging them afterwards silently discards that.
#
# Re-run this whenever this tree's harness changes: the driver skips a staged arm without checking
# how old it is, and a stale binary writes perfectly well-formed records.
#
#   bash sa-benchmarks/stage.sh <session-dir>
set -euo pipefail

# This tree, wherever it has been put. No path in this script is machine-specific.
WORKTREE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SESSION="${1:?usage: stage.sh <session-dir>}"

# fbc9328 has no rust-toolchain.toml, so it takes whatever the default toolchain is. CARGO=... to
# override (e.g. `CARGO="rustup run 1.83.0 cargo"`) if the box's default cannot build it.
CARGO="${CARGO:-cargo}"

echo "building the fbc9328 harness in $WORKTREE ..."
( cd "$WORKTREE" && $CARGO build --release -q -p sa-benchmarks )

BIN="$SESSION/bin"
mkdir -p "$BIN"
cp "$WORKTREE/target/release/sa-benchmarks" "$BIN/fbc9328"
chmod 755 "$BIN/fbc9328"
# The manifest the driver compares against the arm's feature string. fbc9328 has no features and no
# storage selection to make, so this is an empty line — matching `features = []` in the suite files.
printf '\n' > "$BIN/fbc9328.features"

echo "staged -> $BIN/fbc9328"
echo "  from commit : $(cd "$WORKTREE" && git rev-parse HEAD)"
echo "  the driver will now report 'skip build fbc9328' and run it as an ordinary arm."
