#!/usr/bin/env bash
# Correctness gate for per-structure storage selection: every one of the nine configurations must
# return byte-identical answers from one index. Run this BEFORE any storage comparison — a fast arm
# is worthless if it answers differently, and none of the suites would notice.
#
# This drives `sa-server` rather than the benchmark harness, which is why it is a script and not a
# suite: it is about answers, not timing.
#
# Paths come from the same machine profile the suites use, so this cannot drift from what they
# measure. Override the profile with PROFILE=..., or the paths directly with IDX= / PEP=.
#
#   bash sa-benchmarks/check_answers.sh
set -uo pipefail

# No `-e`: a configuration that fails to build or to answer is reported and skipped, so the other
# eight still get compared. Each failure path below is checked explicitly instead.

# The server runs in the background, so it has to be killed on every exit — not only the normal
# one. Left alive it keeps $PORT bound, and the next configuration's readiness probe then succeeds
# against the STALE server: every remaining config would record the previous build's answers and
# the gate would print PASS. A leaked process is not the risk; a false PASS is.
pid=""
cleanup() { [ -n "$pid" ] && kill "$pid" 2>/dev/null; return 0; }
trap 'cleanup; exit 130' INT TERM
trap cleanup EXIT

# ============================================================================
TREE="$(git rev-parse --show-toplevel)"
PROFILE="${PROFILE:-local}"

# Read index_dir and the "mixed" peptide file out of the profile, so there is one place per machine
# where these live. Falls back to IDX= / PEP= if the profile is missing.
read -r PROFILE_IDX PROFILE_PEP <<EOF
$(python3 - "$TREE" "$PROFILE" <<'PY'
import sys
sys.path.insert(0, sys.argv[1] + "/sa-benchmarks")
try:
    from bench.profile import load
    profile = load(sys.argv[2], __import__("pathlib").Path(sys.argv[1]))
    print(profile.index_dir, profile.peptides.get("mixed", ""))
except Exception:
    print(" ")
PY
)
EOF

IDX="${IDX:-$PROFILE_IDX}"
PEP="${PEP:-$PROFILE_PEP}"
OUT="${OUT:-/tmp/storage-answers}"
PORT="${PORT:-3999}"
# ============================================================================

[ -n "$IDX" ] || { echo "ERROR: no index dir — set IDX= or create sa-benchmarks/profiles/$PROFILE.toml"; exit 1; }
# `-s`, not just `-n`: the profile loader checks the file exists, but PEP= bypasses it, and an
# empty or missing file leaves `$PEPTIDES` as `[]`. Every configuration then answers the empty
# query identically and the gate reports PASS having compared nothing.
[ -n "$PEP" ] || { echo "ERROR: no peptide file — set PEP= or add [peptides].mixed to the profile"; exit 1; }
[ -s "$PEP" ] || { echo "ERROR: peptide file is missing or empty: $PEP"; exit 1; }

mkdir -p "$OUT"
# Cleared, not merely written into. A configuration whose build fails is skipped below, and its
# answer file from a PREVIOUS run would otherwise take part in the comparison — matching the others
# and reporting PASS for a configuration this invocation never built, never started and never
# queried. The count check at the end is the other half of that: nine files, or the gate says so.
rm -f "$OUT"/*.json
[ -d "$IDX" ] || { echo "ERROR: index dir not found: $IDX"; exit 1; }
for f in sa.bin proteins.bin mapping.bin; do
  [ -s "$IDX/$f" ] || { echo "ERROR: missing $IDX/$f"; exit 1; }
done

# A fixed query: five real peptides from the index's own peptide file.
PEPTIDES=$(head -5 "$PEP" | python3 -c 'import sys,json; print(json.dumps([l.strip() for l in sys.stdin if l.strip()]))')
# Checked on the QUERY, not just on the file. `-s` above rejects a missing or zero-byte file, but a
# file holding only blank lines is non-empty and still yields no peptides — and an empty query is
# answered identically by all nine configurations, so the gate would compare nothing and pass.
[ "$PEPTIDES" != "[]" ] || { echo "ERROR: no usable peptides in the first five lines of $PEP"; exit 1; }
BODY="{\"peptides\": $PEPTIDES, \"equate_il\": true, \"cutoff\": 1000}"

# All nine configurations: everything preloaded, everything mapped, and the seven mixtures.
CONFIGS=(
  ""
  "mmap"
  "mmap,preloaded-text"
  "mmap,preloaded-proteins"
  "mmap,preloaded-mapping"
  "mmap,preloaded-text,preloaded-proteins"
  "mmap,preloaded-text,preloaded-mapping"
  "mmap,preloaded-proteins,preloaded-mapping"
  "mmap,preloaded-text,preloaded-proteins,preloaded-mapping"
)

for f in "${CONFIGS[@]}"; do
  tag="${f:-preloaded}"; tag="${tag//,/+}"
  echo "===== $tag ====="
  ( cd "$TREE" && if [ -z "$f" ]; then
      cargo build --release -q -p sa-server --no-default-features
    else
      cargo build --release -q -p sa-server --no-default-features --features "$f"
    fi ) || { echo "  BUILD FAILED"; continue; }

  "$TREE/target/release/sa-server" \
      -d "$IDX/proteins.bin" -i "$IDX/sa.bin" --mapping-file "$IDX/mapping.bin" \
      --address "127.0.0.1:$PORT" > "$OUT/$tag.log" 2>&1 &
  pid=$!

  # Index loads are slow on the full DB; wait until the endpoint answers rather than guessing.
  ready=0
  for _ in $(seq 1 1800); do
    sleep 2
    if curl -s -o /dev/null "http://127.0.0.1:$PORT/search" -X POST \
      -H 'content-type: application/json' -d '{"peptides":[]}'; then ready=1; break; fi
    kill -0 $pid 2>/dev/null || { echo "  SERVER DIED — see $OUT/$tag.log"; break; }
  done
  # Querying a server that never came up leaves an empty file, which the comparison at the end
  # reports as "the configurations disagree" — the wrong diagnosis for a server that died on
  # startup. Leave no file instead: the count check names this configuration as the missing one.
  if [ "$ready" != "1" ]; then
    echo "  NOT READY — no answer recorded for $tag"
    kill $pid 2>/dev/null; wait $pid 2>/dev/null
    pid=""
    continue
  fi

  # Written to a temporary and moved only once it parses, for the reason the NOT READY branch
  # above gives: `>` creates the file before the pipeline runs, so a failed request leaves a
  # 0-byte answer behind and the comparison at the end reports "the configurations disagree" —
  # blaming the storage backends for a request that never succeeded. `curl -f` also rejects a
  # non-2xx body, which would otherwise be recorded as this configuration's answer.
  if curl -sf "http://127.0.0.1:$PORT/search" -H 'content-type: application/json' -d "$BODY" \
      | python3 -c 'import sys,json; print(json.dumps(json.load(sys.stdin), sort_keys=True, indent=1))' \
      > "$OUT/$tag.json.tmp"; then
    mv "$OUT/$tag.json.tmp" "$OUT/$tag.json"
    grep -i "Storage backends" "$OUT/$tag.log" || echo "  (no backend line — did it start?)"
    echo "  answer bytes: $(wc -c < "$OUT/$tag.json")"
  else
    echo "  QUERY FAILED — no answer recorded for $tag (see $OUT/$tag.log)"
    rm -f "$OUT/$tag.json.tmp"
  fi
  kill $pid 2>/dev/null; wait $pid 2>/dev/null
  pid=""
done

echo "===== comparison ====="
sha1sum "$OUT"/*.json 2>/dev/null || shasum "$OUT"/*.json

# Every configuration must have produced an answer THIS run. Comparing only the files that happen
# to be there would report PASS on a subset — a build that failed or a server that never came up
# leaves no file, and a gate that is silent about the arm it could not test is not a gate.
answers=$(ls -1 "$OUT"/*.json 2>/dev/null | wc -l | tr -d ' ')
echo "answers recorded: $answers of ${#CONFIGS[@]}"
[ "$answers" = "${#CONFIGS[@]}" ] || {
  echo "FAIL: ${#CONFIGS[@]} configurations were asked for, $answers answered — see the log above"
  exit 1
}

n=$( { sha1sum "$OUT"/*.json 2>/dev/null || shasum "$OUT"/*.json; } | awk '{print $1}' | sort -u | wc -l | tr -d ' ')
echo "distinct answer hashes: $n  (must be 1)"
[ "$n" = "1" ] || { echo "FAIL: configurations disagree"; exit 1; }
echo "PASS: every configuration returns identical answers"
