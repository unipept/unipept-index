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
[ -n "$PEP" ] || { echo "ERROR: no peptide file — set PEP= or add [peptides].mixed to the profile"; exit 1; }

mkdir -p "$OUT"
[ -d "$IDX" ] || { echo "ERROR: index dir not found: $IDX"; exit 1; }
for f in sa.bin proteins.bin mapping.bin; do
  [ -s "$IDX/$f" ] || { echo "ERROR: missing $IDX/$f"; exit 1; }
done

# A fixed query: five real peptides from the index's own peptide file.
PEPTIDES=$(head -5 "$PEP" | python3 -c 'import sys,json; print(json.dumps([l.strip() for l in sys.stdin if l.strip()]))')
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
  for _ in $(seq 1 1800); do
    sleep 2
    curl -s -o /dev/null "http://127.0.0.1:$PORT/search" -X POST \
      -H 'content-type: application/json' -d '{"peptides":[]}' && break
    kill -0 $pid 2>/dev/null || { echo "  SERVER DIED — see $OUT/$tag.log"; break; }
  done

  curl -s "http://127.0.0.1:$PORT/search" -H 'content-type: application/json' -d "$BODY" \
    | python3 -c 'import sys,json; print(json.dumps(json.load(sys.stdin), sort_keys=True, indent=1))' \
    > "$OUT/$tag.json"

  grep -i "Storage backends" "$OUT/$tag.log" || echo "  (no backend line — did it start?)"
  echo "  answer bytes: $(wc -c < "$OUT/$tag.json")"
  kill $pid 2>/dev/null; wait $pid 2>/dev/null
done

echo "===== comparison ====="
sha1sum "$OUT"/*.json 2>/dev/null || shasum "$OUT"/*.json
n=$( { sha1sum "$OUT"/*.json 2>/dev/null || shasum "$OUT"/*.json; } | awk '{print $1}' | sort -u | wc -l | tr -d ' ')
echo "distinct answer hashes: $n  (must be 1)"
[ "$n" = "1" ] || { echo "FAIL: configurations disagree"; exit 1; }
echo "PASS: every configuration returns identical answers"
