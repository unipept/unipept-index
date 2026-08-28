#!/usr/bin/env bash
# Answer-equivalence gate between fbc9328 and the branch. RUN THIS BEFORE ANY TIMING COMPARISON.
#
# Runs BOTH trees' own `sa-server` against ONE index and compares what `/search` returns, peptide by
# peptide. This is the check that makes the timing comparison mean anything: a baseline that answers
# differently is not a baseline, it is a different program, and none of the benchmark suites would
# notice. A faster tree that answers differently is not faster.
#
# Each side goes through its OWN production path rather than through anything written for this
# comparison. Peptide normalisation, the min-length filter and the cutoff all differ between the
# trees; reimplementing them here would be testing the reimplementation.
#
# The index and peptide file come from the BRANCH tree's machine profile, so this cannot drift from
# what the suites measure. Override with IDX= / PEP=.
#
#   bash sa-benchmarks/check_answers_cross.sh <branch-tree> [n_peptides]
set -uo pipefail

OLD_TREE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BRANCH_TREE="${1:?usage: check_answers_cross.sh <branch-tree> [n_peptides]}"
N="${2:-2000}"

PROFILE="${PROFILE:-local}"
CUTOFF="${CUTOFF:-10000}"
OUT="${OUT:-${TMPDIR:-/tmp}/cross-answers}"
# Two free ports. fbc9328's server has no --address flag and binds 0.0.0.0:3000, so only the
# branch's is configurable.
OLD_PORT=3000
NEW_PORT="${NEW_PORT:-3998}"
CARGO="${CARGO:-cargo}"

read -r PROFILE_IDX PROFILE_PEP <<EOF2
$(python3 - "$BRANCH_TREE" "$PROFILE" <<'PY'
import sys, pathlib
sys.path.insert(0, sys.argv[1] + "/sa-benchmarks")
try:
    from bench.profile import load
    profile = load(sys.argv[2], pathlib.Path(sys.argv[1]))
    print(profile.index_dir, profile.peptides.get("mixed", ""))
except Exception:
    print(" ")
PY
)
EOF2
IDX="${IDX:-$PROFILE_IDX}"
PEP="${PEP:-$PROFILE_PEP}"

[ -n "$IDX" ] || { echo "ERROR: no index dir — set IDX= or create $BRANCH_TREE/sa-benchmarks/profiles/$PROFILE.toml"; exit 1; }
[ -n "$PEP" ] || { echo "ERROR: no peptide file — set PEP= or add [peptides].mixed to the profile"; exit 1; }
# fbc9328 reads the database TSV, not proteins.bin. Checked up front: it is the one file the branch
# no longer needs at runtime, so it is the one an index directory can plausibly be missing.
for f in sa.bin proteins.bin mapping.bin proteins.tsv; do
  [ -s "$IDX/$f" ] || { echo "ERROR: missing $IDX/$f"; exit 1; }
done

mkdir -p "$OUT"; rm -f "$OUT"/*.json
echo "index    : $IDX"
echo "peptides : $PEP (first $N)"

PEPTIDES=$(head -"$N" "$PEP" | python3 -c 'import sys,json; print(json.dumps([l.strip() for l in sys.stdin if l.strip()]))')
for il in true false; do
  echo "{\"peptides\": $PEPTIDES, \"equate_il\": $il, \"cutoff\": $CUTOFF}" > "$OUT/body_il$il.json"
done

wait_ready () {  # port pid
  for _ in $(seq 1 1800); do
    sleep 2
    curl -s -o /dev/null "http://127.0.0.1:$1/search" -X POST -H 'content-type: application/json' \
      -d '{"peptides":[]}' && return 0
    kill -0 "$2" 2>/dev/null || return 1
  done
  return 1
}

query () {  # port tag
  for il in true false; do
    curl -s --max-time 1800 "http://127.0.0.1:$1/search" -H 'content-type: application/json' \
      --data-binary "@$OUT/body_il$il.json" > "$OUT/$2_il$il.json"
    echo "  il=$il -> $(wc -c < "$OUT/$2_il$il.json") bytes"
  done
}

echo "===== branch ====="
( cd "$BRANCH_TREE" && $CARGO build --release -q -p sa-server --no-default-features ) || exit 1
"$BRANCH_TREE/target/release/sa-server" -d "$IDX/proteins.bin" -i "$IDX/sa.bin" \
  --mapping-file "$IDX/mapping.bin" --address "127.0.0.1:$NEW_PORT" > "$OUT/branch.log" 2>&1 &
PID=$!
wait_ready "$NEW_PORT" $PID || { echo "branch server never came up — see $OUT/branch.log"; kill $PID 2>/dev/null; exit 1; }
query "$NEW_PORT" branch
kill $PID 2>/dev/null; wait $PID 2>/dev/null

# `-d` is the database TSV here, not proteins.bin: at this commit the server builds the protein store
# and the mapping from it at startup. Same index, different entry point — which is the thing being
# checked. This build pulls in sa-builder -> libsais64-rs, whose build script clones and cmake-builds
# libsais-packed; that needs cmake, make and network on first run.
echo "===== fbc9328 ====="
( cd "$OLD_TREE" && $CARGO build --release -q -p sa-server ) || exit 1
"$OLD_TREE/target/release/sa-server" -d "$IDX/proteins.tsv" -i "$IDX/sa.bin" > "$OUT/old.log" 2>&1 &
PID=$!
wait_ready "$OLD_PORT" $PID || { echo "fbc9328 server never came up — see $OUT/old.log"; kill $PID 2>/dev/null; exit 1; }
query "$OLD_PORT" old
kill $PID 2>/dev/null; wait $PID 2>/dev/null

python3 - "$OUT" <<'PY'
import json, sys, pathlib
out = pathlib.Path(sys.argv[1])

def canon(path):
    """peptide -> (cutoff_used, sorted accessions).

    Sorted because hit ORDER follows suffix order, which is not part of the answer — two trees may
    walk a match range differently and still agree on what is in it."""
    data = json.loads(path.read_text())
    return {r["sequence"]: (r["cutoff_used"],
                            tuple(sorted(p["uniprot_accession"] for p in r["proteins"])))
            for r in data}

status = 0
for il in ("true", "false"):
    new, old = canon(out / f"branch_il{il}.json"), canon(out / f"old_il{il}.json")
    both = new.keys() & old.keys()
    differ = [p for p in both if new[p] != old[p]]
    only_new, only_old = new.keys() - old.keys(), old.keys() - new.keys()
    print(f"\n=== equate_il={il} ===")
    print(f"  answered by both : {len(both):,}   identical: {len(both) - len(differ):,}   DIFFER: {len(differ):,}")
    print(f"  only branch      : {len(only_new):,}")
    print(f"  only fbc9328     : {len(only_old):,}")
    for p in sorted(differ)[:5]:
        n_, o_ = new[p], old[p]
        print(f"    {p}: branch cutoff={n_[0]} {len(n_[1])} hits | fbc9328 cutoff={o_[0]} {len(o_[1])} hits")
        print(f"      only branch : {sorted(set(n_[1]) - set(o_[1]))[:6]}")
        print(f"      only fbc9328: {sorted(set(o_[1]) - set(n_[1]))[:6]}")
    for p in sorted(only_new)[:5]:
        print(f"    only branch answered : {p!r} ({len(new[p][1])} hits)")
    for p in sorted(only_old)[:5]:
        print(f"    only fbc9328 answered: {p!r} ({len(old[p][1])} hits)")
    if differ or only_new or only_old:
        status = 1

print()
print("PASS: the two trees return identical answers" if status == 0
      else "DIFFERENCES FOUND — read them above before trusting any timing comparison")
sys.exit(status)
PY
