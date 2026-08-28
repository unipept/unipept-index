#!/usr/bin/env bash
# Answer-equivalence gate between fbc9328 and the branch. RUN THIS BEFORE ANY TIMING COMPARISON.
#
# A baseline that answers differently is not a baseline, it is a different program, and none of the
# benchmark suites would notice. A faster tree that answers differently is not faster.
#
# Each side goes through ITS OWN production entry point rather than anything written for this
# comparison — peptide normalisation, the sparseness-factor length filter and the cutoff all differ
# between the trees, and reimplementing them here would be testing the reimplementation:
#
#   branch   `sa-server`, queried over HTTP, built with --features mmap. The mmap arm because it is
#            by far the cheapest to load and `check_answers.sh` on the branch already establishes
#            that all nine of its storage configurations answer identically; this gate is about the
#            gap between COMMITS, not within one.
#   fbc9328  the ported harness's `--answers`, which calls `peptide_search::search_all_peptides` —
#            the exact function this commit's own `sa-server` handler calls.
#
# Why not fbc9328's sa-server: it loads the 160 GiB suffix array BEFORE building the protein store,
# and that store is built via a 69 GiB intermediate String. On the full database that is a 352 GiB
# peak and the process dies. The harness loads in the opposite order and can also stream the text
# directly into its packed form; see `stream_proteins_from_database_file`. Same library call either
# way, so the same answers — with a load that survives.
#
# The index and peptide file come from the BRANCH tree's machine profile, so this cannot drift from
# what the suites measure. Override with IDX= / PEP=.
#
#   bash sa-benchmarks/check_answers_cross.sh <branch-tree> [n_peptides]
#
#   IL="true"            only equate_il=true (default is both; each costs one fbc9328 index load)
#   STREAM=1             use --stream-proteins on the fbc9328 side (needed on a box that cannot
#                        hold fbc9328's own loader — the answers are identical either way, which
#                        this gate is itself the check for)
set -uo pipefail

OLD_TREE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BRANCH_TREE="${1:?usage: check_answers_cross.sh <branch-tree> [n_peptides]}"
N="${2:-2000}"

PROFILE="${PROFILE:-local}"
CUTOFF="${CUTOFF:-10000}"
OUT="${OUT:-${TMPDIR:-/tmp}/cross-answers}"
PORT="${PORT:-3998}"
CARGO="${CARGO:-cargo}"
IL="${IL:-true false}"
STREAM_ARG=""; [ "${STREAM:-0}" = "1" ] && STREAM_ARG="--stream-proteins"

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
# proteins.tsv is checked because it is the one file the branch no longer needs at runtime, and so
# the one an index directory can plausibly be missing. fbc9328 cannot start without it.
for f in sa.bin proteins.bin mapping.bin proteins.tsv; do
  [ -s "$IDX/$f" ] || { echo "ERROR: missing $IDX/$f"; exit 1; }
done

mkdir -p "$OUT"; rm -f "$OUT"/*.json
head -"$N" "$PEP" > "$OUT/peptides.txt"
ACTUAL=$(wc -l < "$OUT/peptides.txt" | tr -d ' ')
echo "index    : $IDX"
echo "peptides : $PEP (first $ACTUAL)"
echo "equate_il: $IL${STREAM_ARG:+   (fbc9328 loader: streaming)}"

# ---- branch: its own sa-server ------------------------------------------
echo
echo "===== branch ====="
( cd "$BRANCH_TREE" && $CARGO build --release -q -p sa-server --no-default-features --features mmap ) || exit 1
"$BRANCH_TREE/target/release/sa-server" -d "$IDX/proteins.bin" -i "$IDX/sa.bin" \
  --mapping-file "$IDX/mapping.bin" --address "127.0.0.1:$PORT" > "$OUT/branch.log" 2>&1 &
PID=$!
ready=0
for _ in $(seq 1 1800); do
  sleep 2
  curl -s -o /dev/null "http://127.0.0.1:$PORT/search" -X POST -H 'content-type: application/json' \
    -d '{"peptides":[]}' && { ready=1; break; }
  kill -0 $PID 2>/dev/null || break
done
[ "$ready" = 1 ] || { echo "branch server never came up — see $OUT/branch.log"; kill $PID 2>/dev/null; exit 1; }
for il in $IL; do
  python3 -c 'import sys,json; print(json.dumps({"peptides":[l.strip() for l in open(sys.argv[1]) if l.strip()],"equate_il":sys.argv[2]=="true","cutoff":int(sys.argv[3])}))' \
    "$OUT/peptides.txt" "$il" "$CUTOFF" > "$OUT/body_$il.json"
  curl -s --max-time 3600 "http://127.0.0.1:$PORT/search" -H 'content-type: application/json' \
    --data-binary "@$OUT/body_$il.json" > "$OUT/branch_il$il.json"
  echo "  il=$il -> $(wc -c < "$OUT/branch_il$il.json") bytes"
done
kill $PID 2>/dev/null; wait $PID 2>/dev/null

# ---- fbc9328: the ported harness ----------------------------------------
echo
echo "===== fbc9328 ====="
( cd "$OLD_TREE" && $CARGO build --release -q -p sa-benchmarks ) || exit 1
for il in $IL; do
  "$OLD_TREE/target/release/sa-benchmarks" \
    --index-dir "$IDX" --output "$OUT/unused" --label answers \
    --peptide-file "$OUT/peptides.txt" --amount-of-peptides "$ACTUAL" \
    --equate-il "$il" --tryptic false --max-matches "$CUTOFF" \
    --answers "$OUT/old_il$il.json" $STREAM_ARG 2>&1 | grep -E "proteins,|items,|Wrote|!!" || exit 1
done

# ---- compare -------------------------------------------------------------
python3 - "$OUT" "$IL" <<'PY'
import json, sys, pathlib
out, ils = pathlib.Path(sys.argv[1]), sys.argv[2].split()

def canon(path):
    """peptide -> (cutoff_used, sorted accessions).

    Sorted because hit ORDER follows suffix order, which is not part of the answer — two trees may
    walk a match range differently and still agree on what is in it."""
    data = json.loads(path.read_text())
    return {r["sequence"]: (r["cutoff_used"],
                            tuple(sorted(p["uniprot_accession"] for p in r["proteins"])))
            for r in data}

status = 0
for il in ils:
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
