#!/usr/bin/env bash
# Sweep the MLP cross-query batch size (SA_MLP_BATCH) for one backend and print
# throughput vs the scalar (B=1) baseline. Batching interleaves B independent
# peptide searches per rayon task to hide random-access DRAM latency.
#
# Usage (from the repo root):
#   bash scripts/mlp_sweep.sh [mmap|preloaded]
#
# Override any of these via env vars:
#   MLP_INDEX  = index dir with sa.bin/proteins.bin/mapping.bin/warmup.txt
#   MLP_PEP    = peptide file (needs >= AMT*RUNS lines)
#   MLP_RUNS   = timed runs per B (default 100)
#   MLP_AMT    = peptides per run (default 10000)
#   MLP_BATCHES= space-separated B values (default "1 4 8 16 32 64 128")
#   MLP_WORK   = scratch/results dir (default /tmp/mlp-sweep)
#   MLP_KMER   = k for an in-memory k-mer bounds table (empty = off). Big win: fewer
#                random probes -> fewer page walks. Try 5. (In production, prefer a
#                pre-built table via sa-builder --output-kmer-table + sa-server.)
#   MLP_HUGEPAGE = 1 to set SA_MADV_HUGEPAGE (preloaded only; huge pages for its Vecs).
set -euo pipefail

BACKEND="${1:-mmap}"
REPO="$(git rev-parse --show-toplevel)"
IDX="${MLP_INDEX:-$REPO/uniprot-2025-04/suffix-array}"
PEP="${MLP_PEP:-$REPO/uniprot-2025-04/peptides/peptides_5_50.txt}"
RUNS="${MLP_RUNS:-100}"
AMT="${MLP_AMT:-10000}"
read -r -a BATCHES <<< "${MLP_BATCHES:-1 4 8 16 32 64 128}"
WORK="${MLP_WORK:-/tmp/mlp-sweep}"
KMER="${MLP_KMER:-}"
KMER_ARG=(); [ -n "$KMER" ] && KMER_ARG=(--build-kmer-table "$KMER")

[ -f "$IDX/sa.bin" ] || { echo "ERROR: no sa.bin in $IDX (set MLP_INDEX)"; exit 1; }
[ -f "$PEP" ]        || { echo "ERROR: peptide file not found: $PEP (set MLP_PEP)"; exit 1; }

# Tag the results dir by config so k-mer/huge-page variants don't collide with the
# per-B result cache (files are labeled b_<B>).
TAG="$BACKEND"; [ -n "$KMER" ] && TAG="$TAG-k$KMER"; [ "${MLP_HUGEPAGE:-}" = 1 ] && TAG="$TAG-hp"
OUT="$WORK/$TAG"; mkdir -p "$OUT"
if [ "$BACKEND" = mmap ]; then FEAT="--features mmap"; WARMUP="all:$((AMT*10))"; else FEAT=""; WARMUP="$((AMT*10))"; fi
echo "== Config: backend=$BACKEND kmer=${KMER:-off} hugepage=${MLP_HUGEPAGE:-off} -> $OUT =="

echo "== Build $BACKEND =="
(cd "$REPO" && cargo build --release -q -p sa-benchmarks --no-default-features $FEAT)
BIN="$REPO/target/release/sa-benchmarks"

for b in "${BATCHES[@]}"; do
  lbl="b_$b"
  [ -s "$OUT/$lbl.jsonl" ] && { echo "  skip B=$b (exists in $OUT)"; continue; }
  echo "[$(date +%H:%M:%S)] B=$b"
  if [ "$b" = 1 ]; then unset SA_MLP_BATCH; else export SA_MLP_BATCH="$b"; fi
  if [ "${MLP_HUGEPAGE:-}" = 1 ]; then export SA_MADV_HUGEPAGE=1; else unset SA_MADV_HUGEPAGE; fi
  "$BIN" --index-dir "$IDX" --output "$OUT" --label "$lbl" \
    --peptide-file "$PEP" --amount-of-peptides "$AMT" --runs "$RUNS" --warmup "$WARMUP" \
    "${KMER_ARG[@]}" >/dev/null 2>&1
  unset SA_MLP_BATCH SA_MADV_HUGEPAGE
done

echo "== Results ($BACKEND) =="
python3 - "$OUT" "${BATCHES[@]}" <<'PY'
import json, os, sys, statistics as st
OUT=sys.argv[1]; batches=[int(x) for x in sys.argv[2:]]
def med(b):
    p=os.path.join(OUT,f"b_{b}.jsonl")
    return st.median(json.loads(l)["result"]["throughput_qps"] for l in open(p) if l.strip()) if os.path.exists(p) else None
base=med(1)
print(f"\n{'B':>5s} {'median qps':>13s} {'vs B=1':>9s}")
for b in batches:
    m=med(b)
    if m is None: continue
    d=f"{(m-base)/base*100:+.0f}%" if (base and b!=1) else ""
    print(f"{b:>5d} {m:>13,.0f} {d:>9s}")
print(f"\n(raw jsonl in {OUT})")
PY
echo "== DONE =="
