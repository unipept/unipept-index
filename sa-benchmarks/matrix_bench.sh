#!/usr/bin/env bash
# Master benchmark matrix: compare the preloaded and mmap backends across the full
# parameter grid, for the small / medium / large peptide files.
#
#   equate_il {true,false} x tryptic {true,false}   (4)
#   x searcher {scalar, batched}                    (2)
#   x k-mer {none, 5-mer, 6-mer}                    (3)   = 24 configs
#
# 24 configs x {preloaded, mmap} x {small, medium, large}. max_matches is fixed at 10000.
#
# Uses the benchmark's in-process matrix mode: ONE process per backend loads the index
# once and sweeps all configs (the 6-mer table is built/loaded once, not per config). So
# the whole run is just 2 index loads instead of 144.
#
# Usage (from the repo root):  bash sa-benchmarks/matrix_bench.sh
#
# Env overrides:
#   MATRIX_INDEX      index dir (sa.bin/proteins.bin/mapping.bin)
#   MATRIX_PEP_DIR    dir holding <file>.txt for each file in MATRIX_FILES
#   MATRIX_FILES      peptide files, no extension (default "small medium large")
#   MATRIX_BACKENDS   default "preloaded mmap"
#   MATRIX_RUNS       timed runs per config (default 20)
#   MATRIX_AMT        peptides per run (default 10000)
#   MATRIX_BATCHES    comma-separated MLP batch sizes to sweep, 1=scalar (default "1,16")
#   MATRIX_KMER5_FILE / MATRIX_KMER6_FILE   pre-built table files; if absent they are built
#                     in-process (once each)
#   MATRIX_WORK       scratch + results dir (default /tmp/matrix-bench)
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
IDX="${MATRIX_INDEX:-$REPO/uniprot-2025-04/suffix-array}"
PEP_DIR="${MATRIX_PEP_DIR:-$REPO/uniprot-2025-04/peptides}"
read -r -a FILES    <<< "${MATRIX_FILES:-small medium large}"
read -r -a BACKENDS <<< "${MATRIX_BACKENDS:-preloaded mmap}"
RUNS="${MATRIX_RUNS:-20}"
AMT="${MATRIX_AMT:-10000}"
BATCHES="${MATRIX_BATCHES:-1,16}"
IFS=',' read -r -a _BATCH_ARR <<< "$BATCHES"; NBATCH=${#_BATCH_ARR[@]}
KMER5_FILE="${MATRIX_KMER5_FILE:-$IDX/kmer-tables/5mer_table.bin}"
KMER6_FILE="${MATRIX_KMER6_FILE:-$IDX/kmer-tables/6mer_table.bin}"
WORK="${MATRIX_WORK:-/tmp/matrix-bench}"
BIN_DIR="$WORK/bin"; OUT="$WORK/results"
mkdir -p "$BIN_DIR" "$OUT"

[ -f "$IDX/sa.bin" ] || { echo "ERROR: no sa.bin in $IDX (set MATRIX_INDEX)"; exit 1; }

# Comma-separated matrix peptide files.
csv=""
for f in "${FILES[@]}"; do
  p="$PEP_DIR/$f.txt"
  [ -f "$p" ] || { echo "ERROR: missing peptide file $p"; exit 1; }
  csv+="$p,"
done
csv="${csv%,}"

# Optional pre-built k-mer tables (built in-process if absent).
KARGS=()
[ -f "$KMER5_FILE" ] && KARGS+=(--kmer5-file "$KMER5_FILE")
[ -f "$KMER6_FILE" ] && KARGS+=(--kmer6-file "$KMER6_FILE")

expected=$(( ${#FILES[@]} * 4 * NBATCH * 3 ))   # one aggregated line per config: files x (eq*tr) x batches x kmer
for be in "${BACKENDS[@]}"; do
  feat=""; [ "$be" = mmap ] && feat="--features mmap"
  echo "== build $be =="
  (cd "$REPO" && cargo build --release -q -p sa-benchmarks --no-default-features $feat)
  bin="$REPO/target/release/sa-benchmarks"; cp "$bin" "$BIN_DIR/$be"

  # Resumable at backend granularity: skip if this backend already has a full result set.
  if [ -f "$OUT/$be.jsonl" ] && [ "$(wc -l < "$OUT/$be.jsonl")" -ge "$expected" ]; then
    echo "== $be already complete, skipping =="; continue
  fi
  echo "== run $be matrix ($((${#FILES[@]} * 4 * NBATCH * 3)) configs x $RUNS runs) =="
  "$BIN_DIR/$be" --matrix --index-dir "$IDX" --matrix-files "$csv" "${KARGS[@]}" \
    --matrix-batches "$BATCHES" --amount-of-peptides "$AMT" --runs "$RUNS" --max-matches 10000 \
    --output "$OUT" --label "$be"
done

echo ""
echo "== Results: preloaded vs mmap per config (median qps) =="
python3 - "$OUT" <<'PY'
import json, os, sys, statistics as st, glob
OUT = sys.argv[1]

def pct(sorted_vals, p):
    n = len(sorted_vals)
    if n == 0: return None
    if n == 1: return sorted_vals[0]
    r = p * (n - 1); lo = int(r); frac = r - lo
    return sorted_vals[lo] + (sorted_vals[min(lo + 1, n - 1)] - sorted_vals[lo]) * frac

# data[(file, equate_il, tryptic, batch, kmer)][backend] = {"p50","p10","p90"}
# v2 records carry stats (one line per config); older per-run records fall back to a list.
data, raw = {}, {}
for path in glob.glob(os.path.join(OUT, "*.jsonl")):
    backend = os.path.basename(path)[:-6]
    for line in open(path):
        if not line.strip():
            continue
        r = json.loads(line); c = r["config"]
        key = (c["peptide_source"], c["equate_il"], c["tryptic"], c["batch_size"],
               {0: "none", 5: "5-mer", 6: "6-mer"}.get(c["kmer_k"], str(c["kmer_k"])))
        s = r.get("stats")
        if s:
            data.setdefault(key, {})[backend] = {"p50": s["qps_p50"], "p10": s["qps_p10"], "p90": s["qps_p90"]}
        else:
            raw.setdefault(key, {}).setdefault(backend, []).append(r["result"]["throughput_qps"])
# fold any legacy per-run records into the same shape
for key, bes in raw.items():
    for be, vals in bes.items():
        if data.get(key, {}).get(be):
            continue
        v = sorted(vals)
        data.setdefault(key, {})[be] = {"p50": st.median(v), "p10": pct(v, 0.10), "p90": pct(v, 0.90)}

def band(s):   # half the p10..p90 spread as a percent of the median
    return (s["p90"] - s["p10"]) / 2 / s["p50"] * 100 if s and s["p50"] else 0.0

for file in sorted({k[0] for k in data}):
    print(f"\n### {file}")
    print(f"{'equate_il':<9} {'tryptic':<7} {'batch':<7} {'kmer':<6} {'preloaded':>12} {'mmap':>12} {'mmap/pre':>9} {'noise':>7}")
    for k in sorted((k for k in data if k[0] == file), key=lambda k: (k[1], k[2], k[3], k[4])):
        d = data[k]; sp, sm = d.get("preloaded"), d.get("mmap")
        pre = sp["p50"] if sp else None; mm = sm["p50"] if sm else None
        ratio = f"{mm/pre:.2f}x" if (pre and mm) else ""
        ps = f"{pre:,.0f}" if pre else "-"; ms = f"{mm:,.0f}" if mm else "-"
        # widest band across the backends present — flags cells to distrust
        nz = max(band(sp), band(sm))
        noise = f"±{nz:.1f}%" if (sp or sm) else ""
        batch = "scalar" if k[3] == 1 else str(k[3])
        print(f"{str(k[1]):<9} {str(k[2]):<7} {batch:<7} {k[4]:<6} {ps:>12} {ms:>12} {ratio:>9} {noise:>7}")
print(f"\nmedian of {os.getenv('MATRIX_RUNS','?')} reps; noise = half the p10..p90 spread. raw jsonl in {OUT}")
PY
echo "== DONE =="
