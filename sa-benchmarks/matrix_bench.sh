#!/usr/bin/env bash
# Master benchmark matrix: compare the preloaded and mmap backends across the full
# parameter grid, for the small / medium / large peptide files.
#
#   equate_il {true,false} x tryptic {true,false}   (4)
#   x searcher {scalar, batched}                    (2)
#   x k-mer {none, 5-mer, 6-mer}                    (3)   = 24 configs
#
# 24 configs x {preloaded, mmap} x {small, medium, large}. max_matches is fixed at 10000.
# Idempotent: any config whose .jsonl already exists is skipped, so it is safe to resume.
#
# Usage (from the repo root):  bash sa-benchmarks/matrix_bench.sh
#
# Env overrides:
#   MATRIX_INDEX      index dir with sa.bin/proteins.bin/mapping.bin/warmup.txt
#   MATRIX_PEP_DIR    dir holding <file>.txt for each file in MATRIX_FILES
#   MATRIX_FILES      peptide files, no extension (default "small medium large")
#   MATRIX_BACKENDS   default "preloaded mmap"
#   MATRIX_RUNS       timed runs per config (default 20)
#   MATRIX_AMT        peptides per run  (default 10000)
#   MATRIX_BATCH      batch size for the batched searcher (default 16)
#   MATRIX_KMER5_FILE / MATRIX_KMER6_FILE   pre-built table files; if absent the script
#                     falls back to --build-kmer-table 5 / 6 (slower, rebuilt per run)
#   MATRIX_WARM_PRE / MATRIX_WARM_MMAP      warmup args (default "5000" / "all:5000")
#   MATRIX_WORK       scratch + results dir (default /tmp/matrix-bench)
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
IDX="${MATRIX_INDEX:-$REPO/uniprot-2025-04/suffix-array}"
PEP_DIR="${MATRIX_PEP_DIR:-$REPO/uniprot-2025-04/peptides}"
read -r -a FILES    <<< "${MATRIX_FILES:-small medium large}"
read -r -a BACKENDS <<< "${MATRIX_BACKENDS:-preloaded mmap}"
read -r -a KMERS    <<< "${MATRIX_KMERS:-none k5 k6}"   # subset the k-mer dimension if desired
RUNS="${MATRIX_RUNS:-20}"
AMT="${MATRIX_AMT:-10000}"
MAX_MATCHES=10000
BATCH="${MATRIX_BATCH:-16}"
KMER5_FILE="${MATRIX_KMER5_FILE:-$IDX/kmer-tables/5mer_table.bin}"
KMER6_FILE="${MATRIX_KMER6_FILE:-$IDX/kmer-tables/6mer_table.bin}"
WARM_PRE="${MATRIX_WARM_PRE:-5000}"
WARM_MMAP="${MATRIX_WARM_MMAP:-all:5000}"
WORK="${MATRIX_WORK:-/tmp/matrix-bench}"
BIN_DIR="$WORK/bin"; OUT="$WORK/results"
mkdir -p "$BIN_DIR" "$OUT"

[ -f "$IDX/sa.bin" ] || { echo "ERROR: no sa.bin in $IDX (set MATRIX_INDEX)"; exit 1; }

# Build the requested backend binaries once each.
for be in "${BACKENDS[@]}"; do
  feat=""; [ "$be" = mmap ] && feat="--features mmap"
  echo "== build $be =="
  (cd "$REPO" && cargo build --release -q -p sa-benchmarks --no-default-features $feat)
  cp "$REPO/target/release/sa-benchmarks" "$BIN_DIR/$be"
done

# KMER_ARG is set as an array (handles paths with spaces / empty).
set_kmer_arg() {
  KMER_ARG=()
  case "$1" in
    k5) if [ -f "$KMER5_FILE" ]; then KMER_ARG=(--kmer-table-file "$KMER5_FILE"); else KMER_ARG=(--build-kmer-table 5); fi ;;
    k6) if [ -f "$KMER6_FILE" ]; then KMER_ARG=(--kmer-table-file "$KMER6_FILE"); else KMER_ARG=(--build-kmer-table 6); fi ;;
  esac
}

total=$(( ${#BACKENDS[@]} * ${#FILES[@]} * 4 * 2 * ${#KMERS[@]} )); n=0
echo "== running up to $total configs (resumable) =="
for be in "${BACKENDS[@]}"; do
  bin="$BIN_DIR/$be"
  [ "$be" = mmap ] && warm="$WARM_MMAP" || warm="$WARM_PRE"
  for file in "${FILES[@]}"; do
    pep="$PEP_DIR/$file.txt"
    [ -f "$pep" ] || { echo "  MISSING peptide file, skipping: $pep"; continue; }
    for eq in true false; do for tr in true false; do
      for searcher in scalar batched; do
        for kmer in "${KMERS[@]}"; do
          n=$((n+1))
          label="$be~$file~$eq~$tr~$searcher~$kmer"
          [ -s "$OUT/$label.jsonl" ] && { echo "  [$n/$total] skip $label"; continue; }
          set_kmer_arg "$kmer"
          [ "$searcher" = batched ] && export SA_MLP_BATCH="$BATCH" || unset SA_MLP_BATCH
          echo "  [$n/$total] $(date +%H:%M:%S) $label"
          "$bin" --index-dir "$IDX" --output "$OUT" --label "$label" \
            --peptide-file "$pep" --amount-of-peptides "$AMT" --runs "$RUNS" --warmup "$warm" \
            --equate-il "$eq" --tryptic "$tr" --max-matches "$MAX_MATCHES" \
            "${KMER_ARG[@]}" >/dev/null 2>&1
          unset SA_MLP_BATCH
        done
      done
    done; done
  done
done

echo ""
echo "== Results: preloaded vs mmap per config (median qps) =="
python3 - "$OUT" <<'PY'
import json, os, sys, statistics as st, glob
OUT = sys.argv[1]
def med(path):
    v = [json.loads(l)["result"]["throughput_qps"] for l in open(path) if l.strip()]
    return st.median(v) if v else None
data = {}
for p in glob.glob(os.path.join(OUT, "*.jsonl")):
    parts = os.path.basename(p)[:-6].split("~")
    if len(parts) != 6:
        continue
    be, file, eq, tr, searcher, kmer = parts
    data.setdefault((file, eq, tr, searcher, kmer), {})[be] = med(p)
for file in sorted({k[0] for k in data}):
    print(f"\n### {file}")
    print(f"{'equate_il':<9} {'tryptic':<7} {'searcher':<8} {'kmer':<5} {'preloaded':>12} {'mmap':>12} {'mmap/pre':>9}")
    keys = sorted((k for k in data if k[0] == file), key=lambda k: (k[1], k[2], k[3], k[4]))
    for k in keys:
        d = data[k]; pre = d.get("preloaded"); mm = d.get("mmap")
        ratio = f"{mm/pre:.2f}x" if (pre and mm) else ""
        ps = f"{pre:,.0f}" if pre else "-"; ms = f"{mm:,.0f}" if mm else "-"
        print(f"{k[1]:<9} {k[2]:<7} {k[3]:<8} {k[4]:<5} {ps:>12} {ms:>12} {ratio:>9}")
print(f"\n(raw jsonl in {OUT})")
PY
echo "== DONE =="
