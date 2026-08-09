#!/usr/bin/env bash
# Master benchmark matrix: compare the preloaded and mmap backends across the small / medium /
# large peptide files. Uses the benchmark's in-process matrix mode: ONE process per backend
# loads the index once (index loads are extremely expensive — 37e9 suffixes on the full DB)
# and sweeps every config in-process, so the whole run is just 2 index loads total, not one
# per config.
#
# The grid (see `expand_cells` in main.rs for the exact cell list) is
# kmer {none,5-mer} x equate_il x tryptic x mlp_batch {1,16}, with tryptic collapsed to one
# representative cell on the small/medium files (all 30 small/tryptic cells in the last full
# run landed at 653-684 qps — a flat line, not worth a grid). 34 configs/backend by default
# (9 small + 9 medium + 16 large); 50 with MATRIX_KMER6=1.
#
# The one-factor-at-a-time knob sweep and the combined-tuning confirm sweep this script used
# to drive are gone: they had settled which SearchTuning knobs matter, and the two knobs they
# measured dead (retrieval_batch, +1.7% median; scalar_kmer_prefetch, +0.3% — both inside the
# 3.9% noise floor) have since been deleted from the searcher. Set MATRIX_TUNING_* below to
# put the remaining three knobs under test across the same grid.
#
# Usage (from the repo root):  bash sa-benchmarks/matrix_bench.sh
#
# Env overrides:
#   INDEX_DIR         dataset root holding suffix-array/ and peptides/ (default
#                     $REPO/uniprot-2025-04 — set this, the checked-in default is just a
#                     placeholder and no dataset is committed to the repo)
#   MATRIX_INDEX      index dir (sa.bin/proteins.bin/mapping.bin); default $INDEX_DIR/suffix-array
#   MATRIX_PEP_DIR    dir holding <file>.txt for each file in MATRIX_FILES;
#                     default $INDEX_DIR/peptides
#   MATRIX_FILES      peptide files, no extension (default "small medium large") — stems must
#                     be small/medium/large for the grid's tryptic-collapse rule to apply
#   MATRIX_BACKENDS   default "preloaded mmap"
#   MATRIX_RUNS       timed runs per config (default 20 — matches the measured noise floor of
#                     p90 = 3.9%; don't lower this without re-checking that floor)
#   MATRIX_AMT        peptides per run (default 10000)
#   MATRIX_BATCHES    comma-separated MLP batch sizes for the grid, 1=scalar (default "1,16")
#   MATRIX_KMER6      1 = include the 6-mer table in the grid (default 0 — costs 3.06 GB vs
#                     127 MB for a sub-noise-floor difference; also gates whether the 6-mer
#                     table is built/loaded at all)
#   MATRIX_METRICS    1 = build with the `metrics` feature for this run (adds per-candidate
#                     counters and internal timing breakdowns). Default 0: timing runs build
#                     WITHOUT `metrics`, because the instrumentation perturbs what it measures
#                     (extra atomics on the hot path). Use MATRIX_METRICS=1 for a separate
#                     diagnostic run when you need the examined/accepted/accept-rate numbers
#                     (e.g. to settle whether tryptic's slowdown is a low acceptance rate or
#                     exhaustive scanning) — don't read its qps numbers as the real throughput.
#   MATRIX_TUNING_VALIDATE_BATCH / _VALIDATE_PREFETCH_THRESHOLD / _RETRIEVAL_PREFETCH_DISTANCE
#                     SearchTuning knobs applied to every cell. Unset knobs fall back to
#                     SearchTuning::default(), which is what the reference numbers were taken
#                     at — so leave them unset unless you are testing a candidate tuning.
#   MATRIX_KMER5_FILE / MATRIX_KMER6_FILE   pre-built table files; if absent they are built
#                     in-process (once each)
#   MATRIX_WORK       scratch + results dir (default /tmp/matrix-bench)
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
DATASET="${INDEX_DIR:-$REPO/uniprot-2025-04}"
IDX="${MATRIX_INDEX:-$DATASET/suffix-array}"
PEP_DIR="${MATRIX_PEP_DIR:-$DATASET/peptides}"
read -r -a FILES    <<< "${MATRIX_FILES:-small medium large}"
read -r -a BACKENDS <<< "${MATRIX_BACKENDS:-preloaded mmap}"
RUNS="${MATRIX_RUNS:-20}"
AMT="${MATRIX_AMT:-10000}"
BATCHES="${MATRIX_BATCHES:-1,16}"
KMER6="${MATRIX_KMER6:-0}"
METRICS="${MATRIX_METRICS:-0}"
KMER5_FILE="${MATRIX_KMER5_FILE:-$IDX/kmer-tables/5mer_table.bin}"
KMER6_FILE="${MATRIX_KMER6_FILE:-$IDX/kmer-tables/6mer_table.bin}"
WORK="${MATRIX_WORK:-/tmp/matrix-bench}"
BIN_DIR="$WORK/bin"; OUT="$WORK/results"
mkdir -p "$BIN_DIR" "$OUT"

[ -f "$IDX/sa.bin" ] || { echo "ERROR: no sa.bin in $IDX (set INDEX_DIR or MATRIX_INDEX)"; exit 1; }

# Comma-separated matrix peptide files.
csv=""
for f in "${FILES[@]}"; do
  p="$PEP_DIR/$f.txt"
  [ -f "$p" ] || { echo "ERROR: missing peptide file $p (set INDEX_DIR or MATRIX_PEP_DIR)"; exit 1; }
  csv+="$p,"
done
csv="${csv%,}"

# Optional pre-built k-mer tables (built in-process if absent). --kmer6-file is harmless to
# pass even when the 6-mer cells aren't selected: the binary only loads/builds table6 at all
# when --matrix-kmer6 is also set.
KARGS=()
[ -f "$KMER5_FILE" ] && KARGS+=(--kmer5-file "$KMER5_FILE")
[ -f "$KMER6_FILE" ] && KARGS+=(--kmer6-file "$KMER6_FILE")
[ "$KMER6" = "1" ] && KARGS+=(--matrix-kmer6)

# SearchTuning under test; any knob left unset falls back to SearchTuning::default().
TUNING_ARGS=()
[ -n "${MATRIX_TUNING_VALIDATE_BATCH:-}" ] && TUNING_ARGS+=(--validate-batch "$MATRIX_TUNING_VALIDATE_BATCH")
[ -n "${MATRIX_TUNING_VALIDATE_PREFETCH_THRESHOLD:-}" ] && TUNING_ARGS+=(--validate-prefetch-threshold "$MATRIX_TUNING_VALIDATE_PREFETCH_THRESHOLD")
[ -n "${MATRIX_TUNING_RETRIEVAL_PREFETCH_DISTANCE:-}" ] && TUNING_ARGS+=(--retrieval-prefetch-distance "$MATRIX_TUNING_RETRIEVAL_PREFETCH_DISTANCE")

for be in "${BACKENDS[@]}"; do
  feat=""
  if [ "$be" = mmap ] && [ "$METRICS" = "1" ]; then feat="--features mmap,metrics"
  elif [ "$be" = mmap ]; then feat="--features mmap"
  elif [ "$METRICS" = "1" ]; then feat="--features metrics"
  fi
  echo "== build $be (features: ${feat:-none}) =="
  (cd "$REPO" && cargo build --release -q -p sa-benchmarks --no-default-features $feat)
  bin="$REPO/target/release/sa-benchmarks"; cp "$bin" "$BIN_DIR/$be"

  # Ask the binary itself for the config count (--dry-run needs no index and can't drift from
  # the actual grid generator, unlike a hand-maintained arithmetic formula here).
  expected=$("$BIN_DIR/$be" --matrix --dry-run --index-dir "$IDX" --output "$OUT" --matrix-files "$csv" \
    --matrix-batches "$BATCHES" --runs "$RUNS" "${KARGS[@]}" "${TUNING_ARGS[@]}" \
    | grep -oE 'TOTAL this backend: [0-9]+' | grep -oE '[0-9]+')

  # Resumable at backend granularity: skip if this backend already has a full result set.
  if [ -f "$OUT/$be.jsonl" ] && [ "$(wc -l < "$OUT/$be.jsonl")" -ge "$expected" ]; then
    echo "== $be already complete ($expected configs), skipping =="; continue
  fi
  echo "== run $be matrix ($expected configs x $RUNS runs) =="
  "$BIN_DIR/$be" --matrix --index-dir "$IDX" --matrix-files "$csv" "${KARGS[@]}" "${TUNING_ARGS[@]}" \
    --matrix-batches "$BATCHES" --amount-of-peptides "$AMT" --runs "$RUNS" --max-matches 10000 \
    --output "$OUT" --label "$be"
done

echo ""
echo "== Results =="
python3 - "$OUT" "$RUNS" <<'PY'
import json, os, sys, glob
from collections import defaultdict
OUT, RUNS = sys.argv[1], sys.argv[2]

# The measured run-to-run noise floor (see matrix_bench.sh header). Any delta smaller than
# this — or smaller than the wider of the two cells' own measured bands — is noise, not
# signal, and must be flagged rather than read as a real effect.
NOISE_FLOOR_PCT = 3.9

def band_pct(s):  # half the p10..p90 spread as a percent of the median
    return (s["p90"] - s["p10"]) / 2 / s["p50"] * 100 if s and s["p50"] else 0.0

records = []
for path in glob.glob(os.path.join(OUT, "*.jsonl")):
    backend = os.path.basename(path)[:-6]
    for line in open(path):
        if not line.strip():
            continue
        r = json.loads(line)
        r["_backend"] = backend
        records.append(r)

def stats_of(r):
    s = r.get("stats")
    if s:
        return {"p50": s["qps_p50"], "p10": s["qps_p10"], "p90": s["qps_p90"]}
    v = [r["result"]["throughput_qps"]]
    return {"p50": v[0], "p10": v[0], "p90": v[0]}

# ---------------------------------------------------------------------------
# preloaded vs mmap per config (median qps)
# ---------------------------------------------------------------------------
grid = [r for r in records if r["config"].get("phase") == "grid"]
if grid:
    data = defaultdict(dict)
    for r in grid:
        c = r["config"]
        key = (c["peptide_source"], c["equate_il"], c["tryptic"], c["batch_size"],
               {0: "none", 5: "5-mer", 6: "6-mer"}.get(c["kmer_k"], str(c["kmer_k"])))
        data[key][r["_backend"]] = stats_of(r)

    print("\n### grid: preloaded vs mmap per config (median qps) ###")
    for file in sorted({k[0] for k in data}):
        print(f"\n-- {file} --")
        print(f"{'equate_il':<9} {'tryptic':<7} {'batch':<7} {'kmer':<6} {'preloaded':>12} {'mmap':>12} {'mmap/pre':>9} {'noise':>7}")
        for k in sorted((k for k in data if k[0] == file), key=lambda k: (k[1], k[2], k[3], k[4])):
            d = data[k]; sp, sm = d.get("preloaded"), d.get("mmap")
            pre = sp["p50"] if sp else None; mm = sm["p50"] if sm else None
            ratio = f"{mm/pre:.2f}x" if (pre and mm) else ""
            ps = f"{pre:,.0f}" if pre else "-"; ms = f"{mm:,.0f}" if mm else "-"
            nz = max(band_pct(sp), band_pct(sm))
            noise = f"±{nz:.1f}%" if (sp or sm) else ""
            batch = "scalar" if k[3] == 1 else str(k[3])
            print(f"{str(k[1]):<9} {str(k[2]):<7} {batch:<7} {k[4]:<6} {ps:>12} {ms:>12} {ratio:>9} {noise:>7}")

print(f"\nmedian of {RUNS} reps; noise = half the p10..p90 spread. "
      f"deltas below {NOISE_FLOOR_PCT}% (or below a cell's own band) are noise. raw jsonl in {OUT}")
PY
echo "== DONE =="
