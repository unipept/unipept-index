#!/usr/bin/env bash
# Master benchmark matrix: compare the preloaded and mmap backends across three phases, for
# the small / medium / large peptide files. Uses the benchmark's in-process matrix mode: ONE
# process per backend loads the index once (index loads are extremely expensive — 37e9
# suffixes on the full DB) and sweeps every phase/config in-process, so the whole run is just
# 2 index loads total, not one per config.
#
# Phases (see --matrix-phases in `sa-benchmarks --help` and grid_for_file / build_ofat_cells
# in main.rs for the exact cell lists):
#
#   grid    - Task 1 trimmed default grid: kmer {none,5-mer} x equate_il x tryptic x
#             mlp_batch {1,16}, with tryptic collapsed to one representative cell on the
#             small/medium files (all 30 small/tryptic cells in the last full run landed at
#             653-684 qps — a flat line, not worth a grid). 34 configs/backend by default
#             (9 small + 9 medium + 16 large); 50 with MATRIX_KMER6=1.
#   ofat    - Task 2 one-factor-at-a-time sweep of the 5 SearchTuning knobs (validate_batch,
#             validate_prefetch_threshold, retrieval_prefetch_distance, retrieval_batch,
#             scalar_kmer_prefetch) around two fixed baselines instead of the 46,080-config
#             full cross-product. 26 configs/(file,backend), 78/backend.
#   confirm - Task 3: the MATRIX_CONFIRM_* combo tuning (all 5 knobs at once — the "best of
#             each OFAT knob" candidate) run across the same grid as "grid", to check whether
#             the knobs are separable (validate_batch and retrieval_prefetch_distance both
#             consume line-fill buffers, so they plausibly interact). Opt-in — only useful
#             once OFAT has identified candidate values.
#
# Default MATRIX_PHASES is "grid,ofat" (comparable in size to the old undifferentiated
# 180-config grid, but pointed at what's actually worth measuring). Add "confirm" once you
# have OFAT results to combine.
#
# Usage (from the repo root):  bash sa-benchmarks/matrix_bench.sh
#
# Env overrides:
#   MATRIX_INDEX      index dir (sa.bin/proteins.bin/mapping.bin)
#   MATRIX_PEP_DIR    dir holding <file>.txt for each file in MATRIX_FILES
#   MATRIX_FILES      peptide files, no extension (default "small medium large") — stems must
#                     be small/medium/large for the grid's tryptic-collapse rule to apply
#   MATRIX_BACKENDS   default "preloaded mmap"
#   MATRIX_RUNS       timed runs per config (default 20 — matches the measured noise floor of
#                     p90 = 3.9%; don't lower this without re-checking that floor)
#   MATRIX_AMT        peptides per run (default 10000)
#   MATRIX_BATCHES    comma-separated MLP batch sizes for the "grid"/"confirm" phases,
#                     1=scalar (default "1,16")
#   MATRIX_PHASES     comma-separated phases to run (default "grid,ofat")
#   MATRIX_KMER6      1 = include the 6-mer table in "grid"/"confirm" (default 0 — costs
#                     3.06 GB vs 127 MB for a sub-noise-floor difference; also gates whether
#                     the 6-mer table is built/loaded at all)
#   MATRIX_METRICS    1 = build with the `metrics` feature for this run (adds per-candidate
#                     counters and internal timing breakdowns). Default 0: timing runs build
#                     WITHOUT `metrics`, because the instrumentation perturbs what it measures
#                     (extra atomics on the hot path). Use MATRIX_METRICS=1 for a separate
#                     diagnostic run when you need the examined/accepted/accept-rate numbers
#                     (e.g. to settle whether tryptic's slowdown is a low acceptance rate or
#                     exhaustive scanning) — don't read its qps numbers as the real throughput.
#   MATRIX_CONFIRM_VALIDATE_BATCH / _VALIDATE_PREFETCH_THRESHOLD / _RETRIEVAL_PREFETCH_DISTANCE
#   MATRIX_CONFIRM_RETRIEVAL_BATCH / _SCALAR_KMER_PREFETCH
#                     Only used when MATRIX_PHASES includes "confirm" — the combined tuning to
#                     test. Unset knobs fall back to SearchTuning::default() (a no-op combo),
#                     so set the ones you actually want to combine.
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
PHASES="${MATRIX_PHASES:-grid,ofat}"
KMER6="${MATRIX_KMER6:-0}"
METRICS="${MATRIX_METRICS:-0}"
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

# Optional pre-built k-mer tables (built in-process if absent). --kmer6-file is harmless to
# pass even when the 6-mer phase isn't selected: the binary only loads/builds table6 at all
# when --matrix-kmer6 is also set.
KARGS=()
[ -f "$KMER5_FILE" ] && KARGS+=(--kmer5-file "$KMER5_FILE")
[ -f "$KMER6_FILE" ] && KARGS+=(--kmer6-file "$KMER6_FILE")
[ "$KMER6" = "1" ] && KARGS+=(--matrix-kmer6)

# Task 3 confirm-phase combo tuning — only meaningful (and only passed) when "confirm" is
# selected; any knob left unset falls back to SearchTuning::default() in the binary.
CONFIRM_ARGS=()
if [[ ",$PHASES," == *",confirm,"* ]]; then
  [ -n "${MATRIX_CONFIRM_VALIDATE_BATCH:-}" ] && CONFIRM_ARGS+=(--validate-batch "$MATRIX_CONFIRM_VALIDATE_BATCH")
  [ -n "${MATRIX_CONFIRM_VALIDATE_PREFETCH_THRESHOLD:-}" ] && CONFIRM_ARGS+=(--validate-prefetch-threshold "$MATRIX_CONFIRM_VALIDATE_PREFETCH_THRESHOLD")
  [ -n "${MATRIX_CONFIRM_RETRIEVAL_PREFETCH_DISTANCE:-}" ] && CONFIRM_ARGS+=(--retrieval-prefetch-distance "$MATRIX_CONFIRM_RETRIEVAL_PREFETCH_DISTANCE")
  [ -n "${MATRIX_CONFIRM_RETRIEVAL_BATCH:-}" ] && CONFIRM_ARGS+=(--retrieval-batch "$MATRIX_CONFIRM_RETRIEVAL_BATCH")
  [ -n "${MATRIX_CONFIRM_SCALAR_KMER_PREFETCH:-}" ] && CONFIRM_ARGS+=(--scalar-kmer-prefetch "$MATRIX_CONFIRM_SCALAR_KMER_PREFETCH")
fi

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
  # the actual grid/ofat/confirm generators, unlike a hand-maintained arithmetic formula here).
  expected=$("$BIN_DIR/$be" --matrix --dry-run --index-dir "$IDX" --output "$OUT" --matrix-files "$csv" \
    --matrix-batches "$BATCHES" --matrix-phases "$PHASES" --runs "$RUNS" "${KARGS[@]}" "${CONFIRM_ARGS[@]}" \
    | grep -oE 'TOTAL this backend: [0-9]+' | grep -oE '[0-9]+')

  # Resumable at backend granularity: skip if this backend already has a full result set.
  if [ -f "$OUT/$be.jsonl" ] && [ "$(wc -l < "$OUT/$be.jsonl")" -ge "$expected" ]; then
    echo "== $be already complete ($expected configs), skipping =="; continue
  fi
  echo "== run $be matrix ($expected configs x $RUNS runs, phases: $PHASES) =="
  "$BIN_DIR/$be" --matrix --index-dir "$IDX" --matrix-files "$csv" "${KARGS[@]}" "${CONFIRM_ARGS[@]}" \
    --matrix-batches "$BATCHES" --matrix-phases "$PHASES" --amount-of-peptides "$AMT" --runs "$RUNS" --max-matches 10000 \
    --output "$OUT" --label "$be"
done

echo ""
echo "== Results =="
python3 - "$OUT" "$RUNS" <<'PY'
import json, os, sys, statistics as st, glob
from collections import defaultdict
OUT, RUNS = sys.argv[1], sys.argv[2]

# The measured run-to-run noise floor (see matrix_bench.sh header / Task 4). Any delta smaller
# than this — or smaller than the wider of the two cells' own measured bands — is noise, not
# signal, and must be flagged rather than read as a knob effect.
NOISE_FLOOR_PCT = 3.9

def pct(sorted_vals, p):
    n = len(sorted_vals)
    if n == 0: return None
    if n == 1: return sorted_vals[0]
    r = p * (n - 1); lo = int(r); frac = r - lo
    return sorted_vals[lo] + (sorted_vals[min(lo + 1, n - 1)] - sorted_vals[lo]) * frac

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
# Phase "grid": preloaded vs mmap per config (median qps)
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

# ---------------------------------------------------------------------------
# Phase "ofat": one knob at a time vs its baseline, with an explicit noise flag
# ---------------------------------------------------------------------------
ofat = [r for r in records if r["config"].get("phase") == "ofat"]
if ofat:
    print("\n\n### ofat: knob sweep (median qps); '(noise)' = delta smaller than the wider measured "
          f"band or the {NOISE_FLOOR_PCT}% guardrail floor ###")
    knob_field = {
        "validate_batch": "validate_batch",
        "validate_prefetch_threshold": "validate_prefetch_threshold",
        "retrieval_prefetch_distance": "retrieval_prefetch_distance",
        "retrieval_batch": "retrieval_batch",
        "scalar_kmer_prefetch": "scalar_kmer_prefetch",
    }
    groups = defaultdict(list)
    for r in ofat:
        c = r["config"]
        key = (c["peptide_source"], r["_backend"], c["ofat_baseline"], c["ofat_knob"])
        groups[key].append(r)

    for file in sorted({k[0] for k in groups}):
        print(f"\n-- {file} --")
        for backend in sorted({k[1] for k in groups if k[0] == file}):
            for baseline in sorted({k[2] for k in groups if k[0] == file and k[1] == backend}):
                knobs = sorted({k[3] for k in groups if k[0] == file and k[1] == backend and k[2] == baseline})
                for knob in knobs:
                    rows = groups[(file, backend, baseline, knob)]
                    rows.sort(key=lambda r: r["config"][knob_field[knob]])
                    print(f"  {backend:<10} {baseline} {knob}")
                    prev = None
                    for r in rows:
                        s = stats_of(r)
                        v = r["config"][knob_field[knob]]
                        band = band_pct(s)
                        line = f"    {knob}={v!s:<6} -> {s['p50']:>10,.0f} qps  (±{band:.1f}%)"
                        if prev is not None:
                            prev_val, prev_s = prev
                            delta = (s["p50"] - prev_s["p50"]) / prev_s["p50"] * 100 if prev_s["p50"] else 0.0
                            threshold = max(band, band_pct(prev_s), NOISE_FLOOR_PCT)
                            flag = "  <- within noise" if abs(delta) < threshold else ""
                            line += f"   [{delta:+.1f}% vs {knob}={prev_val!s}{flag}]"
                        print(line)
                        prev = (v, s)

# ---------------------------------------------------------------------------
# Phase "confirm": combo tuning vs the same cell under grid's default tuning
# ---------------------------------------------------------------------------
confirm = [r for r in records if r["config"].get("phase") == "confirm"]
if confirm:
    print("\n\n### confirm: combo tuning vs default-tuning grid cell (median qps) ###")
    grid_by_cell = {}
    for r in grid:
        c = r["config"]
        key = (c["peptide_source"], r["_backend"], c["equate_il"], c["tryptic"], c["batch_size"], c["kmer_k"])
        grid_by_cell[key] = stats_of(r)

    for r in confirm:
        c = r["config"]
        key = (c["peptide_source"], r["_backend"], c["equate_il"], c["tryptic"], c["batch_size"], c["kmer_k"])
        base = grid_by_cell.get(key)
        s = stats_of(r)
        if base:
            delta = (s["p50"] - base["p50"]) / base["p50"] * 100 if base["p50"] else 0.0
            threshold = max(band_pct(s), band_pct(base), NOISE_FLOOR_PCT)
            flag = "  <- within noise" if abs(delta) < threshold else ""
            print(f"  {c['peptide_source']:<8} {r['_backend']:<10} il={c['equate_il']!s:<5} tr={c['tryptic']!s:<5} "
                  f"batch={c['batch_size']:<3} kmer={c['kmer_k']}  confirm {s['p50']:>10,.0f} qps  vs default {base['p50']:>10,.0f} qps  "
                  f"[{delta:+.1f}%{flag}]")
        else:
            print(f"  {c['peptide_source']:<8} {r['_backend']:<10} il={c['equate_il']!s:<5} tr={c['tryptic']!s:<5} "
                  f"batch={c['batch_size']:<3} kmer={c['kmer_k']}  confirm {s['p50']:>10,.0f} qps  (no matching grid cell — run phase 'grid' too)")

print(f"\nmedian of {RUNS} reps; noise = half the p10..p90 spread. raw jsonl in {OUT}")
PY
echo "== DONE =="
