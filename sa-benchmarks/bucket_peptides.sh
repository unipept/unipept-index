#!/usr/bin/env bash
# Split a peptide file into length-bucketed files for per-regime benchmarking
# (short peptides are retrieval-bound; long peptides are search-bound).
#
# Usage:
#   bash sa-benchmarks/bucket_peptides.sh <source.txt> <outdir> ["lo-hi lo-hi ..."]
# Default buckets: 5-7 15-25 35-50
#
# Prints how many peptides landed in each bucket so you can size AMT*RUNS accordingly.
set -euo pipefail

SRC="${1:?usage: bucket_peptides.sh <source.txt> <outdir> [\"lo-hi ...\"]}"
OUTDIR="${2:?output dir required}"
read -r -a BUCKETS <<< "${3:-5-7 15-25 35-50}"

[ -f "$SRC" ] || { echo "ERROR: source not found: $SRC"; exit 1; }
mkdir -p "$OUTDIR"

echo "Source: $SRC ($(wc -l < "$SRC") peptides)"
printf '%-26s %12s\n' "bucket file" "peptides"
for b in "${BUCKETS[@]}"; do
  lo="${b%-*}"; hi="${b#*-}"
  out="$OUTDIR/peptides_${lo}_${hi}.txt"
  awk -v lo="$lo" -v hi="$hi" 'length($0)>=lo && length($0)<=hi' "$SRC" > "$out"
  printf '%-26s %12d\n' "$(basename "$out")" "$(wc -l < "$out")"
done
echo "Wrote to $OUTDIR/"
