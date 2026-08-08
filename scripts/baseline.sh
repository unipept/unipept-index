#!/usr/bin/env bash
set -euo pipefail

# Captures (or re-captures) the behavioural golden set used to prove that a
# refactor changed nothing observable.
#
# Sweeps the full config matrix -- {preloaded, mmap} x {kmer off, kmer on}
# x {equate_il} x {tryptic} -- and writes one normalised JSON file per cell.
#
# Usage:
#   scripts/baseline.sh <index-dir> <peptides-file> <output-dir> [kmer-table]
#
# Compare two capture directories with:
#   scripts/baseline.sh --compare <dir-a> <dir-b>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

normalize() {
    jq --sort-keys '[.[] | .proteins |= sort_by(.uniprot_accession)] | sort_by(.sequence)' "$1"
}

if [ "${1:-}" = "--compare" ]; then
    DIR_A="$2"
    DIR_B="$3"
    status=0
    for f in "$DIR_A"/*.json; do
        name="$(basename "$f")"
        if [ ! -f "$DIR_B/$name" ]; then
            echo "MISSING in $DIR_B: $name"
            status=1
            continue
        fi
        if diff -q "$f" "$DIR_B/$name" > /dev/null 2>&1; then
            echo "  ok    $name"
        else
            echo "  DIFF  $name"
            diff <(normalize "$f") <(normalize "$DIR_B/$name") | head -40
            status=1
        fi
    done
    [ $status -eq 0 ] && echo "All cells identical."
    exit $status
fi

INDEX_DIR="$(cd "$1" && pwd)"
PEPTIDES="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
OUT_DIR="$3"
KMER_TABLE_PATH="${4:-}"

mkdir -p "$OUT_DIR"

for features in "" "mmap"; do
    backend_name="preloaded"
    [ -n "$features" ] && backend_name="mmap"

    kmer_settings=("")
    [ -n "$KMER_TABLE_PATH" ] && kmer_settings=("" "$KMER_TABLE_PATH")

    for kmer in "${kmer_settings[@]}"; do
        kmer_name="nokmer"
        [ -n "$kmer" ] && kmer_name="kmer"

        for equate_il in false true; do
            for tryptic in false true; do
                cell="${backend_name}-${kmer_name}-il${equate_il}-tryp${tryptic}"
                raw="$(mktemp)"
                echo "=== $cell ==="
                SA_FEATURES="$features" \
                INDEX_DIR="$INDEX_DIR" \
                KMER_TABLE="$kmer" \
                PEPTIDES_FILE="$PEPTIDES" \
                EQUATE_IL="$equate_il" \
                TRYPTIC="$tryptic" \
                    "$SCRIPT_DIR/run_index.sh" "$raw" > /dev/null
                normalize "$raw" > "$OUT_DIR/$cell.json"
                rm -f "$raw"
                echo "    -> $OUT_DIR/$cell.json ($(jq 'length' "$OUT_DIR/$cell.json") results)"
            done
        done
    done
done

echo
echo "Captured $(ls -1 "$OUT_DIR"/*.json | wc -l | tr -d ' ') cells into $OUT_DIR"
