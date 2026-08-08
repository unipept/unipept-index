#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROTEINS_TSV="$SCRIPT_DIR/../uniprot-2026-01/suffix-array/proteins.tsv"
OUTPUT_FILE="$SCRIPT_DIR/peptides.txt"

if [ ! -f "$PROTEINS_TSV" ]; then
    echo "Error: proteins.tsv not found at $PROTEINS_TSV" >&2
    exit 1
fi

echo "Generating peptides from $PROTEINS_TSV..."

python3 - "$PROTEINS_TSV" "$OUTPUT_FILE" <<'EOF'
import sys
import random

random.seed(42)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f:
    lines = [line for line in f if len(line.split('\t')) >= 3 and len(line.split('\t')[2]) >= 20]

selected = random.sample(lines, min(20000, len(lines)))

count = 0
with open(output_file, 'w') as out:
    for line in selected:
        sequence = line.rstrip('\n').split('\t')[2]
        max_start = len(sequence) - 20
        start = random.randint(0, max_start)
        remaining = len(sequence) - start
        length = random.randint(20, min(35, remaining))
        peptide = sequence[start:start + length]
        out.write(peptide + '\n')
        count += 1

print(f"Written {count} peptides to {output_file}")
EOF
