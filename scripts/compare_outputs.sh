#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>" >&2
    exit 1
fi

FILE1="$1"
FILE2="$2"

if [ ! -f "$FILE1" ]; then
    echo "Error: file not found: $FILE1" >&2
    exit 1
fi

if [ ! -f "$FILE2" ]; then
    echo "Error: file not found: $FILE2" >&2
    exit 1
fi

normalize() {
    jq --sort-keys '[.[] | .proteins |= sort_by(.uniprot_accession)] | sort_by(.sequence)' "$1"
}

if diff <(normalize "$FILE1") <(normalize "$FILE2") > /dev/null 2>&1; then
    echo "Files are identical"
else
    diff <(normalize "$FILE1") <(normalize "$FILE2")
    echo "Files differ"
    exit 1
fi
