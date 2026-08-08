#!/usr/bin/env bash
set -euo pipefail

# Boots sa-server against a built index, replays a fixed peptide file through
# /search, and writes the raw JSON response.
#
# The backend (mmap vs preloaded) is a COMPILE-TIME choice, so it is selected
# here via SA_FEATURES rather than by a server flag.
#
# Environment:
#   SA_FEATURES  cargo features for sa-server. "" (default) builds the
#                preloaded backend; "mmap" builds the memory-mapped one.
#   INDEX_DIR    directory holding sa.bin / proteins.bin / mapping.bin
#   KMER_TABLE   optional path to a pre-built k-mer bounds table
#   EQUATE_IL    "true"/"false" (default false) - equate I and L during search
#   TRYPTIC      "true"/"false" (default false) - restrict to tryptic peptides
#   CUTOFF       max matches per peptide (default: server default)
#
# Usage: [env ...] scripts/run_index.sh [output.json]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR/.."
INDEX_DIR="${INDEX_DIR:-$REPO_ROOT/uniprot-2026-01/suffix-array}"
PEPTIDES_FILE="${PEPTIDES_FILE:-$SCRIPT_DIR/peptides.txt}"
OUTPUT_FILE="${1:-$SCRIPT_DIR/output.json}"
SA_FEATURES="${SA_FEATURES:-}"
KMER_TABLE="${KMER_TABLE:-}"
EQUATE_IL="${EQUATE_IL:-false}"
TRYPTIC="${TRYPTIC:-false}"
SERVER_LOG="$(mktemp)"

if [ ! -f "$PEPTIDES_FILE" ]; then
    echo "Error: peptides.txt not found at $PEPTIDES_FILE" >&2
    echo "Run generate_peptides.sh first." >&2
    exit 1
fi

# Build sa-server to ensure we're using the current version.
# Each feature set gets its own target dir so that flipping SA_FEATURES back and
# forth does not force a full rebuild every time.
if [ -n "$SA_FEATURES" ]; then
    echo "Building sa-server (features: $SA_FEATURES)..."
    TARGET_DIR="$REPO_ROOT/target/features-$SA_FEATURES"
    cargo build --release -p sa-server --features "$SA_FEATURES" \
        --manifest-path "$REPO_ROOT/Cargo.toml" --target-dir "$TARGET_DIR"
else
    echo "Building sa-server (features: none -> preloaded backend)..."
    TARGET_DIR="$REPO_ROOT/target"
    cargo build --release -p sa-server --manifest-path "$REPO_ROOT/Cargo.toml"
fi

SA_SERVER="$TARGET_DIR/release/sa-server"

# Kill any process already using port 3000
if lsof -ti :3000 > /dev/null 2>&1; then
    echo "Port 3000 already in use, killing existing process..."
    kill "$(lsof -ti :3000)" 2>/dev/null || true
    sleep 1
fi

SERVER_ARGS=(
    -d "$INDEX_DIR/proteins.bin"
    -i "$INDEX_DIR/sa.bin"
    --mapping-file "$INDEX_DIR/mapping.bin"
)
if [ -n "$KMER_TABLE" ]; then
    SERVER_ARGS+=(--kmer-table-file "$KMER_TABLE")
fi

# Start the server in the background
echo "Starting sa-server..."
"$SA_SERVER" "${SERVER_ARGS[@]}" > "$SERVER_LOG" 2>&1 &
SERVER_PID=$!

PAYLOAD_FILE="$(mktemp)"
cleanup() {
    kill "$SERVER_PID" 2>/dev/null || true
    wait "$SERVER_PID" 2>/dev/null || true
    rm -f "$SERVER_LOG" "$PAYLOAD_FILE"
}
trap cleanup EXIT

# Wait for the server to be ready
echo "Waiting for server to be ready..."
while ! grep -q "Server is ready" "$SERVER_LOG" 2>/dev/null; do
    sleep 1
    if ! kill -0 "$SERVER_PID" 2>/dev/null; then
        echo "Error: server process exited unexpectedly" >&2
        cat "$SERVER_LOG" >&2
        exit 1
    fi
done
echo "Server is ready!"

# Build JSON payload from peptides.txt and query the server
echo "Querying index with $(wc -l < "$PEPTIDES_FILE") peptides (equate_il=$EQUATE_IL tryptic=$TRYPTIC)..."
jq -R . "$PEPTIDES_FILE" \
    | jq -s --argjson equate_il "$EQUATE_IL" --argjson tryptic "$TRYPTIC" \
        '{peptides: ., equate_il: $equate_il, tryptic: $tryptic}' \
    > "$PAYLOAD_FILE"
if [ -n "${CUTOFF:-}" ]; then
    jq --argjson cutoff "$CUTOFF" '. + {cutoff: $cutoff}' "$PAYLOAD_FILE" > "$PAYLOAD_FILE.tmp"
    mv "$PAYLOAD_FILE.tmp" "$PAYLOAD_FILE"
fi

curl -s -X POST http://localhost:3000/search \
    -H "Content-Type: application/json" \
    -d "@$PAYLOAD_FILE" > "$OUTPUT_FILE"

echo "Output saved to $OUTPUT_FILE"
