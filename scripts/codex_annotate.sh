#!/usr/bin/env bash
# CellTypePilot — Quick annotation script for Codex CLI
# Usage: bash scripts/codex_annotate.sh <input.h5ad> <cluster_key> [output_dir] [tissue]
#
# This script wraps the full CellTypePilot annotation pipeline into a single
# command that Codex can execute without multi-step orchestration.

set -euo pipefail

INPUT="${1:?Usage: codex_annotate.sh <input.h5ad> <cluster_key> [output_dir] [tissue]}"
CLUSTER_KEY="${2:?Usage: codex_annotate.sh <input.h5ad> <cluster_key> [output_dir] [tissue]}"
OUTPUT_DIR="${3:-./ctp_output}"
TISSUE="${4:-}"

echo "=== CellTypePilot — Codex Quick Annotate ==="
echo ""

# Step 1: Environment check
echo "[1/4] Checking environment..."
celltypepilot doctor 2>/dev/null || {
    echo "ERROR: celltypepilot not found. Run: pip install -e ."
    exit 1
}

# Step 2: Inspect data
echo "[2/4] Inspecting data..."
INSPECT_JSON=$(celltypepilot inspect --input "$INPUT" --cluster-key "$CLUSTER_KEY" --json 2>/dev/null)
SPECIES=$(echo "$INSPECT_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('species','human'))" 2>/dev/null || echo "human")

if [ -z "$TISSUE" ]; then
    TISSUE=$(echo "$INSPECT_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('tissue') or 'general')" 2>/dev/null || echo "general")
fi

EMBEDDING=$(echo "$INSPECT_JSON" | python3 -c "
import sys, json
data = json.load(sys.stdin)
keys = data.get('embedding_keys', [])
for k in keys:
    if 'umap' in k.lower():
        print(k)
        sys.exit(0)
if keys:
    print(keys[0])
" 2>/dev/null || echo "")

echo "  Species: $SPECIES"
echo "  Tissue:  $TISSUE"
echo "  Embedding: ${EMBEDDING:-auto}"
echo ""

# Step 3: Run annotation
echo "[3/4] Running annotation pipeline..."
ANNOTATE_ARGS="--input $INPUT --cluster-key $CLUSTER_KEY --output $OUTPUT_DIR --species $SPECIES --tissue $TISSUE"
if [ -n "$EMBEDDING" ]; then
    ANNOTATE_ARGS="$ANNOTATE_ARGS --embedding-key $EMBEDDING"
fi

celltypepilot annotate $ANNOTATE_ARGS

# Step 4: Summary
echo ""
echo "[4/4] Done! Output files:"
echo "  $OUTPUT_DIR/data.annotated.h5ad"
echo "  $OUTPUT_DIR/evidence_table.csv"
echo "  $OUTPUT_DIR/report_draft.html"
echo "  $OUTPUT_DIR/manifest.json"
echo "  $OUTPUT_DIR/methodology_draft.txt"
if [ -d "$OUTPUT_DIR/figures" ]; then
    echo "  $OUTPUT_DIR/figures/ ($(ls "$OUTPUT_DIR/figures/" | wc -l) files)"
fi
echo ""
echo "=== Annotation complete ==="
