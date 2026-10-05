#!/bin/bash
# Run the full pipeline (preprocess -> machima -> evaluate -> plot).
echo "================================================"
echo "AIRE-experiments: full pipeline (rule all)"
echo "================================================"
echo ""

source "$(dirname "$0")/_run_helper.sh"
aire_run all "$@"

echo ""
echo "================================================"
echo "Full pipeline complete."
echo "================================================"
