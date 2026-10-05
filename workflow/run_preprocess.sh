#!/bin/bash
# Run preprocess phase only (data download / peaks / mudata / mcool / bins
# / cluster->celltype map).
echo "================================================"
echo "AIRE-experiments: preprocess phase"
echo "================================================"
echo ""

source "$(dirname "$0")/_run_helper.sh"
aire_run preprocess "$@"

echo ""
echo "================================================"
echo "Preprocess complete."
echo "================================================"
