#!/bin/bash
# Run machima phase (Stage A NMF + Stage B Machima2 for every stage x T).
# Snakemake auto-builds any missing preprocess outputs.
echo "================================================"
echo "AIRE-experiments: machima phase"
echo "================================================"
echo ""

source "$(dirname "$0")/_run_helper.sh"
aire_run machima "$@"

echo ""
echo "================================================"
echo "Machima phase complete."
echo "================================================"
