#!/bin/bash
# Run evaluate phase (self-fit metrics + held-out per-celltype validation).
echo "================================================"
echo "AIRE-experiments: evaluate phase"
echo "================================================"
echo ""

source "$(dirname "$0")/_run_helper.sh"
aire_run evaluate "$@"

echo ""
echo "================================================"
echo "Evaluate phase complete."
echo "================================================"
