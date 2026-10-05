#!/bin/bash
# Run plot phase (summary PDF).
echo "================================================"
echo "AIRE-experiments: plot phase"
echo "================================================"
echo ""

source "$(dirname "$0")/_run_helper.sh"
aire_run plot "$@"

echo ""
echo "================================================"
echo "Plot phase complete."
echo "================================================"
