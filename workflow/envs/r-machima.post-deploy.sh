#!/usr/bin/env bash
# Post-deploy script for the r-machima conda env.
# Snakemake automatically runs this script after creating the env from
# workflow/envs/r-machima.yaml (matched by basename).
#
# Installs Seurat, Signac, MuData, nnTensor (from CRAN/Bioc) and Machima
# (from GitHub). Kept out of the conda yaml because:
#   - Seurat / Signac stacks are heavy and unstable to resolve via conda
#   - MuData (Bioconductor) is not always mirrored on bioconda
#   - Machima is GitHub-only

set -euo pipefail
Rscript "$(dirname "$0")/../../setup_r_packages.R"
