#!/bin/bash

#################################
# Run mTEC SRA Hi-C: SRA fastq -> nf-core/hic -> bulk mcool
#################################
# Heavy: ~30 GB SRA download, ~200 GB intermediate, 12-24 h on 32 cores.
# Run separately from the main preprocess pipeline.

set -e

echo "================================================"
echo "AIRE-experiments: mTEC SRA Hi-C (nf-core)"
echo "================================================"
echo ""

# See workflow/run_preprocess.sh for the rationale on the conda symlink shim.
SNAKEMAKE_ENV="${SNAKEMAKE_ENV:-/home/koki/anaconda3/envs/snakemake}"
SNAKEMAKE_BIN="${SNAKEMAKE_ENV}/bin/snakemake"

CONDA_SHIM=$(mktemp -d -t aire_conda_shim_XXXXXX)
trap "rm -rf '$CONDA_SHIM'" EXIT
ln -s "${SNAKEMAKE_ENV}/bin/conda" "${CONDA_SHIM}/conda"
[ -e "${SNAKEMAKE_ENV}/bin/mamba" ] && ln -s "${SNAKEMAKE_ENV}/bin/mamba" "${CONDA_SHIM}/mamba"
export PATH="${CONDA_SHIM}:${PATH}"

"$SNAKEMAKE_BIN" \
  --snakefile workflow/sra_hic.smk \
  --cores 16 \
  --printshellcmds \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  "$@"

echo ""
echo "================================================"
echo "mTEC SRA Hi-C complete."
echo "================================================"
echo ""
echo "Outputs:"
echo "  - Per-sample cools  : data/mtec/processed/nfcore_hic/contact_maps/cool/"
echo "  - Bulk mTEC mcool   : data/mtec/processed/mtec_hic.mcool"
echo ""
echo "Next: re-run workflow/run_preprocess.sh; bin_hic_by_chr (tissue=mtec)"
echo "      will pick up the new mcool."
echo ""
