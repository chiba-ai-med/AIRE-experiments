#!/bin/bash
# Shared snakemake wrapper used by run.sh / run_<phase>.sh.
# Source this from a wrapper, then call `aire_run <target> [extra-args...]`.
#
# Snakemake env handling -- two competing constraints:
#   (a) Snakemake (8.x) requires conda >= 24.7.1, which lives in the dedicated
#       snakemake env, not the system base. It locates conda via
#       `shutil.which("conda")`, so `conda` MUST be on PATH.
#   (b) Putting the entire snakemake env's bin in PATH makes `conda activate
#       <rule_env>` (run by snakemake for each rule) misbehave -- the rule
#       env lands at PATH position 3 instead of 1, so `python`/`Rscript`
#       resolve to the snakemake env instead of the rule env. Every rule
#       fails with ModuleNotFoundError.
# Fix: a temp dir holding ONLY a symlink to conda (and mamba), prepended to
# PATH. Snakemake finds conda; rule envs activate cleanly.

set -e

SNAKEMAKE_ENV="${SNAKEMAKE_ENV:-/home/koki/anaconda3/envs/snakemake}"
SNAKEMAKE_BIN="${SNAKEMAKE_ENV}/bin/snakemake"

CONDA_SHIM=$(mktemp -d -t aire_conda_shim_XXXXXX)
trap "rm -rf '$CONDA_SHIM'" EXIT
ln -s "${SNAKEMAKE_ENV}/bin/conda" "${CONDA_SHIM}/conda"
[ -e "${SNAKEMAKE_ENV}/bin/mamba" ] && ln -s "${SNAKEMAKE_ENV}/bin/mamba" "${CONDA_SHIM}/mamba"
export PATH="${CONDA_SHIM}:${PATH}"

aire_run() {
  local target="$1"
  shift
  "$SNAKEMAKE_BIN" \
    --snakefile workflow/Snakefile \
    --cores 4 \
    --printshellcmds \
    --keep-going \
    --rerun-incomplete \
    --use-conda \
    "$target" \
    "$@"
}
