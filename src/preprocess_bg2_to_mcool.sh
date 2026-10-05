#!/usr/bin/env bash
# Aggregate per-GSM BedGraph2D files into a single bulk multi-resolution mcool.
#
# Usage:
#   preprocess_bg2_to_mcool.sh <chrom_sizes> <base_resolution> <zoom_resolutions> <output_mcool> <bg2_1> [<bg2_2> ...]
#
# Steps:
#   1. For each input bg2.gz: `cooler load -f bg2 chrom_sizes:base_res in.bg2 out.cool`
#   2. `cooler merge merged.cool <per-gsm cools>` (sums counts at common bins)
#   3. `cooler zoomify --balance --resolutions <zoom_resolutions> -o output.mcool merged.cool`
#
# Per-GSM bg2 files are produced by preprocess_extract_bonev_track.R from
# the Bonev/Tanay misha 2D tracks (GSE96107 tarballs).

set -euo pipefail

chrom_sizes="${1:?chrom_sizes required}"
base_resolution="${2:?base_resolution required}"
zoom_resolutions="${3:?zoom_resolutions (comma-separated) required}"
output_mcool="${4:?output_mcool required}"
shift 4
inputs=( "$@" )
if [ ${#inputs[@]} -eq 0 ]; then
  echo "ERROR: at least one bg2 input required" >&2
  exit 1
fi

# Per-output work dir so parallel rule firings (e.g. bg2_to_mcool_brain_celltype
# for npc + cn) don't collide on the same hdf5 temp files.
work_basename=$(basename "$output_mcool" .mcool)
work_dir="$(dirname "$output_mcool")/_bg2_work_${work_basename}"
mkdir -p "$work_dir"

#################################
# Per-input: bg2 -> .cool at base resolution
#################################
per_gsm_cools=()
for bg2 in "${inputs[@]}"; do
  base=$(basename "$bg2" .bg2.gz)
  cool="$work_dir/${base}_${base_resolution}.cool"
  echo "[$(date -Is)] cooler load: $bg2 -> $cool"
  # Bonev/Tanay misha tracks:
  #   * store float-valued normalised contact intensities (--count-as-float)
  #   * emit BOTH upper- and lower-triangle for intra-chromosomal pairs
  #     (--input-copy-status duplex drops the lower triangle without halving
  #     the value -- correct for already-normalised intensities)
  cooler load -f bg2 --count-as-float --input-copy-status duplex \
      "${chrom_sizes}:${base_resolution}" "$bg2" "$cool"
  per_gsm_cools+=( "$cool" )
done

#################################
# Merge per-GSM cools (cooler merge sums counts at coincident bins)
#################################
merged_cool="$work_dir/bulk_${base_resolution}.cool"
echo "[$(date -Is)] cooler merge: ${#per_gsm_cools[@]} cools -> $merged_cool"
cooler merge "$merged_cool" "${per_gsm_cools[@]}"

#################################
# Zoomify + balance
#################################
echo "[$(date -Is)] cooler zoomify --balance --resolutions $zoom_resolutions"
cooler zoomify --balance --resolutions "$zoom_resolutions" -o "$output_mcool" "$merged_cool"

echo "[$(date -Is)] Done: $output_mcool"
ls -lh "$output_mcool"
