#!/usr/bin/env bash
# Merge a subset of nf-core/hic per-sample cools into a single bulk mcool.
#
# Usage:
#   sra_merge_bulk_mcool.sh <nfcore_outdir> "<sample1 sample2 ...>" \
#                           <base_resolution> <zoom_resolutions> <output_mcool>
#
# nf-core/hic v2.1.0 writes per-sample contact maps to:
#   <nfcore_outdir>/contact_maps/cool/<SAMPLE>.<RESOLUTION>_balanced.cool
# (period before resolution, "_balanced" suffix, weights stored alongside
# raw counts; cooler merge sums the count column ignoring the weights).
#
# We pick the {base_resolution} cool for each requested sample, sum them
# with `cooler merge` to a single bulk cool, then `cooler zoomify --balance`
# to produce the multi-resolution mcool (re-balanced for the merged total).

set -euo pipefail

nfcore_outdir="${1:?nfcore_outdir required}"
samples_str="${2:?samples space-separated required}"
base_resolution="${3:?base_resolution required}"
zoom_resolutions="${4:?zoom_resolutions required}"
output_mcool="${5:?output_mcool required}"

work_dir="$(dirname "$output_mcool")/_merge_work"
mkdir -p "$work_dir"

# Locate the per-sample cools at the base resolution.
inputs=()
for sample in $samples_str; do
  candidate="$nfcore_outdir/contact_maps/cool/${sample}.${base_resolution}_balanced.cool"
  if [ ! -f "$candidate" ]; then
    echo "ERROR: missing $candidate" >&2
    echo "Available cools:" >&2
    ls "$nfcore_outdir/contact_maps/cool/" >&2 || true
    exit 1
  fi
  inputs+=( "$candidate" )
done
echo "[$(date -Is)] Merging ${#inputs[@]} cools at ${base_resolution} bp"
printf '  %s\n' "${inputs[@]}"

merged_cool="$work_dir/bulk_${base_resolution}.cool"
cooler merge "$merged_cool" "${inputs[@]}"

echo "[$(date -Is)] zoomify --balance to $zoom_resolutions"
cooler zoomify --balance --resolutions "$zoom_resolutions" -o "$output_mcool" "$merged_cool"

echo "[$(date -Is)] Done: $output_mcool"
ls -lh "$output_mcool"
