#!/usr/bin/env bash
# MACS2 peak calling on a 10x-style fragments.tsv.gz file.
#
# Fragments format (tab-separated, gzipped):
#   chrom  start  end  cell_barcode  count
#
# We feed it to MACS2 with -f BED (only first 3 cols used) and the standard
# ATAC shift/extend recipe. Output is a narrowPeak BED that downstream
# rules use to build the cell x peak matrix.
#
# Usage:
#   preprocess_call_peaks.sh <fragments.tsv.gz> <name> <outdir> [<gsize=mm>]

set -euo pipefail

fragments="${1:?fragments path required}"
name="${2:?sample name required}"
outdir="${3:?outdir required}"
gsize="${4:-mm}"

mkdir -p "$outdir"

echo "[$(date -Is)] MACS2 callpeak on $fragments -> $outdir/${name}_peaks.narrowPeak"

# --shift -75 --extsize 150 centres a 150 bp window on each Tn5 cut site.
# --keep-dup all preserves PCR-deduped fragments (cellranger-arc already deduped).
# --call-summits enables sub-peak resolution for densely co-occurring peaks.
macs3 callpeak \
  --treatment "$fragments" \
  --format BED \
  --gsize "$gsize" \
  --nomodel \
  --shift -75 \
  --extsize 150 \
  --keep-dup all \
  --call-summits \
  --name "$name" \
  --outdir "$outdir"

echo "[$(date -Is)] Peaks: $(wc -l < "$outdir/${name}_peaks.narrowPeak") narrowPeak records"
