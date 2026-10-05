#!/usr/bin/env bash
# Download GEO supplementary files via NCBI FTP (over HTTPS).
#
# Usage:
#   preprocess_download_geo.sh series  GSE210747               <outdir>
#   preprocess_download_geo.sh samples "GSM2533835 GSM2533836" <outdir>
#
# - series mode: downloads everything in
#     ftp://ftp.ncbi.nlm.nih.gov/geo/series/<GSE-prefix>nnn/<GSE>/suppl/
#   into <outdir>/.
# - samples mode: takes a space-separated list of GSM IDs and downloads each
#     ftp://ftp.ncbi.nlm.nih.gov/geo/samples/<GSM-prefix>nnn/<GSM>/suppl/
#   into <outdir>/<GSM>/.
#
# Notes:
# - `--execute robots=off` is required: NCBI's robots.txt disallows recursion,
#   which silently makes wget exit 0 with only the index.html fetched.
# - `--timestamping` makes re-runs idempotent (skips files already on disk
#   with the same size and mtime).
# - After each download we validate that >=1 payload file was actually
#   fetched, otherwise wget's exit-0-on-empty-recursion would let us touch
#   sentinels with no data.

set -euo pipefail

WGET_ARGS=(
  --no-verbose
  --recursive
  --no-parent
  --no-host-directories
  --cut-dirs=5
  --execute robots=off
  --timestamping
  --reject "index.html*,robots.txt"
)

# Count "real" payload files (anything not index.html or robots.txt).
count_payload() {
  local dir="$1"
  find "$dir" -type f ! -name 'index.html*' ! -name 'robots.txt' ! -name '.downloaded' \
    | wc -l
}

mode="${1:-}"
shift || { echo "missing mode" >&2; exit 1; }

case "$mode" in
  series)
    accession="${1:?accession required}"
    outdir="${2:?outdir required}"
    mkdir -p "$outdir"
    # GSE prefix dir replaces the last 3 chars with "nnn"
    prefix="${accession:0:${#accession}-3}nnn"
    url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${accession}/suppl/"
    echo "[$(date -Is)] Downloading $url -> $outdir"
    wget "${WGET_ARGS[@]}" --directory-prefix="$outdir" "$url"
    n=$(count_payload "$outdir")
    echo "[$(date -Is)] $accession: downloaded $n payload file(s)"
    if [ "$n" -eq 0 ]; then
      echo "ERROR: no payload files downloaded for $accession" >&2
      exit 1
    fi
    ;;

  samples)
    gsm_list="${1:?gsm list required}"
    outdir="${2:?outdir required}"
    mkdir -p "$outdir"
    for gsm in $gsm_list; do
      prefix="${gsm:0:${#gsm}-3}nnn"
      url="https://ftp.ncbi.nlm.nih.gov/geo/samples/${prefix}/${gsm}/suppl/"
      sample_outdir="$outdir/$gsm"
      mkdir -p "$sample_outdir"
      echo "[$(date -Is)] Downloading $url -> $sample_outdir"
      wget "${WGET_ARGS[@]}" --directory-prefix="$sample_outdir" "$url"
      n=$(count_payload "$sample_outdir")
      echo "[$(date -Is)] $gsm: downloaded $n payload file(s)"
      if [ "$n" -eq 0 ]; then
        echo "ERROR: no payload files downloaded for $gsm" >&2
        exit 1
      fi
    done
    ;;

  *)
    echo "Unknown mode: $mode (expected: series | samples)" >&2
    exit 1
    ;;
esac
