#!/usr/bin/env python
"""Re-bin the ATAC modality of a labeled multiome MuData onto fixed genomic
bins, split by chromosome, and write one MatrixMarket file per chromosome.

Output layout (under <output_dir>):
    chr1.mtx, chr2.mtx, ..., chrX.mtx, chrY.mtx  -- cell x bin sparse counts
    chroms.txt                                    -- chromosome names, one per line
    bins.tsv.gz                                   -- BED-like (chrom, start, end) for ALL bins
    cells.tsv                                     -- cell barcodes
    labels.tsv                                    -- per-cell labels (mdata.obs['cell_type'])

Per-chromosome files are sized n_cells x n_bins(chr). Loading on the R side:
    library(Matrix)
    chroms <- readLines(file.path(d, "chroms.txt"))
    X_RNA <- lapply(chroms, function(c) t(readMM(file.path(d, paste0(c, ".mtx")))))
    # transpose: file is cell x bin, Machima2 expects feature x cell (n_k x m).
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import mudata as md
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import snapatac2 as snap


def parse_peak_coords(var: pd.DataFrame) -> pd.DataFrame:
    """Extract (chrom, start, end) from peak feature names like 'chr1:12345-23456'."""
    if {"chrom", "start", "end"}.issubset(var.columns):
        return var[["chrom", "start", "end"]].copy()
    # snapatac2 typically stores peak coords in var.index as 'chrom:start-end'.
    parts = var.index.to_series().str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    if parts.isna().any().any():
        raise RuntimeError("Could not parse peak coordinates from var.index")
    return parts.assign(start=lambda df: df["start"].astype(int),
                        end=lambda df: df["end"].astype(int))


def make_bins(chrom_sizes: dict[str, int], resolution: int) -> pd.DataFrame:
    rows = []
    for chrom, length in chrom_sizes.items():
        starts = np.arange(0, length, resolution)
        ends = np.minimum(starts + resolution, length)
        rows.append(pd.DataFrame({"chrom": chrom, "start": starts, "end": ends}))
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True, help="Labeled .h5mu")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--genome", default="mm10")
    p.add_argument("--resolution", type=int, required=True, help="Bin size in bp")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input}")
    mdata = md.read(str(args.input))
    atac = mdata["atac"]
    print(f"[atac] shape: {atac.shape} (cells x peaks)")

    chrom_sizes = dict(getattr(snap.genome, args.genome).chrom_sizes.items())
    # Restrict to canonical autosomes + X + Y; drop alt/random contigs.
    canonical = [c for c in chrom_sizes if c.startswith("chr") and "_" not in c]
    chrom_sizes = {c: chrom_sizes[c] for c in canonical}

    bins = make_bins(chrom_sizes, args.resolution)
    bins["bin_idx"] = np.arange(len(bins))
    print(f"[bins] {len(bins)} total bins at {args.resolution} bp across {len(chrom_sizes)} chromosomes")

    peaks = parse_peak_coords(atac.var)
    peaks = peaks.reset_index(drop=True)
    peaks["peak_idx"] = np.arange(len(peaks))

    # Map each peak to its bin (assignment via floor(start / resolution) within chrom).
    peaks_in = peaks[peaks["chrom"].isin(chrom_sizes)].copy()
    peaks_in["bin_local"] = (peaks_in["start"].astype(int) // args.resolution).astype(int)

    chrom_first_bin = bins.groupby("chrom", sort=False)["bin_idx"].first().to_dict()
    peaks_in["bin_global"] = peaks_in["chrom"].map(chrom_first_bin) + peaks_in["bin_local"]

    # Build a (peak -> bin) sparse mapping matrix M (n_peaks x n_bins_total).
    n_peaks = atac.n_vars
    n_bins = len(bins)
    M = sp.csr_matrix(
        (np.ones(len(peaks_in)),
         (peaks_in["peak_idx"].values, peaks_in["bin_global"].values)),
        shape=(n_peaks, n_bins),
    )

    X = atac.X.tocsr() if sp.issparse(atac.X) else sp.csr_matrix(atac.X)
    cell_by_bin = (X @ M).tocsr()
    print(f"[rebin] cell x bin: {cell_by_bin.shape}, nnz={cell_by_bin.nnz}")

    # Write per-chromosome submatrices.
    chroms_out = []
    for chrom, group in bins.groupby("chrom", sort=False):
        idx = group["bin_idx"].values
        sub = cell_by_bin[:, idx]
        out = args.output_dir / f"{chrom}.mtx"
        scipy.io.mmwrite(str(out), sub)
        chroms_out.append(chrom)
        print(f"  {chrom}: {sub.shape} -> {out.name}")

    (args.output_dir / "chroms.txt").write_text("\n".join(chroms_out) + "\n")
    bins[["chrom", "start", "end"]].to_csv(
        args.output_dir / "bins.tsv.gz", sep="\t", index=False, header=False, compression="gzip")
    pd.Series(atac.obs_names, name="barcode").to_csv(
        args.output_dir / "cells.tsv", sep="\t", index=False, header=False)
    if "cell_type" in mdata.obs:
        mdata.obs["cell_type"].to_csv(
            args.output_dir / "labels.tsv", sep="\t", header=False)
    print(f"[write] {args.output_dir}")


if __name__ == "__main__":
    main()
