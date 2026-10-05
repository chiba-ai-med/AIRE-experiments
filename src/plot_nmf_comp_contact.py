#!/usr/bin/env python3
"""Visualize NMF component W columns as contact maps."""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb_oe"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
resolution = 100000

print("[info] loading...", flush=True)
W = np.loadtxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), delimiter="\t")
feat = pd.read_csv(os.path.join(data_dir, "features.tsv"), sep="\t",
                   header=None, names=["chrom", "bin_i", "bin_j"])
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["sample", "label"])
H = np.loadtxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), delimiter="\t")
print(f"  W: {W.shape}, features: {len(feat)}", flush=True)

sample_names = cells["sample"].tolist()
H_frac = H / H.sum(axis=0, keepdims=True)

regions = [
    ("chr10", 170, 270, "Aire locus (chr10:17-27Mb)"),
    ("chr1", 500, 700, "chr1:50-70Mb"),
    ("chr6", 500, 650, "chr6:50-65Mb"),
]


def w_to_matrix(w_col, chrom, bin_start, bin_end):
    """Extract W column values for a region and fill into a symmetric matrix."""
    mask = feat["chrom"] == chrom
    feat_chr = feat[mask]
    w_chr = w_col[mask.values]
    size = bin_end - bin_start
    mat = np.zeros((size, size))
    for idx in range(len(feat_chr)):
        bi = feat_chr.iloc[idx]["bin_i"]
        bj = feat_chr.iloc[idx]["bin_j"]
        if bi >= bin_start and bi < bin_end and bj >= bin_start and bj < bin_end:
            li = bi - bin_start
            lj = bj - bin_start
            mat[li, lj] = w_chr[idx]
            mat[lj, li] = w_chr[idx]
    return mat


fig, axes = plt.subplots(len(regions), K, figsize=(3.5 * K, 3.5 * len(regions)))

for ri, (chrom, bstart, bend, title) in enumerate(regions):
    for k in range(K):
        ax = axes[ri, k]
        mat = w_to_matrix(W[:, k], chrom, bstart, bend)
        start_mb = bstart * resolution / 1e6
        end_mb = bend * resolution / 1e6
        extent = [start_mb, end_mb, end_mb, start_mb]
        vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax,
                  extent=extent, interpolation="none")
        if ri == 0:
            top_sample = sample_names[np.argmax(H_frac[k, :])]
            ax.set_title(f"comp{k+1}\n({top_sample})", fontsize=8)
        if k == 0:
            ax.set_ylabel(f"{title}\n(Mb)", fontsize=8)
        ax.tick_params(labelsize=6)

plt.suptitle(f"NMF W columns as contact maps (O/E, K={K})", fontsize=13, y=1.01)
plt.tight_layout()
out_png = "plot/nmf_comp_contact_maps.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\n[saved] {out_png}", flush=True)
