#!/usr/bin/env python3
"""Distance-band NMF with shared H: X_d ≈ W_d @ H for each distance d."""
import sys, os
import numpy as np
import pandas as pd

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
n_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 500

print(f"[info] loading {data_dir}...", flush=True)
X = np.load(os.path.join(data_dir, "X_bulk.npy"))
feat = pd.read_csv(os.path.join(data_dir, "features.tsv"), sep="\t",
                   header=None, names=["chrom", "bin_i", "bin_j"])
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["sample", "label"])
N, M = X.shape
feat_dist = (feat["bin_j"] - feat["bin_i"]).values
max_d = feat_dist.max()
print(f"  X: {N} x {M}, K={K}, max_dist={max_d}", flush=True)

# Filter out sex chromosomes
autosome_mask = ~feat["chrom"].isin(["chrX", "chrY"])
n_auto = autosome_mask.sum()
print(f"[info] filtering to autosomes: {n_auto}/{N} features ({100*n_auto/N:.1f}%)", flush=True)
X = X[autosome_mask.values, :]
feat_dist = feat_dist[autosome_mask.values]
feat = feat[autosome_mask].reset_index(drop=True)
N = X.shape[0]

# Group features by distance
print("[info] grouping features by distance band...", flush=True)
band_indices = {}
for d in range(max_d + 1):
    idx = np.where(feat_dist == d)[0]
    if len(idx) > 0:
        band_indices[d] = idx

n_bands = len(band_indices)
print(f"  {n_bands} distance bands, sizes: d=0 -> {len(band_indices[0])}, "
      f"d=1 -> {len(band_indices[1])}, ..., d={max_d} -> {len(band_indices[max_d])}")

# Per-band per-sample L1 normalization
print("[info] per-band per-sample normalization...", flush=True)
X_bands = {}
for d, idx in band_indices.items():
    Xd = X[idx, :].copy()
    col_sums = Xd.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    X_bands[d] = Xd / col_sums[np.newaxis, :]

# Coupled NMF: X_d ≈ W_d @ H, H shared, with L2 regularization on W
lambda_w = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
print(f"\n[info] Coupled NMF (distance-band, shared H), K={K}, {n_iter} iters, lambda_W={lambda_w}...", flush=True)
np.random.seed(42)
eps = 1e-10

H = np.random.rand(K, M).astype(np.float64) + 1e-6
W_bands = {}
for d, idx in band_indices.items():
    nd = len(idx)
    W_bands[d] = np.random.rand(nd, K).astype(np.float64) + 1e-6

HHt = H @ H.T

for it in range(n_iter):
    # W_d update (per band, with L2 regularization)
    for d in band_indices:
        Xd = X_bands[d]
        Wd = W_bands[d]
        XdHt = Xd @ H.T
        Wd *= XdHt / (Wd @ HHt + lambda_w * Wd + eps)
        W_bands[d] = Wd

    # H update (aggregate across all bands)
    num = np.zeros((K, M), dtype=np.float64)
    den = np.zeros((K, K), dtype=np.float64)
    for d in band_indices:
        Wd = W_bands[d]
        Xd = X_bands[d]
        num += Wd.T @ Xd
        den += Wd.T @ Wd
    H *= num / (den @ H + eps)
    HHt = H @ H.T

    if (it + 1) % 50 == 0 or it == 0:
        total_err = 0.0
        for d in band_indices:
            diff = X_bands[d] - W_bands[d] @ H
            total_err += np.sum(diff ** 2)
        total_err = np.sqrt(total_err)
        print(f"  iter {it+1}: total ||X - WH||_F = {total_err:.6f}", flush=True)

# Save H
out_dir = data_dir.replace("bulk_nmf", "bulk_nmf_distband")
os.makedirs(out_dir, exist_ok=True)

np.savetxt(os.path.join(out_dir, f"H_nmf_k{K}.tsv"), H, delimiter="\t")

# Save W_bands as single concatenated W (same order as features)
W_full = np.zeros((N, K), dtype=np.float64)
for d, idx in band_indices.items():
    W_full[idx, :] = W_bands[d]
np.savetxt(os.path.join(out_dir, f"W_nmf_k{K}.tsv"), W_full, delimiter="\t")

# Save filtered features and cells
feat.to_csv(os.path.join(out_dir, "features.tsv"), sep="\t", header=False, index=False)
import shutil
src = os.path.join(data_dir, "cells.tsv")
dst = os.path.join(out_dir, "cells.tsv")
if os.path.exists(src):
    shutil.copy2(src, dst)

print(f"\n[saved] {out_dir}/W_nmf_k{K}.tsv ({W_full.shape})", flush=True)
print(f"[saved] {out_dir}/H_nmf_k{K}.tsv ({H.shape})", flush=True)

# Show H matrix
H_frac = H / H.sum(axis=0, keepdims=True)
df = pd.DataFrame(H_frac, index=[f"comp{k+1}" for k in range(K)], columns=cells["sample"])
print(f"\nH matrix (K={K} x {M} samples):", flush=True)
print(df.round(3).to_string(), flush=True)
