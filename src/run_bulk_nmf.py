#!/usr/bin/env python3
"""MU-NMF on bulk Hi-C feature matrix (dense, small M)."""
import sys, os
import numpy as np
import pandas as pd

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
n_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 500

print(f"[info] loading {data_dir}...", flush=True)
X = np.load(os.path.join(data_dir, "X_bulk.npy"))
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["sample", "label"])
N, M = X.shape
print(f"  X: {N} x {M}", flush=True)

# Per-sample L1 normalization (sum to 1) to handle scale differences
col_sums = X.sum(axis=0)
print(f"  Column sums before normalization: {col_sums.round(1)}", flush=True)
X_norm = X / col_sums[np.newaxis, :]
print(f"  Column sums after normalization: {X_norm.sum(axis=0).round(4)}", flush=True)

# MU-NMF
print(f"\n[info] MU-NMF, K={K}, {n_iter} iterations...", flush=True)
np.random.seed(42)
W = np.random.rand(N, K).astype(np.float64) + 1e-6
H = np.random.rand(K, M).astype(np.float64) + 1e-6
eps = 1e-10

x_sq = np.sum(X_norm ** 2)
for it in range(n_iter):
    # H update
    WtX = W.T @ X_norm
    WtW = W.T @ W
    H = H * WtX / (WtW @ H + eps)

    # W update
    XHt = X_norm @ H.T
    HHt = H @ H.T
    W = W * XHt / (W @ HHt + eps)

    if (it + 1) % 50 == 0 or it == 0:
        wh_sq = np.trace((W.T @ W) @ (H @ H.T))
        cross = np.sum((X_norm.T @ W) * H.T)
        err = np.sqrt(max(x_sq + wh_sq - 2 * cross, 0))
        print(f"  iter {it+1}: ||X - WH||_F = {err:.6f}", flush=True)

# Save
np.savetxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), W, delimiter="\t")
np.savetxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), H, delimiter="\t")
print(f"\n[saved] W_nmf_k{K}.tsv ({W.shape}), H_nmf_k{K}.tsv ({H.shape})", flush=True)

# Show H matrix
print(f"\nH matrix (K={K} x {M} samples):", flush=True)
H_frac = H / H.sum(axis=0, keepdims=True)
df = pd.DataFrame(H_frac, index=[f"comp{k+1}" for k in range(K)], columns=cells["sample"])
print(df.round(3).to_string(), flush=True)
