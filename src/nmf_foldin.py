#!/usr/bin/env python3
"""NMF folding-in: fix W (learned from scHi-C), estimate H for new bulk Hi-C samples."""
import sys, os
import numpy as np
import pandas as pd
import cooler
from scipy import sparse

data_dir = sys.argv[1]  # NMF results dir (W, features)
K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
resolution = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
n_iter = int(sys.argv[4]) if len(sys.argv) > 4 else 200

print("[info] loading W and feature index...", flush=True)
W = np.loadtxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), delimiter="\t")
feat = pd.read_csv(os.path.join(data_dir, "features.tsv"), sep="\t",
                   header=None, names=["chrom", "bin_i", "bin_j"])
N = W.shape[0]
print(f"  W: {W.shape}, features: {N}", flush=True)

chrom_sizes = {
    "chr1": 195471971, "chr2": 182113224, "chr3": 160039680,
    "chr4": 156508116, "chr5": 151834684, "chr6": 149736546,
    "chr7": 145441459, "chr8": 129401213, "chr9": 124595110,
    "chr10": 130694993, "chr11": 122082543, "chr12": 120129022,
    "chr13": 120421639, "chr14": 124902244, "chr15": 104043685,
    "chr16": 98207768, "chr17": 94987271, "chr18": 90702639,
    "chr19": 61431566, "chrX": 171031299, "chrY": 91744698,
}

def chrom_sort_key(c):
    n = c.replace("chr", "")
    return 100 if n == "X" else 101 if n == "Y" else int(n)

chroms = sorted(chrom_sizes.keys(), key=chrom_sort_key)
chrom_nbins = {c: (chrom_sizes[c] + resolution - 1) // resolution for c in chroms}

max_dist_bins = 50  # 5Mb / 100kb

chrom_row_offsets = {}
for c in chroms:
    nb = chrom_nbins[c]
    offsets = [0] * (nb + 1)
    for i in range(nb):
        offsets[i+1] = offsets[i] + min(max_dist_bins, nb - 1 - i) + 1
    chrom_row_offsets[c] = offsets

chrom_feat_offset = {}
total = 0
for c in chroms:
    chrom_feat_offset[c] = total
    total += chrom_row_offsets[c][chrom_nbins[c]]
assert total == N, f"Feature count mismatch: {total} vs {N}"


def cool_to_feature_vector(cool_uri):
    """Extract upper-triangle feature vector from a .cool file matching NMF feature index."""
    clr = cooler.Cooler(cool_uri)
    x = np.zeros(N, dtype=np.float64)
    for chrom in chroms:
        nb = chrom_nbins[chrom]
        try:
            mat = clr.matrix(balance=True).fetch(chrom)
        except Exception:
            mat = clr.matrix(balance=False).fetch(chrom)
        mat = np.nan_to_num(mat, 0.0)
        off = chrom_feat_offset[chrom]
        row_off = chrom_row_offsets[chrom]
        for i in range(min(nb, mat.shape[0])):
            j_max = min(i + max_dist_bins, nb - 1, mat.shape[1] - 1)
            for j in range(i, j_max + 1):
                idx = off + row_off[i] + (j - i)
                if idx < N:
                    x[idx] = mat[i, j]
    return x


def foldin_mu(W, X_new, n_iter=200, eps=1e-10):
    """MU update for H only, with W fixed. X_new: N x M_new, W: N x K -> H: K x M_new."""
    K = W.shape[1]
    M = X_new.shape[1]
    np.random.seed(42)
    H = np.random.rand(K, M) + 1e-4
    WtW = W.T @ W  # K x K (precompute, constant)
    WtX = W.T @ X_new  # K x M
    for it in range(n_iter):
        H = H * WtX / (WtW @ H + eps)
        if (it + 1) % 50 == 0 or it == 0:
            err = np.sqrt(np.sum((X_new - W @ H) ** 2))
            print(f"  iter {it+1}: ||X - WH||_F = {err:.1f}", flush=True)
    return H


# --- Extract feature vectors from Bonev bulk Hi-C ---
bulk_cools = {
    "ES": "data/brain/processed/brain_hic_es.mcool::resolutions/100000",
    "NPC": "data/brain/processed/brain_hic_npc.mcool::resolutions/100000",
    "CN": "data/brain/processed/brain_hic_cn.mcool::resolutions/100000",
}

print("\n[step1] Extracting feature vectors from bulk Hi-C...", flush=True)
X_bulk = np.zeros((N, len(bulk_cools)))
sample_names = list(bulk_cools.keys())
for si, (name, uri) in enumerate(bulk_cools.items()):
    print(f"  {name}...", flush=True)
    X_bulk[:, si] = cool_to_feature_vector(uri)
    nnz = np.sum(X_bulk[:, si] > 0)
    print(f"    nnz={nnz}/{N} ({100*nnz/N:.1f}%), mean={X_bulk[:, si].mean():.4f}", flush=True)

# --- Folding-in ---
print(f"\n[step2] Folding-in (MU, H only, {n_iter} iterations)...", flush=True)
H_bulk = foldin_mu(W, X_bulk, n_iter=n_iter)

print(f"\n[result] H_bulk (K={K} x {len(sample_names)} samples):")
df_h = pd.DataFrame(H_bulk, index=[f"comp{k+1}" for k in range(K)], columns=sample_names)
print(df_h.round(4).to_string())

# Normalize per sample (fraction)
H_frac = H_bulk / H_bulk.sum(axis=0, keepdims=True)
print(f"\nH_bulk (fraction):")
df_frac = pd.DataFrame(H_frac, index=[f"comp{k+1}" for k in range(K)], columns=sample_names)
print(df_frac.round(3).to_string())

# Save
os.makedirs("output/tad_calling", exist_ok=True)
df_h.to_csv("output/tad_calling/H_foldin_bonev.tsv", sep="\t")
df_frac.to_csv("output/tad_calling/H_foldin_bonev_frac.tsv", sep="\t")

# --- Also fold-in HiRES per-celltype pseudo-bulk ---
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["cell", "label"])
H_orig = np.loadtxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), delimiter="\t")
ct_names = sorted(cells["label"].unique())

print(f"\n[step3] HiRES per-celltype H (original NMF, mean):")
H_ct = np.zeros((K, len(ct_names)))
for ci, ct in enumerate(ct_names):
    H_ct[:, ci] = H_orig[:, cells["label"].values == ct].mean(axis=1)
df_ct = pd.DataFrame(H_ct / H_ct.sum(axis=0, keepdims=True),
                     index=[f"comp{k+1}" for k in range(K)], columns=ct_names)
print(df_ct.round(3).to_string())

# --- Plot ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
x = np.arange(K)
width = 0.25
for si, name in enumerate(sample_names):
    ax.bar(x + si * width, H_frac[:, si], width, label=name)
ax.set_xticks(x + width)
ax.set_xticklabels([f"comp{k+1}" for k in range(K)])
ax.set_ylabel("H fraction")
ax.set_title("Bonev bulk Hi-C folded-in")
ax.legend()

ax = axes[1]
ct_frac = H_ct / H_ct.sum(axis=0, keepdims=True)
width = 0.15
for ci, ct in enumerate(ct_names):
    ax.bar(x + ci * width, ct_frac[:, ci], width, label=ct)
ax.set_xticks(x + width * 2)
ax.set_xticklabels([f"comp{k+1}" for k in range(K)])
ax.set_ylabel("H fraction")
ax.set_title("HiRES scHi-C (original NMF)")
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("plot/diagnostic_foldin_bonev.png", dpi=150)
print("\n[saved] plot/diagnostic_foldin_bonev.png")
