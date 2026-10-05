#!/usr/bin/env python3
"""Vectorize bulk Hi-C as ICE-balanced values for NMF."""
import sys, os, glob
import numpy as np
import pandas as pd
import cooler

out_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb_balanced"
resolution = 100000
max_dist_bp = 5000000
max_dist_bins = max_dist_bp // resolution
os.makedirs(out_dir, exist_ok=True)

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

chrom_feat_offset = {}
chrom_row_offsets = {}
total_features = 0
for c in chroms:
    chrom_feat_offset[c] = total_features
    nb = chrom_nbins[c]
    offsets = [0] * (nb + 1)
    for i in range(nb):
        offsets[i+1] = offsets[i] + min(max_dist_bins, nb - 1 - i) + 1
    chrom_row_offsets[c] = offsets
    total_features += offsets[nb]

feat_dist = np.zeros(total_features, dtype=np.int32)
idx = 0
for c in chroms:
    nb = chrom_nbins[c]
    for i in range(nb):
        j_max = min(i + max_dist_bins, nb - 1)
        for j in range(i, j_max + 1):
            feat_dist[idx] = j - i
            idx += 1

print(f"[info] {total_features} features, max_dist_bins={max_dist_bins}", flush=True)

with open(os.path.join(out_dir, "features.tsv"), "w") as f:
    for c in chroms:
        nb = chrom_nbins[c]
        for i in range(nb):
            j_max = min(i + max_dist_bins, nb - 1)
            for j in range(i, j_max + 1):
                f.write(f"{c}\t{i}\t{j}\n")


def cool_to_feature_vector_balanced(cool_uri):
    clr = cooler.Cooler(cool_uri)
    x = np.zeros(total_features, dtype=np.float64)
    has_balance = "weight" in clr.bins().columns
    for chrom in chroms:
        nb = chrom_nbins[chrom]
        if has_balance:
            mat = clr.matrix(balance=True).fetch(chrom)
        else:
            mat = clr.matrix(balance=False).fetch(chrom)
        mat = np.nan_to_num(mat, 0.0)
        mat = np.maximum(mat, 0.0)
        off = chrom_feat_offset[chrom]
        row_off = chrom_row_offsets[chrom]
        for i in range(min(nb, mat.shape[0])):
            j_max = min(i + max_dist_bins, nb - 1, mat.shape[1] - 1)
            for j in range(i, j_max + 1):
                feat_idx = off + row_off[i] + (j - i)
                x[feat_idx] = mat[i, j]
    return x


samples = []
for f in sorted(glob.glob("data/4dn_thymocyte/*.mcool")):
    name = os.path.basename(f).replace(".mcool", "")
    samples.append((name, f + "::resolutions/100000", "thymocyte"))
samples.append(("mTEC_WT", "data/mtec/processed/mtec_hic.mcool::resolutions/100000", "mTEC"))
samples.append(("Bonev_ES", "data/brain/processed/brain_hic_es.mcool::resolutions/100000", "brain"))
samples.append(("Bonev_NPC", "data/brain/processed/brain_hic_npc.mcool::resolutions/100000", "brain"))
samples.append(("Bonev_CN", "data/brain/processed/brain_hic_cn.mcool::resolutions/100000", "brain"))

print(f"[info] {len(samples)} samples", flush=True)

X = np.zeros((total_features, len(samples)), dtype=np.float64)
sample_names, sample_labels = [], []
for si, (name, uri, label) in enumerate(samples):
    print(f"  [{si+1}/{len(samples)}] {name}...", end="", flush=True)
    x = cool_to_feature_vector_balanced(uri)
    nnz = np.sum(x > 0)
    print(f" nnz={nnz}, mean={x[x>0].mean():.6f}", flush=True)
    X[:, si] = x
    sample_names.append(name)
    sample_labels.append(label)

with open(os.path.join(out_dir, "cells.tsv"), "w") as f:
    for n, l in zip(sample_names, sample_labels):
        f.write(f"{n}\t{l}\n")

np.save(os.path.join(out_dir, "X_bulk.npy"), X)
print(f"\n[saved] X_bulk.npy ({X.shape})", flush=True)
