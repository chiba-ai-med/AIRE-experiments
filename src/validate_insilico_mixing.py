#!/usr/bin/env python3
"""In-silico mixing validation: create synthetic bulk Hi-C mixtures, fold-in, check recovery."""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb_oe"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8

print("[info] loading...", flush=True)
X = np.load(os.path.join(data_dir, "X_bulk_oe.npy"))
W = np.loadtxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), delimiter="\t")
H = np.loadtxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), delimiter="\t")
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["sample", "label"])
N, M = X.shape
sample_names = cells["sample"].tolist()
print(f"  X: {N}x{M}, W: {W.shape}, K={K}", flush=True)

# Define mixing experiments — exclude mTEC (production data)
representatives = {
    "DN": "DN_thymocyte",
    "DP": "DP_thymocyte",
    "CD4_SP": "CD4_SP_thymocyte",
    "CD8_SP": "CD8_SP_thymocyte",
    "Treg": "Mature_Treg_GFP",
    "TregPre": "Treg_precursor",
    "Brain_ES": "Bonev_ES",
    "Brain_NPC": "Bonev_NPC",
    "Brain_CN": "Bonev_CN",
}

rep_idx = {name: sample_names.index(sname) for name, sname in representatives.items()}
rep_vectors = {name: X[:, idx] for name, idx in rep_idx.items()}


def foldin_mu(W, x_new, n_iter=300, eps=1e-10):
    """Fold-in: fix W, estimate h for a single sample."""
    K = W.shape[1]
    np.random.seed(0)
    h = np.random.rand(K) + 1e-6
    WtW = W.T @ W
    Wtx = W.T @ x_new
    for _ in range(n_iter):
        h = h * Wtx / (WtW @ h + eps)
    return h / h.sum()


# Experiment 1: two-component mixtures at various ratios
print("\n[exp1] Two-component mixtures", flush=True)
mix_pairs = [
    ("DN", "Brain_ES"),       # cross-lineage: thymocyte vs brain
    ("CD4_SP", "Brain_CN"),   # cross-lineage: thymocyte vs brain
    ("DN", "DP"),             # within thymocyte: early stages
    ("CD4_SP", "Treg"),       # within thymocyte: close lineages
    ("Brain_ES", "Brain_CN"), # within brain: ES vs CN
    ("CD8_SP", "TregPre"),    # within thymocyte: distant stages
]
ratios = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]

results_2comp = []
for c1, c2 in mix_pairs:
    for r in ratios:
        x_mix = r * rep_vectors[c1] + (1 - r) * rep_vectors[c2]
        h_est = foldin_mu(W, x_mix)
        # Find which NMF components correspond to c1, c2
        h_c1 = h_est[np.argmax(H[:, rep_idx[c1]])]
        h_c2 = h_est[np.argmax(H[:, rep_idx[c2]])]
        results_2comp.append({
            "pair": f"{c1} vs {c2}",
            "true_frac_c1": r,
            "est_frac_c1": h_c1,
            "est_frac_c2": h_c2,
            "c1": c1, "c2": c2,
        })

df_2comp = pd.DataFrame(results_2comp)

# Experiment 2: three-component mixtures (DN + CD4_SP + Brain_ES)
print("[exp2] Three-component mixtures (DN + CD4_SP + Brain_ES)", flush=True)
results_3comp = []
tri_names = ["DN", "CD4_SP", "Brain_ES"]
for r0 in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    for r1 in np.arange(0, 1.01 - r0, 0.2):
        r2 = 1.0 - r0 - r1
        if r2 < -0.01:
            continue
        r2 = max(r2, 0)
        x_mix = r0 * rep_vectors[tri_names[0]] + r1 * rep_vectors[tri_names[1]] + r2 * rep_vectors[tri_names[2]]
        h_est = foldin_mu(W, x_mix)
        comps = [np.argmax(H[:, rep_idx[n]]) for n in tri_names]
        results_3comp.append({
            "true_DN": r0, "true_CD4SP": r1, "true_Brain": r2,
            "est_DN": h_est[comps[0]], "est_CD4SP": h_est[comps[1]], "est_Brain": h_est[comps[2]],
        })

df_3comp = pd.DataFrame(results_3comp)

# Print results
print("\n=== Two-component mixing ===", flush=True)
for pair in df_2comp["pair"].unique():
    sub = df_2comp[df_2comp["pair"] == pair]
    corr = np.corrcoef(sub["true_frac_c1"], sub["est_frac_c1"])[0, 1]
    rmse = np.sqrt(np.mean((sub["true_frac_c1"] - sub["est_frac_c1"]) ** 2))
    print(f"  {pair:25s}: r={corr:.3f}, RMSE={rmse:.3f}", flush=True)

print("\n=== Three-component mixing (DN + CD4_SP + Brain_ES) ===", flush=True)
for col in ["DN", "CD4SP", "Brain"]:
    corr = np.corrcoef(df_3comp[f"true_{col}"], df_3comp[f"est_{col}"])[0, 1]
    rmse = np.sqrt(np.mean((df_3comp[f"true_{col}"] - df_3comp[f"est_{col}"]) ** 2))
    print(f"  {col:10s}: r={corr:.3f}, RMSE={rmse:.3f}", flush=True)

# Plot
n_pairs = len(mix_pairs)
n_cols = 3
n_rows_2comp = (n_pairs + n_cols - 1) // n_cols
fig, axes = plt.subplots(n_rows_2comp + 1, n_cols, figsize=(18, 6 * (n_rows_2comp + 1)))

# Rows for two-component scatter
for pi, (c1, c2) in enumerate(mix_pairs):
    ax = axes[pi // n_cols, pi % n_cols]
    sub = df_2comp[df_2comp["pair"] == f"{c1} vs {c2}"]
    ax.scatter(sub["true_frac_c1"], sub["est_frac_c1"], s=80, color="steelblue", zorder=5)
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    corr = np.corrcoef(sub["true_frac_c1"], sub["est_frac_c1"])[0, 1]
    rmse = np.sqrt(np.mean((sub["true_frac_c1"] - sub["est_frac_c1"]) ** 2))
    ax.set_title(f"{c1} vs {c2}\nr={corr:.3f}, RMSE={rmse:.3f}")
    ax.set_xlabel(f"True {c1} fraction")
    ax.set_ylabel(f"Estimated {c1} fraction")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")

# Bottom row: three-component
for pi, col in enumerate(["DN", "CD4SP", "Brain"]):
    ax = axes[n_rows_2comp, pi]
    ax.scatter(df_3comp[f"true_{col}"], df_3comp[f"est_{col}"], s=80, color="coral", zorder=5)
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    corr = np.corrcoef(df_3comp[f"true_{col}"], df_3comp[f"est_{col}"])[0, 1]
    rmse = np.sqrt(np.mean((df_3comp[f"true_{col}"] - df_3comp[f"est_{col}"]) ** 2))
    ax.set_title(f"3-comp: {col}\nr={corr:.3f}, RMSE={rmse:.3f}")
    ax.set_xlabel(f"True {col} fraction")
    ax.set_ylabel(f"Estimated {col} fraction")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")

plt.suptitle("In-silico mixing validation (O/E, K=8, no mTEC)", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig("plot/diagnostic_insilico_mixing.png", dpi=150)
print("\n[saved] plot/diagnostic_insilico_mixing.png")
