#!/usr/bin/env python3
"""TAD calling on bulk Hi-C NMF components. Compare comp4 (mTEC) vs actual mTEC bulk Hi-C."""
import sys, os
import numpy as np
import pandas as pd
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data/bulk_nmf_100kb_oe"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
resolution = 100000
out_dir = "output/tad_calling_bulk"
os.makedirs(out_dir, exist_ok=True)

print("[info] loading NMF results...", flush=True)
W = np.loadtxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), delimiter="\t")
H = np.loadtxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), delimiter="\t")
feat = pd.read_csv(os.path.join(data_dir, "features.tsv"), sep="\t",
                   header=None, names=["chrom", "bin_i", "bin_j"])
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["sample", "label"])
print(f"  W: {W.shape}, H: {H.shape}, features: {len(feat)}", flush=True)

sample_names = cells["sample"].tolist()
mtec_idx = sample_names.index("mTEC_WT")
mtec_comp = np.argmax(H[:, mtec_idx])
print(f"  mTEC dominant component: comp{mtec_comp+1} (loading={H[mtec_comp, mtec_idx]:.3f})")

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
    if n == "X": return 100
    if n == "Y": return 101
    return int(n)

chroms = sorted(chrom_sizes.keys(), key=chrom_sort_key)

bins_list = []
for c in chroms:
    n_bins = (chrom_sizes[c] + resolution - 1) // resolution
    for b in range(n_bins):
        start = b * resolution
        end = min(start + resolution, chrom_sizes[c])
        bins_list.append((c, start, end))
bins_df = pd.DataFrame(bins_list, columns=["chrom", "start", "end"])

chrom_bin_offset = {}
offset = 0
for c in chroms:
    chrom_bin_offset[c] = offset
    offset += (chrom_sizes[c] + resolution - 1) // resolution


def w_column_to_cool(w_col, cool_path, label):
    """Convert W column to .cool file."""
    bin1_ids, bin2_ids, counts = [], [], []
    for chrom in chroms:
        mask = feat["chrom"] == chrom
        if not mask.any():
            continue
        feat_chr = feat[mask]
        w_chr = w_col[mask.values]
        off = chrom_bin_offset[chrom]
        pos = w_chr > 0
        if not pos.any():
            continue
        feat_pos = feat_chr[pos]
        w_pos = w_chr[pos]
        bi = (feat_pos["bin_i"].values + off).astype(int)
        bj = (feat_pos["bin_j"].values + off).astype(int)
        bin1_ids.extend(bi)
        bin2_ids.extend(bj)
        counts.extend(w_pos)
    pixels = pd.DataFrame({"bin1_id": bin1_ids, "bin2_id": bin2_ids, "count": counts})
    cooler.create_cooler(cool_path, bins_df, pixels, ordered=True, columns=["count"])
    print(f"  [saved] {cool_path} ({len(pixels)} pixels)", flush=True)
    return cool_path


def run_insulation(cool_path, label):
    """Run cooltools insulation and return DataFrame."""
    clr = cooler.Cooler(cool_path)
    window_sizes = [500000]
    ins = cooltools.insulation(clr, window_sizes, ignore_diags=2, clr_weight_name=None)
    bcol = "is_boundary_500000"
    n_bound = ins[bcol].sum() if bcol in ins.columns else 0
    print(f"  {label}: {n_bound} boundaries @ 500kb window", flush=True)
    return ins


# Step 1: Create .cool for each NMF component
print("\n[step1] NMF components -> .cool files", flush=True)
for k in range(K):
    cool_path = os.path.join(out_dir, f"nmf_comp{k+1}.cool")
    w_column_to_cool(W[:, k], cool_path, f"comp{k+1}")

# Step 2: Insulation score for all components
print("\n[step2] Insulation scores for NMF components", flush=True)
ins_nmf = {}
for k in range(K):
    cool_path = os.path.join(out_dir, f"nmf_comp{k+1}.cool")
    ins_nmf[k] = run_insulation(cool_path, f"comp{k+1}")

# Step 3: Reference insulation from actual mTEC bulk Hi-C
print("\n[step3] Reference insulation from mTEC bulk Hi-C", flush=True)
mtec_cool = "data/mtec/processed/mtec_hic.mcool::resolutions/100000"
clr_ref = cooler.Cooler(mtec_cool)
ins_ref = cooltools.insulation(clr_ref, [500000], ignore_diags=2)
bcol = "is_boundary_500000"
n_ref = ins_ref[bcol].sum() if bcol in ins_ref.columns else 0
print(f"  mTEC reference: {n_ref} boundaries @ 500kb window", flush=True)

# Step 4: Compare comp4 vs mTEC reference per chromosome
print("\n[step4] Boundary comparison: NMF comp vs mTEC reference", flush=True)
ins_col = "log2_insulation_score_500000"

ref_boundaries = set()
if bcol in ins_ref.columns:
    ref_b = ins_ref[ins_ref[bcol] == True]
    for _, row in ref_b.iterrows():
        ref_boundaries.add((row["chrom"], int(row["start"])))

for k in range(K):
    ins_k = ins_nmf[k]
    if bcol not in ins_k.columns:
        continue
    nmf_b = ins_k[ins_k[bcol] == True]
    nmf_boundaries = set()
    for _, row in nmf_b.iterrows():
        nmf_boundaries.add((row["chrom"], int(row["start"])))

    exact = len(ref_boundaries & nmf_boundaries)
    within_1bin = 0
    for rc, rp in ref_boundaries:
        for delta in [-resolution, 0, resolution]:
            if (rc, rp + delta) in nmf_boundaries:
                within_1bin += 1
                break
    jaccard_exact = exact / max(len(ref_boundaries | nmf_boundaries), 1)
    recall_1bin = within_1bin / max(len(ref_boundaries), 1)
    print(f"  comp{k+1}: {len(nmf_boundaries)} NMF boundaries, "
          f"exact overlap={exact}, Jaccard={jaccard_exact:.3f}, "
          f"recall@1bin={recall_1bin:.3f}", flush=True)

# Step 5: Insulation profile correlation (comp4 vs reference)
print(f"\n[step5] Insulation profile correlation per chrom (comp{mtec_comp+1} vs ref)", flush=True)
ins_comp = ins_nmf[mtec_comp]
corrs = []
for chrom in chroms:
    ref_chr = ins_ref[ins_ref["chrom"] == chrom][ins_col].values
    nmf_chr = ins_comp[ins_comp["chrom"] == chrom][ins_col].values
    n = min(len(ref_chr), len(nmf_chr))
    if n < 10:
        continue
    ref_v = ref_chr[:n]
    nmf_v = nmf_chr[:n]
    valid = np.isfinite(ref_v) & np.isfinite(nmf_v)
    if valid.sum() < 10:
        continue
    r = np.corrcoef(ref_v[valid], nmf_v[valid])[0, 1]
    corrs.append({"chrom": chrom, "r": r, "n_bins": int(valid.sum())})
    print(f"  {chrom:6s}: r={r:.3f} (n={valid.sum()})", flush=True)

df_corr = pd.DataFrame(corrs)
mean_r = df_corr["r"].mean()
print(f"\n  Mean correlation across chroms: r={mean_r:.3f}", flush=True)

# Step 6: Plot insulation profiles for selected chroms
print("\n[step6] Plotting insulation profiles...", flush=True)
plot_chroms = ["chr1", "chr10", "chr19"]
fig, axes = plt.subplots(len(plot_chroms), 1, figsize=(16, 4 * len(plot_chroms)))
for pi, chrom in enumerate(plot_chroms):
    ax = axes[pi]
    ref_chr = ins_ref[ins_ref["chrom"] == chrom]
    nmf_chr = ins_comp[ins_comp["chrom"] == chrom]
    n = min(len(ref_chr), len(nmf_chr))
    x_mb = ref_chr["start"].values[:n] / 1e6
    ref_v = ref_chr[ins_col].values[:n]
    nmf_v = nmf_chr[ins_col].values[:n]
    ax.plot(x_mb, ref_v, color="black", lw=0.8, alpha=0.8, label="mTEC bulk (ref)")
    ax.plot(x_mb, nmf_v, color="tomato", lw=0.8, alpha=0.8, label=f"NMF comp{mtec_comp+1}")
    valid = np.isfinite(ref_v) & np.isfinite(nmf_v)
    if valid.sum() > 10:
        r = np.corrcoef(ref_v[valid], nmf_v[valid])[0, 1]
        ax.set_title(f"{chrom} — insulation score (r={r:.3f})")
    else:
        ax.set_title(f"{chrom} — insulation score")
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("log2 insulation score")
    ax.legend(loc="upper right", fontsize=8)

plt.tight_layout()
out_png = "plot/diagnostic_bulk_nmf_tad.png"
plt.savefig(out_png, dpi=150)
print(f"\n[saved] {out_png}", flush=True)

# Step 7: Heatmap of correlation across all components vs reference
fig2, ax2 = plt.subplots(figsize=(10, 5))
corr_matrix = np.zeros((K, len(chroms)))
for k in range(K):
    ins_k = ins_nmf[k]
    for ci, chrom in enumerate(chroms):
        ref_chr = ins_ref[ins_ref["chrom"] == chrom][ins_col].values
        nmf_chr = ins_k[ins_k["chrom"] == chrom][ins_col].values
        n = min(len(ref_chr), len(nmf_chr))
        if n < 10:
            corr_matrix[k, ci] = np.nan
            continue
        valid = np.isfinite(ref_chr[:n]) & np.isfinite(nmf_chr[:n])
        if valid.sum() < 10:
            corr_matrix[k, ci] = np.nan
            continue
        corr_matrix[k, ci] = np.corrcoef(ref_chr[:n][valid], nmf_chr[:n][valid])[0, 1]

im = ax2.imshow(corr_matrix, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
ax2.set_yticks(range(K))
ax2.set_yticklabels([f"comp{k+1}" for k in range(K)])
ax2.set_xticks(range(len(chroms)))
ax2.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=8)
ax2.set_xlabel("Chromosome")
ax2.set_ylabel("NMF component")
ax2.set_title("Insulation profile correlation: NMF components vs mTEC bulk Hi-C reference")
plt.colorbar(im, ax=ax2, label="Pearson r")
plt.tight_layout()
out_png2 = "plot/diagnostic_bulk_nmf_tad_heatmap.png"
plt.savefig(out_png2, dpi=150)
print(f"[saved] {out_png2}", flush=True)

print("\n[done]", flush=True)
