#!/usr/bin/env python3
"""Compare NMF component TAD boundaries against Bonev et al. reference (ES/NPC/CN)."""
import sys, os
import numpy as np
import pandas as pd
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out_dir = sys.argv[1] if len(sys.argv) > 1 else "output/tad_calling"
resolution = 100000

# --- Step 1: Compute reference insulation scores ---
ref_cools = {
    "ES": "data/brain/processed/brain_hic_es.mcool::resolutions/100000",
    "NPC": "data/brain/processed/brain_hic_npc.mcool::resolutions/100000",
    "CN": "data/brain/processed/brain_hic_cn.mcool::resolutions/100000",
}
window_sizes = [3 * resolution, 5 * resolution, 10 * resolution]

ref_boundaries = {}
print("[step1] Computing reference insulation scores...", flush=True)
for name, uri in ref_cools.items():
    clr = cooler.Cooler(uri)
    ins = cooltools.insulation(clr, window_sizes, ignore_diags=2)
    ins.to_csv(os.path.join(out_dir, f"insulation_ref_{name}.tsv"), sep="\t", index=False)
    for ws in window_sizes:
        ws_kb = ws // 1000
        bcol = f"is_boundary_{ws}"
        if bcol in ins.columns:
            bnd = ins[ins[bcol] == True][["chrom", "start", "end"]].copy()
            key = f"{name}_{ws_kb}kb"
            ref_boundaries[key] = bnd
            print(f"  {key}: {len(bnd)} boundaries", flush=True)

# --- Step 2: Load NMF component boundaries ---
nmf_boundaries = {}
print("\n[step2] Loading NMF component boundaries...", flush=True)
for k in range(1, 6):
    f = os.path.join(out_dir, f"insulation_nmf_comp{k}.tsv")
    if not os.path.exists(f):
        continue
    ins = pd.read_csv(f, sep="\t")
    for ws in window_sizes:
        ws_kb = ws // 1000
        bcol = f"is_boundary_{ws}"
        if bcol in ins.columns:
            bnd = ins[ins[bcol] == True][["chrom", "start", "end"]].copy()
            key = f"comp{k}_{ws_kb}kb"
            nmf_boundaries[key] = bnd
            print(f"  {key}: {len(bnd)} boundaries", flush=True)

# --- Step 3: Compute Jaccard overlap ---
def boundary_overlap(bnd1, bnd2, slop=1):
    """Fraction of bnd1 boundaries within `slop` bins of any bnd2 boundary."""
    if len(bnd1) == 0 or len(bnd2) == 0:
        return 0.0
    hits = 0
    for _, row in bnd1.iterrows():
        chrom, start = row["chrom"], row["start"]
        bnd2_chr = bnd2[bnd2["chrom"] == chrom]
        if len(bnd2_chr) == 0:
            continue
        min_dist = np.abs(bnd2_chr["start"].values - start).min()
        if min_dist <= slop * resolution:
            hits += 1
    return hits / len(bnd1)

ws_focus = 500  # 500kb window
print(f"\n[step3] Boundary overlap (500kb window, slop=1 bin)...", flush=True)
overlap_matrix = np.zeros((5, 3))
ref_names = ["ES", "NPC", "CN"]
for ki in range(5):
    nmf_key = f"comp{ki+1}_{ws_focus}kb"
    if nmf_key not in nmf_boundaries:
        continue
    for ri, rname in enumerate(ref_names):
        ref_key = f"{rname}_{ws_focus}kb"
        if ref_key not in ref_boundaries:
            continue
        frac = boundary_overlap(nmf_boundaries[nmf_key], ref_boundaries[ref_key])
        overlap_matrix[ki, ri] = frac
        print(f"  comp{ki+1} vs {rname}: {frac:.3f} ({int(frac*len(nmf_boundaries[nmf_key]))}/{len(nmf_boundaries[nmf_key])})", flush=True)

# Reverse: what fraction of reference boundaries are captured by each component
print(f"\n  Reverse (ref captured by NMF comp):", flush=True)
overlap_rev = np.zeros((5, 3))
for ki in range(5):
    nmf_key = f"comp{ki+1}_{ws_focus}kb"
    if nmf_key not in nmf_boundaries:
        continue
    for ri, rname in enumerate(ref_names):
        ref_key = f"{rname}_{ws_focus}kb"
        if ref_key not in ref_boundaries:
            continue
        frac = boundary_overlap(ref_boundaries[ref_key], nmf_boundaries[nmf_key])
        overlap_rev[ki, ri] = frac
        print(f"  {rname} captured by comp{ki+1}: {frac:.3f}", flush=True)

# --- Step 4: Insulation profile correlation ---
print(f"\n[step4] Insulation profile correlation (500kb window)...", flush=True)
ws_col = f"log2_insulation_score_{ws_focus * 1000}"
ref_profiles = {}
for rname in ref_names:
    ins = pd.read_csv(os.path.join(out_dir, f"insulation_ref_{rname}.tsv"), sep="\t")
    ref_profiles[rname] = ins.set_index(["chrom", "start"])[ws_col]

nmf_profiles = {}
for k in range(1, 6):
    ins = pd.read_csv(os.path.join(out_dir, f"insulation_nmf_comp{k}.tsv"), sep="\t")
    nmf_profiles[f"comp{k}"] = ins.set_index(["chrom", "start"])[ws_col]

corr_matrix = np.zeros((5, 3))
for ki in range(5):
    for ri, rname in enumerate(ref_names):
        merged = pd.DataFrame({
            "nmf": nmf_profiles[f"comp{ki+1}"],
            "ref": ref_profiles[rname]
        }).dropna()
        if len(merged) > 10:
            corr_matrix[ki, ri] = merged["nmf"].corr(merged["ref"])

print("Pearson correlation of insulation profiles:")
df_corr = pd.DataFrame(corr_matrix, index=[f"comp{k+1}" for k in range(5)], columns=ref_names)
print(df_corr.round(3).to_string())

# --- Step 5: Plot ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Overlap heatmap
ax = axes[0]
im = ax.imshow(overlap_matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(3)); ax.set_xticklabels(ref_names)
ax.set_yticks(range(5)); ax.set_yticklabels([f"comp{k+1}" for k in range(5)])
for i in range(5):
    for j in range(3):
        ax.text(j, i, f"{overlap_matrix[i,j]:.2f}", ha="center", va="center", fontsize=10)
ax.set_title("NMF boundary overlap with ref")
plt.colorbar(im, ax=ax, fraction=0.046)

# Reverse overlap
ax = axes[1]
im = ax.imshow(overlap_rev, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(3)); ax.set_xticklabels(ref_names)
ax.set_yticks(range(5)); ax.set_yticklabels([f"comp{k+1}" for k in range(5)])
for i in range(5):
    for j in range(3):
        ax.text(j, i, f"{overlap_rev[i,j]:.2f}", ha="center", va="center", fontsize=10)
ax.set_title("Ref boundaries captured by NMF comp")
plt.colorbar(im, ax=ax, fraction=0.046)

# Insulation correlation
ax = axes[2]
im = ax.imshow(corr_matrix, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)
ax.set_xticks(range(3)); ax.set_xticklabels(ref_names)
ax.set_yticks(range(5)); ax.set_yticklabels([f"comp{k+1}" for k in range(5)])
for i in range(5):
    for j in range(3):
        ax.text(j, i, f"{corr_matrix[i,j]:.3f}", ha="center", va="center", fontsize=10)
ax.set_title("Insulation profile correlation")
plt.colorbar(im, ax=ax, fraction=0.046)

plt.tight_layout()
plt.savefig("plot/diagnostic_tad_overlap_ref.png", dpi=150)
print("\n[saved] plot/diagnostic_tad_overlap_ref.png")
