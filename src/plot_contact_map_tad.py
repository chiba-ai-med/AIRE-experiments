#!/usr/bin/env python3
"""Visualize Hi-C contact maps with TAD boundaries highlighted."""
import sys, os
import numpy as np
import pandas as pd
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm

out_dir = "plot"
os.makedirs(out_dir, exist_ok=True)

resolution = 100000
window_bp = 500000

samples = {
    "mTEC": "data/mtec/processed/mtec_hic.mcool::resolutions/100000",
    "Bonev_ES": "data/brain/processed/brain_hic_es.mcool::resolutions/100000",
    "DN_thymocyte": "data/4dn_thymocyte/DN_thymocyte.mcool::resolutions/100000",
}

regions = [
    ("chr10", 17_000_000, 27_000_000, "Aire locus (chr10:17-27Mb)"),
    ("chr1", 50_000_000, 70_000_000, "chr1:50-70Mb"),
    ("chr6", 50_000_000, 65_000_000, "chr6:50-65Mb"),
]

print("[info] Computing insulation scores...", flush=True)
ins_data = {}
tad_blocks = {}
for name, uri in samples.items():
    clr = cooler.Cooler(uri)
    try:
        ins = cooltools.insulation(clr, [window_bp], ignore_diags=2)
    except Exception:
        ins = cooltools.insulation(clr, [window_bp], ignore_diags=2, clr_weight_name=None)
    ins_data[name] = ins
    bcol = "is_boundary_500000"
    nb = ins[bcol].sum() if bcol in ins.columns else 0
    print(f"  {name}: {nb} boundaries", flush=True)

print("\n[info] Plotting contact maps with TAD boundaries...", flush=True)
n_samples = len(samples)
n_regions = len(regions)
fig, axes = plt.subplots(n_regions, n_samples, figsize=(7 * n_samples, 6 * n_regions))

for ri, (chrom, start, end, title) in enumerate(regions):
    region_str = f"{chrom}:{start}-{end}"
    for si, (name, uri) in enumerate(samples.items()):
        ax = axes[ri, si]
        clr = cooler.Cooler(uri)
        try:
            mat = clr.matrix(balance=True).fetch(region_str)
        except Exception:
            mat = clr.matrix(balance=False).fetch(region_str)
        mat = np.nan_to_num(mat, 0.0)
        np.fill_diagonal(mat, 0)

        n_bins = mat.shape[0]
        extent = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]
        vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        vmin = max(vmax * 1e-4, 1e-10)
        im = ax.imshow(mat, cmap="YlOrRd", norm=LogNorm(vmin=vmin, vmax=vmax),
                       extent=extent, interpolation="none")

        ins = ins_data[name]
        bcol = "is_boundary_500000"
        if bcol in ins.columns:
            bounds = ins[(ins["chrom"] == chrom) & (ins[bcol] == True)]
            bounds = bounds[(bounds["start"] >= start) & (bounds["end"] <= end)]
            for _, row in bounds.iterrows():
                pos_mb = row["start"] / 1e6
                ax.axvline(pos_mb, color="blue", lw=0.5, alpha=0.6)
                ax.axhline(pos_mb, color="blue", lw=0.5, alpha=0.6)

            bound_starts = sorted(bounds["start"].values.tolist())
            all_starts = [start] + bound_starts + [end]
            for ti in range(len(all_starts) - 1):
                t_start = all_starts[ti] / 1e6
                t_end = all_starts[ti + 1] / 1e6
                t_size = t_end - t_start
                rect = Rectangle((t_start, t_start), t_size, t_size,
                                 linewidth=1.5, edgecolor="blue", facecolor="none",
                                 linestyle="-", alpha=0.7)
                ax.add_patch(rect)

        ax.set_xlim(start / 1e6, end / 1e6)
        ax.set_ylim(end / 1e6, start / 1e6)
        if ri == 0:
            ax.set_title(f"{name}\n{title}", fontsize=10)
        else:
            ax.set_title(title, fontsize=10)
        ax.set_xlabel("Position (Mb)")
        ax.set_ylabel("Position (Mb)")

plt.suptitle("Hi-C contact maps with TAD boundaries (blue)", fontsize=14, y=1.01)
plt.tight_layout()
out_png = os.path.join(out_dir, "contact_map_tad_boundaries.png")
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\n[saved] {out_png}", flush=True)

# --- Insulation profile comparison plot ---
print("\n[info] Plotting insulation profiles...", flush=True)
ins_col = "log2_insulation_score_500000"
fig2, axes2 = plt.subplots(n_regions, 1, figsize=(16, 4 * n_regions))

colors = {"mTEC": "tomato", "Bonev_ES": "steelblue", "DN_thymocyte": "seagreen"}
for ri, (chrom, start, end, title) in enumerate(regions):
    ax = axes2[ri]
    for name in samples:
        ins = ins_data[name]
        sub = ins[(ins["chrom"] == chrom) & (ins["start"] >= start) & (ins["start"] < end)]
        x_mb = sub["start"].values / 1e6
        y = sub[ins_col].values
        ax.plot(x_mb, y, color=colors[name], lw=1, alpha=0.8, label=name)

        bcol = "is_boundary_500000"
        if bcol in ins.columns:
            bounds = sub[sub[bcol] == True]
            ax.scatter(bounds["start"].values / 1e6, bounds[ins_col].values,
                       color=colors[name], s=20, zorder=5, marker="v")

    ax.set_title(title)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("log2 insulation score")
    ax.legend(fontsize=8)
    ax.axhline(0, color="grey", lw=0.5, ls="--")

plt.tight_layout()
out_png2 = os.path.join(out_dir, "insulation_profile_comparison.png")
plt.savefig(out_png2, dpi=150)
print(f"[saved] {out_png2}", flush=True)

print("\n[done]", flush=True)
