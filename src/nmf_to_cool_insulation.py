#!/usr/bin/env python3
"""Convert NMF components to .cool files and run cooltools insulation score.
Threshold: left-tail trimming (cut where density drops to 10% of mode on log10 scale).
"""
import sys, os
import numpy as np
import pandas as pd
import cooler
import cooltools

data_dir = sys.argv[1]
out_dir = sys.argv[2] if len(sys.argv) > 2 else "output/tad_calling"
K = int(sys.argv[3]) if len(sys.argv) > 3 else 5
resolution = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
os.makedirs(out_dir, exist_ok=True)

print(f"[info] loading NMF results from {data_dir}...", flush=True)
W = np.loadtxt(os.path.join(data_dir, f"W_nmf_k{K}.tsv"), delimiter="\t")
H = np.loadtxt(os.path.join(data_dir, f"H_nmf_k{K}.tsv"), delimiter="\t")
feat = pd.read_csv(os.path.join(data_dir, "features.tsv"), sep="\t",
                   header=None, names=["chrom", "bin_i", "bin_j"])
cells = pd.read_csv(os.path.join(data_dir, "cells.tsv"), sep="\t",
                    header=None, names=["cell", "label"])
print(f"  W: {W.shape}, H: {H.shape}, features: {len(feat)}", flush=True)

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


def left_tail_threshold(values, frac=0.10, n_bins=300):
    """Find threshold where histogram density drops to frac of peak, on left side of mode."""
    pos = values[values > 1e-10]
    if len(pos) < 100:
        return 0.0
    log_vals = np.log10(pos)
    counts, edges = np.histogram(log_vals, bins=n_bins)
    centers = (edges[:-1] + edges[1:]) / 2
    mode_idx = np.argmax(counts)
    peak_count = counts[mode_idx]
    cutoff = centers[0]
    for i in range(mode_idx, -1, -1):
        if counts[i] < peak_count * frac:
            cutoff = centers[i]
            break
    return 10 ** cutoff


def w_column_to_pixels(w_col, threshold):
    """Convert a W column to sparse pixel DataFrame with thresholding."""
    bin1_ids, bin2_ids, counts = [], [], []
    for chrom in chroms:
        mask = feat["chrom"] == chrom
        if not mask.any():
            continue
        feat_chr = feat[mask]
        w_chr = w_col[mask.values]
        off = chrom_bin_offset[chrom]
        above = w_chr >= threshold
        if not above.any():
            continue
        feat_above = feat_chr[above]
        w_above = w_chr[above]
        bi = (feat_above["bin_i"].values + off).astype(int)
        bj = (feat_above["bin_j"].values + off).astype(int)
        bin1_ids.extend(bi)
        bin2_ids.extend(bj)
        counts.extend(w_above)
    return pd.DataFrame({"bin1_id": bin1_ids, "bin2_id": bin2_ids, "count": counts})


def create_cool(pixels_df, cool_path, label):
    """Create a .cool file from pixels DataFrame."""
    if len(pixels_df) == 0:
        print(f"  WARNING: no pixels for {label}, skipping", flush=True)
        return None
    cooler.create_cooler(cool_path, bins_df, pixels_df, ordered=True, columns=["count"])
    print(f"  [saved] {cool_path} ({len(pixels_df)} pixels)", flush=True)
    return cool_path


# --- Per-component cool files ---
print(f"\n[step1] Left-tail thresholding + .cool for {K} components...", flush=True)
for k in range(K):
    w_col = W[:, k]
    thresh = left_tail_threshold(w_col)
    n_above = np.sum(w_col >= thresh)
    print(f"  comp{k+1}: threshold={thresh:.6f}, "
          f"pixels={n_above}/{len(w_col)} ({100*n_above/len(w_col):.1f}%)", flush=True)
    pixels = w_column_to_pixels(w_col, thresh)
    cool_path = os.path.join(out_dir, f"nmf_comp{k+1}.cool")
    create_cool(pixels, cool_path, f"comp{k+1}")

# --- Per-celltype reconstruction ---
ct_names = sorted(cells["label"].unique())
print(f"\n[step2] Per-celltype reconstruction (left-tail threshold)...", flush=True)
for ct in ct_names:
    h_mean = H[:, cells["label"].values == ct].mean(axis=1)
    w_recon = W @ h_mean
    thresh = left_tail_threshold(w_recon)
    n_above = np.sum(w_recon >= thresh)
    ct_safe = ct.replace(" ", "_")
    print(f"  {ct}: threshold={thresh:.6f}, pixels={n_above}/{len(w_recon)} ({100*n_above/len(w_recon):.1f}%)", flush=True)
    pixels = w_column_to_pixels(w_recon, thresh)
    cool_path = os.path.join(out_dir, f"nmf_recon_{ct_safe}.cool")
    create_cool(pixels, cool_path, ct)

# --- Insulation score ---
print(f"\n[step3] Running cooltools insulation score...", flush=True)
cool_files = sorted([f for f in os.listdir(out_dir) if f.endswith(".cool")])
window_sizes = [3 * resolution, 5 * resolution, 10 * resolution]

all_insulation = {}
for cf in cool_files:
    cool_path = os.path.join(out_dir, cf)
    label = cf.replace(".cool", "")
    try:
        clr = cooler.Cooler(cool_path)
        ins = cooltools.insulation(clr, window_sizes, ignore_diags=2, clr_weight_name=None)
        ins.to_csv(os.path.join(out_dir, f"insulation_{label}.tsv"), sep="\t", index=False)
        all_insulation[label] = ins

        for ws in window_sizes:
            ws_kb = ws // 1000
            bcol = f"is_boundary_{ws}"
            if bcol in ins.columns:
                boundaries = ins[ins[bcol] == True]
                print(f"  {label} @ {ws_kb}kb window: {len(boundaries)} boundaries", flush=True)
    except Exception as e:
        print(f"  ERROR on {label}: {e}", flush=True)

# --- Save boundary summary ---
print(f"\n[step4] Boundary summary...", flush=True)
summary_rows = []
for label, ins in all_insulation.items():
    for ws in window_sizes:
        ws_kb = ws // 1000
        bcol = f"is_boundary_{ws}"
        if bcol in ins.columns:
            boundaries = ins[ins[bcol] == True]
            for _, row in boundaries.iterrows():
                summary_rows.append({
                    "label": label,
                    "window_kb": ws_kb,
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "boundary_strength": row.get(f"boundary_strength_{ws}", np.nan),
                })

if summary_rows:
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(out_dir, "tad_boundaries_all.tsv"), sep="\t", index=False)
    print(f"[saved] tad_boundaries_all.tsv ({len(summary)} boundaries total)", flush=True)

    for label in sorted(summary["label"].unique()):
        sub = summary[summary["label"] == label]
        for ws_kb in sorted(sub["window_kb"].unique()):
            n = len(sub[sub["window_kb"] == ws_kb])
            print(f"  {label} @ {ws_kb}kb: {n} boundaries", flush=True)
else:
    print("  No boundaries detected.", flush=True)

print("\n[done]", flush=True)
