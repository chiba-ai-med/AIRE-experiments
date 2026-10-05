#!/usr/bin/env python
"""Annotate leiden clusters via marker-gene scoring across N celltypes.

Generalises annotate_cluster_celltype.py (brain, hardcoded NPC/CN) to
N celltypes read from a TSV. Used for mTEC where mimetic + Aire+ subtypes
total ~10 celltypes.

Input marker TSV format (header line + N rows):
    Cell type<TAB>Marker genes
    Proliferate<TAB>Mki67, Top2
    Aire<TAB>Aire
    tuft<TAB>Pou2f3, Dclk1, Trpm5, Avil
    ...

For each leiden cluster, scanpy.tl.score_genes scores cells against each
celltype's marker set; per-cluster mean scores are tabulated. Assignment
takes the argmax celltype, falling back to 'other' if (top - second) < min_margin.

Output TSV columns:
    cluster celltype n_cells top_score second_score margin score_<celltype1> ...
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True, help="multiome_labeled.h5mu")
    p.add_argument("--markers", type=Path, required=True,
                   help="TSV with header 'Cell type<TAB>Marker genes' and one row per celltype")
    p.add_argument("--output", type=Path, required=True, help="Output TSV")
    p.add_argument("--min-margin", type=float, default=0.10,
                   help="Min (top_score - second_score) to assign; below → 'other'")
    return p.parse_args()


def load_markers(path: Path) -> dict[str, list[str]]:
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] < 2:
        sys.exit(f"ERROR: {path} must have >=2 columns (got {df.shape[1]})")
    name_col, marker_col = df.columns[0], df.columns[1]
    out: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        name = str(row[name_col]).strip()
        markers = [g.strip() for g in str(row[marker_col]).split(",") if g.strip()]
        out[name] = markers
    return out


def main():
    args = parse_args()
    print(f"[load] {args.input}")
    mdata = md.read(str(args.input))
    rna = mdata["rna"].copy()
    if "cell_type" not in mdata.obs.columns:
        sys.exit("ERROR: mdata.obs['cell_type'] missing -- run cluster_label first")
    rna.obs["cell_type"] = mdata.obs.loc[rna.obs_names, "cell_type"].values

    markers = load_markers(args.markers)
    var_names = set(rna.var_names)
    print(f"[markers] loaded {len(markers)} celltypes from {args.markers}")
    present_per_celltype: dict[str, list[str]] = {}
    for ct, gs in markers.items():
        present = [g for g in gs if g in var_names]
        present_per_celltype[ct] = present
        print(f"  {ct:<15} {len(present)}/{len(gs)} present: {','.join(present)}")
    usable = {ct: gs for ct, gs in present_per_celltype.items() if gs}
    if not usable:
        sys.exit("ERROR: no marker overlap with var_names for any celltype")
    if len(usable) < len(markers):
        skipped = sorted(set(markers) - set(usable))
        print(f"[warn] skipping celltypes with zero markers present: {skipped}")

    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)

    score_cols: list[str] = []
    for ct, gs in usable.items():
        col = f"score_{ct}"
        sc.tl.score_genes(rna, gene_list=gs, score_name=col)
        score_cols.append(col)

    df = rna.obs[["cell_type", *score_cols]].copy()
    summary = df.groupby("cell_type", observed=True).agg(
        **{c: (c, "mean") for c in score_cols},
        n_cells=(score_cols[0], "size"),
    ).reset_index()

    score_arr = summary[score_cols].to_numpy()
    order = np.argsort(-score_arr, axis=1)
    top_idx = order[:, 0]
    second_idx = order[:, 1] if score_arr.shape[1] > 1 else order[:, 0]
    top_celltype_full = [score_cols[i].removeprefix("score_") for i in top_idx]
    top_score = score_arr[np.arange(len(top_idx)), top_idx]
    second_score = score_arr[np.arange(len(second_idx)), second_idx]
    margin = top_score - second_score
    assigned = [
        ct if m >= args.min_margin else "other"
        for ct, m in zip(top_celltype_full, margin)
    ]

    summary.insert(1, "celltype", assigned)
    summary.insert(2, "top_score", top_score)
    summary.insert(3, "second_score", second_score)
    summary.insert(4, "margin", margin)
    summary = summary.rename(columns={"cell_type": "cluster"})
    summary = summary.sort_values("cluster")
    print()
    print(summary.to_string(index=False))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)
    print(f"[done] wrote {args.output}")


if __name__ == "__main__":
    main()
