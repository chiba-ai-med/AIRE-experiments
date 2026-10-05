#!/usr/bin/env python
"""Annotate leiden clusters as NPC / CN / other via marker-gene scoring.

Reads multiome_labeled.h5mu and uses scanpy.tl.score_genes to score each
cell on the NPC and CN marker signatures. Per-cluster mean scores are then
compared and each cluster is labeled with the higher-scoring celltype, or
'other' when the absolute score gap is below `--min-diff`.

Output: TSV with columns
    cluster cell_type score_npc score_cn score_diff n_cells
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
    p.add_argument("--input", type=Path, required=True, help="Input multiome_labeled.h5mu")
    p.add_argument("--output", type=Path, required=True, help="Output TSV")
    p.add_argument("--npc-markers", type=str, required=True, help="Comma-separated mouse symbols")
    p.add_argument("--cn-markers", type=str, required=True, help="Comma-separated mouse symbols")
    p.add_argument("--min-diff", type=float, default=0.15,
                   help="Min |score_npc - score_cn| to assign npc/cn (else 'other')")
    return p.parse_args()


def main():
    args = parse_args()
    print(f"[load] {args.input}")
    mdata = md.read(str(args.input))
    rna = mdata["rna"].copy()
    cell_type_col = mdata.obs.get("cell_type")
    if cell_type_col is None:
        sys.exit("ERROR: mdata.obs['cell_type'] missing -- run cluster_label first")

    # Move cluster labels onto the RNA modality so scanpy ops can find them.
    # bin_atac_by_chr already reorders cells; we re-attach via barcode.
    rna.obs["cell_type"] = mdata.obs.loc[rna.obs_names, "cell_type"].values

    npc_markers = [g.strip() for g in args.npc_markers.split(",") if g.strip()]
    cn_markers  = [g.strip() for g in args.cn_markers.split(",")  if g.strip()]
    var_names = set(rna.var_names)
    npc_present = [g for g in npc_markers if g in var_names]
    cn_present  = [g for g in cn_markers  if g in var_names]
    print(f"[markers] npc: {len(npc_present)}/{len(npc_markers)} present "
          f"({','.join(npc_present)})")
    print(f"[markers] cn : {len(cn_present)}/{len(cn_markers)}  present "
          f"({','.join(cn_present)})")
    if not npc_present or not cn_present:
        sys.exit("ERROR: no markers from one of the two signatures present in var_names")

    # score_genes wants log-normalised data; do that on a copy if not already.
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)

    sc.tl.score_genes(rna, gene_list=npc_present, score_name="score_npc")
    sc.tl.score_genes(rna, gene_list=cn_present,  score_name="score_cn")

    df = rna.obs[["cell_type", "score_npc", "score_cn"]].copy()
    grouped = df.groupby("cell_type")
    summary = grouped.agg(score_npc=("score_npc", "mean"),
                           score_cn=("score_cn", "mean"),
                           n_cells=("score_npc", "size")).reset_index()
    summary["score_diff"] = summary["score_npc"] - summary["score_cn"]

    def assign(row):
        if abs(row.score_diff) < args.min_diff:
            return "other"
        return "npc" if row.score_diff > 0 else "cn"

    summary["celltype"] = summary.apply(assign, axis=1)
    summary = summary[["cell_type", "celltype", "score_npc", "score_cn",
                        "score_diff", "n_cells"]]
    summary = summary.rename(columns={"cell_type": "cluster"})
    summary = summary.sort_values("cluster")
    print(summary.to_string(index=False))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)
    print(f"[done] wrote {args.output}")


if __name__ == "__main__":
    main()
