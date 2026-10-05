#!/usr/bin/env python
"""Cluster the RNA modality of a multiome MuData and write per-cell labels.

For multiome data the RNA and ATAC barcodes are identical, so the RNA
clustering result becomes the cell label for both modalities -- no
explicit label transfer is needed. We store the result as obs['cell_type']
on the parent MuData.

Pipeline:
  1. Load .h5mu, take rna modality
  2. Normalize, log1p, HVG, PCA
  3. Neighbors + leiden
  4. Write cluster IDs (as 'cluster_<k>' strings) into mdata.obs['cell_type']
  5. Save updated .h5mu

Cell-type biological annotation (e.g. by marker genes or reference
mapping) is left to a follow-up step. Machima2's `label` argument only
needs human-readable per-cell strings; cluster IDs satisfy that.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mudata as md
import scanpy as sc


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True, help="Input .h5mu")
    p.add_argument("--output", type=Path, required=True, help="Output .h5mu (with cell_type)")
    p.add_argument("--n-hvg", type=int, default=3000)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--leiden-resolution", type=float, default=0.5)
    args = p.parse_args()

    print(f"[load] {args.input}")
    mdata = md.read(str(args.input))
    rna = mdata["rna"].copy()
    print(f"[rna] shape: {rna.shape}")

    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    # flavor="seurat" uses log-normalised data and avoids the scikit-misc
    # dependency required by "seurat_v3".
    sc.pp.highly_variable_genes(rna, n_top_genes=args.n_hvg, flavor="seurat")
    sc.pp.scale(rna, max_value=10)
    sc.tl.pca(rna, n_comps=args.n_pcs)
    sc.pp.neighbors(rna, n_neighbors=args.n_neighbors)
    sc.tl.leiden(rna, resolution=args.leiden_resolution, key_added="leiden")

    labels = rna.obs["leiden"].astype(str).map(lambda x: f"cluster_{x}")
    print(f"[cluster] {labels.nunique()} leiden clusters at resolution={args.leiden_resolution}")

    # mdata["rna"] = rna does not work in mudata >=0.3 (no item assignment).
    # Use mdata.mod (the underlying dict).
    mdata.mod["rna"] = rna
    mdata.obs["cell_type"] = labels.values

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"[write] {args.output}")
    mdata.write(str(args.output))


if __name__ == "__main__":
    main()
