#!/usr/bin/env python
"""Slice an mcool at a chosen resolution into per-chromosome symmetric
bin x bin contact matrices (MatrixMarket).

Output layout (under <output_dir>):
    chr1.mtx, chr2.mtx, ..., chrX.mtx, chrY.mtx   -- bin x bin (symmetric, balanced if available)
    chroms.txt                                     -- chromosome names

The bin grid here MUST match the ATAC bin grid produced by
preprocess_bin_atac_by_chr.py at the same resolution. Both use a
left-aligned grid over chrom_sizes from snapatac2.genome (mm10), so they
align by construction at every resolution.

If `weight` is present in the cool's bins (cooler balance applied), the
balanced matrix is exported; otherwise raw counts.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True, help="Input .mcool")
    p.add_argument("--resolution", type=int, required=True, help="Bin size in bp")
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cool_uri = f"{args.input}::resolutions/{args.resolution}"
    print(f"[open] {cool_uri}")
    c = cooler.Cooler(cool_uri)
    has_weight = "weight" in c.bins().columns
    print(f"[cool] {c.binsize} bp, {c.chromnames}, balanced={has_weight}")

    chroms_out = []
    for chrom in c.chromnames:
        if "_" in chrom or not chrom.startswith("chr"):
            continue
        mat = c.matrix(balance=has_weight, sparse=True).fetch(chrom)
        # Replace NaN (from balance weights) with 0 to keep .mtx round-trip safe.
        if has_weight:
            mat = mat.tocoo()
            valid = ~np.isnan(mat.data)
            mat = sp.coo_matrix((mat.data[valid], (mat.row[valid], mat.col[valid])),
                                shape=mat.shape).tocsr()
        else:
            mat = mat.tocsr()
        # Ensure symmetry (cooler returns upper triangle for some fetches).
        mat = mat.maximum(mat.T)
        out = args.output_dir / f"{chrom}.mtx"
        scipy.io.mmwrite(str(out), mat, symmetry="symmetric")
        chroms_out.append(chrom)
        print(f"  {chrom}: {mat.shape}, nnz={mat.nnz} -> {out.name}")

    (args.output_dir / "chroms.txt").write_text("\n".join(chroms_out) + "\n")
    print(f"[write] {args.output_dir}")


if __name__ == "__main__":
    main()
