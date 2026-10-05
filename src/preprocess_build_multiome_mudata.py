#!/usr/bin/env python
"""Build a multiome MuData from RNA gene-cell counts + ATAC fragments + peaks.

Pipeline:
  1. Load RNA gene-cell counts -> AnnData
       --rna-format dgem-tsv : tab-separated genes x cells matrix
       --rna-format 10x-mtx  : directory containing matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz
  2. Build cell x peak count matrix from fragments + peaks (snapatac2)
  3. Intersect cell barcodes between RNA and ATAC
  4. Compose MuData {rna, atac} and write .h5mu

This is the unified loader for both brain (DGEM-TSV + bare fragments) and
mTEC (cellranger-arc 10x mtx + atac_fragments). See CLAUDE.md and
project_multiome_input_unification memory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap


def load_rna(rna_format: str, rna_input: Path, primary_sample: str | None = None) -> ad.AnnData:
    if rna_format == "dgem-tsv":
        # SCENIC+ DGEM: rows = genes, columns = cells. Read full into memory --
        # ~46 MB compressed for GSE210747 brain.
        df = pd.read_csv(rna_input, sep="\t", index_col=0)
        # Brain DGEM cell barcodes are "<16bp>-1___<sample>" combining all 5
        # protocol conditions. We must filter to cells whose ___<sample> tag
        # matches `primary_sample` BEFORE stripping the suffix -- otherwise
        # different conditions can collide on plain barcodes and break
        # downstream indexing.
        cols = df.columns.astype(str)
        if cols.str.contains("___").any():
            if primary_sample is None:
                raise ValueError(
                    "DGEM has '___<sample>' tags (multi-sample combined matrix); "
                    "pass --primary-sample to pick one condition."
                )
            keep = cols.str.endswith(f"___{primary_sample}")
            n_total, n_keep = len(cols), int(keep.sum())
            if n_keep == 0:
                raise RuntimeError(
                    f"primary_sample={primary_sample!r} matched 0 of {n_total} cells"
                )
            df = df.loc[:, keep]
            df.columns = df.columns.str.replace(r"___.*$", "", regex=True)
            print(f"[load] RNA filtered to primary_sample={primary_sample}: "
                  f"{n_keep}/{n_total} cells kept")
        adata = ad.AnnData(
            X=df.T.values.astype(np.float32),
            obs=pd.DataFrame(index=df.columns.astype(str)),
            var=pd.DataFrame(index=df.index.astype(str)),
        )
    elif rna_format == "10x-mtx":
        # cellranger-arc filtered_feature_bc_matrix/ directory.
        adata = sc.read_10x_mtx(str(rna_input), var_names="gene_symbols", cache=False)
        # Some cellranger-arc outputs include both Gene Expression and Peaks
        # feature types; subset to GEX only if the feature_type column exists.
        if "feature_types" in adata.var.columns:
            adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    else:
        raise ValueError(f"Unknown rna-format: {rna_format}")
    adata.var_names_make_unique()
    return adata


def build_atac(fragments: Path, peaks_bed: Path, genome: str, tmp_h5ad: Path) -> ad.AnnData:
    chrom_sizes = getattr(snap.genome, genome)
    # snapatac2 2.9.x: import_data was renamed to import_fragments.
    atac = snap.pp.import_fragments(
        fragment_file=str(fragments),
        chrom_sizes=chrom_sizes,
        file=str(tmp_h5ad),
        sorted_by_barcode=False,
    )
    # make_peak_matrix takes the peak set via peak_file= (a BED-like path) and
    # returns an in-memory AnnData when no `file=` is given.
    peak_mat = snap.pp.make_peak_matrix(atac, peak_file=str(peaks_bed))
    if hasattr(peak_mat, "to_memory"):
        peak_mat = peak_mat.to_memory()
    return peak_mat


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rna-format", choices=["dgem-tsv", "10x-mtx"], required=True)
    p.add_argument("--rna-input", type=Path, required=True)
    p.add_argument("--primary-sample", default=None,
                   help="Required for combined dgem-tsv with '___<sample>' barcode suffix")
    p.add_argument("--atac-fragments", type=Path, required=True)
    p.add_argument("--atac-peaks", type=Path, required=True,
                   help="MACS2 narrowPeak BED produced by preprocess_call_peaks.sh")
    p.add_argument("--genome", default="mm10")
    p.add_argument("--output", type=Path, required=True, help="Output .h5mu")
    p.add_argument("--tmp-dir", type=Path, default=None,
                   help="Where snapatac2 stores its backed h5ad (default: <output>.snap.h5ad)")
    args = p.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = args.tmp_dir or args.output.with_suffix(".snap.h5ad")

    print(f"[load] RNA  ({args.rna_format}) <- {args.rna_input}")
    rna = load_rna(args.rna_format, args.rna_input, args.primary_sample)
    print(f"[load] RNA shape: {rna.shape} (cells x genes)")

    print(f"[build] ATAC peak matrix from {args.atac_fragments} + {args.atac_peaks}")
    atac = build_atac(args.atac_fragments, args.atac_peaks, args.genome, tmp_h5ad)
    print(f"[load] ATAC shape: {atac.shape} (cells x peaks)")

    # Intersect barcodes. snapatac2 may strip cellranger -1 suffix; harmonise.
    def strip_suffix(idx: pd.Index) -> pd.Index:
        return idx.str.replace(r"-\d+$", "", regex=True)

    rna.obs_names = strip_suffix(rna.obs_names)
    atac.obs_names = strip_suffix(atac.obs_names)
    common = sorted(set(rna.obs_names) & set(atac.obs_names))
    if not common:
        raise RuntimeError("No shared cell barcodes between RNA and ATAC.")
    print(f"[intersect] Common cells: {len(common)} (rna={rna.n_obs}, atac={atac.n_obs})")
    rna = rna[common, :].copy()
    atac = atac[common, :].copy()

    mdata = md.MuData({"rna": rna, "atac": atac})
    mdata.update()
    print(f"[write] {args.output}")
    mdata.write(str(args.output))


if __name__ == "__main__":
    main()
