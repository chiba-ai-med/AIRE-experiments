#!/usr/bin/env python
"""Bin scRNA UMI counts to genomic intervals using TES (3' end) anchor.

10X Chromium 3' sequencing piles reads near the transcript 3' end. To align
RNA-derived bin counts with ATAC bin counts (peak coordinates) and Hi-C bin
counts (contact bin pairs), we anchor each gene's UMI counts at its **TES**
(not TSS):

    + strand gene: TES = gene.end   (largest genomic coordinate)
    - strand gene: TES = gene.start (smallest genomic coordinate)

Reads from cellranger filtered_feature_bc_matrix (raw UMI integer counts,
gene x cell), filters to Gene Expression features, bins by TES, then applies
library-size normalization (target_sum=10000) and log1p transform per cell.

Output schema:
    <out_dir>/chr<N>.mtx     # cells x bins, float (log1p of lib-norm counts)
    <out_dir>/cells.tsv      # one barcode per line
    <out_dir>/labels.tsv     # barcode<TAB>cellType (filtered to keep list)
    <out_dir>/chroms.txt     # one chrom per line, matches atac_bins ordering
    <out_dir>/bins.tsv.gz    # chrom<TAB>start<TAB>end, full bin index

Usage:
    python preprocess_bin_rna_by_chr.py \
        --matrix-dir data/mtec/reference/filtered_feature_bc_matrix \
        --gtf  /home/godayuki/OP/ref_genome/refdata-gex-mm10-2020-A/genes/genes.gtf \
        --chroms-template data/mtec/processed/atac_bins_100000_celltype/chroms.txt \
        --bins-template data/mtec/processed/atac_bins_100000_celltype/bins.tsv.gz \
        --labels-template data/mtec/processed/atac_bins_100000_celltype/labels.tsv \
        --cells-template data/mtec/processed/atac_bins_100000_celltype/cells.tsv \
        --resolution 100000 \
        --target-sum 10000 \
        --output-dir data/mtec/processed/rna_bins_100000_celltype
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser()
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--matrix-dir", type=Path,
                       help="cellranger filtered_feature_bc_matrix dir "
                            "(matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)")
    group.add_argument("--dgem-tsv", type=Path,
                       help="GSE-style DGEM TSV (gene x cell, gene_name in col 1, "
                            "cell barcodes in header row). Optionally gzipped.")
    p.add_argument("--dgem-cell-suffix", default=None,
                   help="When using --dgem-tsv, keep only cells whose header "
                        "ends with this suffix; suffix is then stripped to "
                        "match cells-template. Example: '-1___10x_complex'.")
    p.add_argument("--gtf", type=Path, required=True,
                   help="mm10 GTF (cellranger reference)")
    p.add_argument("--chroms-template", type=Path, required=True)
    p.add_argument("--bins-template", type=Path, required=True)
    p.add_argument("--labels-template", type=Path, required=True)
    p.add_argument("--cells-template", type=Path, required=True)
    p.add_argument("--resolution", type=int, required=True)
    p.add_argument("--target-sum", type=float, default=10000.0,
                   help="Library size normalization target (default 10000)")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def parse_gtf_genes(gtf_path: Path) -> pd.DataFrame:
    """Extract gene records with TES (strand-aware 3' end)."""
    records = []
    open_fn = gzip.open if gtf_path.suffix == ".gz" else open
    with open_fn(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attrs = {}
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, _, val = attr.partition(" ")
                attrs[key] = val.strip('"')
            tes = end if strand == "+" else start
            records.append((attrs.get("gene_id", ""), attrs.get("gene_name", ""),
                            chrom, start, end, strand, tes))
    return pd.DataFrame(records, columns=[
        "gene_id", "gene_name", "chrom", "start", "end", "strand", "tes"])


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Parse GTF ---
    print(f"[gtf] reading {args.gtf}")
    df_genes = parse_gtf_genes(args.gtf)
    print(f"[gtf] {len(df_genes)} gene records")

    # --- 2. Read raw UMI matrix ---
    if args.matrix_dir is not None:
        # cellranger filtered_feature_bc_matrix format
        print(f"[cellranger] reading {args.matrix_dir}")
        mtx = sio.mmread(str(args.matrix_dir / "matrix.mtx.gz")).tocsc()
        print(f"[cellranger] matrix.mtx.gz shape: {mtx.shape} (features x cells)")
        features = pd.read_csv(args.matrix_dir / "features.tsv.gz", sep="\t",
                               header=None, names=["gene_id", "gene_name", "type",
                                                   "chrom", "start", "end"])
        with gzip.open(args.matrix_dir / "barcodes.tsv.gz", "rt") as f:
            barcodes = [b.strip() for b in f if b.strip()]
        print(f"[cellranger] {len(features)} features, {len(barcodes)} cells")
        barcodes_stripped = [b.split("-")[0] for b in barcodes]
        gene_mask = (features["type"] == "Gene Expression").values
        gene_features = features.loc[gene_mask].reset_index(drop=True)
        print(f"[cellranger] {gene_mask.sum()} Gene Expression features")
        X_gene = mtx[gene_mask, :].tocsr()  # genes x cells
    else:
        # DGEM TSV format: rows = genes, columns = cells, first col = gene_name,
        # header row = cell barcodes.
        print(f"[dgem] reading {args.dgem_tsv}")
        open_fn = gzip.open if str(args.dgem_tsv).endswith(".gz") else open
        # Read header first
        with open_fn(args.dgem_tsv, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
        # The header may or may not have a leading empty column for the index.
        # If it does, the gene_name column index is 0; else gene_name is implicit
        # in row index. Use pd.read_csv with index_col=0 either way.
        # First pass: get all barcodes + filter by suffix
        all_barcodes = header  # column headers, gene_name column may or may not be present
        # If first cell is "" or a non-barcode label like "gene"/"gene_id"/"GENE",
        # it's the index column header; drop it so that all_barcodes length matches
        # the cell-column count read by pd.read_csv(..., index_col=0).
        if all_barcodes and (all_barcodes[0] == ""
                             or all_barcodes[0].lower() in ("gene", "gene_id", "gene_name", "geneid", "genename", "index")):
            all_barcodes = all_barcodes[1:]
        if args.dgem_cell_suffix:
            keep_mask = np.array([b.endswith(args.dgem_cell_suffix) for b in all_barcodes])
            barcodes = [b for b, k in zip(all_barcodes, keep_mask) if k]
            barcodes_stripped = [b[:-len(args.dgem_cell_suffix)] for b in barcodes]
            print(f"[dgem] {len(all_barcodes)} total cells, "
                  f"{len(barcodes)} after suffix filter '{args.dgem_cell_suffix}'")
        else:
            barcodes = list(all_barcodes)
            barcodes_stripped = list(barcodes)
            keep_mask = np.ones(len(all_barcodes), dtype=bool)
        # Read full DGEM as DataFrame (chunked if huge, but 19k x 16k is OK)
        df = pd.read_csv(args.dgem_tsv, sep="\t", index_col=0,
                         compression="infer", header=0)
        # df is gene_name (index) x cell_barcode (columns), values are counts
        df = df.loc[:, keep_mask]  # filter cells by suffix
        # Sparse conversion
        X_gene = sp.csr_matrix(df.values.astype(np.int32))  # genes x cells
        gene_features = pd.DataFrame({
            "gene_id": df.index,  # DGEM has gene_name only; use as gene_id for matching
            "gene_name": df.index,
            "type": "Gene Expression",
        })
        print(f"[dgem] {len(gene_features)} genes, {X_gene.shape[1]} cells; "
              f"sparse density {X_gene.nnz / (X_gene.shape[0]*X_gene.shape[1]):.4f}")
    print(f"[matrix] X_gene shape: {X_gene.shape}, "
          f"sample values: min={X_gene.data.min() if X_gene.nnz>0 else 0}, "
          f"max={X_gene.data.max() if X_gene.nnz>0 else 0}, "
          f"dtype={X_gene.dtype}")
    assert X_gene.data.min() >= 0, "Negative UMI counts -- not raw!"

    # --- 3. Read chroms list and bin grid ---
    chroms = args.chroms_template.read_text().splitlines()
    bins_df = pd.read_csv(args.bins_template, sep="\t", header=None,
                          names=["chrom", "start", "end"])
    bins_per_chrom = bins_df.groupby("chrom").size().to_dict()
    print(f"[chroms] {len(chroms)} chromosomes, bins/chrom example: chr1={bins_per_chrom.get('chr1', 0)}")

    # --- 4. Match gene_features to GTF gene records ---
    # cellranger has both ENSMUSG gene_id and gene_name; DGEM has only gene_name.
    # Auto-detect: if first few gene_ids start with "ENSMUSG", match by gene_id;
    # else match by gene_name.
    sample_ids = gene_features["gene_id"].head(10).astype(str).tolist()
    lookup_col = "gene_id" if any(s.startswith("ENSMUSG") for s in sample_ids) else "gene_name"
    gene_map = df_genes.drop_duplicates(lookup_col, keep="first").set_index(lookup_col)
    print(f"[gtf] matching by '{lookup_col}', {len(gene_map)} unique entries")

    # For each gene-row in cellranger, find TES and bin
    chrom_gene_rows = {c: {"row_idx": [], "bin_idx": []} for c in chroms}
    n_matched = 0; n_unmatched = 0; n_off_chrom = 0
    for row_idx, gene_id in enumerate(gene_features[lookup_col].values):
        if gene_id not in gene_map.index:
            n_unmatched += 1
            continue
        rec = gene_map.loc[gene_id]
        if rec["chrom"] not in bins_per_chrom:
            n_off_chrom += 1
            continue
        bin_local = (rec["tes"] - 1) // args.resolution
        if bin_local >= bins_per_chrom.get(rec["chrom"], 0):
            n_off_chrom += 1
            continue
        chrom_gene_rows[rec["chrom"]]["row_idx"].append(row_idx)
        chrom_gene_rows[rec["chrom"]]["bin_idx"].append(bin_local)
        n_matched += 1
    print(f"[map] matched={n_matched}, unmatched={n_unmatched}, off_chrom={n_off_chrom}")

    # --- 5. Filter cells to keep list ---
    keep_bc = args.cells_template.read_text().splitlines()
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes_stripped)}
    keep_idx = []
    keep_kept_bc = []
    for bc in keep_bc:
        if bc in bc_to_idx:
            keep_idx.append(bc_to_idx[bc])
            keep_kept_bc.append(bc)
    print(f"[cells] keeping {len(keep_idx)} / {len(keep_bc)} cells")
    keep_arr = np.array(keep_idx)
    X_keep = X_gene[:, keep_arr].tocsc()  # genes x cells (kept)
    print(f"[cells] X_keep shape: {X_keep.shape}")

    # --- 6. Compute library size per cell BEFORE binning (using all genes, not just binned) ---
    # Note: we use the SUM OF GENES THAT WILL BE BINNED (not all 32k genes) as library size
    # so that bin counts and norm factor are consistent.
    binned_row_idxs = sorted(set(
        idx for c in chroms for idx in chrom_gene_rows[c]["row_idx"]))
    X_binnable = X_keep[binned_row_idxs, :]
    lib_size = np.asarray(X_binnable.sum(axis=0)).flatten()  # length = n_cells_kept
    lib_size = np.maximum(lib_size, 1.0)  # avoid div by zero
    print(f"[norm] library size per cell: median={np.median(lib_size):.0f}, "
          f"min={lib_size.min():.0f}, max={lib_size.max():.0f}")

    # --- 7. Per chrom: aggregate gene rows to bin rows, normalize, log1p, transpose to (cells, bins) ---
    for c in chroms:
        rows = chrom_gene_rows[c]
        n_bins_c = bins_per_chrom.get(c, 0)
        n_cells_kept = X_keep.shape[1]
        if not rows["row_idx"] or n_bins_c == 0:
            M = sp.csr_matrix((n_cells_kept, n_bins_c))
            sio.mmwrite(str(args.output_dir / f"{c}.mtx"), M, field="real", symmetry="general")
            print(f"  {c}: 0 genes -> empty {n_cells_kept} x {n_bins_c}")
            continue
        # Build (bin_idx x gene_row_idx) mapping matrix P: P[b, g] = 1 if gene g maps to bin b
        gene_idxs = np.array(rows["row_idx"])
        bin_idxs = np.array(rows["bin_idx"])
        P = sp.csr_matrix(
            (np.ones(len(gene_idxs)), (bin_idxs, gene_idxs)),
            shape=(n_bins_c, X_keep.shape[0]),
        )
        # P @ X_keep = (n_bins_c x n_cells_kept) of summed raw counts per bin per cell
        M_bin_cell_raw = P @ X_keep  # sparse
        # Normalize per cell: divide column j by lib_size[j], multiply by target_sum
        scale = (args.target_sum / lib_size).reshape(1, -1)  # 1 x n_cells
        # multiply: each column j scaled by scale[0, j]
        D = sp.diags(scale.flatten())  # diagonal n_cells x n_cells
        M_norm = M_bin_cell_raw @ D
        # log1p
        M_lognorm = M_norm.copy()
        M_lognorm.data = np.log1p(M_lognorm.data)
        # transpose to (cells, bins) for ATAC convention
        M_out = M_lognorm.T.tocsr()
        sio.mmwrite(str(args.output_dir / f"{c}.mtx"), M_out, field="real",
                    symmetry="general")
        nnz_raw = M_bin_cell_raw.nnz
        print(f"  {c}: {len(gene_idxs)} genes -> ({n_cells_kept}, {n_bins_c}) "
              f"nnz={M_out.nnz}, value range=[{M_out.data.min():.4g}, {M_out.data.max():.4g}]")

    # --- 8. Write cells.tsv, labels.tsv, chroms.txt, bins.tsv.gz ---
    (args.output_dir / "cells.tsv").write_text("\n".join(keep_kept_bc) + "\n")
    lab_df = pd.read_csv(args.labels_template, sep="\t", header=None,
                         names=["barcode", "celltype"])
    lab_map = dict(zip(lab_df["barcode"], lab_df["celltype"]))
    with (args.output_dir / "labels.tsv").open("w") as f:
        for bc in keep_kept_bc:
            f.write(f"{bc}\t{lab_map.get(bc, 'other')}\n")
    shutil.copy(args.chroms_template, args.output_dir / "chroms.txt")
    shutil.copy(args.bins_template, args.output_dir / "bins.tsv.gz")
    print(f"[done] wrote outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
