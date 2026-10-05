# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project purpose

Validation experiments for **Machima2** (R package: https://github.com/kokitsuyuzaki/Machima), a cell-type deconvolution method for Hi-C using paired scRNA+scATAC multiome as the cell-type basis. The model is:

- `X_RNA  ≈ W_RNA · H_RNA`  (scRNA-side, per cell)
- `X_Epi  ≈ (T · W_RNA) · H_Sym · (T · W_RNA)ᵀ`  (Hi-C contact matrix, symmetric)

`T` (`l × n`, Hi-C bins × ATAC features) can be **identity (frozen with `fixT=TRUE`)** or **dense (learned)**. Both variants are evaluated.

## Two-tier dataset design

| tier | dataset | scRNA+ATAC multiome | bulk Hi-C | scHi-C |
|---|---|---|---|---|
| **practice** | mouse brain cortex | GSE210749 (10x), GSE152020 (Paired-Tag, H3K27ac as ATAC proxy) | GSE96107 (pairs) | GSE238001 (GAGE-seq), GSE253407 (Droplet Hi-C) |
| **production** | mTEC / AIRE | `data/mtec/reference/` (10x Multiome, already present) | SRA SRP330308 (fastq, no contact matrix on GEO) | none |

Production has **no ground truth** — this is by design. Pipeline must close without scHi-C so it transfers from brain → mTEC unchanged. scHi-C tooling lives in a separate module and never blocks the main DAG.

**mTEC Hi-C source caveat:** GSE180937 supplementary contains only PC1 compartment bedGraphs and ChIP-seq bigWigs — no contact matrices. Contact matrices must be regenerated from raw fastq in SRA project SRP330308 (~30 GB compressed, ~12-24 h alignment on 32 cores). This is handled in a separate Snakefile `workflow/sra_hic.smk` that wraps **nf-core/fetchngs + nf-core/hic** so resolution / bin grid match the ATAC side exactly. The output `data/mtec/processed/mtec_hic.mcool` then feeds into the main `bin_hic_by_chr` rule.

## Validation strategy

- **Primary**: `W_RNA` column-label match against scRNA clusters; bulk Hi-C reconstruction error; Identity-T vs Dense-T comparison.
- **Secondary**: TAD/compartment plausibility of reconstructed cell-type contact maps; comparison with sorted-population bulk Hi-C if available.
- **Reference only (not gold)**: scHi-C aggregate. scHi-C is too noisy/sparse to falsify Machima2 outputs.

## Bin-resolution staging

ATAC peaks are re-binned to match Hi-C bins (required for Identity T where `n = l`).

- **100kb** — default; full genome × both datasets × both T variants. Used for development iteration and TAD-level views.
- **25kb** — secondary pass after 100kb shake-down; sub-TAD / boundaries.
- **10kb** — selective (e.g., chr10 around *Aire* for mTEC); loop-level. Identity T only due to cost.

## Machima2 input shape (list mode)

Per-chromosome list inputs are required:

- `X_RNA[[k]]`: `n_k × m`  (m = cells, **identical across all k**)
- `X_Epi[[k]]`: `l_k × l_k` symmetric
- `T[[k]]`:    `l_k × n_k` (or `diag(l_k)` with `fixT=TRUE` for identity variant)

`label` argument receives RNA-derived per-cell labels (multiome → direct barcode copy, no diagonal-integration assumption).

Output: `W_RNA` becomes a list (per chr), `H_RNA` and `H_Sym` are shared single matrices, `T` is a list.

## Repository layout

```
data/{brain,mtec}/
  reference/    # scRNA + scATAC multiome (mtec/reference/ already populated with cellranger-arc output)
  query/        # Hi-C (brain only -- mTEC Hi-C lives under data/mtec/sra/ + data/mtec/processed/)
  processed/    # peaks, .h5mu, mcool, per-chr bin .mtx dirs
  sra/          # mTEC only: nf-core/fetchngs SRA workspace
workflow/
  preprocess.smk    # main pipeline: download -> MuData -> binning -> Machima2
  sra_hic.smk       # mTEC-only: SRA fastq -> nf-core/hic -> mtec_hic.mcool
  envs/             # conda env yamls (r-machima, py-hic, nextflow)
  config.yaml
  run_preprocess.sh / run_sra_hic.sh
src/                # R/Python/shell scripts called by Snakemake rules
output/             # Machima2 results (.rds), keyed by {tissue, T-variant, resolution}
plot/  logs/  benchmarks/
```

Output naming convention: `output/machima2_{tissue}_T{identity,dense}_{resolution}.rds` (resolution in bp).

## Environment & execution

Snakemake-driven. Style follows `/home/koki/dev/eCCI-experiments/workflow/phase0_loss_detection.smk`:

- `min_version("8.10.0")`
- Top-of-file `DOCKER_IMAGE` and `CONDA_ENV` constants
- **Every rule has both `conda:` and `container:` directives** — pipeline runs identically under `--use-conda` (current) and `--use-singularity` (future Dockerization). The two modes differ only in the snakemake flag.
- Every rule: triple-quoted docstring, `log: logs/<phase>/<rule>.log`, `benchmark: benchmarks/<phase>/<rule>.txt`
- Rules invoke single `Rscript src/<phase>_<task>.R` with positional args ending in `>& {log}`
- Wrapper: `workflow/run_<phase>.sh` with `snakemake --snakefile ... --cores N --printshellcmds --keep-going --rerun-incomplete "$@"`

Stack: R (Seurat, Signac, nnTensor, Machima) + Python (scanpy, snapatac2, cooler) + Snakemake.

**`Machima` is GitHub-only** (not on CRAN/Bioconductor) — env setup needs a post-install step in `setup_r_packages.R` or in the Dockerfile:

```sh
Rscript -e 'remotes::install_github("kokitsuyuzaki/Machima")'
```

## Data acquisition

Prefer **GEO supplementary processed matrices** over raw fastq → realignment. Use cellranger-arc output structure (already present for mTEC) as the canonical multiome shape.
