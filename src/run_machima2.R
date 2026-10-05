#!/usr/bin/env Rscript
# Run Machima2 on per-chromosome multiome ATAC + Hi-C lists.
#
# Inputs are two directories of MatrixMarket files (one .mtx per chromosome):
#   atac_dir/  : cell x bin counts  (from preprocess_bin_atac_by_chr.py)
#   hic_dir/   : bin x bin symmetric (from preprocess_bin_hic_by_chr.py)
# plus optional per-cell labels (atac_dir/labels.tsv).
#
# Both directories must agree on chroms.txt (same chromosomes in same order)
# and on the per-chromosome bin count for each chromosome.
#
# Machima2 list-mode shapes (see CLAUDE.md):
#   X_RNA[[k]] : n_k x m   (transpose of file)
#   X_Epi[[k]] : l_k x l_k symmetric
#   T[[k]]     : l_k x n_k (l_k == n_k for identity, since same bin grid)
#
# Stages:
#   joint           : original Machima2 joint factorisation (X_RNA + X_Epi
#                     fitted simultaneously; W_RNA shared)
#   transferFlog    : Stage A pure NMF on log1p(X_RNA) with Frobenius
#                     -> Stage B Machima2 with frozen W_RNA, H_RNA
#   transferKraw    : Stage A pure NMF on raw X_RNA with KL divergence
#                     -> Stage B as above
#   supervisedHfix  : Stage A NMF with H_RNA pinned to leiden one-hot
#                     (W = per-cluster centroid, ARI=1 by construction)
#                     -> Stage B as above
#   supervisedGNMF  : Stage A graph-regularised NMF with cell-side Laplacian
#                     built from leiden labels (rikenbit/nnTensor only)
#                     -> Stage B as above
#   supervisedWinit : Stage A NMF initialised from per-cluster mean W,
#                     refined unconstrained -> Stage B as above
#   jointLowRank    : like joint, but T_regularization="low_rank" with
#                     T[k]=U[k]*t(V[k]). Tdense only (incompatible with
#                     fixT=TRUE in upstream Machima >= 1.1.0).

suppressPackageStartupMessages({
  library(Matrix)
  library(nnTensor)
  library(Machima)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: run_machima2.R <atac_dir> <hic_dir> ",
       "<stage: joint|transferFlog|transferKraw|",
       "supervisedHfix|supervisedGNMF|supervisedWinit> ",
       "<T_variant: identity|dense> <out.rds> ",
       "[J=7] [num_iter=30] [transfer_num_iter=100] ",
       "[nmf_num_iter=200] [nmf_n_restart=5] [gnmf_lambda_V=1.0] ",
       "[T_regularization=frobenius_unit] [lambda_T=0] [T_rank=14] ",
       "[H_Sym_structure=diagonal] [lambda_balance=0.5]")
}
atac_dir          <- args[[1]]
hic_dir           <- args[[2]]
stage             <- args[[3]]
T_variant         <- args[[4]]
out_rds           <- args[[5]]
J                 <- if (length(args) >=  6) as.integer(args[[6]])  else 7L
num_iter          <- if (length(args) >=  7) as.integer(args[[7]])  else 30L
transfer_num_iter <- if (length(args) >=  8) as.integer(args[[8]])  else 100L
nmf_num_iter      <- if (length(args) >=  9) as.integer(args[[9]])  else 200L
nmf_n_restart     <- if (length(args) >= 10) as.integer(args[[10]]) else 5L
gnmf_lambda_V     <- if (length(args) >= 11) as.numeric(args[[11]]) else 1.0
T_regularization  <- if (length(args) >= 12) args[[12]]              else "frobenius_unit"
lambda_T          <- if (length(args) >= 13) as.numeric(args[[13]]) else 0
T_rank            <- if (length(args) >= 14) as.integer(args[[14]]) else 14L
H_Sym_structure   <- if (length(args) >= 15) args[[15]]              else "diagonal"
lambda_balance    <- if (length(args) >= 16) as.numeric(args[[16]])  else 0.5
stopifnot(T_regularization %in% c("none", "frobenius_unit", "l2", "low_rank"))
stopifnot(H_Sym_structure %in% c("symmetric", "diagonal"))
stopifnot(lambda_balance >= 0, lambda_balance <= 1)

KNOWN_STAGES <- c("joint", "transferFlog", "transferKraw",
                  "supervisedHfix", "supervisedGNMF", "supervisedWinit",
                  "jointLowRank")
stopifnot(stage %in% KNOWN_STAGES)
stopifnot(T_variant %in% c("identity", "dense"))

# jointLowRank requires a learned T (upstream errors at .checkMachima2 if
# fixT=TRUE && T_regularization="low_rank"). Force the regularization mode
# regardless of the global config -- this stage IS the low-rank evaluation.
if (stage == "jointLowRank") {
  if (T_variant != "dense") {
    stop(sprintf("stage=jointLowRank requires T_variant=dense (got %s)", T_variant))
  }
  T_regularization <- "low_rank"
  message(sprintf("[stage] jointLowRank -> forcing T_regularization=low_rank, T_rank=%d",
                  T_rank))
}

is_supervised <- startsWith(stage, "supervised")
is_joint_like <- stage %in% c("joint", "jointLowRank")

#################################
# Read per-chrom matrices as base dense (Machima2 historical S3-dispatch
# / which() / t.default issues with S4 sparse classes; per-chr dense at
# 100 kb tops out at ~30 MB, fine).
#################################
read_mtx_list <- function(dir, chroms, transpose = FALSE) {
  setNames(lapply(chroms, function(chr) {
    m <- readMM(file.path(dir, paste0(chr, ".mtx")))
    if (transpose) m <- t(m)
    as.matrix(m)
  }), chroms)
}

atac_chroms <- readLines(file.path(atac_dir, "chroms.txt"))
hic_chroms  <- readLines(file.path(hic_dir,  "chroms.txt"))
chroms <- intersect(atac_chroms, hic_chroms)
message(sprintf("[chroms] %d shared (atac=%d, hic=%d)",
                length(chroms), length(atac_chroms), length(hic_chroms)))

# Files are cell x bin; Machima2 expects feature x cell (n_k x m), so transpose.
X_RNA <- read_mtx_list(atac_dir, chroms, transpose = TRUE)
X_Epi <- read_mtx_list(hic_dir,  chroms, transpose = FALSE)

# Sanity: bin counts must match per chr (l_k == n_k for identity T).
for (chr in chroms) {
  if (nrow(X_RNA[[chr]]) != nrow(X_Epi[[chr]])) {
    stop(sprintf("%s: ATAC n_k=%d != HiC l_k=%d (bin grids do not match)",
                 chr, nrow(X_RNA[[chr]]), nrow(X_Epi[[chr]])))
  }
}

# Drop chroms whose ATAC or Hi-C matrix is all zero (e.g. chrM in mm10 since
# Bonev tracks don't cover it). Machima2's .checkMachima2 stopifnots otherwise.
nz <- vapply(chroms, function(chr) {
  any(X_RNA[[chr]] != 0) && any(X_Epi[[chr]] != 0)
}, logical(1))
if (any(!nz)) {
  message(sprintf("[drop] all-zero chroms: %s",
                  paste(chroms[!nz], collapse = ", ")))
  chroms <- chroms[nz]
  X_RNA  <- X_RNA[chroms]
  X_Epi  <- X_Epi[chroms]
}
if (length(chroms) == 0) stop("no chroms remain after all-zero filter")

# Drop chroms with too few bins for rank J (Machima2's .checkMachima2 requires
# J <= min(dim(X_RNA), nrow(X_Epi)). chrM at 100 kb is 1 bin in mm10 -- when
# X_RNA has non-zero mitochondrial signal (RNA bins) but X_Epi is empty there
# (Hi-C tracks don't cover chrM), the all-zero filter above might still keep
# it via a stricter input, but here we add the J-dim guard explicitly).
dim_ok <- vapply(chroms, function(chr) {
  min(dim(X_RNA[[chr]]), nrow(X_Epi[[chr]])) >= J
}, logical(1))
if (any(!dim_ok)) {
  message(sprintf("[drop] chroms with min(dim) < J=%d: %s",
                  J, paste(chroms[!dim_ok], collapse = ", ")))
  chroms <- chroms[dim_ok]
  X_RNA  <- X_RNA[chroms]
  X_Epi  <- X_Epi[chroms]
}
if (length(chroms) == 0) stop("no chroms remain after J-dim filter")

#################################
# Cell labels (required for supervised* stages, optional otherwise).
#################################
labels_path <- file.path(atac_dir, "labels.tsv")
labels <- if (file.exists(labels_path)) {
  lab <- read.table(labels_path, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("barcode", "cell_type"))
  setNames(lab$cell_type, lab$barcode)
} else NULL

if (is_supervised && is.null(labels)) {
  stop(sprintf("stage=%s requires labels.tsv in %s", stage, atac_dir))
}

#################################
# Build aligned label vector matching X_RNA cells (m).
# bin_atac_by_chr writes cells.tsv in the same column order as the .mtx, so
# match by barcode against labels (which is name-keyed by barcode).
#################################
cells <- readLines(file.path(atac_dir, "cells.tsv"))
labels_aligned <- NULL
clusters <- NULL
if (!is.null(labels)) {
  labels_aligned <- unname(labels[cells])
  if (any(is.na(labels_aligned))) {
    stop(sprintf("labels.tsv missing entries for %d / %d cells",
                 sum(is.na(labels_aligned)), length(cells)))
  }
  clusters <- sort(unique(labels_aligned))
  message(sprintf("[labels] %d unique clusters: %s",
                  length(clusters), paste(clusters, collapse = ", ")))
  if (is_supervised && length(clusters) != J) {
    stop(sprintf("supervised stages require J (=%d) == #clusters (=%d)",
                 J, length(clusters)))
  }
}

#################################
# T construction. Identity = base dense diag (Matrix::Diagonal triggers
# S4 dispatch failures inside Machima2's .BetaDivergence -> base::t.default).
#################################
Ts <- if (T_variant == "identity") {
  lapply(X_RNA, function(m) diag(nrow(m)))
} else {
  NULL  # let Machima2 randomly initialise dense T
}
fixT <- (T_variant == "identity")

#################################
# Stage A: pure RNA NMF (every stage except `joint`).
#
# Stack per-chrom X_RNA vertically into one (sum(n_k) x m) matrix, then
# either:
#   transferFlog   : log1p + Frobenius
#   transferKraw   : raw + KL
#   supervisedHfix : raw + Frobenius + initV=one-hot, fixV=TRUE
#   supervisedGNMF : raw + Frobenius + L_graph_V=label-Laplacian
#   supervisedWinit: raw + Frobenius + initU=cluster-mean W, fixU=FALSE
# split U back into a per-chrom list via Mat2List.
#################################

build_label_laplacian <- function(labels_aligned) {
  # Same-cluster adjacency (block-diagonal up to permutation), then L = D - A.
  # Returns a base dense matrix (m x m). At m=2086 this is 35 MB -- fine.
  m <- length(labels_aligned)
  lab_mat <- outer(labels_aligned, labels_aligned, FUN = "==")
  diag(lab_mat) <- FALSE
  A <- 1 * lab_mat
  D <- diag(rowSums(A))
  D - A
}

build_onehot_V <- function(labels_aligned, clusters) {
  # nnTensor::NMF convention: X = U V^T with V being (m x J), so initV must
  # be (m x J). Machima2's init_H_RNA expects (J x m); we transpose later.
  J <- length(clusters)
  m <- length(labels_aligned)
  V <- matrix(0, nrow = m, ncol = J)
  for (j in seq_along(clusters)) {
    V[labels_aligned == clusters[j], j] <- 1
  }
  V
}

build_cluster_mean_W <- function(X_stacked, labels_aligned, clusters) {
  W <- matrix(0, nrow = nrow(X_stacked), ncol = length(clusters))
  for (j in seq_along(clusters)) {
    cells_j <- labels_aligned == clusters[j]
    W[, j]  <- rowMeans(X_stacked[, cells_j, drop = FALSE])
  }
  W
}

init_W_RNA   <- NULL
init_H_RNA   <- NULL
nmf_recerror <- NULL

if (!is_joint_like) {
  prep_X    <- if (stage == "transferFlog") function(M) log1p(M) else identity
  X_RNA_for_nmf <- lapply(X_RNA, prep_X)
  X_stacked     <- do.call(rbind, X_RNA_for_nmf)

  base_args <- list(
    X         = X_stacked,
    J         = J,
    num.iter  = nmf_num_iter,
    verbose   = FALSE
  )

  stage_args <- switch(stage,
    transferFlog = list(algorithm = "Frobenius"),
    transferKraw = list(algorithm = "KL"),
    supervisedHfix = list(
      algorithm = "Frobenius",
      initV     = build_onehot_V(labels_aligned, clusters),
      fixV      = TRUE
    ),
    supervisedGNMF = list(
      algorithm      = "Frobenius",
      L_graph_V      = build_label_laplacian(labels_aligned),
      lambda_graph_V = gnmf_lambda_V
    ),
    supervisedWinit = list(
      algorithm = "Frobenius",
      initU     = build_cluster_mean_W(X_stacked, labels_aligned, clusters)
    )
  )

  message(sprintf(
    "[stageA] nnTensor::NMF: stage=%s algorithm=%s J=%d num.iter=%d ",
    stage, stage_args$algorithm, J, nmf_num_iter),
    sprintf("n_restart=%d X_stacked=%dx%d",
            nmf_n_restart, nrow(X_stacked), ncol(X_stacked)))

  set.seed(1L)
  # supervisedHfix is deterministic (V is pinned). Just 1 restart.
  effective_restart <- if (stage == "supervisedHfix") 1L else nmf_n_restart

  nmf_runs <- lapply(seq_len(effective_restart), function(i) {
    do.call(nnTensor::NMF, c(base_args, stage_args))
  })
  final_err <- vapply(nmf_runs, function(r) tail(r$RecError, 1), numeric(1))
  best_idx  <- which.min(final_err)
  nmf_best  <- nmf_runs[[best_idx]]
  nmf_recerror <- nmf_best$RecError
  message(sprintf("[stageA] best restart %d/%d, final RecError=%.6g",
                  best_idx, effective_restart, final_err[best_idx]))

  init_W_RNA <- Mat2List(X_RNA_for_nmf, nmf_best$U)
  # nnTensor convention: V is (m x J). Machima2 expects init_H_RNA as (J x m).
  init_H_RNA <- t(nmf_best$V)
}

#################################
# Stage B (or joint): Machima2.
#
# joint  : no init; standard joint factorisation
# others : init_W/H from Stage A; fixW_RNA + fixH_RNA = TRUE so only
#          H_Sym (and T if dense) is fit against Hi-C
#################################
if (is_joint_like) {
  fixW_RNA <- FALSE
  fixH_RNA <- FALSE
  iter     <- num_iter
} else {
  fixW_RNA <- TRUE
  fixH_RNA <- TRUE
  iter     <- transfer_num_iter
}

message(sprintf("[stageB] Machima2: stage=%s J=%d num.iter=%d T=%s fixT=%s ",
                stage, J, iter, T_variant, fixT),
        sprintf("fixW_RNA=%s fixH_RNA=%s ", fixW_RNA, fixH_RNA),
        sprintf("T_regularization=%s lambda_T=%g T_rank=%s H_Sym_structure=%s lambda_balance=%g",
                T_regularization, lambda_T,
                if (T_regularization == "low_rank") T_rank else "NA",
                H_Sym_structure, lambda_balance))

machima_args <- list(
  X_RNA            = X_RNA,
  X_Epi            = X_Epi,
  label            = labels_aligned,
  T                = Ts,
  fixT             = fixT,
  fixW_RNA         = fixW_RNA,
  fixH_RNA         = fixH_RNA,
  init_W_RNA       = init_W_RNA,
  init_H_RNA       = init_H_RNA,
  J                = J,
  num.iter         = iter,
  T_regularization = T_regularization,
  lambda_T         = lambda_T,
  H_Sym_structure  = H_Sym_structure,
  lambda_balance   = lambda_balance,
  verbose          = TRUE
)
if (T_regularization == "low_rank") {
  machima_args$T_rank <- T_rank
}
res <- do.call(Machima2, machima_args)

#################################
# Stamp metadata + chrom names so downstream eval can index by name.
#################################
if (is.list(res$W_RNA) && is.null(names(res$W_RNA)) &&
    length(res$W_RNA) == length(chroms)) {
  names(res$W_RNA) <- chroms
}
if (is.list(res$T) && is.null(names(res$T)) &&
    length(res$T) == length(chroms)) {
  names(res$T) <- chroms
}
res$.meta <- list(stage = stage, T_variant = T_variant, J = J,
                  chroms = chroms, nmf_recerror = nmf_recerror,
                  gnmf_lambda_V = if (stage == "supervisedGNMF") gnmf_lambda_V else NULL,
                  T_regularization = T_regularization,
                  lambda_T = lambda_T,
                  T_rank = if (T_regularization == "low_rank") T_rank else NULL,
                  H_Sym_structure = H_Sym_structure,
                  lambda_balance = lambda_balance)

dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(res, out_rds)
message(sprintf("[done] saved %s", out_rds))
