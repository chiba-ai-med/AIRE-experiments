#!/usr/bin/env Rscript
# Evaluate a Machima2 result against the scRNA-derived cluster labels and
# the bulk Hi-C input.
#
# Inputs:
#   <result.rds>  : Machima2 output (list with W_RNA, H_RNA, H_Sym, T, RecError, RelChange)
#   <atac_dir>    : directory of per-chr ATAC bin .mtx files (also has labels.tsv, cells.tsv)
#   <hic_dir>     : directory of per-chr Hi-C bin .mtx files
#   <out_dir>     : where to write metrics.csv + plots
#
# Metrics:
#   - convergence: RecError / RelChange traces (PDF)
#   - cell-type alignment: argmax(H_RNA, by component) vs scRNA leiden labels
#       reported as a confusion matrix CSV + ARI / NMI single numbers
#   - reconstruction: per-chr Frobenius error between X_Epi[[k]] and
#       reconstructed (T W) H_Sym (T W)^T  (PDF, one panel per chr)
#
# scHi-C is intentionally NOT compared here -- it is qualitative-only per
# project_machima2_validation_strategy memory.

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: evaluate_machima2.R <result.rds> <atac_dir> <hic_dir> <out_dir>")
}
result_rds <- args[[1]]
atac_dir   <- args[[2]]
hic_dir    <- args[[3]]
out_dir    <- args[[4]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

res <- readRDS(result_rds)
chroms <- readLines(file.path(atac_dir, "chroms.txt"))

# Machima2 may drop chroms (e.g. chrM via the all-zero filter in run_machima2)
# and may not carry names on W_RNA / T (Dense T variant: input T is NULL so
# Machima creates the list internally without names). Recover the chrom
# ordering: prefer T's names, then W_RNA's names, else re-derive by dropping
# all-zero chroms from chroms.txt the same way run_machima2 does.
infer_chrom_names <- function(res, chroms, hic_dir) {
  # Newer run_machima2 stamps res$.meta$chroms -- authoritative.
  if (!is.null(res$.meta$chroms)) return(res$.meta$chroms)
  if (is.list(res$T)     && !is.null(names(res$T)))     return(names(res$T))
  if (is.list(res$W_RNA) && !is.null(names(res$W_RNA))) return(names(res$W_RNA))
  if (length(res$W_RNA) == length(chroms)) return(chroms)
  # Drop chroms whose Hi-C bin matrix is all zero -- mirrors run_machima2.R.
  nz <- vapply(chroms, function(chr) {
    m <- Matrix::readMM(file.path(hic_dir, paste0(chr, ".mtx")))
    Matrix::nnzero(m) > 0
  }, logical(1))
  kept <- chroms[nz]
  if (length(res$W_RNA) != length(kept)) {
    stop(sprintf("cannot infer chrom names: W_RNA len=%d, kept chroms=%d",
                 length(res$W_RNA), length(kept)))
  }
  kept
}
chrom_names <- infer_chrom_names(res, chroms, hic_dir)
if (is.list(res$W_RNA) && is.null(names(res$W_RNA))) names(res$W_RNA) <- chrom_names
if (is.list(res$T)     && is.null(names(res$T)))     names(res$T)     <- chrom_names
chroms <- intersect(chroms, chrom_names)
message(sprintf("[chroms] %d for evaluation", length(chroms)))

#################################
# 1. Convergence traces
#################################
trace_df <- data.frame(
  iter      = seq_along(res$RecError),
  RecError  = as.numeric(res$RecError),
  RelChange = as.numeric(res$RelChange)
)
write.csv(trace_df, file.path(out_dir, "convergence.csv"), row.names = FALSE)

p_err <- ggplot(trace_df, aes(iter, RecError)) +
  geom_line() + scale_y_log10() +
  labs(title = "Machima2 reconstruction error", x = "iteration", y = "RecError (log10)") +
  theme_bw()
ggsave(file.path(out_dir, "convergence_RecError.pdf"), p_err, width = 5, height = 3.5)

p_rc <- ggplot(trace_df, aes(iter, RelChange)) +
  geom_line() + scale_y_log10() +
  labs(title = "Machima2 relative change", x = "iteration", y = "RelChange (log10)") +
  theme_bw()
ggsave(file.path(out_dir, "convergence_RelChange.pdf"), p_rc, width = 5, height = 3.5)

#################################
# 2. Cell-type alignment (H_RNA argmax vs scRNA leiden)
#################################
labels_path <- file.path(atac_dir, "labels.tsv")
ari <- NA_real_
if (file.exists(labels_path) && !is.null(res$H_RNA)) {
  H <- res$H_RNA   # J x m
  cells <- readLines(file.path(atac_dir, "cells.tsv"))
  lab <- read.table(labels_path, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("barcode", "cell_type"))
  # Order labels to match cells in H.
  lab <- lab[match(cells, lab$barcode), ]
  predicted <- apply(H, 2, which.max)
  truth     <- lab$cell_type

  cm <- table(predicted = paste0("comp_", predicted), truth = truth)
  write.csv(as.matrix(cm), file.path(out_dir, "confusion_matrix.csv"))

  # ARI without depending on mclust: compute via pair-counting.
  ari <- {
    pairs <- function(x) choose(table(x), 2) |> sum()
    cm_choose <- sum(choose(cm, 2))
    a <- pairs(predicted); b <- pairs(truth)
    n2 <- choose(length(predicted), 2)
    expected <- a * b / n2
    max_idx  <- (a + b) / 2
    if (max_idx == expected) NA_real_ else (cm_choose - expected) / (max_idx - expected)
  }
  message(sprintf("[align] ARI(predicted, leiden) = %.4f", ari))
}

#################################
# 3. Per-chromosome Hi-C reconstruction error
#################################
read_mtx <- function(path) as(Matrix::readMM(path), "CsparseMatrix")

reconstruct_chr <- function(W, H_Sym, T_k = NULL) {
  G <- if (is.null(T_k)) W else T_k %*% W
  G %*% H_Sym %*% t(G)
}

frob <- function(M) sqrt(sum(M^2))

frob_rows <- list()
for (chr in chroms) {
  X_Epi <- read_mtx(file.path(hic_dir, paste0(chr, ".mtx")))
  W_chr <- if (is.list(res$W_RNA)) res$W_RNA[[chr]] else res$W_RNA
  T_chr <- if (is.list(res$T))     res$T[[chr]]     else res$T
  if (is.null(W_chr)) next
  X_hat <- reconstruct_chr(W_chr, res$H_Sym, T_chr)
  X_Epi <- as(X_Epi, "denseMatrix")
  X_hat <- as.matrix(X_hat)
  err <- frob(X_Epi - X_hat) / max(frob(X_Epi), .Machine$double.eps)
  frob_rows[[chr]] <- data.frame(chrom = chr, rel_frob = err,
                                 nnz_epi = Matrix::nnzero(read_mtx(file.path(hic_dir, paste0(chr, ".mtx")))))
  message(sprintf("[recon] %s: rel_frob = %.4f", chr, err))
}
recon_df <- do.call(rbind, frob_rows)
write.csv(recon_df, file.path(out_dir, "reconstruction_per_chr.csv"), row.names = FALSE)

p_rec <- ggplot(recon_df, aes(reorder(chrom, rel_frob), rel_frob)) +
  geom_col() + coord_flip() +
  labs(title = "Per-chromosome relative Frobenius error",
       x = "chromosome", y = "||X - X_hat|| / ||X||") +
  theme_bw()
ggsave(file.path(out_dir, "reconstruction_per_chr.pdf"), p_rec, width = 5, height = 7)

#################################
# 4. Summary metrics CSV
#################################
summary_df <- data.frame(
  metric = c("ari_predicted_vs_leiden",
             "final_RecError",
             "final_RelChange",
             "n_iter",
             "mean_rel_frob_per_chr"),
  value  = c(ari,
             tail(trace_df$RecError, 1),
             tail(trace_df$RelChange, 1),
             nrow(trace_df),
             mean(recon_df$rel_frob, na.rm = TRUE))
)
write.csv(summary_df, file.path(out_dir, "metrics.csv"), row.names = FALSE)
message(sprintf("[done] wrote evaluation to %s", out_dir))
