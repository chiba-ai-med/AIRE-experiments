#!/usr/bin/env Rscript
# Diagnostics for the differential-Pearson-near-zero finding.
#
# Part 1: H_Sym inspection
#   For each Machima2 .rds:
#     - assign each component to NPC / CN / other (same majority-vote rule as
#       validate_against_sorted_bulk.R)
#     - print H_Sym (J x J) and summary stats:
#         diag mean per group, within-group off-diag mean, cross-group mean
#     - flag if NPC and CN components are not separable in H_Sym
#
# Part 2: Bonev NPC vs CN intrinsic difference
#   For each chromosome, compute against the two bulks themselves
#   (no Machima2 involvement):
#     - cor(X_NPC, X_CN)                 ~= how similar the two bulks are
#     - ||X_NPC - X_CN||_F / ||X_NPC||_F ~= relative diff magnitude
#     - sd(diff) / sd(X_NPC)             ~= signal-to-shared ratio
#
# Usage:
#   diagnose_hsym_and_bonev_diff.R <output_dir> <atac_dir> <hic_npc_dir>
#       <hic_cn_dir> <cluster_celltype_tsv> <out_prefix>

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: diagnose_hsym_and_bonev_diff.R <output_dir> <atac_dir> ",
       "<hic_npc_dir> <hic_cn_dir> <cluster_celltype_tsv> <out_prefix>")
}
output_dir   <- args[[1]]
atac_dir     <- args[[2]]
npc_dir      <- args[[3]]
cn_dir       <- args[[4]]
cluster_tsv  <- args[[5]]
out_prefix   <- args[[6]]

cells   <- readLines(file.path(atac_dir, "cells.tsv"))
labels  <- read.table(file.path(atac_dir, "labels.tsv"), sep = "\t",
                      header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("barcode", "cell_type"))
labels  <- labels[match(cells, labels$barcode), ]
cluster_per_cell <- labels$cell_type
cluster_map <- read.table(cluster_tsv, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
cluster_to_celltype <- setNames(cluster_map$celltype, cluster_map$cluster)

assign_components <- function(H_RNA) {
  J <- nrow(H_RNA)
  predicted_per_cell <- apply(H_RNA, 2, which.max)
  dominant_clust <- character(J)
  for (j in seq_len(J)) {
    members <- cluster_per_cell[predicted_per_cell == j]
    if (length(members) == 0) {
      fallback_cell <- which.max(H_RNA[j, ])
      dominant_clust[j] <- cluster_per_cell[fallback_cell]
    } else {
      tab <- sort(table(members), decreasing = TRUE)
      dominant_clust[j] <- names(tab)[1]
    }
  }
  ct <- cluster_to_celltype[dominant_clust]
  ct[is.na(ct)] <- "other"
  list(celltype = ct, dominant_clust = dominant_clust)
}

#################################
# Part 1: H_Sym inspection
#################################
rds_files <- Sys.glob(file.path(output_dir, "machima2_brain_*_T*_100000.rds"))
# Drop the legacy "machima2_brain_T{dense,identity}_..." entries (no stage).
rds_files <- rds_files[!grepl("machima2_brain_T(dense|identity)_", rds_files)]

hsym_rows <- list()
detail_log <- character()

for (rds in rds_files) {
  base <- sub("\\.rds$", "", sub(".*machima2_brain_", "", basename(rds)))
  m <- regmatches(base, regexec("^(.*)_T(identity|dense)_\\d+$", base))[[1]]
  if (length(m) != 3) next
  stage <- m[2]; Tvar <- m[3]
  res <- readRDS(rds)
  H_Sym <- as.matrix(res$H_Sym)
  H_RNA <- as.matrix(res$H_RNA)
  J <- nrow(H_Sym)
  asg <- assign_components(H_RNA)
  ct  <- asg$celltype
  is_npc <- ct == "npc"
  is_cn  <- ct == "cn"
  is_other <- ct == "other"

  # Per-cell-group means in H_Sym.
  diag_vec <- diag(H_Sym)
  off_mat  <- H_Sym; diag(off_mat) <- NA

  mean_diag_npc <- mean(diag_vec[is_npc])
  mean_diag_cn  <- mean(diag_vec[is_cn])
  # within-group off-diag
  within_npc <- if (sum(is_npc) >= 2)
                  mean(off_mat[is_npc, is_npc], na.rm = TRUE) else NA_real_
  within_cn  <- if (sum(is_cn) >= 2)
                  mean(off_mat[is_cn, is_cn], na.rm = TRUE) else NA_real_
  cross_np_cn <- if (any(is_npc) && any(is_cn))
                   mean(off_mat[is_npc, is_cn], na.rm = TRUE) else NA_real_
  diag_overall <- mean(diag_vec)
  off_overall  <- mean(off_mat, na.rm = TRUE)

  hsym_rows[[length(hsym_rows) + 1]] <- data.frame(
    stage = stage, T_variant = Tvar, J = J,
    n_npc = sum(is_npc), n_cn = sum(is_cn), n_other = sum(is_other),
    diag_mean_npc = mean_diag_npc, diag_mean_cn = mean_diag_cn,
    within_off_npc = within_npc, within_off_cn = within_cn,
    cross_off_npc_cn = cross_np_cn,
    diag_overall = diag_overall, off_overall = off_overall,
    npc_minus_cn_diag = mean_diag_npc - mean_diag_cn,
    cross_minus_within = cross_np_cn - mean(c(within_npc, within_cn),
                                             na.rm = TRUE)
  )

  # Override Machima2's internal dimnames with comp_<j>+celltype labels so
  # the printed matrix is unambiguous. (Machima2 sometimes assigns NA to a
  # "background absorbing" component; see notes/machima2_hsym_dimnames_na_prompt.md)
  H_print <- H_Sym
  comp_labels <- sprintf("c%d_%s", seq_len(J), ct)
  dimnames(H_print) <- list(comp_labels, comp_labels)

  detail_log <- c(detail_log,
    sprintf("===== %s / T%s   J=%d  npc=%d cn=%d other=%d =====",
            stage, Tvar, J, sum(is_npc), sum(is_cn), sum(is_other)),
    sprintf("Machima2 dimnames (info only): %s",
            paste(replace(dimnames(H_Sym)[[1]], is.na(dimnames(H_Sym)[[1]]),
                          "<NA>"), collapse = ", ")),
    sprintf("our assign per comp: %s", paste(ct, collapse = ", ")),
    "H_Sym (rounded 3 sig, rows/cols = c<j>_<celltype>):",
    paste(capture.output(print(round(H_print, 3))), collapse = "\n"),
    sprintf("diag mean: npc=%.3g  cn=%.3g  overall=%.3g",
            mean_diag_npc, mean_diag_cn, diag_overall),
    sprintf("off-diag mean: within_npc=%.3g  within_cn=%.3g  cross=%.3g",
            within_npc, within_cn, cross_np_cn),
    sprintf("interp: cross-within=%.3g (>0 means NPC-CN are CLOSER than within-group)",
            cross_np_cn - mean(c(within_npc, within_cn), na.rm = TRUE)),
    ""
  )
}

hsym_df <- do.call(rbind, hsym_rows)
hsym_csv <- paste0(out_prefix, "_hsym_summary.csv")
write.csv(hsym_df, hsym_csv, row.names = FALSE)
hsym_log <- paste0(out_prefix, "_hsym_detail.txt")
writeLines(detail_log, hsym_log)
message(sprintf("[hsym] %s  +  %s", hsym_csv, hsym_log))

#################################
# Part 2: Bonev intrinsic NPC vs CN difference
#################################
chroms <- list.files(npc_dir, pattern = "\\.mtx$")
chroms <- sub("\\.mtx$", "", chroms)

upper_tri_vec <- function(M) M[upper.tri(M, diag = TRUE)]

bonev_rows <- list()
for (chr in chroms) {
  np_path <- file.path(npc_dir, paste0(chr, ".mtx"))
  cn_path <- file.path(cn_dir,  paste0(chr, ".mtx"))
  if (!file.exists(np_path) || !file.exists(cn_path)) next
  X_np <- as.matrix(Matrix::readMM(np_path))
  X_cn <- as.matrix(Matrix::readMM(cn_path))
  if (any(dim(X_np) != dim(X_cn))) next
  v_np <- upper_tri_vec(X_np)
  v_cn <- upper_tri_vec(X_cn)
  sd_np <- sd(v_np); sd_cn <- sd(v_cn)
  if (!is.finite(sd_np) || !is.finite(sd_cn) || sd_np == 0 || sd_cn == 0) {
    message(sprintf("[skip bonev/%s] sd_np=%s sd_cn=%s (zero or NA)",
                    chr, format(sd_np), format(sd_cn)))
    next
  }
  fr_np   <- sqrt(sum(X_np^2))
  fr_cn   <- sqrt(sum(X_cn^2))
  fr_diff <- sqrt(sum((X_np - X_cn)^2))
  cor_np_cn <- suppressWarnings(cor(v_np, v_cn, method = "pearson"))
  rel_diff_to_np <- fr_diff / fr_np
  rel_diff_to_cn <- fr_diff / fr_cn
  sd_ratio <- sd(v_np - v_cn) / sd(v_np)
  bonev_rows[[length(bonev_rows) + 1]] <- data.frame(
    chrom = chr,
    cor_npc_cn = cor_np_cn,
    frob_npc = fr_np, frob_cn = fr_cn, frob_diff = fr_diff,
    rel_diff_to_npc = rel_diff_to_np,
    rel_diff_to_cn  = rel_diff_to_cn,
    sd_ratio_diff_to_npc = sd_ratio
  )
  message(sprintf(
    "[bonev/%s] cor=%.4f rel_diff=%.4f sd_ratio=%.4f",
    chr, cor_np_cn, rel_diff_to_np, sd_ratio))
}

bonev_df <- do.call(rbind, bonev_rows)
# Mean row.
mean_row <- data.frame(
  chrom = "MEAN",
  cor_npc_cn       = mean(bonev_df$cor_npc_cn,           na.rm = TRUE),
  frob_npc         = mean(bonev_df$frob_npc,             na.rm = TRUE),
  frob_cn          = mean(bonev_df$frob_cn,              na.rm = TRUE),
  frob_diff        = mean(bonev_df$frob_diff,            na.rm = TRUE),
  rel_diff_to_npc  = mean(bonev_df$rel_diff_to_npc,      na.rm = TRUE),
  rel_diff_to_cn   = mean(bonev_df$rel_diff_to_cn,       na.rm = TRUE),
  sd_ratio_diff_to_npc = mean(bonev_df$sd_ratio_diff_to_npc, na.rm = TRUE)
)
bonev_df <- rbind(bonev_df, mean_row)
bonev_csv <- paste0(out_prefix, "_bonev_intrinsic_diff.csv")
write.csv(bonev_df, bonev_csv, row.names = FALSE)
message(sprintf("[bonev] %s", bonev_csv))

message("[done]")
