#!/usr/bin/env Rscript
# Validate Machima2's per-component reconstruction against held-out
# per-celltype Bonev bulk Hi-C (NPC-only, CN-only).
#
# Idea
# ----
# Machima2 fits to the *combined* NPC+CN bulk; we never gave it the per-
# celltype bulks. After fitting, we can construct a "celltype-g-only"
# reconstruction by zeroing the components that don't correspond to g:
#
#     X_hat_g[k] = (T[k] W[k]) M_g H_Sym M_g^T (T[k] W[k])^T
#
# where M_g is a J x J diagonal mask: M_g[j,j] = 1 if component j corresponds
# to celltype g, else 0. The "correspondence" is inferred from H_RNA: each
# component j's dominant cell is the one with largest H_RNA[j,*]; the
# component's celltype is the celltype of that cell's leiden cluster.
#
# We then compare X_hat_g[k] to the held-out per-celltype bulk X_g_true[k]
# at three angles:
#   pearson  : Pearson correlation of upper-triangle pixels (scale-invariant)
#   cosine   : cosine similarity of upper-triangle pixel vectors
#   rel_frob : ||X_g_true_norm - X_hat_g_norm||_F / ||X_g_true_norm||_F
#              (after each side is divided by its own Frobenius norm)
#
# Differential metric (celltype = "diff")
# ---------------------------------------
# NPC and CN bulks share most of their structure (TADs, compartments,
# distance decay), so absolute-pearson against a single celltype bulk is
# dominated by shared signal even when W has no celltype information. To
# isolate cell-type-*specific* fit we additionally compute, per chrom:
#
#     pearson_diff = cor( upper(X_hat_NPC - X_hat_CN),
#                          upper(X_true_NPC - X_true_CN) )
#
# emitted as a row with celltype = "diff". cosine and rel_frob on the
# difference matrices are reported in the same row.
#
# Outputs a CSV with one row per (chrom, celltype in {npc, cn, diff}) and
# overall summary rows with per-celltype means.
#
# Usage:
#   validate_against_sorted_bulk.R <res.rds> <atac_dir> <hic_npc_dir>
#     <hic_cn_dir> <cluster_celltype_tsv> <out.csv>

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: validate_against_sorted_bulk.R <res.rds> <atac_dir> ",
       "<hic_npc_dir> <hic_cn_dir> <cluster_celltype_tsv> <out.csv>")
}
rds_path     <- args[[1]]
atac_dir     <- args[[2]]
npc_dir      <- args[[3]]
cn_dir       <- args[[4]]
cluster_tsv  <- args[[5]]
out_csv      <- args[[6]]

#################################
# Load Machima output + cell labels.
#################################
res <- readRDS(rds_path)
# R parses standalone `else` on a new line as a separate statement; chain
# everything inline OR wrap in {}.
chroms_all <- {
  if (!is.null(res$.meta$chroms))         res$.meta$chroms
  else if (!is.null(names(res$T)))        names(res$T)
  else if (!is.null(names(res$W_RNA)))    names(res$W_RNA)
  else                                    readLines(file.path(atac_dir, "chroms.txt"))
}
J <- ncol(res$H_RNA)
if (is.null(J)) J <- nrow(res$H_RNA)  # safety; H_RNA is J x m so nrow

H_RNA <- res$H_RNA
H_Sym <- res$H_Sym
W_RNA <- res$W_RNA
T_lst <- res$T
fixT  <- isTRUE(res$.meta$T_variant == "identity")  # informational only

cells   <- readLines(file.path(atac_dir, "cells.tsv"))
labels  <- read.table(file.path(atac_dir, "labels.tsv"), sep = "\t",
                      header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("barcode", "cell_type"))
labels  <- labels[match(cells, labels$barcode), ]
cluster_per_cell <- labels$cell_type

cluster_map <- read.table(cluster_tsv, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
cluster_to_celltype <- setNames(cluster_map$celltype, cluster_map$cluster)

#################################
# Component -> celltype mapping (majority vote)
# For each component j, look at the cells that pick j as their top component
# (predicted = argmax_j H_RNA[j, i] for each cell i). The component's
# dominant cluster is the modal leiden cluster among those cells. Falls
# back to argmax-row when no cell picks the component as top.
#################################
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
component_celltype <- cluster_to_celltype[dominant_clust]
component_celltype[is.na(component_celltype)] <- "other"
message(sprintf("[components] J=%d  celltype assignment (majority vote):", J))
for (j in seq_len(J)) {
  n_members <- sum(predicted_per_cell == j)
  message(sprintf("  comp_%d -> cluster=%s celltype=%s (n_top=%d)",
                  j, dominant_clust[j], component_celltype[j], n_members))
}

build_mask <- function(g) {
  m <- diag(0, J)
  diag(m) <- as.numeric(component_celltype == g)
  m
}
M_npc <- build_mask("npc")
M_cn  <- build_mask("cn")

#################################
# Helpers.
#################################
read_mtx_dense <- function(dir, chr) {
  as.matrix(Matrix::readMM(file.path(dir, paste0(chr, ".mtx"))))
}

upper_tri_vec <- function(M) M[upper.tri(M, diag = TRUE)]

cmp <- function(true, hat) {
  # both are matrices of same shape
  vt <- upper_tri_vec(true); vh <- upper_tri_vec(hat)
  list(
    pearson  = if (sd(vt) == 0 || sd(vh) == 0) NA_real_
               else suppressWarnings(cor(vt, vh, method = "pearson")),
    cosine   = sum(vt * vh) / max(sqrt(sum(vt^2) * sum(vh^2)), .Machine$double.eps),
    rel_frob = {
      tn <- sqrt(sum(true^2)); hn <- sqrt(sum(hat^2))
      if (tn < .Machine$double.eps || hn < .Machine$double.eps) NA_real_
      else sqrt(sum((true / tn - hat / hn)^2)) / 1.0  # normalised both
    }
  )
}

reconstruct_g <- function(W_chr, T_chr, H_Sym, M_g) {
  TW <- if (is.null(T_chr)) W_chr else T_chr %*% W_chr
  G  <- TW %*% M_g
  G %*% H_Sym %*% t(G)
}

#################################
# Per-chrom validation.
#################################
hic_dirs <- list(npc = npc_dir, cn = cn_dir)
rows <- list()
for (chr in chroms_all) {
  W_chr <- if (is.list(W_RNA)) W_RNA[[chr]] else W_RNA
  T_chr <- if (is.list(T_lst)) T_lst[[chr]]  else T_lst
  if (is.null(W_chr)) next

  # First pass: build per-celltype reconstruction and truth for this chrom.
  X_true_list <- list()
  X_hat_list  <- list()
  for (g in c("npc", "cn")) {
    Mg <- if (g == "npc") M_npc else M_cn
    if (sum(diag(Mg)) == 0) {
      message(sprintf("[skip] %s/%s: no components assigned to %s", chr, g, g))
      next
    }
    true_path <- file.path(hic_dirs[[g]], paste0(chr, ".mtx"))
    if (!file.exists(true_path)) next
    X_true <- read_mtx_dense(hic_dirs[[g]], chr)
    if (nrow(X_true) != ncol(X_true) || nrow(X_true) != nrow(W_chr)) next
    if (sum(X_true^2) < .Machine$double.eps) {
      # truly empty target (e.g. chrY with ~0 contacts in single-rep bulk)
      next
    }
    X_true_list[[g]] <- X_true
    X_hat_list[[g]]  <- reconstruct_g(W_chr, T_chr, H_Sym, Mg)
  }

  # Per-celltype rows.
  for (g in names(X_true_list)) {
    Mg <- if (g == "npc") M_npc else M_cn
    metrics <- cmp(X_true_list[[g]], X_hat_list[[g]])
    rows[[length(rows) + 1]] <- data.frame(
      chrom = chr, celltype = g,
      n_components = sum(diag(Mg)),
      pearson  = metrics$pearson,
      cosine   = metrics$cosine,
      rel_frob = metrics$rel_frob
    )
    message(sprintf("[chr=%s g=%s] pearson=%.4f cosine=%.4f rel_frob=%.4f n_comp=%d",
                    chr, g, metrics$pearson, metrics$cosine, metrics$rel_frob,
                    sum(diag(Mg))))
  }

  # Differential row: cancel shared TAD/compartment by subtracting CN from NPC
  # on both sides. Tests cell-type-specific signal in isolation.
  if (all(c("npc", "cn") %in% names(X_true_list))) {
    diff_true <- X_true_list[["npc"]] - X_true_list[["cn"]]
    diff_hat  <- X_hat_list[["npc"]]  - X_hat_list[["cn"]]
    metrics_d <- cmp(diff_true, diff_hat)
    rows[[length(rows) + 1]] <- data.frame(
      chrom = chr, celltype = "diff",
      n_components = sum(diag(M_npc)) + sum(diag(M_cn)),
      pearson  = metrics_d$pearson,
      cosine   = metrics_d$cosine,
      rel_frob = metrics_d$rel_frob
    )
    message(sprintf("[chr=%s g=diff] pearson=%.4f cosine=%.4f rel_frob=%.4f",
                    chr, metrics_d$pearson, metrics_d$cosine, metrics_d$rel_frob))
  }
}

per_chr_df <- do.call(rbind, rows)
if (is.null(per_chr_df) || nrow(per_chr_df) == 0) {
  stop("no per-chrom validation rows -- check component->celltype assignment")
}

#################################
# Mean per celltype + overall.
#################################
mean_per_g <- aggregate(per_chr_df[, c("pearson", "cosine", "rel_frob")],
                         by = list(celltype = per_chr_df$celltype),
                         FUN = mean, na.rm = TRUE)
mean_per_g$chrom <- "MEAN"
mean_per_g$n_components <- NA_integer_
mean_per_g <- mean_per_g[, c("chrom", "celltype", "n_components",
                              "pearson", "cosine", "rel_frob")]

matched_only <- per_chr_df[per_chr_df$celltype %in% c("npc", "cn"), ]
overall <- data.frame(chrom = "OVERALL", celltype = "matched",
                      n_components = sum(component_celltype != "other"),
                      pearson  = mean(matched_only$pearson,  na.rm = TRUE),
                      cosine   = mean(matched_only$cosine,   na.rm = TRUE),
                      rel_frob = mean(matched_only$rel_frob, na.rm = TRUE))

out_df <- rbind(per_chr_df, mean_per_g, overall)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(out_df, out_csv, row.names = FALSE)
message(sprintf("[done] %s", out_csv))
