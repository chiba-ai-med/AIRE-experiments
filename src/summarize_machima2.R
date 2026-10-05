#!/usr/bin/env Rscript
# Generate a single multi-page PDF that walks the reader through
# (a) the Machima2 model + experimental setup, then
# (b) the brain 100kb shake-down results across stages x T variants.
#
# Output text is in Japanese; technical terms (Pearson, Frobenius, ARI,
# H_Sym, Identity T etc.) are kept in English to preserve precision.
# Device: cairo_pdf with Noto Sans CJK JP.
#
# Usage:
#   summarize_machima2.R <eval_root_glob> <out_pdf> [<sorted_glob>]

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
eval_glob   <- if (length(args) >= 1) args[[1]] else "output/eval_brain_*_100000"
out_pdf     <- if (length(args) >= 2) args[[2]] else "output/summary_brain_100000.pdf"
sorted_glob <- if (length(args) >= 3) args[[3]] else "output/sorted_validation_brain_*_100000.csv"

JP_FONT <- "Noto Sans CJK JP"
theme_jp <- function(base = 11) {
  theme_bw(base_size = base, base_family = JP_FONT) +
    theme(text = element_text(family = JP_FONT))
}

#################################
# Collect metrics from each eval directory.
#################################
parse_eval_dir <- function(d) {
  base <- sub("_100000$", "", sub("^eval_brain_", "", basename(d)))
  m <- regmatches(base, regexec("^(.*)_T(identity|dense)$", base))[[1]]
  if (length(m) != 3) return(NULL)
  stage <- m[2]; Tvar <- m[3]

  metrics_path <- file.path(d, "metrics.csv")
  recon_path   <- file.path(d, "reconstruction_per_chr.csv")
  conv_path    <- file.path(d, "convergence.csv")
  if (!file.exists(metrics_path)) return(NULL)

  metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
  v <- function(k) {
    x <- metrics$value[metrics$metric == k]
    if (length(x) == 0) NA_real_ else x[1]
  }
  recon <- if (file.exists(recon_path)) read.csv(recon_path, stringsAsFactors = FALSE) else NULL
  conv  <- if (file.exists(conv_path))  read.csv(conv_path,  stringsAsFactors = FALSE) else NULL

  list(
    stage = stage, T_variant = Tvar,
    ARI = v("ari_predicted_vs_leiden"),
    rel_frob = v("mean_rel_frob_per_chr"),
    RecError = v("final_RecError"),
    n_iter = v("n_iter"),
    recon = recon,
    conv  = conv
  )
}

eval_dirs <- Sys.glob(eval_glob)
results   <- Filter(Negate(is.null), lapply(eval_dirs, parse_eval_dir))
if (length(results) == 0) stop("no eval directories found at ", eval_glob)

summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(stage = r$stage, T_variant = r$T_variant,
             ARI = r$ARI, rel_frob = r$rel_frob,
             RecError = r$RecError, n_iter = r$n_iter,
             stringsAsFactors = FALSE)
}))
summary_df <- summary_df[!is.na(summary_df$ARI), ]
stage_order <- c("joint", "transferFlog", "transferKraw",
                 "supervisedHfix", "supervisedGNMF", "supervisedWinit")
summary_df$stage <- factor(summary_df$stage,
                            levels = intersect(stage_order, summary_df$stage))
summary_df$T_variant <- factor(summary_df$T_variant, levels = c("identity", "dense"))
summary_df <- summary_df[order(summary_df$stage, summary_df$T_variant), ]
message(sprintf("[summary] %d eval dirs parsed", nrow(summary_df)))

#################################
# Page helpers.
#################################
A4_W <- 8.27; A4_H <- 11.69

text_page <- function(title, body, footer = "AIRE-experiments / brain 100kb shake-down") {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.92, height = 0.94))
  grid.text(title, x = 0, y = 1, hjust = 0, vjust = 1,
            gp = gpar(fontsize = 20, fontface = "bold", fontfamily = JP_FONT))
  grid.text(body, x = 0, y = 0.92, hjust = 0, vjust = 1,
            gp = gpar(fontsize = 10, lineheight = 1.45,
                      fontfamily = JP_FONT))
  grid.text(footer, x = 1, y = -0.02, hjust = 1, vjust = 0,
            gp = gpar(fontsize = 8, col = "grey50", fontface = "italic",
                      fontfamily = JP_FONT))
  popViewport()
}

table_page <- function(title, df, body_above = NULL) {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.92, height = 0.94))
  grid.text(title, x = 0, y = 1, hjust = 0, vjust = 1,
            gp = gpar(fontsize = 20, fontface = "bold", fontfamily = JP_FONT))
  if (!is.null(body_above)) {
    grid.text(body_above, x = 0, y = 0.92, hjust = 0, vjust = 1,
              gp = gpar(fontsize = 10, lineheight = 1.4, fontfamily = JP_FONT))
    tbl_y <- 0.55
  } else {
    tbl_y <- 0.7
  }
  tg <- tableGrob(df, rows = NULL,
                  theme = ttheme_minimal(
                    base_size = 9, base_family = JP_FONT,
                    core = list(fg_params = list(hjust = 0, x = 0.02))))
  pushViewport(viewport(x = 0.5, y = tbl_y, width = 1, height = 0.6))
  grid.draw(tg)
  popViewport()
  popViewport()
}

#################################
# Open PDF (cairo_pdf to render Japanese).
#################################
cairo_pdf(out_pdf, width = A4_W, height = A4_H, onefile = TRUE,
          family = JP_FONT)

#################################
# Page 1: Title + scope.
#################################
text_page(
  "Machima2: brain 100 kb 検証 -- モデル + 結果",
  paste(
    "目的",
    "  Machima2 (scRNA + scATAC のマルチオームと bulk Hi-C を結合した",
    "  cell-type deconvolution 手法) を、leiden による cell-type ラベル",
    "  (gold ではないが proxy) が利用可能なデータで検証する。",
    "",
    "本 PDF の構成",
    "  1. 背景: 入力データ、モデル方程式、各因子の意味",
    "  2. 試した 6 つの stage (joint vs transfer vs supervised)",
    "  3. 2 つの training-fit 指標 (ARI / Frobenius reconstruction)",
    "  4. Stage x T variant ごとの結果 (100 kb)",
    "  5. ホールドアウト検証 (sorted bulk + differential Pearson)",
    "  6. 考察と次の論点",
    "",
    "実行設定",
    sprintf("  実行日時      : %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    "  reference     : GSE210747 (10x multiome, mouse cortex E15.5)",
    "  query         : GSE96107 (Bonev/Tanay bulk Hi-C, FACS NPC + CN)",
    "  resolution    : 100 kb (mm10、単一の bin grid)",
    "  J             : 7 (= scRNA leiden res=0.5 のクラスタ数)",
    "  cells         : 2086 multiome バーコード",
    "  bins / chrom  : 168 (chr19) - 1955 (chr1)、22 染色体",
    sep = "\n"
  )
)

#################################
# Page 2: Machima2 model.
#################################
text_page(
  "1. Machima2 モデル",
  paste(
    "1 つの W を共有する 2 つの非負行列分解:",
    "",
    "    X_RNA[k]  ~=  W[k] . H_RNA",
    "    X_Epi[k]  ~=  (T[k] . W[k]) . H_Sym . (T[k] . W[k])^T",
    "",
    "(k は染色体)。各行列の形:",
    "",
    "    X_RNA[k]  : n_k x m   ATAC bin counts (cell x bin、転置)",
    "    X_Epi[k]  : l_k x l_k Hi-C コンタクト行列 (対称)",
    "    W[k]      : n_k x J   bin x コンポーネント (染色体ごと)",
    "    H_RNA     : J x m     コンポーネント x cell (全染色体共通)",
    "    H_Sym     : J x J     対称行列 (全染色体共通)",
    "    T[k]      : l_k x n_k bin grid 変換子。Identity か Dense",
    "",
    "命名上の注意",
    "  本 pipeline では X_RNA スロットに *ATAC* bin counts を入れている",
    "  (Hi-C bin grid に re-bin、Identity T の場合 l_k = n_k 必須)。",
    "  スロット名が API 上 'X_RNA' というだけで中身は ATAC。",
    "",
    "設計の意図",
    "  W は 2 つの方程式の両方に現れる唯一の行列。これが",
    "  cross-modality bridge: 'cell-type コンポーネント j を特徴づける bin'",
    "  が ATAC でも Hi-C でも一貫して説明力を持つよう拘束される。",
    "",
    "  H_RNA は各 cell のコンポーネント load (J x m、ソフトクラスタ割当)。",
    "",
    "  H_Sym は J x J の対称行列で、コンポーネント間の接触強度を表す",
    "  (cell-type i x cell-type j のコンタクト)。",
    sep = "\n"
  )
)

#################################
# Page 3: Six stages.
#################################
stages_tbl <- data.frame(
  stage = c("joint", "transferFlog", "transferKraw",
            "supervisedHfix", "supervisedGNMF", "supervisedWinit"),
  description = c(
    "Machima2 を joint で W, H_RNA, H_Sym, T を同時更新 (オリジナル)",
    "Stage A: nnTensor::NMF on log1p(X_RNA), Frobenius. Stage B: W,H_RNA freeze",
    "Stage A: nnTensor::NMF on raw X_RNA, KL divergence. Stage B: 同上",
    "Stage A: H_RNA を leiden one-hot 固定 (W = cluster mean). Stage B: 同上",
    "Stage A: graph regularised NMF (same-cluster Laplacian). Stage B: 同上",
    "Stage A: cluster mean を init として NMF を unconstrained で精緻化. Stage B: 同上"
  ),
  supervised = c("no", "no", "no", "hard", "soft", "init only"),
  stringsAsFactors = FALSE
)

table_page(
  "2. 6 つの stage (x 2 種の T)",
  stages_tbl,
  body_above = paste(
    "各 stage は (W, H_RNA, H_Sym, T) を生成する。違いは",
    "*Hi-C の H_Sym フィット前に W と H_RNA をどう得るか*。",
    "",
    "Stage = cross-modality basis である W の取得法:",
    sep = "\n"
  )
)

text_page(
  "2 (続). T variant -- bin grid 変換子",
  paste(
    "T[k] : l_k x n_k。各 stage に対して 2 通りを試した:",
    "",
    "  identity : T[k] = diag(l_k)。Hi-C bin grid と ATAC bin grid が",
    "             一致 (l_k = n_k) しているので 1-to-1 対応で十分。",
    "             T の自由度は 0 で、Hi-C の再構築は H_Sym (J x J) と",
    "             W だけで行うことになる。",
    "",
    "  dense    : T[k] を l_k x n_k の密行列として学習する",
    "             (joint stage では Machima2 の RandomEpi init で初期化)。",
    "             染色体ごとに数百万の自由 parameter を追加し、これが",
    "             W と Hi-C 構造との 'ズレ' を吸収できる。",
    "",
    "結果から得られる含意",
    "  Identity T は 'W が両 modality を固定 bridge 経由で同時に満たす'",
    "  必要があり、ARI と rel_frob が鋭く trade-off する。",
    "",
    "  Dense T は 2 modality を実質的に decouple する -- T が (T.W) を",
    "  Hi-C 向きの基底に reshape できるので、W が強く制約されても",
    "  Hi-C 側の fit (rel_frob) はよく出る。ただし W の 'bridge' 解釈",
    "  は弱まる。",
    sep = "\n"
  )
)

#################################
# Page 4: Validation metrics.
#################################
text_page(
  "3. 2 つの training-fit 指標",
  paste(
    "指標 A. ARI (Adjusted Rand Index) -- cell-type 整合性",
    "",
    "  predicted = argmax_j H_RNA[j, i]   (cell i ごと、J クラス)",
    "  truth     = cell i の leiden cluster (7 クラス)",
    "  contingency n[ij] = |{predicted = i かつ truth = j}|",
    "",
    "         ARI = (Index - Expected) / (Max - Expected)",
    "",
    "  Index = sum C(n[ij], 2)、Expected は marginal から、Max = (a+b)/2。",
    "  ARI = 1 で完全一致、0 で chance level、< 0 で chance より悪い。",
    "",
    "  注: ここでの 'truth' 自体が unsupervised (scRNA leiden res=0.5)。",
    "  proxy alignment であって生物学的 ground truth ではない。",
    "",
    "指標 B. 平均 per-chrom relative Frobenius -- Hi-C 再構築",
    "",
    "  X_hat[k]    = (T[k] . W[k]) . H_Sym . (T[k] . W[k])^T",
    "  rel_frob[k] = || X_Epi[k] - X_hat[k] ||_F / || X_Epi[k] ||_F",
    "  報告値は染色体間の平均",
    "",
    "  rel_frob = 0 で完全 fit、~= 1 で X_hat が X_Epi の信号をほぼ",
    "  含まない。これは *training data 上の fit 品質* であり、",
    "  汎化測定でも held-out 比較でもないことに注意。",
    "",
    "2 指標を合わせて読むと",
    "  ARI 高 / rel_frob 低 : sweet spot. RNA も Hi-C も基底が説明",
    "  ARI 高 / rel_frob 高 : RNA に忠実だが Hi-C を説明できない",
    "  ARI 低 / rel_frob 低 : 両 modality fit するが label と乖離",
    "  ARI 低 / rel_frob 高 : 最悪。どちらも説明できていない",
    sep = "\n"
  )
)

#################################
# Page 5: Results table.
#################################
fmt_tbl <- summary_df
fmt_tbl$ARI       <- sprintf("%+.3f", fmt_tbl$ARI)
fmt_tbl$rel_frob  <- sprintf("%.3f",  fmt_tbl$rel_frob)
fmt_tbl$RecError  <- sprintf("%.2f",  fmt_tbl$RecError)
fmt_tbl$n_iter    <- as.integer(fmt_tbl$n_iter)

table_page(
  "4. 結果一覧",
  fmt_tbl,
  body_above = paste(
    "stage x T variant 順。ARI は leiden との一致 (1 が best)、",
    "rel_frob は染色体平均 Hi-C 再構築誤差 (0 が best、1 で no signal)。",
    "RecError は Machima2 の最終 beta-divergence 値 (低いほど良いが",
    "loss 形が stage で異なるので stage 間の直接比較は不可)。",
    sep = "\n"
  )
)

#################################
# Page 6: Pareto scatter.
#################################
p_pareto <- ggplot(summary_df,
                    aes(x = ARI, y = rel_frob,
                        colour = stage, shape = T_variant)) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", colour = "grey80") +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey80") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(aes(label = paste0(stage, "/", T_variant)),
            size = 2.6, vjust = -1.2, hjust = 0.5,
            family = JP_FONT, show.legend = FALSE) +
  scale_x_continuous(limits = c(-0.1, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_shape_manual(values = c(identity = 16, dense = 17)) +
  labs(title = "5. Pareto plot -- ARI vs Hi-C 再構築",
       subtitle = "右上ほど cell-type ラベル一致、下ほど Hi-C fit 良。右下が sweet spot。",
       x = "ARI (predicted argmax(H_RNA) vs leiden)",
       y = "平均 per-chrom rel_frob (低いほど Hi-C fit 良)",
       colour = "Stage", shape = "T variant",
       caption = "Identity T (○) は trade-off が鋭い、Dense T (△) は W 制約に関わらず Hi-C を回復する") +
  theme_jp(11) +
  theme(legend.position = "right")
print(p_pareto)

#################################
# Page 7: per-chrom rel_frob heatmap.
#################################
recon_long <- do.call(rbind, lapply(results, function(r) {
  if (is.null(r$recon)) return(NULL)
  data.frame(stage = r$stage, T_variant = r$T_variant,
             chrom = r$recon$chrom, rel_frob = r$recon$rel_frob,
             stringsAsFactors = FALSE)
}))
if (!is.null(recon_long) && nrow(recon_long) > 0) {
  recon_long$stage <- factor(recon_long$stage,
                              levels = intersect(stage_order, unique(recon_long$stage)))
  recon_long$T_variant <- factor(recon_long$T_variant, levels = c("identity", "dense"))
  recon_long$run <- paste(recon_long$stage, recon_long$T_variant, sep = "/")
  chrom_order <- c(paste0("chr", c(1:19)), "chrX", "chrY")
  recon_long$chrom <- factor(recon_long$chrom,
                              levels = intersect(chrom_order, unique(recon_long$chrom)))

  p_heat <- ggplot(recon_long, aes(chrom, run, fill = pmin(rel_frob, 1))) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", rel_frob)),
              size = 2.3, colour = "white", family = JP_FONT) +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                         name = "rel_frob\n(1 で頭打ち)") +
    labs(title = "6. 染色体ごとの Hi-C 再構築誤差",
         subtitle = "各セル = ||X - X_hat||_F / ||X||_F (0 が best)。chrY ~= 1 は Bonev cortex で chrY コンタクトが ~7 とほぼ皆無 (input ~= 0) のため。",
         x = NULL, y = NULL) +
    theme_jp(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  print(p_heat)
}

#################################
# Page 8: convergence (RecError) -- facet by run.
#################################
conv_long <- do.call(rbind, lapply(results, function(r) {
  if (is.null(r$conv)) return(NULL)
  data.frame(stage = r$stage, T_variant = r$T_variant,
             iter = r$conv$iter, RecError = r$conv$RecError,
             stringsAsFactors = FALSE)
}))
if (!is.null(conv_long) && nrow(conv_long) > 0) {
  conv_long$run <- paste(conv_long$stage, conv_long$T_variant, sep = "/")
  conv_long$stage <- factor(conv_long$stage,
                             levels = intersect(stage_order, unique(conv_long$stage)))

  p_conv <- ggplot(conv_long, aes(iter, RecError, colour = stage,
                                    linetype = T_variant, group = run)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    scale_y_log10() +
    scale_linetype_manual(values = c(identity = "dashed", dense = "solid")) +
    labs(title = "7. Stage B 収束 (Machima2 RecError、log scale)",
         subtitle = "joint = 30 iter、transfer/supervised = 100 iter (W, H_RNA を freeze、H_Sym + T のみ更新)",
         x = "イテレーション",
         y = "RecError (log10)") +
    theme_jp(11) +
    theme(legend.position = "right")
  print(p_conv)
}

#################################
# Page 9: Held-out per-celltype validation (sorted-bulk + differential).
#################################
sorted_files <- Sys.glob(sorted_glob)
sorted_long <- NULL
for (f in sorted_files) {
  base <- sub("\\.csv$", "", sub("^sorted_validation_brain_", "", basename(f)))
  m <- regmatches(base, regexec("^(.*)_T(identity|dense)_\\d+$", base))[[1]]
  if (length(m) != 3) next
  stage <- m[2]; Tvar <- m[3]
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) next
  per_chr <- df[!df$chrom %in% c("MEAN", "OVERALL"), ]
  if (nrow(per_chr) == 0) next
  per_chr$stage <- stage
  per_chr$T_variant <- Tvar
  sorted_long <- rbind(sorted_long, per_chr)
}

if (!is.null(sorted_long) && nrow(sorted_long) > 0) {
  sorted_long$stage <- factor(sorted_long$stage,
                               levels = intersect(stage_order, unique(sorted_long$stage)))
  sorted_long$T_variant <- factor(sorted_long$T_variant, levels = c("identity", "dense"))
  sorted_long$celltype  <- factor(sorted_long$celltype,  levels = c("npc", "cn", "diff"))
  sorted_long$run <- paste(sorted_long$stage, sorted_long$T_variant, sep = "/")

  mean_per_run <- aggregate(pearson ~ stage + T_variant + celltype,
                             data = sorted_long, FUN = mean, na.rm = TRUE)
  mean_per_run$run <- paste(mean_per_run$stage, mean_per_run$T_variant, sep = "/")

  p_sorted_pearson <- ggplot(mean_per_run,
                              aes(x = run, y = pearson, fill = celltype)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_hline(yintercept = 0, colour = "grey60") +
    scale_fill_manual(values = c(npc = "#377eb8", cn = "#e41a1c", diff = "#4daf4a"),
                      labels = c(npc = "NPC (matched)",
                                 cn  = "CN (matched)",
                                 diff = "NPC - CN (差分)")) +
    scale_y_continuous(limits = c(min(c(0, mean_per_run$pearson), na.rm = TRUE) - 0.05,
                                    max(mean_per_run$pearson, na.rm = TRUE) + 0.05)) +
    labs(title = "8. ホールドアウト検証 (matched + 差分)",
         subtitle = paste(
           "matched: g = NPC, CN について Pearson(X_hat_g, X_g_true)。NPC と CN bulk が共有する TAD/compartment が支配的なので、",
           "高くても 'Hi-C らしい再構築' を意味するだけで、cell-type 特異性とは別物。",
           "差分: Pearson(X_hat_NPC - X_hat_CN, X_NPC - X_CN)。共有構造を打ち消すので、cell-type 特異信号だけが寄与する。",
           "判断は差分バーで。",
           sep = "\n"),
         x = NULL, y = "Pearson 相関 (染色体平均)",
         fill = "比較",
         caption = "training input は combined NPC+CN bulk、per-celltype bulk はホールドアウト。コンポーネント -> NPC/CN は cluster_celltype_map.tsv 経由。") +
    theme_jp(10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          plot.subtitle = element_text(size = 8, lineheight = 1.1, family = JP_FONT))
  print(p_sorted_pearson)
} else {
  text_page(
    "8. ホールドアウト検証",
    paste(
      "sorted_validation_brain_*_*.csv が見つからない:",
      sprintf("  %s", sorted_glob),
      "",
      "validate_against_sorted_bulk が未実行 (または glob にマッチしない)。",
      "次を実行してから PDF を再生成:",
      "",
      "  bash workflow/run_evaluate.sh",
      sep = "\n"
    )
  )
}

#################################
# Page 10: Discussion.
#################################
disc <- paste(
  "Pareto plot (training-fit 指標) からわかること",
  "",
  "  Identity T 系は rel_frob ~= 0.93-0.97 の上帯に散らばる。W が",
  "  制約された stage (supervisedHfix は cluster centroid に固定、",
  "  supervisedWinit と transferKraw もそこから出発) では H_Sym (28",
  "  free params) だけでは百万 pixel の Hi-C を fit しきれない。",
  "",
  "  Dense T 系は stage を問わず rel_frob ~= 0.12 に集まる。T (染色体",
  "  ごと l_k x n_k) が W 制約を吸収して (T.W) を Hi-C 向き低 rank 基底",
  "  に整形できるため。Hi-C fit には良いが W の 'bridge' 解釈は弱まる。",
  "",
  "ホールドアウト検証 (matched + 差分) からわかること",
  "",
  "  training に使ったのは combined (NPC+CN) bulk。per-celltype bulk",
  "  (NPC-only / CN-only) はホールドアウトなので、matched Pearson は",
  "  overfitting ではなく真の汎化測定。",
  "",
  "  ただし NPC bulk と CN bulk は TAD/compartment/distance decay の",
  "  大部分を共有している。matched Pearson はこの共有信号で支配される",
  "  ので、'NPC 再構築が NPC bulk と一致' '同 CN' が ~0.9 だけでは",
  "  cell-type 区別ができていることにはならない。",
  "",
  "  差分 Pearson は共有部分を打ち消す:",
  "      cor( X_hat_NPC - X_hat_CN ,  X_NPC^true - X_CN^true ).",
  "  cell-type 特異な W を生成できた stage はここでスコアが上がる。",
  "  本実験の結果: 12 stage 全てで差分 ~= 0 (range -0.11 - +0.12)。",
  "  どの stage も cell-type 区別はできていない。",
  "",
  "原因の候補 (診断進行中)",
  "",
  "  (a) Bonev NPC/CN bulk 自体の cell-type 特異性が小さい",
  "      - 100 kb 解像度では TAD 境界の差は粗い",
  "      - (確認中: src/diagnose_hsym_and_bonev_diff.R で",
  "         cor(X_NPC, X_CN), ||X_NPC - X_CN|| 等を測定)",
  "",
  "  (b) Machima2 の H_Sym が NPC コンポーネント間と CN コンポーネント",
  "      間で似た値に収束してしまっている",
  "      - 確認中: 同じ診断 script で H_Sym の within/cross 比較",
  "",
  "  (c) コンポーネント -> celltype の majority-vote 割当てがノイジー",
  "      - J=7 で 'other' に落ちるコンポーネントがある場合に問題化",
  "",
  "Caveats",
  "",
  "  1. 'training truth' = leiden(scRNA, res=0.5)。unsupervised proxy。",
  "     marker score での annotation 化や res sweep で硬化させる余地。",
  "  2. J = 7 は leiden cluster 数固定。J を sweep すれば ARI 上限と",
  "     stage 間 geometry が変わる可能性。",
  "  3. Stage A NMF: 5 restart x 200 iter (HfixV のみ 1 restart)。",
  "     supervisedHfix は決定的、それ以外は init の確率性が残る。",
  "  4. コンポーネント -> celltype 割当ては majority vote (J=7 では",
  "     通常 clean だが境界 stage で flip しうる)。",
  "",
  "次の検討事項",
  "",
  "  - 診断結果 (Bonev intrinsic diff + H_Sym 構造) が出てから、",
  "    モデル側 / データ側のどちらが原因か切り分け",
  "  - 解像度を上げて (25 kb / 10 kb) 再評価。100 kb で cell-type 信号が",
  "    元々小さい場合にこれで救われる可能性",
  "  - gnmf_lambda_V (現状 1.0) の sweep で soft-supervised の sub-Pareto",
  "  - leiden 'truth' を marker-based annotation に置換 (E15.5 cortex 用",
  "    の panel または SingleR)",
  "  - Bonev replicate 群 (GSM2533836-8 NPC, GSM2533840-2 CN) を取り込み、",
  "    per-celltype target の single-rep ノイズを下げる",
  sep = "\n"
)
text_page("9. 考察と次の検討事項", disc)

dev.off()
message(sprintf("[done] %s", out_pdf))
