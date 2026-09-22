#!/usr/bin/env Rscript
# ============================================================================
# fig_scale_uniform.R — 年龄维度面板 + Fig3 轨迹面板（统一标尺）
#
# 统一规则（全仓库通用）：
#   R1. 骨系 UMAP（umap.bone）统一 extent + coord_equal
#   R2. 所有 DriftIndex 小提琴 y 轴统一 c(-0.80, 0.45)
#   R3. 所有 "% Hybrid" 柱/线 y 轴统一 c(0, 100)
#   R4. 同一样本同一颜色：S1 红 / S2 橙 / S3 紫 / S4 绿
#   R5. cell-state 名称统一为 Osteo-biased / Adipo-biased
#
# 本脚本产出（每个面板仅此一个产家；其余面板见下）：
#   02_DriftIndex_age_dimension.png   Fig2C  年龄维度（S3 66M vs S4 31M）
#   02_hybrid_pct_age.png             Fig2C  年龄维度 hybrid %
#   04_trajectory_umap.png            Fig3A  轨迹 UMAP（DriftIndex + pseudotime）
#   04_genes_along_trajectory.png     Fig3C  基因沿 pseudotime（统一 y 轴）
#
# 不在本脚本内（各有专责脚本，避免重复写同一文件）：
#   02_DriftIndex_by_group.png / 02_cell_state_proportion.png /
#   05_cross_species_driftindex.png / 05_cross_species_hybrid_pct.png
#       → fig_strata_rebuild.R（按密度分层 / 供体重制的版本）
#   02b_DriftIndex_by_sample.png / 02b_hybrid_pct_by_sample.png
#       → figS5_axis_uniform.R
#
# 输出：覆盖 results/figures/PNG_hires/ 下对应 PNG（300 dpi）
# 用法：Rscript scripts/fig_scale_uniform.R
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})

# ---- 自动定位项目根目录（脚本位于 <root>/scripts/；可用环境变量 DRIFT_PROJECT_DIR 覆盖）----
base_dir <- local({
  env <- Sys.getenv("DRIFT_PROJECT_DIR", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) return(normalizePath(file.path(dirname(a[1]), ".."), mustWork = FALSE))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), ".."), mustWork = FALSE))
  }
  normalizePath(".", mustWork = FALSE)
})
proc_dir <- file.path(base_dir, "data/processed")
out_dir  <- file.path(base_dir, "results/figures/PNG_hires")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 统一常量 -------------------------------------------------------------
DRIFT_YLIM <- c(-0.80, 0.45)
PCT_YLIM   <- c(0, 100)
PCT_LABEL  <- "% Hybrid (AD+/OS+) cells"

col_sample <- c("S1" = "#E41A1C", "S2" = "#FF7F00", "S3" = "#984EA3", "S4" = "#4DAF4A")

theme_set(theme_classic(base_size = 12))

## ===========================================================================
## PART 1 — 载入人类骨系对象，计算统一 UMAP extent（供 Fig1/Fig3A 共用）
## ===========================================================================
cat("===== [1] 载入人类骨系对象 =====\n")
seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))

um <- as.data.frame(Embeddings(seu, reduction = "umap.bone"))
xlim0 <- range(um[, 1]); ylim0 <- range(um[, 2])
span  <- max(diff(xlim0), diff(ylim0))
pad   <- 0.06 * span
XLIM  <- c(min(xlim0) - pad, max(xlim0) + pad)
YLIM  <- c(min(ylim0) - pad, max(ylim0) + pad)
cat(sprintf("统一 UMAP extent: x[%.2f, %.2f] y[%.2f, %.2f]\n", XLIM[1], XLIM[2], YLIM[1], YLIM[2]))

## ===========================================================================
## PART 2 — Fig2C：年龄维度（S3 紫 / S4 绿），统一 y 轴 + 统一标签
## ===========================================================================
cat("\n===== [2] Fig2C 年龄维度 =====\n")
seu_age <- subset(seu, sample %in% c("S3", "S4"))
seu_age$age_label <- factor(ifelse(seu_age$sample == "S3", "S3 (66M)", "S4 (31M)"),
                            levels = c("S3 (66M)", "S4 (31M)"))
col_age <- c("S3 (66M)" = col_sample[["S3"]], "S4 (31M)" = col_sample[["S4"]])

p <- ggplot(seu_age@meta.data, aes(x = age_label, y = DriftIndex, fill = age_label)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = col_age) +
  scale_y_continuous(limits = DRIFT_YLIM) +
  labs(title = "Age dimension: S3 (66M) vs S4 (31M)",
       x = "Sample", y = "Cell-level DriftIndex") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "02_DriftIndex_age_dimension.png"), p,
       width = 7, height = 6, dpi = 300)
cat("  02_DriftIndex_age_dimension.png\n")

age_hyb <- seu_age@meta.data %>%
  group_by(age_label) %>%
  summarise(Hybrid = 100 * mean(cell_state == "Hybrid"), .groups = "drop")
p <- ggplot(age_hyb, aes(x = age_label, y = Hybrid, fill = age_label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = col_age) +
  scale_y_continuous(limits = PCT_YLIM, expand = c(0, 0)) +
  labs(title = "Hybrid cell % by age", x = "Sample", y = PCT_LABEL) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "02_hybrid_pct_age.png"), p,
       width = 7, height = 6, dpi = 300)
cat("  02_hybrid_pct_age.png\n")

## ===========================================================================
## PART 3 — Fig3 面板（轨迹 UMAP 与 Fig1 同标尺；基因曲线统一 y 轴）
## ===========================================================================
cat("\n===== [3] Fig3 面板 =====\n")
seu_t <- readRDS(file.path(proc_dir, "seurat_gse169396_traj.rds"))
pdf_df <- data.frame(
  umap1      = seu_t[["umap.bone"]]@cell.embeddings[, 1],
  umap2      = seu_t[["umap.bone"]]@cell.embeddings[, 2],
  DriftIndex = seu_t$DriftIndex,
  pseudotime = seu_t$pseudotime
)
pdf_df <- pdf_df[!is.na(pdf_df$pseudotime), ]
cat(sprintf("  轨迹细胞数: %d\n", nrow(pdf_df)))

# Fig3A：与 Fig1 同 extent + coord_equal
p1 <- ggplot(pdf_df, aes(umap1, umap2, color = DriftIndex)) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_color_gradient2(low = "#377EB8", mid = "#FFFFFF", high = "#E41A1C", midpoint = 0) +
  coord_equal(xlim = XLIM, ylim = YLIM) +
  labs(title = "Cell-level DriftIndex (AD \u2212 OS)", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")
p2 <- ggplot(pdf_df, aes(umap1, umap2, color = pseudotime)) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_color_viridis_c(option = "C") +
  coord_equal(xlim = XLIM, ylim = YLIM) +
  labs(title = "slingshot pseudotime", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "04_trajectory_umap.png"), p1 | p2,
       width = 12, height = 6, dpi = 300)
cat("  04_trajectory_umap.png (与 Fig1 同 extent + coord_equal)\n")

# Fig3C：基因沿 pseudotime，统一 y 轴
genes_show <- c("RUNX2", "SP7", "ALPL", "BGLAP", "SOST",
                "PPARG", "CEBPA", "ADIPOQ", "FABP4", "VIM")
genes_show <- intersect(genes_show, rownames(seu_t))
expr <- GetAssayData(seu_t, assay = "SCT", layer = "data")[genes_show, ]
ord  <- order(seu_t$pseudotime, na.last = NA)
expr_ord <- as.matrix(expr)[, ord]
n_win <- 20
win <- floor(seq(1, length(ord), length.out = n_win + 1))
mat <- sapply(genes_show, function(g)
  sapply(1:n_win, function(i) mean(expr_ord[g, win[i]:win[i + 1]])))
cat(sprintf("  基因曲线范围: [%.3f, %.3f]\n", min(mat), max(mat)))
YMAX <- signif(max(mat) * 1.05, 2)
plot_list <- lapply(genes_show, function(g)
  ggplot(data.frame(x = 1:n_win, y = mat[, g]), aes(x, y)) +
    geom_line(color = "#185FA5", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, YMAX)) +
    labs(title = g, x = "Pseudotime window", y = "Expression") +
    theme_minimal(base_size = 9))
ggsave(file.path(out_dir, "04_genes_along_trajectory.png"),
       wrap_plots(plot_list, ncol = 5), width = 14, height = 6, dpi = 300)
cat(sprintf("  04_genes_along_trajectory.png (统一 y 轴 0-%.2f)\n", YMAX))

cat("\n===== 完成 =====\n")
