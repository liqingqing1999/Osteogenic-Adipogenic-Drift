#!/usr/bin/env Rscript
# ============================================================================
# figS5_axis_uniform.R — 重制 FigS5 两个面板，统一 x 轴标签
#   问题：(A) x 轴 "Donor (S1/S2 densitometrically stratified; S3/S4 not
#   stratified)"，(B) x 轴只有 "Donor"，两面板不一致。
#   本脚本用同一 x 轴标签重制两面板，其余标尺沿用 fig_scale_uniform.R /
#   fig_strata_rebuild.R（DRIFT_YLIM c(-0.80,0.45)、PCT_YLIM c(0,100)、
#   S1 红 / S2 橙 / S3 紫 / S4 绿）。
# 输出（覆盖）：results/figures/PNG_hires/02b_DriftIndex_by_sample.png
#                results/figures/PNG_hires/02b_hybrid_pct_by_sample.png
# 用法：Rscript scripts/figS5_axis_uniform.R
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
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

DRIFT_YLIM <- c(-0.80, 0.45)
PCT_YLIM   <- c(0, 100)
PCT_LABEL  <- "% Hybrid (AD+/OS+) cells"
XLAB       <- "Donor (S1, S2 densitometrically stratified; S3, S4 not stratified)"
col_sample <- c("S1" = "#E41A1C", "S2" = "#FF7F00", "S3" = "#984EA3", "S4" = "#4DAF4A")
theme_set(theme_classic(base_size = 12))

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
seu$sample_f <- factor(seu$sample, levels = c("S1", "S2", "S3", "S4"))

# --- FigS5A：DriftIndex 小提琴（四例，描述性）
p <- ggplot(seu@meta.data, aes(x = sample_f, y = DriftIndex, fill = sample_f)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = col_sample) +
  scale_y_continuous(limits = DRIFT_YLIM) +
  labs(title = "DriftIndex by donor (all four samples)",
       x = XLAB, y = "Cell-level DriftIndex") +
  theme(legend.position = "none", axis.title.x = element_text(size = 9))
ggsave(file.path(out_dir, "02b_DriftIndex_by_sample.png"), p, width = 8, height = 6, dpi = 300)

# --- FigS5B：Hybrid 比例（四例，描述性）
hyb_s <- seu@meta.data %>%
  group_by(sample_f) %>%
  summarise(Hybrid = 100 * mean(cell_state == "Hybrid"), .groups = "drop")
p <- ggplot(hyb_s, aes(x = sample_f, y = Hybrid, fill = sample_f)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = col_sample) +
  scale_y_continuous(limits = PCT_YLIM, expand = c(0, 0)) +
  labs(title = "Hybrid proportion by donor (all four samples)",
       x = XLAB, y = PCT_LABEL) +
  theme(legend.position = "none", axis.title.x = element_text(size = 9))
ggsave(file.path(out_dir, "02b_hybrid_pct_by_sample.png"), p, width = 8, height = 6, dpi = 300)

cat("FigS5 panels re-rendered with unified x-axis label.\n")
