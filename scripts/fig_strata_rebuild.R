#!/usr/bin/env Rscript
# ============================================================================
# fig_strata_rebuild.R — 图件改版：把 S3/S4 移出疾病分层后的重制
#
# 触发：S4（31 y M）的 T 值按 ISCD 不用于 <50 岁男性诊断 → 疾病分层收敛为
#       S1（骨质疏松范围）vs S2（骨量减少范围）；S3/S4 仅作描述性参照。
#
# 输出（覆盖 results/figures/PNG_hires/；每个面板仅此一个产家）：
#   02_DriftIndex_by_group.png       Fig2A  ← 按密度分层（S2 骨量减少 vs S1 骨质疏松）
#   02_cell_state_proportion.png     Fig2B  ← 同上
#   05_cross_species_driftindex.png  Fig4A  ← 人样本按 S1–S4 标注
#   05_cross_species_hybrid_pct.png  Fig4B  ← 同上
#
# 注：FigS5 面板（02b_DriftIndex_by_sample / 02b_hybrid_pct_by_sample）由
#     figS5_axis_uniform.R 产出（该脚本负责统一两面板的 x 轴标签）；
#     本脚本不再重复渲染，以免运行顺序影响成品。
#
# 标尺沿用 fig_scale_uniform.R：DRIFT_YLIM c(-0.80,0.45)、PCT_YLIM c(0,100)、
# 同一样本同一颜色 S1 红 / S2 橙 / S3 紫 / S4 绿。
# 用法：Rscript scripts/fig_strata_rebuild.R
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

LOG <- file.path(base_dir, "logs", "fig_strata_rebuild.log")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, "w"); wr <- function(...) { writeLines(paste0(...), con); flush(con) }
wr("=== fig_strata_rebuild start ===")

DRIFT_YLIM <- c(-0.80, 0.45)
PCT_YLIM   <- c(0, 100)
PCT_LABEL  <- "% Hybrid (AD+/OS+) cells"

col_stratum <- c("Osteopenic range (S2)" = "#FF7F00", "Osteoporotic range (S1)" = "#E41A1C")
col_sample  <- c("S1" = "#E41A1C", "S2" = "#FF7F00", "S3" = "#984EA3", "S4" = "#4DAF4A")
col_state   <- c("Uncommitted" = "#BDBDBD", "Osteo-biased" = "#377EB8",
                 "Adipo-biased" = "#984EA3", "Hybrid" = "#E41A1C")
theme_set(theme_classic(base_size = 12))

## ===========================================================================
## PART 1 — Fig2A / Fig2B：两个可分层供体（S2 骨量减少范围 vs S1 骨质疏松范围）
## ===========================================================================
seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
wr("human bone-lineage cells: ", ncol(seu))

STRATUM_OF <- c(S1 = "Osteoporotic range (S1)", S2 = "Osteopenic range (S2)")
seu$stratum <- unname(STRATUM_OF[as.character(seu$sample)])
seu_s <- subset(seu, cells = colnames(seu)[!is.na(seu$stratum)])
seu_s$stratum <- factor(seu_s$stratum, levels = c("Osteopenic range (S2)", "Osteoporotic range (S1)"))
wr("stratified cells: ", ncol(seu_s), "  (", paste(names(table(seu_s$stratum)), table(seu_s$stratum), sep = "=", collapse = ", "), ")")

# --- Fig2A：DriftIndex 小提琴（两档）
p <- ggplot(seu_s@meta.data, aes(x = stratum, y = DriftIndex, fill = stratum)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = col_stratum) +
  scale_y_continuous(limits = DRIFT_YLIM) +
  labs(title = "DriftIndex by densitometric stratum",
       x = "Densitometric stratum (lumbar T-score)", y = "Cell-level DriftIndex") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "02_DriftIndex_by_group.png"), p, width = 8, height = 6, dpi = 300)
wr("  02_DriftIndex_by_group.png  (2 strata)")

# --- Fig2B：cell-state 堆叠条（两档）
hs <- seu_s@meta.data %>%
  group_by(stratum) %>%
  summarise(`Osteo-biased` = 100 * sum(cell_state == "Osteo-biased") / n(),
            `Adipo-biased` = 100 * sum(cell_state == "Adipo-biased") / n(),
            Uncommitted    = 100 * sum(cell_state == "Uncommitted") / n(),
            Hybrid         = 100 * sum(cell_state == "Hybrid") / n(),
            .groups = "drop") %>%
  tidyr::pivot_longer(-stratum, names_to = "state", values_to = "pct")
hs$state <- factor(hs$state, levels = c("Uncommitted", "Osteo-biased", "Adipo-biased", "Hybrid"))
writeLines(capture.output(print(as.data.frame(hs))), con); flush(con)
p <- ggplot(hs, aes(x = stratum, y = pct, fill = state)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_state) +
  scale_y_continuous(limits = PCT_YLIM, expand = c(0, 0)) +
  labs(title = "Cell state proportions by densitometric stratum",
       x = "Densitometric stratum (lumbar T-score)", y = "% of bone-lineage cells") +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "02_cell_state_proportion.png"), p, width = 8, height = 6, dpi = 300)
wr("  02_cell_state_proportion.png  (2 strata)")

## ===========================================================================
## PART 2 — Fig4：跨物种对照，人样本按 S1–S4 标注（不再用疾病组名）
##   （FigS5 的 02b_* 面板由 figS5_axis_uniform.R 负责，此处不再渲染）
## ===========================================================================
seu$sample_f <- factor(seu$sample, levels = c("S1", "S2", "S3", "S4"))  # 供体重排（本 PART 用）
mb <- readRDS(file.path(proc_dir, "seurat_mouse_metaphysis_bone.rds"))
mouse_hyb <- 100 * mean(mb$cell_state_mouse == "Hybrid")
wr("mouse hybrid %: ", round(mouse_hyb, 2))

LV <- c("S1", "S2", "S3", "S4", "Mouse (healthy)")
col_sp <- c(col_sample, "Mouse (healthy)" = "#377EB8")

plot_df <- data.frame(
  cat        = c(as.character(seu$sample), rep("Mouse (healthy)", ncol(mb))),
  DriftIndex = c(seu$DriftIndex, mb$DriftIndex_mouse)
)
plot_df$cat <- factor(plot_df$cat, levels = LV)
p <- ggplot(plot_df, aes(x = cat, y = DriftIndex, fill = cat)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  scale_fill_manual(values = col_sp) +
  scale_y_continuous(limits = DRIFT_YLIM) +
  labs(title = "Cross-species comparison: DriftIndex",
       x = "Human donor / species", y = "Cell-level DriftIndex") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))
ggsave(file.path(out_dir, "05_cross_species_driftindex.png"), p, width = 8, height = 6, dpi = 300)
wr("  05_cross_species_driftindex.png  (S1-S4 + mouse)")

hyb_by_sample <- seu@meta.data %>%
  group_by(sample_f) %>%
  summarise(pct = 100 * mean(cell_state == "Hybrid"), .groups = "drop")
bar_df <- data.frame(
  label = factor(c("S1", "S2", "S3", "S4", "Mouse (healthy)"), levels = LV),
  hybrid_pct = c(hyb_by_sample$pct[match(c("S1", "S2", "S3", "S4"), hyb_by_sample$sample_f)], mouse_hyb)
)
writeLines(capture.output(print(bar_df)), con); flush(con)
p <- ggplot(bar_df, aes(x = label, y = hybrid_pct, fill = label)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = mouse_hyb, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = col_sp) +
  scale_y_continuous(limits = PCT_YLIM, expand = c(0, 0)) +
  labs(title = "Hybrid cell proportion: cross-species",
       x = "Human donor / species", y = PCT_LABEL) +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))
ggsave(file.path(out_dir, "05_cross_species_hybrid_pct.png"), p, width = 8, height = 6, dpi = 300)
wr(sprintf("  05_cross_species_hybrid_pct.png  (mouse dashed line %.2f%%)", mouse_hyb))

wr("=== done ===")
close(con)
