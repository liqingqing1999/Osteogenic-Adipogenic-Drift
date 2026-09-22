#!/usr/bin/env Rscript
# ============================================================================
# figS3_scale_uniform.R — FigS3（阈值敏感性折线）统一 y 标尺与轴标签
#   与 Fig2C-right / Fig4B 一致：y 轴 "% Hybrid (AD+/OS+) cells"，范围 0-100
# 输出：results/figures/PNG_hires/03_sensitivity_trend.png
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

# 对外口径：按供体着色（S1 红 / S2 橙 / S3 紫 / S4 绿），不再使用疾病组名
col_sample <- c("S1" = "#E41A1C", "S2" = "#FF7F00",
                "S3" = "#984EA3", "S4" = "#4DAF4A")

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
thrs <- c(-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06)
sens <- bind_rows(lapply(thrs, function(thr) {
  seu@meta.data %>%
    mutate(hybrid = AD_drift_ucell > thr & OS_identity_ucell > thr) %>%
    group_by(sample, group) %>%
    summarise(hybrid_pct = 100 * mean(hybrid), .groups = "drop") %>%
    mutate(threshold = thr)
}))
sens$sample_f <- factor(sens$sample, levels = c("S1", "S2", "S3", "S4"))

p <- ggplot(sens, aes(x = threshold, y = hybrid_pct, color = sample_f, group = sample_f)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = col_sample,
                     labels = c("S1 (osteoporotic range)", "S2 (osteopenic range)",
                                "S3 (not stratified)", "S4 (not stratified)")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Hybrid cell % across thresholds (sensitivity)",
       x = "Hybrid threshold (AD & OS > thr)",
       y = "% Hybrid (AD+/OS+) cells",
       color = "Donor") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

ggsave(file.path(out_dir, "03_sensitivity_trend.png"), p, width = 8, height = 6, dpi = 300)
cat("03_sensitivity_trend.png 重渲染完成（y 轴 0-100，标签统一）\n")
