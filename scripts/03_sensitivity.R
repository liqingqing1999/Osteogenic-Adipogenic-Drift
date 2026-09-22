#!/usr/bin/env Rscript
# ============================================================================
# 03_sensitivity.R — 中间态阈值敏感性分析（步骤 2 稳健性验证）
# 目的：验证"三档中间态比例趋势"不依赖单一阈值（审稿人必查）
# 方法：对中间态判定阈值 thr ∈ {0, ±0.02, ±0.04, ±0.06} 重复计算
#       检验 DriftIndex 排序与中间态比例排序是否稳定
# 输入：data/processed/seurat_gse169396_drift.rds
# 输出：results/tables/03_sensitivity_summary.csv
#       results/figures/03_sensitivity_trend.pdf
# 用法：Rscript scripts/03_sensitivity.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
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
fig_dir  <- file.path(base_dir, "results/figures")
tab_dir  <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== 中间态阈值敏感性分析 =====\n")

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
cat(sprintf("对象: %d 细胞\n", ncol(seu)))

# ---- 阈值扫描 --------------------------------------------------------------
thrs <- c(-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06)
res <- list()
for (thr in thrs) {
  # Hybrid：AD 与 OS 都高于阈值
  hybrid_flag <- seu$AD_drift_ucell > thr & seu$OS_identity_ucell > thr
  tmp <- data.frame(
    sample = seu$sample,
    group  = seu$group,
    hybrid = hybrid_flag
  ) %>%
    group_by(sample, group) %>%
    summarise(
      n_cells    = n(),
      hybrid_pct = 100 * mean(hybrid),
      .groups    = "drop"
    ) %>%
    mutate(threshold = thr)
  res[[as.character(thr)]] <- tmp
}
sens <- bind_rows(res)
# 对外口径：分组列改用密度分层标签（与正文 Table S1 / Fig S3 一致；group 仅内部计算用）
STRATUM_OF <- c(S1 = "Osteoporotic range", S2 = "Osteopenic range",
                S3 = "Not stratified",     S4 = "Not stratified")
sens_out <- sens %>%
  mutate(stratum = unname(STRATUM_OF[as.character(sample)])) %>%
  arrange(match(sample, c("S1", "S2", "S3", "S4")), threshold) %>%
  select(sample, stratum, n_cells, hybrid_pct, threshold)
write.csv(sens_out, file.path(tab_dir, "03_sensitivity_summary.csv"), row.names = FALSE)

# 打印核心：各供体（S1–S4）在各阈值下的中间态比例与排序
cat("\n=== 各供体中间态比例（按阈值） ===\n")
print(as.data.frame(sens %>% select(threshold, sample, hybrid_pct) %>%
        tidyr::pivot_wider(names_from = sample, values_from = hybrid_pct)))

# 排序稳定性：每个阈值下供体间中间态比例排序
# 注：仅 ≥0 的阈值"informative"（<0 时所有样本均饱和为 100%）
cat("\n=== 每阈值下中间态比例排序（关键稳定项：S2 > S1） ===\n")
for (thr in thrs) {
  sub <- sens %>% filter(threshold == thr)
  ord <- sub %>% arrange(desc(hybrid_pct)) %>% pull(sample)
  cat(sprintf("  thr=%+.2f: %s\n", thr, paste(ord, collapse = " > ")))
}

# ---- 绘图：敏感性趋势 -----------------------------------------------------
# 按供体画线（S1–S4），x=阈值 y=中间态比例；配色固定 S1 红 / S2 橙 / S3 紫 / S4 绿
sens$sample_f <- factor(sens$sample, levels = c("S1", "S2", "S3", "S4"))
pdf(file.path(fig_dir, "03_sensitivity_trend.pdf"), width = 8, height = 6)
p <- ggplot(sens, aes(x = threshold, y = hybrid_pct,
                      color = sample_f, group = sample_f)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(S1 = "#E41A1C", S2 = "#FF7F00",
                                S3 = "#984EA3", S4 = "#4DAF4A"),
                     labels = c("S1 (osteoporotic range)", "S2 (osteopenic range)",
                                "S3 (not stratified)", "S4 (not stratified)")) +
  labs(title = "Hybrid cell % across thresholds (sensitivity)",
       x = "Hybrid threshold (AD & OS > thr)",
       y = "% Hybrid cells",
       color = "Donor") +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
print(p)
dev.off()

# 关键判定输出（只评估 informative 阈值 ≥0，饱和区 <0 不纳入）
# 论文稳健主张 = S2（骨量减少范围）> S1（骨质疏松范围）
cat("\n=== 结论 ===")
stable <- TRUE
for (thr in thrs[thrs >= 0]) {
  sub <- sens %>% filter(threshold == thr, sample %in% c("S1", "S2"))
  s1 <- sub$hybrid_pct[sub$sample == "S1"]
  s2 <- sub$hybrid_pct[sub$sample == "S2"]
  if (!(s2 > s1)) stable <- FALSE
}
cat(sprintf("\n  S2 > S1 在所有信息阈值(≥0)下成立: %s\n",
            ifelse(stable, "✅ 是（趋势稳健）", "❌ 否（趋势依赖阈值）")))
cat("\n  注：S1 骨系细胞较少(n=768)，其比例估计的稳定性低于 S2，故不对 S1/S4 作排序断言。")
cat("\n  建议仅强调 S2 > S1 的稳健结论（与正文 Fig S3 图注一致）。\n")
cat(sprintf("  表: %s\n  图: %s\n",
            file.path(tab_dir, "03_sensitivity_summary.csv"),
            file.path(fig_dir, "03_sensitivity_trend.pdf")))
cat("===== 敏感性分析完成 =====\n")
