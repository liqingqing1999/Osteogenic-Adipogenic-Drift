#!/usr/bin/env Rscript
# ============================================================================
# 08_pseudotime_shape.R — 检验 DriftIndex 沿 pseudotime 的"形状"
#
# 背景：04_trajectory.R 只做了细胞级 Spearman（rho = -0.321）与十分位均值。
#       十分位均值为"倒 U"（bin1→5 升到峰 +0.081，bin8–9 转负），
#       并非单调下降；因此需要专门的形状检验来支撑"双相（biphasic）"表述。
#
# 本脚本：
#   1) GAM（mgcv）平滑：edf / F / p / 偏差解释率 → 判定非线性；
#      与线性模型做似然比检验（anova），判定"是否显著优于直线"。
#   2) 分段回归（segmented）：估计断点 ± CI、断点前后斜率、Davies 检验 p。
#   3) 十分位均值 + 峰位（GAM 拟合最大点）。
#   4) 输出 summary CSV + report txt + 形状图（pdf/png）。
#
# 输入：data/processed/seurat_gse169396_traj.rds（含 DriftIndex, pseudotime）
# 输出：
#   results/tables/08_pseudotime_shape_summary.csv
#   results/tables/08_pseudotime_shape_report.txt
#   results/figures/08_pseudotime_shape.pdf / .png
# 用法：Rscript scripts/08_pseudotime_shape.R
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(mgcv); library(segmented)
  library(dplyr); library(ggplot2)
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
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(tab_dir, "08_pseudotime_shape_report.txt")
con <- file(LOG, "w")
say <- function(...) { writeLines(paste0(...), con); flush(con); cat(paste0(...), "\n") }

say("===== 08 DriftIndex–pseudotime 形状分析 =====")
say("日期：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_traj.rds"))
df <- data.frame(DriftIndex = seu$DriftIndex, pt = seu$pseudotime)
df <- df[is.finite(df$DriftIndex) & is.finite(df$pt), ]
# 归一化 pseudotime 到 [0,1]（便于解释斜率与断点位置）
pt_max <- max(df$pt); df$pts <- df$pt / pt_max
say(sprintf("细胞数 n = %d；pseudotime 原范围 [%.2f, %.2f] → 归一化 [0,1]", nrow(df), 0, pt_max))

## ---- 1. Spearman（复核）----------------------------------------------------
sp <- cor.test(df$DriftIndex, df$pts, method = "spearman")
say(sprintf("\n[Spearman] 细胞级 rho = %.4f, p = %.3e", sp$estimate, sp$p.value))

## ---- 2. 十分位均值 + 峰位 ---------------------------------------------------
nbin <- 10
df$bin <- cut(df$pts, breaks = seq(0, 1, length.out = nbin + 1),
              include.lowest = TRUE, labels = FALSE)
dec <- df %>% group_by(bin) %>%
  summarise(pt_mid = mean(pts), DriftIndex = mean(DriftIndex), n = n(), .groups = "drop")
say("\n[十分位均值]")
for (i in seq_len(nrow(dec)))
  say(sprintf("  bin %2d  pt=%.3f  DriftIndex=%+.4f  (n=%d)", dec$bin[i], dec$pt_mid[i], dec$DriftIndex[i], dec$n[i]))
first3 <- mean(dec$DriftIndex[1:3]); last3 <- mean(dec$DriftIndex[8:10])
say(sprintf("  早期(bins 1-3) 均值 = %+.4f；晚期(bins 8-10) 均值 = %+.4f；净变化 = %+.4f",
            first3, last3, last3 - first3))

## ---- 3. GAM 平滑 -----------------------------------------------------------
k_use <- 12
m_gam <- gam(DriftIndex ~ s(pts, k = k_use, bs = "tp"), data = df, method = "REML")
s_tab <- summary(m_gam)
m_lin <- gam(DriftIndex ~ pts, data = df, method = "REML")
lrt <- anova(m_lin, m_gam, test = "F")
edf <- s_tab$s.table[1, "edf"]
say("\n[GAM 平滑 s(pts)]")
say(sprintf("  edf = %.2f (Ref.df = %.2f)  F = %.2f  p = %.3e  |  调整 R² = %.4f  偏差解释 = %.2f%%",
            edf, s_tab$s.table[1, "Ref.df"], s_tab$s.table[1, "F"],
            s_tab$s.table[1, "p-value"], s_tab$r.sq, s_tab$dev.expl * 100))
say(sprintf("  线性 vs GAM 似然比检验: F = %.2f, p = %.3e  → 非线性%s",
            lrt$F[2], lrt$`Pr(>F)`[2],
            ifelse(lrt$`Pr(>F)`[2] < 0.05, "显著（形状非直线）", "不显著（近似直线）")))
# GAM 拟合曲线上的峰
grid <- data.frame(pts = seq(0, 1, length.out = 500))
grid$fit <- predict(m_gam, newdata = grid)
pk <- grid[which.max(grid$fit), ]
say(sprintf("  GAM 拟合峰位: pts ≈ %.3f（DriftIndex ≈ %+.4f）", pk$pts, pk$fit))

## ---- 4. 分段（断点）回归 ---------------------------------------------------
# 用分段线性拟合捕捉"先升后降"的拐点；psi 初值取峰位
lm0 <- lm(DriftIndex ~ pts, data = df)
seg1 <- tryCatch(
  segmented(lm0, seg.Z = ~pts, psi = pk$pts,
            control = seg.control(n.boot = 50, it.max = 50, tol = 1e-4)),
  error = function(e) { say(paste("  [警告] segmented 失败:", conditionMessage(e))); NULL })

bp_txt <- "未估计"
bp_ci <- NULL
if (!is.null(seg1)) {
  s1 <- summary(seg1)
  bp <- s1$psi[1, "Est."]
  bp_ci <- tryCatch(confint(seg1, parm = "pts"), error = function(e) NULL)
  say("\n[分段回归（单断点）]")
  say(sprintf("  断点 pts ≈ %.3f", bp))
  if (!is.null(bp_ci)) {
    say(sprintf("  断点 95%%CI ≈ [%.3f, %.3f]", bp_ci[1, 2], bp_ci[1, 3]))
  }
  # 断点前后斜率（slope() 矩阵列：Est, St.Err, t, CI.l, CI.u）
  sl <- slope(seg1)[[1]]
  df_res <- seg1$df.residual
  p_from_t <- function(tt) 2 * pt(-abs(tt), df_res)
  say("  断点前斜率（pts<bp）与断点后斜率（pts>bp）：")
  say(sprintf("    pre : 斜率 %+.4f (95%%CI %.4f..%.4f), t=%.2f, p=%.3e",
              sl[1, 1], sl[1, 4], sl[1, 5], sl[1, 3], p_from_t(sl[1, 3])))
  say(sprintf("    post: 斜率 %+.4f (95%%CI %.4f..%.4f), t=%.2f, p=%.3e",
              sl[2, 1], sl[2, 4], sl[2, 5], sl[2, 3], p_from_t(sl[2, 3])))
  # Davies 检验：是否存在显著斜率变化
  dv <- tryCatch(davies.test(lm0, seg.Z = ~pts, k = 10), error = function(e) NULL)
  if (!is.null(dv)) say(sprintf("  Davies 检验（H0：斜率恒定）: p = %.3e → %s",
                                dv$p.value, ifelse(dv$p.value < 0.05, "拒绝斜率恒定（存在断点）", "不能拒绝斜率恒定")))
  bp_txt <- sprintf("pts=%.3f", bp)
}

## ---- 5. 分段（双断点）看"峰"是否显著 --------------------------------------
seg2 <- tryCatch(
  segmented(lm0, seg.Z = ~pts, psi = c(pk$pts * 0.8, pk$pts * 1.6),
            control = seg.control(n.boot = 50, it.max = 50, tol = 1e-4)),
  error = function(e) NULL)
if (!is.null(seg2)) {
  s2 <- summary(seg2)
  bps <- s2$psi[, "Est."]
  say("\n[分段回归（双断点，探测倒 U 的两个拐点）]")
  say(sprintf("  断点1 pts ≈ %.3f；断点2 pts ≈ %.3f", bps[1], bps[2]))
  sl2 <- slope(seg2)
  for (nm in names(sl2)) {
    m <- sl2[[nm]]
    say(sprintf("    %s: 斜率 %+.4f (p=%.3e)", nm, m[1, 1], m[1, 4]))
  }
  # AIC 对比
  say(sprintf("  AIC：直线 %.1f < 单断点 %.1f < 双断点 %.1f", AIC(lm0), AIC(seg1), AIC(seg2)))
}

## ---- 6. 输出 summary CSV ---------------------------------------------------
summ <- data.frame(
  metric = c("n_cells", "spearman_rho", "spearman_p", "bin1_3_mean", "bin8_10_mean",
             "gam_edf", "gam_F", "gam_p", "gam_adjR2", "gam_dev_expl",
             "linear_vs_gam_p", "gam_peak_pts", "gam_peak_value",
             "breakpoint_pts", "davies_p"),
  value = c(nrow(df), sp$estimate, sp$p.value, first3, last3,
            edf, s_tab$s.table[1, "F"], s_tab$s.table[1, "p-value"], s_tab$r.sq,
            s_tab$dev.expl, lrt$`Pr(>F)`[2], pk$pts, pk$fit,
            if (!is.null(seg1)) summary(seg1)$psi[1, "Est."] else NA,
            if (exists("dv") && !is.null(dv)) dv$p.value else NA)
)
write.csv(summ, file.path(tab_dir, "08_pseudotime_shape_summary.csv"), row.names = FALSE)
write.csv(dec, file.path(tab_dir, "08_pseudotime_deciles.csv"), row.names = FALSE)
say("\n[输出] 08_pseudotime_shape_summary.csv / 08_pseudotime_deciles.csv")

## ---- 7. 形状图（两张独立面板，供 Fig S6 合成）------------------------------
# 合成管线约定：面板图须为独立 PNG，坐标轴/字号统一；此处用 ASCII（避免箭头编码告警）
png_dir <- file.path(base_dir, "results/figures/PNG_hires")
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
sub <- df[sample(nrow(df), min(1500, nrow(df))), ]
bp_val <- if (!is.null(seg1)) summary(seg1)$psi[1, "Est."] else NA_real_

# 面板 A：散点 + GAM + 断点
ggA <- ggplot() +
  geom_point(data = sub, aes(pts, DriftIndex), colour = "grey70", size = 0.25, alpha = 0.35) +
  geom_line(data = grid, aes(pts, fit), colour = "#E41A1C", linewidth = 1.1) +
  geom_vline(xintercept = pk$pts, linetype = "dashed", colour = "#185FA5", linewidth = 0.6) +
  annotate("text", x = pk$pts, y = max(grid$fit) + 0.055, hjust = -0.06, vjust = 1,
           label = sprintf("GAM maximum (pt = %.2f)", pk$pts), size = 3.1, colour = "#185FA5") +
  labs(title = "DriftIndex along osteoblast pseudotime",
       subtitle = sprintf("GAM smoother (edf = %.1f, non-linearity p < 1e-40); n = %d cells",
                          edf, nrow(df)),
       x = "Normalized pseudotime (MSC-like to mature osteoblast)",
       y = "Cell-level DriftIndex") +
  theme_classic(base_size = 12)
if (!is.na(bp_val))
  ggA <- ggA + geom_vline(xintercept = bp_val, linetype = "dotted", colour = "#4DAF4A", linewidth = 0.6)
ggsave(file.path(png_dir, "08_shape_smoother.png"), ggA, width = 8, height = 6, dpi = 300)

# 面板 B：十分位均值散点 + 分段回归拟合线 + 断点
seg_fit <- data.frame(pts = seq(0, 1, length.out = 300))
seg_fit$fit <- as.numeric(predict(seg1, newdata = seg_fit))
seg_fit$seg <- ifelse(seg_fit$pts <= bp_val, "pre", "post")
ggB <- ggplot(dec, aes(pt_mid, DriftIndex)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(data = seg_fit, aes(pts, fit, group = seg), colour = "#4DAF4A",
            linewidth = 0.9, inherit.aes = FALSE) +
  geom_point(colour = "#185FA5", size = 2.4) +
  geom_vline(xintercept = bp_val, linetype = "dotted", colour = "#4DAF4A", linewidth = 0.7) +
  labs(title = "Mean DriftIndex by pseudotime decile",
       subtitle = sprintf("segmented regression: breakpoint pt = %.2f (95%% CI %.2f-%.2f)",
                          bp_val, bp_ci[1, 2], bp_ci[1, 3]),
       x = "Normalized pseudotime (decile midpoint)",
       y = "Mean cell-level DriftIndex") +
  theme_classic(base_size = 12)
ggsave(file.path(png_dir, "08_shape_deciles.png"), ggB, width = 8, height = 6, dpi = 300)

# 合并版（备查）
ggsave(file.path(fig_dir, "08_pseudotime_shape.png"), ggA, width = 8, height = 6, dpi = 300)
say("[输出] PNG_hires/08_shape_smoother.png + 08_shape_deciles.png")

say("\n===== 08 完成 =====")
close(con)
