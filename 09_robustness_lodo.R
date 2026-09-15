#!/usr/bin/env Rscript
# ============================================================================
# 09_robustness_lodo.R — 处理"细胞级 p 值伪重复"：留一供体 + 供体级置换
#
# 背景：V1.9 §3.3 报告的 DriftIndex–pseudotime 关联（Spearman rho = -0.321,
#       p = 3.3e-93）与双相形状（GAM / segmented, p < 1e-30）均为**细胞级**检验，
#       而 3,850 个细胞仅来自 **4 名供体** —— 细胞彼此不独立，细胞级 p 值存在
#       "伪重复（pseudoreplication）"放大。V1.9 Methods 2.6 原本只是"自陈"该局限。
#
# 本脚本把该局限"处理"掉（而非仅声明）：
#   A. 留一供体（leave-one-donor-out, LODO）：逐一剔除 S1..S4，重算
#      Spearman rho、GAM(edf / 非线性 p / 偏差解释) 与 segmented(断点 / 前后斜率 / Davies p)。
#   B. 供体聚类 bootstrap（B = 2000）：以**供体**为重抽样单位 → rho 的 95% CI。
#   C. 供体块置换检验（B = 2000）：在每个供体**内部**打乱 pseudotime，
#      保留各供体的 DriftIndex 分布与 pseudotime 取值集合 → 有效零分布 → p_perm。
#   D. 供体内中心化相关：DriftIndex 与 pseudotime 各自**按供体均值中心化**后再算
#      Spearman → rho_within；检验去除供体层面混杂后关联是否仍成立。
#   E. 供体层面描述：4 个供体均值 DriftIndex ~ 均值 pseudotime（n = 4）。
#   F. S2 > S1 的细胞级置换（**反保守**，仅用于界定下界）+ 有效样本量说明。
#
# 输入：data/processed/seurat_gse169396_traj.rds（含 DriftIndex, pseudotime, sample）
# 输出：
#   results/tables/09_lodo_summary.csv
#   results/tables/09_robustness_summary.csv
#   results/tables/09_robustness_report.txt
#   results/figures/PNG_hires/09_lodo_rho_forest.png
# 用法：Rscript scripts/09_robustness_lodo.R
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
png_dir  <- file.path(fig_dir, "PNG_hires")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(tab_dir, "09_robustness_report.txt")
con <- file(LOG, "w")
say <- function(...) { writeLines(paste0(...), con); flush(con); cat(paste0(...), "\n") }

say("===== 09 稳健性：留一供体 + 供体级置换（处理细胞级 p 值伪重复）=====")
say("日期：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

set.seed(20260913)
B_PERM <- 2000   # 置换 / bootstrap 次数

# ---- 0. 载入与准备 ---------------------------------------------------------
seu <- readRDS(file.path(proc_dir, "seurat_gse169396_traj.rds"))
df <- data.frame(
  sample     = as.character(seu$sample),
  DriftIndex = as.numeric(seu$DriftIndex),
  pt_raw     = as.numeric(seu$pseudotime)
)
df <- df[is.finite(df$DriftIndex) & is.finite(df$pt_raw), ]
df$pt <- df$pt_raw / max(df$pt_raw)          # 归一化到 [0,1]
donors <- sort(unique(df$sample))
say(sprintf("\n细胞数 n = %d；供体数 = %d (%s)", nrow(df), length(donors), paste(donors, collapse = ", ")))
tab_n <- table(df$sample)
for (s in donors) say(sprintf("  %s: n = %d", s, tab_n[[s]]))

rho_of <- function(d) suppressWarnings(cor(d$DriftIndex, d$pt, method = "spearman"))

# ---- 1. 全样本基准 + 形状统计函数 ------------------------------------------
shape_stats <- function(d, tag = "full") {
  out <- list(tag = tag, n_cells = nrow(d), n_donors = length(unique(d$sample)),
              rho = rho_of(d))
  # GAM 非线性
  g <- tryCatch({
    m1 <- gam(DriftIndex ~ s(pt, k = 12, bs = "tp"), data = d, method = "REML")
    m0 <- gam(DriftIndex ~ pt, data = d, method = "REML")
    s  <- summary(m1)
    lr <- anova(m0, m1, test = "F")
    list(edf = unname(s$s.table[1, "edf"]), gam_p = unname(s$s.table[1, "p-value"]),
         dev = unname(s$dev.expl), lr_p = unname(lr$`Pr(>F)`[2]))
  }, error = function(e) list(edf = NA, gam_p = NA, dev = NA, lr_p = NA))
  # segmented 断点
  sg <- tryCatch({
    lm0 <- lm(DriftIndex ~ pt, data = d)
    grid <- data.frame(pt = seq(0, 1, length.out = 400))
    m1 <- gam(DriftIndex ~ s(pt, k = 12, bs = "tp"), data = d, method = "REML")
    psi0 <- grid$pt[which.max(predict(m1, newdata = grid))]
    sg1 <- segmented(lm0, seg.Z = ~pt, psi = psi0,
                     control = seg.control(n.boot = 50, it.max = 50, tol = 1e-4))
    bp  <- summary(sg1)$psi[1, "Est."]
    ci  <- tryCatch(confint(sg1, parm = "pt"), error = function(e) NULL)
    sl  <- slope(sg1)[[1]]
    dv  <- tryCatch(davies.test(lm0, seg.Z = ~pt, k = 10), error = function(e) NULL)
    list(bp = bp,
         bp_lo = if (!is.null(ci)) ci[1, 2] else NA,
         bp_hi = if (!is.null(ci)) ci[1, 3] else NA,
         slope_pre = sl[1, 1], slope_post = sl[2, 1],
         davies_p = if (!is.null(dv)) dv$p.value else NA)
  }, error = function(e) list(bp = NA, bp_lo = NA, bp_hi = NA,
                              slope_pre = NA, slope_post = NA, davies_p = NA))
  c(out, g, sg)
}

full <- shape_stats(df, "full (S1-S4)")
say("\n[全样本基准]")
say(sprintf("  n = %d；rho = %.4f", full$n_cells, full$rho))
say(sprintf("  GAM: edf = %.2f, 非线性 p = %.3e, 偏差解释 = %.2f%%",
            full$edf, full$gam_p, full$dev * 100))
say(sprintf("  segmented: 断点 = %.3f [%.3f, %.3f]; 前斜率 %+.4f; 后斜率 %+.4f; Davies p = %.3e",
            full$bp, full$bp_lo, full$bp_hi, full$slope_pre, full$slope_post, full$davies_p))

# ---- 2. A. 留一供体（LODO）-------------------------------------------------
say("\n[A. 留一供体 LODO] 逐一剔除，看关联与形状是否由单一供体驱动：")
lodo <- lapply(donors, function(s) shape_stats(df[df$sample != s, ], paste0("drop ", s)))
lodo_df <- do.call(rbind, lapply(lodo, function(x) data.frame(
  fold = x$tag, n_cells = x$n_cells, n_donors = x$n_donors, rho = x$rho,
  gam_edf = x$edf, gam_nonlin_p = x$gam_p, dev_expl = x$dev,
  breakpoint = x$bp, bp_lo = x$bp_lo, bp_hi = x$bp_hi,
  slope_pre = x$slope_pre, slope_post = x$slope_post, davies_p = x$davies_p)))
for (i in seq_len(nrow(lodo_df)))
  say(sprintf("  %-10s  n=%4d  rho=%+.4f  edf=%.1f  断点=%.2f  前%+.3f / 后%+.3f",
              lodo_df$fold[i], lodo_df$n_cells[i], lodo_df$rho[i], lodo_df$gam_edf[i],
              lodo_df$breakpoint[i], lodo_df$slope_pre[i], lodo_df$slope_post[i]))
say(sprintf("  → LODO rho 范围 [%.3f, %.3f]，符号全部为负：%s",
            min(lodo_df$rho), max(lodo_df$rho),
            ifelse(all(lodo_df$rho < 0), "是", "否")))
say(sprintf("  → LODO 断点范围 [%.2f, %.2f]；前后斜率符号全部一致：%s",
            min(lodo_df$breakpoint, na.rm = TRUE), max(lodo_df$breakpoint, na.rm = TRUE),
            ifelse(all(lodo_df$slope_pre > 0, na.rm = TRUE) &
                     all(lodo_df$slope_post < 0, na.rm = TRUE), "是", "否")))

# ---- 3. B. 供体聚类 bootstrap ----------------------------------------------
say(sprintf("\n[B. 供体聚类 bootstrap（B = %d，重抽样单位为供体）]", B_PERM))
boot <- numeric(B_PERM)
for (b in seq_len(B_PERM)) {
  take <- sample(donors, length(donors), replace = TRUE)
  db <- do.call(rbind, lapply(seq_along(take), function(i) {
    x <- df[df$sample == take[i], ]; x$sample <- paste0(x$sample, "_", i); x
  }))
  boot[b] <- rho_of(db)
}
boot_ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
say(sprintf("  rho 中位 = %.4f；95%% CI = [%.3f, %.3f]；含 0：%s",
            median(boot, na.rm = TRUE), boot_ci[1], boot_ci[2],
            ifelse(boot_ci[1] < 0 && boot_ci[2] > 0, "是（不确定）", "否")))

# ---- 4. C. 供体块置换检验 ---------------------------------------------------
say(sprintf("\n[C. 供体块置换检验（B = %d，供体内打乱 pseudotime）]", B_PERM))
blocks <- lapply(donors, function(s) which(df$sample == s))
pt_perm_rho <- numeric(B_PERM)
for (b in seq_len(B_PERM)) {
  ptp <- df$pt
  for (idx in blocks) ptp[idx] <- df$pt[sample(idx)]
  pt_perm_rho[b] <- suppressWarnings(cor(df$DriftIndex, ptp, method = "spearman"))
}
p_perm <- (1 + sum(abs(pt_perm_rho) >= abs(full$rho))) / (B_PERM + 1)
say(sprintf("  零分布 rho 中位 = %.4f（即供体构成本身可解释的关联量级）；观测 rho = %.4f",
            median(pt_perm_rho), full$rho))
say(sprintf("  置换 p = %.4f（双侧）；有效零分布最小可达 p = %.4f",
            p_perm, 1 / (B_PERM + 1)))

# ---- 5. D. 供体内中心化相关 ------------------------------------------------
df_c <- df %>% group_by(sample) %>%
  mutate(DriftIndex_c = DriftIndex - mean(DriftIndex),
         pt_c = pt - mean(pt)) %>% ungroup()
rho_within <- suppressWarnings(cor(df_c$DriftIndex_c, df_c$pt_c, method = "spearman"))
# 其置换零分布
wc_perm <- numeric(B_PERM)
for (b in seq_len(B_PERM)) {
  ptp <- df_c$pt_c
  for (idx in blocks) ptp[idx] <- df_c$pt_c[sample(idx)]
  wc_perm[b] <- suppressWarnings(cor(df_c$DriftIndex_c, ptp, method = "spearman"))
}
p_within <- (1 + sum(abs(wc_perm) >= abs(rho_within))) / (B_PERM + 1)
say("\n[D. 供体内中心化（去除供体层面混杂）]")
say(sprintf("  供体内中心化 rho_within = %.4f；置换 p = %.4f", rho_within, p_within))

# ---- 6. E. 供体层面描述（n = 4）--------------------------------------------
dm <- df %>% group_by(sample) %>%
  summarise(mean_DriftIndex = mean(DriftIndex), mean_pt = mean(pt),
            n = n(), .groups = "drop")
rho_donor <- suppressWarnings(cor(dm$mean_DriftIndex, dm$mean_pt, method = "spearman"))
say("\n[E. 供体层面（n = 4）]")
for (i in seq_len(nrow(dm)))
  say(sprintf("  %s: 均值 DriftIndex = %+.4f, 均值 pseudotime = %.3f (n=%d)",
              dm$sample[i], dm$mean_DriftIndex[i], dm$mean_pt[i], dm$n[i]))
say(sprintf("  4 供体均值间 Spearman rho = %+.3f（n = 4，仅描述，p 无意义）", rho_donor))

# ---- 7. F. S2 > S1 的细胞级置换（反保守下界）-------------------------------
s12 <- df[df$sample %in% c("S1", "S2"), ]
obs_diff <- mean(s12$DriftIndex[s12$sample == "S2"]) - mean(s12$DriftIndex[s12$sample == "S1"])
lab <- s12$sample
perm_diff <- numeric(B_PERM)
for (b in seq_len(B_PERM)) {
  L <- sample(lab)
  perm_diff[b] <- mean(s12$DriftIndex[L == "S2"]) - mean(s12$DriftIndex[L == "S1"])
}
p_s12 <- (1 + sum(abs(perm_diff) >= abs(obs_diff))) / (B_PERM + 1)
say("\n[F. S2 vs S1（有效样本量 = 每层 1 名供体）]")
say(sprintf("  观测差值（S2 - S1）= %+.4f；细胞级标签置换 p = %.4f（**反保守**）", obs_diff, p_s12))
say("  注：同一供体内细胞不独立，该 p 值仅给出下界；按供体计有效 n = 1/层，无法做推断，")
say("      故 S2 > S1 仍为描述性结论（与 V1.9 正文一致）。")

# ---- 8. 输出表 -------------------------------------------------------------
lodo_out <- rbind(
  data.frame(fold = full$tag, n_cells = full$n_cells, n_donors = full$n_donors,
             rho = full$rho, gam_edf = full$edf, gam_nonlin_p = full$gam_p,
             dev_expl = full$dev, breakpoint = full$bp, bp_lo = full$bp_lo,
             bp_hi = full$bp_hi, slope_pre = full$slope_pre,
             slope_post = full$slope_post, davies_p = full$davies_p),
  lodo_df)
write.csv(lodo_out, file.path(tab_dir, "09_lodo_summary.csv"), row.names = FALSE)

summ <- data.frame(
  metric = c("n_cells", "n_donors", "rho_cell_level", "rho_lodo_min", "rho_lodo_max",
             "rho_lodo_all_negative", "rho_within_donor", "p_within_donor_perm",
             "boot_ci_lo", "boot_ci_hi", "boot_ci_excludes_zero",
             "p_perm_donor_block", "gam_edf_full", "gam_p_full",
             "breakpoint_full", "bp_lo_full", "bp_hi_full",
             "slope_pre_full", "slope_post_full", "davies_p_full",
             "breakpoint_lodo_min", "breakpoint_lodo_max",
             "donor_level_rho_n4", "s2_vs_s1_cell_perm_p_anticonservative"),
  value = c(nrow(df), length(donors), full$rho, min(lodo_df$rho), max(lodo_df$rho),
            all(lodo_df$rho < 0), rho_within, p_within,
            boot_ci[1], boot_ci[2], !(boot_ci[1] < 0 && boot_ci[2] > 0),
            p_perm, full$edf, full$gam_p,
            full$bp, full$bp_lo, full$bp_hi,
            full$slope_pre, full$slope_post, full$davies_p,
            min(lodo_df$breakpoint, na.rm = TRUE), max(lodo_df$breakpoint, na.rm = TRUE),
            rho_donor, p_s12)
)
write.csv(summ, file.path(tab_dir, "09_robustness_summary.csv"), row.names = FALSE)
say("\n[输出] 09_lodo_summary.csv / 09_robustness_summary.csv / 09_robustness_report.txt")

# ---- 9. 森林图（LODO rho）--------------------------------------------------
fp <- rbind(
  data.frame(label = "All donors (n=3850)", rho = full$rho, lo = NA, hi = NA,
             type = "full"),
  data.frame(label = paste0("Drop ", donors), rho = lodo_df$rho,
             lo = NA, hi = NA, type = "lodo"),
  data.frame(label = "Donor cluster bootstrap", rho = median(boot, na.rm = TRUE),
             lo = boot_ci[1], hi = boot_ci[2], type = "boot"),
  data.frame(label = "Within-donor centred", rho = rho_within, lo = NA, hi = NA,
             type = "within")
)
fp$label <- factor(fp$label, levels = rev(fp$label))
fp$type  <- factor(fp$type, levels = c("full", "lodo", "boot", "within"))
ggF <- ggplot(fp, aes(rho, label, colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(full = "#E41A1C", lodo = "#185FA5",
                                 boot = "#4DAF4A", within = "#984EA3"),
                      labels = c("Full sample", "Leave-one-donor-out (LODO)",
                                 "Donor-cluster bootstrap 95% CI",
                                 "Within-donor centred"),
                      name = NULL) +
  labs(title = "Donor-aware robustness of the DriftIndex-pseudotime association",
       subtitle = sprintf("LODO rho [%.3f, %.3f]; donor-block permutation p < %.4f (B = %d)",
                          min(lodo_df$rho), max(lodo_df$rho), 1 / (B_PERM + 1), B_PERM),
       x = "Spearman's rho (DriftIndex vs pseudotime)", y = NULL) +
  coord_cartesian(xlim = c(-0.40, 0.015)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        plot.subtitle = element_text(size = 9.5)) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(png_dir, "09_lodo_rho_forest.png"), ggF, width = 8, height = 6, dpi = 300)
say("[输出] PNG_hires/09_lodo_rho_forest.png")

say("\n===== 09 完成 =====")
close(con)
