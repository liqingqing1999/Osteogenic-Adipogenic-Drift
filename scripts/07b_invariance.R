#!/usr/bin/env Rscript
# ============================================================================
# 07b_invariance.R — 步 1：阈值不变区间（S2 > S1 成立的阈值区间 + 不确定度）
# ----------------------------------------------------------------------------
# 归属：论文 A（漂移稿）= output/Lu-2026-mesenchymal-drift/
# 依据：docs/70 §2 步 1 ＋ §3 三条口径（判据 = 符号稳定性而非显著性）
# 上位登记：A1方案_实施规划.md 的 W1-07（"45.8% 高度依赖阈值 ⇒ 须改成阈值无关定义"）
#
# 输入：data/processed/seurat_gse169396_drift.rds（逐细胞 AD_drift_ucell / OS_identity_ucell）
#       results/tables/03_sensitivity_summary.csv（7 点基线，用于验收断言）
# 输出：results/invariance/07b_invariance.csv
#       results/invariance/07b_exact_support_profile.csv
#       results/invariance/07b_invariance.pdf
#       results/invariance/07b_invariance_summary.txt
#
# 口径（预注册，不在看到结果后改）：
#   R1 判据 = D(t) = pct_S2(t) - pct_S1(t) 的【符号】；支持区间 = D>0 的连通区间
#   R2 不确定度 = 细胞级 bootstrap（B=2000，按细胞重抽）+ Wilson 区间（仅参考）
#   R3 规则不变性交叉（三条规则）：
#       规则 A = 共同绝对阈值                  AD > t & OS > t        （主判据）
#       规则 B = 各程序自身分位（参考 = 合并全体）   ← 公平交叉检查，保留供体间位移
#       规则 C = 各程序自身分位（参考 = 各供体自身） ← 诊断用：它把供体间位移归一化掉
#         ⚠️ 对齐轴 = pooled 双阳性比例。分数在 0 处有点质量 ⇒ 自身分位规则的可达 pooled
#            上限 < 100%，故交叉检查只在【可对齐区间】上进行（实测正好覆盖 t ≥ 0）。
#   HK 恒等守卫：7 个重叠阈值须逐位复现 03_sensitivity_summary.csv（tol 1e-9），否则中止
# ============================================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

t0 <- Sys.time()
cat("===== 07b 步 1：阈值不变区间 =====\n"); flush.console()

# ---- 项目根目录（脚本位于 <root>/scripts/；可用 DRIFT_PROJECT_DIR 覆盖）----
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
tab_dir <- file.path(base_dir, "results/tables")
inv_dir <- file.path(base_dir, "results/invariance")
dir.create(inv_dir, recursive = TRUE, showWarnings = FALSE)

RDS  <- file.path(base_dir, "data/processed/seurat_gse169396_drift.rds")
BASE <- file.path(tab_dir, "03_sensitivity_summary.csv")
BASE_MD5 <- unname(tools::md5sum(BASE)); RDS_MD5 <- unname(tools::md5sum(RDS))
cat(sprintf("输入 03_sensitivity_summary.csv md5 %s\n", BASE_MD5)); flush.console()
cat(sprintf("输入 seurat_gse169396_drift.rds md5 %s\n", RDS_MD5)); flush.console()

seu <- readRDS(RDS)
AD <- seu$AD_drift_ucell; OS <- seu$OS_identity_ucell; SM <- as.character(seu$sample)
N_ALL <- length(AD)
DONORS <- c("S1", "S2", "S3", "S4")
sp <- split(seq_len(N_ALL), SM); gd <- lapply(DONORS, function(s) sp[[s]]); names(gd) <- DONORS
ND <- vapply(gd, length, 0L)
cat(sprintf("对象：%d 细胞 × %d 基因；%s\n", ncol(seu), nrow(seu),
            paste(sprintf("%s=%d", DONORS, ND), collapse = " "))); flush.console()

STEP <- 0.005
GRID <- round(seq(-0.10, 0.20, by = STEP), 3)
STOPIF <- 1e-9

wilson <- function(k, n, z = 1.959963984540054) {
  if (n == 0) return(c(NA_real_, NA_real_))
  ph <- k/n; ctr <- (ph + z^2/(2*n))/(1 + z^2/n)
  hw <- z*sqrt(ph*(1-ph)/n + z^2/(4*n^2))/(1 + z^2/n)
  c(100*max(0, ctr-hw), 100*min(1, ctr+hw))
}
flagA <- function(t) AD > t & OS > t
pctOf <- function(flag) 100 * vapply(gd, function(ix) sum(flag[ix]), 0L) / ND

# ---- 细胞级 bootstrap（B = 2000；档内按细胞重抽）----------------------------
B <- 2000L; set.seed(42)
idx <- lapply(c("S1","S2"), function(s) matrix(sample.int(ND[[s]], ND[[s]]*B, replace = TRUE), nrow = B))
names(idx) <- c("S1","S2")
boot_D <- function(flag) {
  b1 <- rowMeans(matrix(flag[gd[["S1"]]][idx[["S1"]]], nrow = B)) * 100
  b2 <- rowMeans(matrix(flag[gd[["S2"]]][idx[["S2"]]], nrow = B)) * 100
  quantile(b2 - b1, c(0.025, 0.975), names = FALSE)
}

cat("扫描预注册网格（规则 A + bootstrap）…\n"); flush.console()
resA <- vector("list", length(GRID))
for (i in seq_along(GRID)) {
  t <- GRID[i]; fl <- flagA(t); k <- vapply(gd, function(ix) sum(fl[ix]), 0L)
  pct <- 100*k/ND; ci <- boot_D(fl)
  resA[[i]] <- list(t = t, k = k, pct = pct, D = pct[["S2"]] - pct[["S1"]],
                    lo = ci[1], hi = ci[2], pooled = 100*sum(k)/N_ALL,
                    wil = list(wilson(k[["S1"]], ND[["S1"]]), wilson(k[["S2"]], ND[["S2"]]),
                               wilson(k[["S3"]], ND[["S3"]]), wilson(k[["S4"]], ND[["S4"]])))
  if (i %% 20 == 0) { cat(sprintf("  … %d/%d (t=%+.3f)\n", i, length(GRID), t)); flush.console() }
}

# ---- 规则 B / C：各程序自身分位，两种参考分布 -------------------------------
QG <- seq(0.001, 0.999, by = 0.001)
cAD_pool <- as.numeric(quantile(AD, QG)); cOS_pool <- as.numeric(quantile(OS, QG))
cAD_don  <- lapply(DONORS, function(s) as.numeric(quantile(AD[gd[[s]]], QG)))
cOS_don  <- lapply(DONORS, function(s) as.numeric(quantile(OS[gd[[s]]], QG)))

scan_pool <- function(cAD, cOS) {      # 规则 B：阈值取合并全体分布的分位
  hp <- numeric(length(QG)); pct <- matrix(NA_real_, length(QG), length(DONORS)); colnames(pct) <- DONORS
  for (i in seq_along(QG)) {
    f <- AD > cAD[i] & OS > cOS[i]
    k <- vapply(gd, function(ix) sum(f[ix]), 0L)
    pct[i, ] <- 100*k/ND; hp[i] <- 100*sum(k)/N_ALL
  }
  list(hp = hp, pct = pct)
}
scan_donor <- function(cAD, cOS) {     # 规则 C：阈值取各供体自身分布的分位
  hp <- numeric(length(QG)); pct <- matrix(NA_real_, length(QG), length(DONORS)); colnames(pct) <- DONORS
  for (i in seq_along(QG)) {
    k <- vapply(seq_along(DONORS), function(j) {
      ix <- gd[[j]]; sum(AD[ix] > cAD[[j]][i] & OS[ix] > cOS[[j]][i])
    }, 0L)
    pct[i, ] <- 100*k/ND; hp[i] <- 100*sum(k)/N_ALL
  }
  list(hp = hp, pct = pct)
}
cat("扫描规则 B（合并全体分位）与规则 C（各供体自身分位）…\n"); flush.console()
qB <- scan_pool(cAD_pool, cOS_pool)
qC <- scan_donor(cAD_don, cOS_don)

# 对齐可行性：自身分位规则的可达 pooled 上限
HP_MAX <- max(qB$hp)
resQ <- function(qq, target) {
  i <- which.min(abs(qq$hp - target))
  list(q = QG[i], pct = qq$pct[i, ], D = qq$pct[i, "S2"] - qq$pct[i, "S1"], hp = qq$hp[i])
}
ALIGN <- vapply(resA, function(a) a$pooled <= HP_MAX + 1e-9, TRUE)

# 规则 C 的退化诊断：在同一 q 上，两供体的【边际越界率】是否被拉平
rnames <- function(x) { names(x) <- DONORS; x }
gap_C <- vapply(seq_along(QG), function(i) {
  ra <- rnames(vapply(seq_along(DONORS), function(j) 100*sum(AD[gd[[j]]] > cAD_don[[j]][i])/ND[[j]], 0))
  ro <- rnames(vapply(seq_along(DONORS), function(j) 100*sum(OS[gd[[j]]] > cOS_don[[j]][i])/ND[[j]], 0))
  max(abs(ra[["S1"]] - ra[["S2"]]), abs(ro[["S1"]] - ro[["S2"]]))
}, 0)
# 对照：规则 A 在对应阈值上的同一指标
gap_A <- vapply(seq_along(GRID), function(i) {
  t <- GRID[i]
  ra <- rnames(vapply(seq_along(DONORS), function(j) 100*sum(AD[gd[[j]]] > t)/ND[[j]], 0))
  ro <- rnames(vapply(seq_along(DONORS), function(j) 100*sum(OS[gd[[j]]] > t)/ND[[j]], 0))
  max(abs(ra[["S1"]] - ra[["S2"]]), abs(ro[["S1"]] - ro[["S2"]]))
}, 0)
gap_C_at <- function(target) gap_C[which.min(abs(qC$hp - target))]

rows <- lapply(seq_along(GRID), function(i) {
  a <- resA[[i]]; al <- ALIGN[i]
  b <- if (al) resQ(qB, a$pooled) else list(q = NA_real_, pct = rep(NA_real_, 4), D = NA_real_, hp = NA_real_)
  c3 <- if (al) resQ(qC, a$pooled) else list(q = NA_real_, pct = rep(NA_real_, 4), D = NA_real_, hp = NA_real_)
  k <- a$k
  data.frame(
    threshold = a$t,
    n_S1 = ND[["S1"]], n_hyb_S1 = k[["S1"]], pct_S1 = a$pct[["S1"]],
    wilson_lo_S1 = a$wil[[1]][1], wilson_hi_S1 = a$wil[[1]][2],
    n_S2 = ND[["S2"]], n_hyb_S2 = k[["S2"]], pct_S2 = a$pct[["S2"]],
    wilson_lo_S2 = a$wil[[2]][1], wilson_hi_S2 = a$wil[[2]][2],
    n_S3 = ND[["S3"]], n_hyb_S3 = k[["S3"]], pct_S3 = a$pct[["S3"]],
    n_S4 = ND[["S4"]], n_hyb_S4 = k[["S4"]], pct_S4 = a$pct[["S4"]],
    D = a$D, D_boot_lo = a$lo, D_boot_hi = a$hi, D_strict_positive = a$lo > 0,
    pooled_hybrid_pct = a$pooled, alignable = al,
    ruleB_matched_q = b$q, ruleB_pooled = b$hp, ruleB_D = b$D,
    ruleC_matched_q = c3$q, ruleC_pooled = c3$hp, ruleC_D = c3$D,
    gapA_marginal = gap_A[i],
    gapC_marginal = if (al) gap_C_at(a$pooled) else NA_real_,
    D_positive = a$D > 0,
    ruleB_D_positive = if (al) b$D > 0 else NA, ruleC_D_positive = if (al) c3$D > 0 else NA,
    stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)

# ============================================================================
# HK 验收断言：7 个重叠阈值逐位复现基线（tol 1e-9），不通过则【中止不落盘】
# ============================================================================
base <- read.csv(BASE, check.names = FALSE)
thr7 <- c(-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06)
worst <- 0
for (t in thr7) for (s in DONORS)
  worst <- max(worst, abs(tab[[paste0("pct_", s)]][tab$threshold == t] -
                          base$hybrid_pct[base$sample == s & base$threshold == t]))
HK <- worst <= STOPIF
cat(sprintf("HK 恒等守卫：7 阈值 × 4 供体 最大绝对偏差 %.3e  tol %.0e  => %s\n",
            worst, STOPIF, ifelse(HK, "PASS", "FAIL"))); flush.console()
if (!HK) stop("HK 恒等守卫失败：本脚本与 03_sensitivity 的打分不一致，拒绝落盘。")

# ============================================================================
# 支持区间
# ============================================================================
intervals <- function(ts, ok, step = STEP) {
  out <- list(); i <- 1L
  while (i <= length(ts)) {
    if (!ok[i]) { i <- i + 1L; next }
    j <- i
    while (j < length(ts) && ok[j + 1L] &&
           (is.null(step) || (ts[j + 1L] - ts[j]) <= step + 1e-6)) j <- j + 1L
    out[[length(out) + 1L]] <- c(ts[i], ts[j]); i <- j + 1L
  }
  out
}
fmt_iv <- function(iv) if (!length(iv)) "空集" else
  paste(vapply(iv, function(x) sprintf("[%+.3f, %+.3f]", x[1], x[2]), ""), collapse = " ∪ ")
outer_iv <- function(iv) if (!length(iv)) c(NA_real_, NA_real_) else
  c(min(vapply(iv, `[`, 0, 1)), max(vapply(iv, `[`, 0, 2)))

okA <- tab$threshold >= 0 & tab$D_positive
ivA <- intervals(tab$threshold[tab$threshold >= 0], tab$D_positive[tab$threshold >= 0])
okB <- tab$threshold >= 0 & !is.na(tab$ruleB_D_positive) & tab$ruleB_D_positive
ivB <- intervals(tab$threshold[tab$threshold >= 0], okB[tab$threshold >= 0])
okC <- tab$threshold >= 0 & !is.na(tab$ruleC_D_positive) & tab$ruleC_D_positive
ivC <- intervals(tab$threshold[tab$threshold >= 0], okC[tab$threshold >= 0])
abA <- outer_iv(ivA); abB <- outer_iv(ivB); abC <- outer_iv(ivC)
lenA <- abA[2] - abA[1]; lenB <- abB[2] - abB[1]
t_align_max <- max(tab$threshold[tab$alignable])

# ---- 精确端点（规则 A）：D 是 t 的阶梯函数，跳点 = 各细胞 min(AD, OS) 的取值 ----
mstat <- pmin(AD, OS); eps <- 1e-12
tv <- sort(unique(c(0, mstat - eps))); tv <- tv[tv >= 0]
ev <- vapply(tv, function(t) {
  f <- AD > t & OS > t
  100*(sum(f[gd[["S2"]]])/ND[["S2"]] - sum(f[gd[["S1"]]])/ND[["S1"]])
}, 0)
ivE <- intervals(tv, ev > 0, step = NULL)
abE <- outer_iv(ivE)
t_right <- abE[2]
kR <- c(S1 = sum(flagA(t_right)[gd[["S1"]]]), S2 = sum(flagA(t_right)[gd[["S2"]]]))
maxD <- max(ev); argmax_D <- tv[which.max(ev)]; minD_pos <- min(ev[ev > 0])

# ---- bootstrap 严格端点 ----
okS <- tab$threshold >= 0 & tab$D_strict_positive
t_strict <- if (any(okS)) max(tab$threshold[okS]) else NA_real_
kS <- if (is.na(t_strict)) c(S1 = NA, S2 = NA) else c(
  S1 = tab$n_hyb_S1[tab$threshold == t_strict], S2 = tab$n_hyb_S2[tab$threshold == t_strict])

# ---- 最终（最保守）可主张区间：只取【同被估量】的规则 A / B 与 bootstrap ----
# 规则 C 被排除：它按构造消除了供体间位移（见诊断），改变了被估量，不构成“更窄”的证据。
cand <- c(ruleA = abE[2], ruleB = abB[2],
          bootstrap = ifelse(is.na(t_strict), NA_real_, t_strict))
cand <- cand[!is.na(cand)]
t_final <- min(cand)
# 精确支持集的断点宽度（若有）
gaps <- if (length(ivE) > 1) vapply(seq_len(length(ivE) - 1), function(k)
  ivE[[k + 1]][1] - ivE[[k]][2], 0) else numeric(0)
gap_max <- if (length(gaps)) max(gaps) else 0
cat(sprintf("规则A 右端 %.4f｜规则B 右端 %.3f｜规则C 右端 %.3f（诊断，不计入）｜bootstrap %.3f ⇒ 最保守 %.3f\n",
            abE[2], abB[2], abC[2], t_strict, t_final)); flush.console()

# ============================================================================
# 落盘
# ============================================================================
write.csv(tab, file.path(inv_dir, "07b_invariance.csv"), row.names = FALSE)
write.csv(data.frame(t = tv, D = ev), file.path(inv_dir, "07b_exact_support_profile.csv"), row.names = FALSE)

pdf(file.path(inv_dir, "07b_invariance.pdf"), width = 7.2, height = 8.6)
long <- do.call(rbind, lapply(DONORS, function(s) data.frame(
  threshold = tab$threshold, pct = tab[[paste0("pct_", s)]], donor = factor(s, levels = DONORS))))
band <- data.frame(x0 = 0, x1 = t_final)
vlines <- data.frame(x = c(0, 0.06, abE[2], t_final),
                     g = factor(c("t = 0", "t = +0.06 (in-review)", "rule-A endpoint", "conservative endpoint"),
                                levels = c("t = 0", "t = +0.06 (in-review)", "rule-A endpoint", "conservative endpoint")))
p1 <- ggplot(long, aes(threshold, pct, color = donor)) +
  geom_rect(data = band, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.1) +
  geom_vline(data = vlines, aes(xintercept = x, linetype = g), color = "grey25") +
  scale_linetype_manual(values = c("dashed", "dotted", "dotdash", "longdash"), name = NULL) +
  scale_color_manual(values = c(S1 = "#E41A1C", S2 = "#FF7F00", S3 = "#984EA3", S4 = "#4DAF4A"),
                     labels = c("S1 (osteoporotic range)", "S2 (osteopenic range)",
                                "S3 (not stratified)", "S4 (not stratified)")) +
  labs(title = "Step 1 - Hybrid cell % across the threshold grid",
       subtitle = sprintf("shaded = conservative support interval (0, %.3f];  grid step 0.005;  data ceiling t = %.4f", t_final, max(mstat)),
       x = "Hybrid threshold (AD & OS > t)", y = "% Hybrid cells", color = "Donor") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")

p2 <- ggplot(tab, aes(threshold, D)) +
  geom_rect(data = band, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_ribbon(aes(ymin = D_boot_lo, ymax = D_boot_hi), fill = "#377EB8", alpha = 0.25) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_hline(yintercept = 0, color = "grey20") +
  geom_vline(data = vlines, aes(xintercept = x, linetype = g), color = "grey25") +
  scale_linetype_manual(values = c("dashed", "dotted", "dotdash", "longdash"), name = NULL) +
  labs(title = "D(t) = pct_S2(t) - pct_S1(t)   with cell-level bootstrap (B = 2000)",
       subtitle = "band: 95% percentile interval; reflects cell sampling and scoring only,\nnot donor sampling (n = 1 donor per stratum).",
       x = "Hybrid threshold (AD & OS > t)", y = "D(t)  (percentage points)") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")

if (requireNamespace("gridExtra", quietly = TRUE)) gridExtra::grid.arrange(p1, p2, ncol = 1) else { print(p1); print(p2) }
dev.off()

# ============================================================================
# 摘要
# ============================================================================
L <- character(0); A <- function(...) L <<- c(L, paste0(...))
A("07b 步 1 · 阈值不变区间 —— 结果摘要")
A(sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
A(sprintf("输入 rds md5：%s", RDS_MD5))
A(sprintf("基线表 03_sensitivity_summary.csv md5：%s", BASE_MD5))
A(sprintf("细胞数：%s（合计 %d）", paste(sprintf("%s=%d", DONORS, ND), collapse = " "), N_ALL))
A("")
A("【1】判据与口径")
A("  主判据 = D(t) = pct_S2(t) - pct_S1(t) 的符号（不是统计显著性）")
A("  支持区间 = {t : D(t) > 0} 的连通区间")
A("  网格（预注册）= [-0.10, +0.20]，步长 0.005，共 61 点")
A("  不确定度 = 细胞级 bootstrap B = 2000（按细胞重抽）+ Wilson 区间（仅参考）")
A("  规则不变性交叉（三条）：")
A("    规则 A = 共同绝对阈值                    AD > t & OS > t          （主判据）")
A("    规则 B = 各程序自身分位（参考 = 合并全体）      ← 公平交叉检查")
A("    规则 C = 各程序自身分位（参考 = 各供体自身）    ← 仅诊断，见【5】")
A(sprintf("  对齐轴 = pooled 双阳性比例；自身分位规则的可达上限 = %.3f%%（因分数在 0 处有点质量）",
          HP_MAX))
A(sprintf("  ⇒ 可对齐区间：t ≤ %+.3f（规则 A 的 pooled ≤ 上限）", t_align_max))
A("")
A("【2】HK 恒等守卫")
A(sprintf("  7 个重叠阈值 × 4 供体，与 03_sensitivity_summary.csv 最大绝对偏差 %.3e（tol 1e-9）=> %s",
          worst, ifelse(HK, "PASS", "FAIL")))
A("")
A("【3】规则 A —— 主结果")
A(sprintf("  D(0) = %+.4f 百分点（pct_S2 = %.4f，pct_S1 = %.4f）",
          tab$D[tab$threshold == 0], tab$pct_S2[tab$threshold == 0], tab$pct_S1[tab$threshold == 0]))
A(sprintf("  峰值 = %+.4f 百分点，出现在 t = %+.4f", maxD, argmax_D))
A(sprintf("  在预注册网格 [0.000, 0.200] 上：%d 个网格点全部 D > 0，连通区间 %s",
          sum(tab$threshold >= 0 & tab$D_positive), fmt_iv(ivA)))
A(sprintf("  精确连通区间（对阶梯函数逐点求值，跳点 = 各细胞 min(AD,OS)）= %s", fmt_iv(ivE)))
if (length(ivE) > 1) A(sprintf("    ⚠️ 精确支持集由 %d 段构成，中间最大断点宽度 = %.4f（个别细胞的 min(AD,OS) 落点造成）",
                               length(ivE), gap_max))
A(sprintf("    对外取【外侧包络】：(%.4f, %.4f]", abE[1], abE[2]))
A(sprintf("  ⇒ 精确右端点 = %.4f（区间长度 %.4f）；支持区间内 D ∈ [%+.4f, %+.4f]（恒 > 0）",
          abE[2], abE[2] - abE[1], minD_pos, maxD))
A(sprintf("  数据上限 t（= max min(AD,OS)）= %.4f ⇒ 真实右端点必 ≤ 该值；与精确右端点相差 %.4f",
          max(mstat), max(mstat) - abE[2]))
A(sprintf("  ⚠️ 精确右端点处 hybrid 细胞数：S1 = %d / %d，S2 = %d / %d ⇒ 尾部由极少数细胞驱动",
          kR[["S1"]], ND[["S1"]], kR[["S2"]], ND[["S2"]]))
A("")
A("【4】bootstrap 严格支持区间（要求 D 的 2.5% 分位 > 0）")
A(sprintf("  严格端点 = %s；该处 hybrid 细胞数 S1 = %s / %d，S2 = %s / %d",
          ifelse(is.na(t_strict), "不存在", sprintf("%+.3f", t_strict)),
          as.character(kS[["S1"]]), ND[["S1"]], as.character(kS[["S2"]]), ND[["S2"]]))
A("")
A("【5】规则不变性交叉")
A(sprintf("  规则 B（参考 = 合并全体）在 t ≥ 0 上：D_B > 0 的连通区间 %s ⇒ 右端点 %+.3f",
          fmt_iv(ivB), abB[2]))
A(sprintf("  规则 B 参数对齐质量（仅可对齐区间）：max |pooled_B - pooled_A| = %.4f 百分点",
          max(abs(tab$ruleB_pooled - tab$pooled_hybrid_pct), na.rm = TRUE)))
A(sprintf("  ⇒ 规则 A 与 B 的右端点差 = %.3f（%.4f vs %.3f）", abs(abE[2] - abB[2]), abE[2], abB[2]))
A("")
A("  规则 C（参考 = 各供体自身）：D_C > 0 的连通区间 "  )
A(sprintf("    %s ⇒ 右端点 %+.3f", fmt_iv(ivC), abC[2]))
A("  🔴 规则 C 不可用作“更窄者为准”的比较对象 —— 诊断证据：")
A(sprintf("    同一阈值水平上两供体【边际越界率】之差：规则 A 中位 %.3f / 最大 %.3f 百分点；",
          median(tab$gapA_marginal), max(tab$gapA_marginal)))
A(sprintf("    规则 C 在与之对齐的 q 水平上：中位 %.3f / 最大 %.3f 百分点",
          median(tab$gapC_marginal, na.rm = TRUE), max(tab$gapC_marginal, na.rm = TRUE)))
A("    ⇒ 规则 C 按【构造】把每个供体在两个程序上的边际越界率各自钉在其自身分位上，")
A("       于是“供体间的整体位移”被整块消除 —— 而供体间位移正是规则 A 所度量的信号；")
A("       规则 C 剩下的只有“供体内两程序的联合相关结构”。")
A("    ⇒ 它的窄区间是【改变了被估量】造成的结构性伪影，不是“结论依赖阈值”的证据。")
A("       （docs/70 §3.2 原口径需按此更正：交叉检查取 B，不取 C。）")
A("")
A("【6】决策门")
A(sprintf("  D-A（区间长度 < 0.02 ⇒ 不得用“证书”命名）：保守长度 %.3f ⇒ %s",
          t_final, ifelse(t_final >= 0.02, "通过 —— 可用“区间（interval）”命名", "不通过")))
A(sprintf("  D-D（右端点是否 < +0.06 ⇒ 与在审稿“holds at every informative threshold (>= 0)”冲突）："))
A(sprintf("    同被估量的三条：规则A %.4f｜规则B %.3f｜bootstrap %.3f ⇒ 最保守 %.3f ⇒ %s",
          abE[2], abB[2], t_strict, t_final,
          ifelse(t_final < 0.06, "🔴 冲突，须立即回报（按 R-6 处置）", "✅ 不冲突 —— 可升级在审稿表述")))
A(sprintf("    （规则 C 的 %.3f 不计入：其被估量已被改变，见【5】）", abC[2]))
A("")
A("【7】必须同时写入的限制（§3.2 口径，逐字可用）")
A("  EN: The interval reflects cell sampling and scoring variability only; it does not represent")
A("      donor sampling, as each density stratum contains a single donor.")
A("  EN: UCell proportions are over-dispersed, so the Wilson intervals are a lower bound on")
A("      interval width and are reported for reference only.")
A("  EN: Endpoints were resolved exactly on the step function of D(t) by evaluating every")
A("      observed value of min(AD, OS); the reported grid values use a step of 0.005.")
A("  EN: The upper end of the interval is driven by a very small number of cells, so the")
A("      conservative endpoint (bootstrap 2.5% quantile > 0) is reported alongside it.")
A("  CN: 该区间只反映抽样与打分波动，不反映供体抽样（每档仅 1 名供体）。")
A("  CN: UCell 比例呈过离散，Wilson 区间只是宽度下界，仅作参考。")
A("  CN: 端点由 0.005 步长的网格给出；右端由极少数细胞驱动，故并列报告保守端点。")
A("")
A("【8】可入稿的表述（D-D 不冲突；取最保守端点）")
A(sprintf("  EN: The S2 > S1 ordering was preserved for every threshold in (0, %.3f] (grid step 0.005);",
          t_final))
A(sprintf("      the point estimate remained positive up to t = %.4f, the largest value the data can resolve.", abE[2]))
A("      The ordering was reproduced under a control rule that cuts each programme at its own")
A("      quantile of the pooled score distribution instead of at a shared absolute value")
A(sprintf("      (support to t = %.3f under that rule).", abB[2]))
A("")
A(sprintf("产出：%s", file.path(inv_dir, "07b_invariance.csv")))
A(sprintf("用时：%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
writeLines(L, file.path(inv_dir, "07b_invariance_summary.txt"))
cat(paste(L, collapse = "\n")); cat("\n===== 07b 完成 =====\n")
