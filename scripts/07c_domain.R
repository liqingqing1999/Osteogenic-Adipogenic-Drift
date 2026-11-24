#!/usr/bin/env Rscript
# ============================================================================
# 07c_domain.R — 步 2：三参数成立域（t × r_doublet × c_ambient）
# ----------------------------------------------------------------------------
# 归属：论文 A（漂移稿）= output/Lu-2026-mesenchymal-drift/
# 依据：docs/70 §4（步 2 设计）＋ §8 决策门 D-C
# 上位登记：A1方案_实施规划.md 的 W1-07（"45.8% 高度依赖阈值 ⇒ 须改成阈值无关定义"）
#
# 输入：data/processed/seurat_gse169396_drift.rds   （SCT counts；3850 骨系细胞）
#       data/processed/seurat_gse169396_all.rds     （非骨系细胞 → ambient 背景池）
#       results/invariance/07b_invariance.csv       （步 1 结果 → 恒等守卫基准）
# 输出：(--preflight) results/invariance/07c_preflight.csv / _summary.txt
#       (全量)      results/invariance/07c_domain.csv / 07c_domain.pdf / 07c_domain_summary.txt
#
# 用法：Rscript scripts/07c_domain.R --preflight     # 三点预跑
#       Rscript scripts/07c_domain.R                 # 全量 16 格
#
# ----------------------------------------------------------------------------
# 平台口径（2026-09-30 实测决定，见 docs/75 §2）：
#   UCell 打分落在 Seurat【SCT assay 的 counts 层】。实测 rank(counts) == rank(data)
#   逐细胞 30/30 一致（SCT 的 data 是 counts 的逐细胞单调变换）⇒ 可直接在 counts 上打分。
#   ScoreSignatures_UCell(SCT counts, maxRank = 1500) 与存量 AD_drift_ucell /
#   OS_identity_ucell 的 max|Δ| = 0（逐位复刻）。
#   ⇒ 污染模拟在【SCTransform 校正后的 counts 空间】进行；这不重跑 SCTransform，
#     因此 (r = 0, c = 0) 可做到逐位复现步 1（恒等守卫）。
#   ⚠️ 限制：混合发生在校正空间，不等于液滴层面的原始 counts 相加 —— 须写入限制段。
#
#   随机数口径（2026-10-01 修正，裁定 R-3e）：随机种子由【格点标称值 (r, c) 与重抽号 rep】
#     唯一确定，与网格遍历顺序解耦。旧写法 100000*gi + rep_i 依赖 expand.grid 的行序，
#     使同一标称格点在预跑与全量中抽到不同样本（实测 (0.05,0.05)：+23.635 vs +23.074）。
#     ⇒ 修正后三点预跑与全量在同标称格点上得到【同一批样本】，预跑成为全量的有效预演。
#
#   交付口径（2026-10-01，裁定 R-3f）：交付件只保留 doublet 轴 r（即 c = 0 列）。
#     c（ambient）轴实测为【阶跃伪影】（c = 0.001 即饱和，机制见 docs/76 §4）⇒ 不计入交付结论，
#     仅作记录；交付剖面另出 07c_delivery_r.csv。
#
#   合成方案（零安装，不依赖 scDblFinder／SoupX）：
#     doublet  同供体内随机配对两细胞，按 Dirichlet(1,1) 权重相加，替换同数量槽位
#              ⇒ 总细胞数不变
#     ambient  以【同供体非骨系细胞】SCT counts 的池化比例向量为背景，
#              每细胞 new = (1-c) * x + c * (sum(x) / sum(bg)) * bg  ⇒ 每细胞总 counts 不变
#     顺序：先 ambient 后 doublet（真实场景两者共存）
#
#   ⚠️ 判据沿用步 1 的【规则 A】（共同绝对阈值 AD > t & OS > t），网格 = 步 1 的 t ≥ 0 部分
#      （[0, 0.20]，步长 0.005，共 41 点）。恒等守卫 = (0,0) 时 41 点 pct_S1..S4 逐位复现 07b。
# ============================================================================

suppressPackageStartupMessages({ library(Seurat); library(UCell); library(Matrix); library(ggplot2) })

t0 <- Sys.time()
PREFLIGHT <- "--preflight" %in% commandArgs(TRUE)
cat(sprintf("===== 07c 步 2：三参数成立域%s =====\n",
            if (PREFLIGHT) "（三点预跑）" else "（全量）")); flush.console()

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
proc_dir <- file.path(base_dir, "data/processed")
inv_dir  <- file.path(base_dir, "results/invariance")
dir.create(inv_dir, recursive = TRUE, showWarnings = FALSE)

RDS  <- file.path(proc_dir, "seurat_gse169396_drift.rds")
ALL  <- file.path(proc_dir, "seurat_gse169396_all.rds")
BASE <- file.path(inv_dir, "07b_invariance.csv")
tail_name <- if (PREFLIGHT) "07c_preflight" else "07c_domain"

for (f in c(RDS, ALL, BASE)) {
  cat(sprintf("输入 %s md5 %s\n", basename(f), unname(tools::md5sum(f)))); flush.console()
}

# ---- 0. 参数面 -------------------------------------------------------------
if (PREFLIGHT) {
  # 三点预跑：对角 (0,0) / (0.05,0.05) / (0.10,0.10)
  GRID_GRID <- data.frame(r = c(0, 0.05, 0.10), c = c(0, 0.05, 0.10))
  B <- 10L
} else {
  GRID_GRID <- expand.grid(r = c(0, 0.02, 0.05, 0.10), c = c(0, 0.02, 0.05, 0.10))
  B <- 10L
}
GRID_T <- round(seq(0, 0.20, by = 0.005), 3)
DONORS <- c("S1", "S2", "S3", "S4")
STOPIF <- 1e-9
cat(sprintf("网格：%d 格 × B=%d 次重抽；t 网格 %d 点 [%.3f, %.3f]\n",
            nrow(GRID_GRID), B, length(GRID_T), min(GRID_T), max(GRID_T))); flush.console()

# ---- 1. 数据 ---------------------------------------------------------------
dr  <- readRDS(RDS)
cnt <- LayerData(dr, assay = "SCT", layer = "counts")
SM  <- as.character(dr$sample)
gd  <- lapply(DONORS, function(s) which(SM == s)); names(gd) <- DONORS
ND  <- vapply(gd, length, 0L)
cat(sprintf("骨系对象：%d 基因 × %d 细胞；%s\n", nrow(cnt), ncol(cnt),
            paste(sprintf("%s=%d", DONORS, ND), collapse = " "))); flush.console()

GENE_SETS <- list(
  OS_identity = c("RUNX2","SP7","ALPL","BGLAP","SPP1","IBSP","DMP1",
                  "COL1A1","SOST","PHEX","BMP2","BMP4"),
  AD_drift    = c("PPARG","CEBPA","CEBPB","ADIPOQ","LEP","FABP4",
                  "PLIN1","LPL","GPD1","RETN","ACSL1"))
gs <- lapply(GENE_SETS, function(g) intersect(g, rownames(cnt)))
cat(sprintf("基因集（SCT 基因内）：OS %d／AD %d\n",
            length(gs$OS_identity), length(gs$AD_drift))); flush.console()

# 背景池：同供体非骨系细胞的 SCT counts 池化比例
al <- readRDS(ALL)
cnt_all <- LayerData(al, assay = "SCT", layer = "counts")[rownames(cnt), , drop = FALSE]
SMall   <- setNames(as.character(al$sample), colnames(al))   # 必须保留 names：下面按条形码索引
nb_ids  <- setdiff(colnames(cnt_all), colnames(cnt))
bg <- lapply(DONORS, function(s) {
  ix <- nb_ids[SMall[nb_ids] == s]
  v <- Matrix::rowSums(cnt_all[, ix, drop = FALSE]); v / sum(v)
}); names(bg) <- DONORS
cat(sprintf("背景池（非骨系）：合计 %d 细胞；背景向量长度 = %d 基因\n",
            length(nb_ids), nrow(cnt)))
cat(sprintf("非骨系细胞数按供体：%s\n",
            paste(sprintf("%s=%d", DONORS, vapply(DONORS, function(s) sum(SMall[nb_ids] == s), 0L)),
                  collapse = " "))); flush.console()
rm(cnt_all); invisible(gc())

# ---- 2. 打分与判据 ---------------------------------------------------------
score_all <- function(X) ScoreSignatures_UCell(X, features = gs, maxRank = 1500)

analyse <- function(sc) {
  AD <- sc[, "AD_drift_UCell"]; OS <- sc[, "OS_identity_UCell"]
  pct <- vapply(GRID_T, function(t) {
    fl <- AD > t & OS > t
    vapply(DONORS, function(s) 100 * sum(fl[gd[[s]]]) / ND[[s]], 0)
  }, numeric(length(DONORS)))                     # 4 × nt
  D <- pct["S2", ] - pct["S1", ]
  ok <- D > 0
  # 最大连通区间（网格连续即可）
  run <- rle(ok); a <- b <- NA_real_
  if (any(ok)) {
    ends <- cumsum(run$lengths); starts <- ends - run$lengths + 1
    i <- which(run$values & run$lengths == max(run$lengths[run$values]))[1]
    a <- GRID_T[starts[i]]; b <- GRID_T[ends[i]]
    # 若并列，取最靠左的最长段（已由上面保证）
  }
  list(D = D, pct = pct, n_pos = sum(ok), a = a, b = b,
       D0 = D[GRID_T == 0], Dmax = max(D), t_Dmax = GRID_T[which.max(D)])
}

# ---- 3. 恒等守卫：先跑 (r=0, c=0)，必须逐位复现 07b ------------------------
cat("\n[HK] 恒等守卫：(r=0, c=0) 重算须逐位复现 07b_invariance.csv …\n"); flush.console()
thk <- Sys.time()
sc0 <- score_all(cnt)                              # 不改动 ⇒ 应逐位等存量
a0  <- analyse(sc0)
base <- read.csv(BASE, check.names = FALSE)
worst <- 0
for (s in DONORS) {
  col <- paste0("pct_", s)
  worst <- max(worst, max(abs(a0$pct[s, ] -
                 base[[col]][match(sprintf("%.3f", GRID_T), sprintf("%.3f", base$threshold))])))
}
HK <- worst <= STOPIF
cat(sprintf("[HK] 41 点 × 4 供体 最大绝对偏差 %.3e  tol %.0e  => %s（%.1f s）\n",
            worst, STOPIF, ifelse(HK, "PASS", "FAIL"),
            as.numeric(difftime(Sys.time(), thk, units = "secs")))); flush.console()

# ---- 4. 合成函数 -----------------------------------------------------------
Xd <- as.matrix(cnt)                               # 转稠密：BLAS 向量化，便于列替换
cat(sprintf("稠密矩阵 %.1f MB\n", as.numeric(object.size(Xd)) / 1024^2)); flush.console()

# 随机种子：只取决于格点标称值与重抽号，与网格遍历顺序无关（裁定 R-3e）
seed_of <- function(r, c, rep_i) {
  ri <- as.integer(round(r * 100)); ci <- as.integer(round(c * 100))
  as.integer(ri * 10000L + ci * 100L + rep_i)
}

synthesize <- function(r, c, seed) {
  set.seed(seed)
  X <- Xd
  if (c > 0) {                                     # --- ambient ---
    for (s in DONORS) {
      ix  <- gd[[s]]
      tot <- colSums(X[, ix, drop = FALSE])        # 每细胞总 counts
      X[, ix] <- (1 - c) * X[, ix, drop = FALSE] +
                 c * outer(bg[[s]], tot)           # 保持每细胞总 counts 不变
    }
  }
  if (r > 0) {                                     # --- doublet（同供体，替换）---
    for (s in DONORS) {
      ix <- gd[[s]]; n <- length(ix); k <- round(r * n)
      if (k < 1) next
      pos <- sample(ix, k)                         # 被替换的槽位
      ia  <- sample(ix, k, replace = TRUE)
      ib  <- sample(ix, k, replace = TRUE)
      w   <- matrix(rgamma(2L * k, shape = 1), nrow = 2L)
      w   <- sweep(w, 2L, colSums(w), "/")         # Dirichlet(1,1)
      X[, pos] <- w[1, ] * X[, ia, drop = FALSE] + w[2, ] * X[, ib, drop = FALSE]
    }
  }
  X
}

# ---- 5. 主循环 -------------------------------------------------------------
cat(sprintf("\n开始 %d 格 × B=%d ...\n", nrow(GRID_GRID), B)); flush.console()
rows <- list(); ri <- 0L
for (gi in seq_len(nrow(GRID_GRID))) {
  r <- GRID_GRID$r[gi]; c <- GRID_GRID$c[gi]
  for (rep_i in seq_len(B)) {
    ts <- Sys.time()
    X   <- synthesize(r, c, seed = seed_of(r, c, rep_i))
    sc  <- score_all(X)
    an  <- analyse(sc)
    ri  <- ri + 1L
    rows[[ri]] <- data.frame(
      r = r, c = c, rep = rep_i,
      n_grid_D_pos = an$n_pos, support_a = an$a, support_b = an$b,
      holds_full_grid = an$n_pos == length(GRID_T),
      D_at_0 = an$D0, D_max = an$Dmax, t_at_Dmax = an$t_Dmax,
      pct_S1_t0 = an$pct["S1", GRID_T == 0], pct_S2_t0 = an$pct["S2", GRID_T == 0],
      pct_S3_t0 = an$pct["S3", GRID_T == 0], pct_S4_t0 = an$pct["S4", GRID_T == 0],
      secs = as.numeric(difftime(Sys.time(), ts, units = "secs")),
      stringsAsFactors = FALSE)
    cat(sprintf("  [%2d/%2d] r=%.2f c=%.2f rep=%2d → D(0)=%+.3f Dmax=%+.3f 支持[%s, %s] %s (%.1fs)\n",
                ri, nrow(GRID_GRID) * B, r, c, rep_i, an$D0, an$Dmax,
                ifelse(is.na(an$a), "NA", sprintf("%+.3f", an$a)),
                ifelse(is.na(an$b), "NA", sprintf("%+.3f", an$b)),
                ifelse(an$n_pos == length(GRID_T), "全覆盖", sprintf("%d/%d点", an$n_pos, length(GRID_T))),
                as.numeric(difftime(Sys.time(), ts, units = "secs"))))
    flush.console()
    rm(X); invisible(gc(FALSE))
  }
}
tab <- do.call(rbind, rows)

# 按格点汇总
agg <- aggregate(cbind(D_at_0, D_max, support_b, n_grid_D_pos, secs) ~ r + c, data = tab, FUN = mean)
agg$holds_rate <- aggregate(holds_full_grid ~ r + c, data = tab, FUN = mean)$holds_full_grid
agg$support_b_min <- aggregate(support_b ~ r + c, data = tab, FUN = min)$support_b
agg$support_b_max <- aggregate(support_b ~ r + c, data = tab, FUN = max)$support_b

write.csv(tab, file.path(inv_dir, paste0(tail_name, ".csv")), row.names = FALSE)

# ---- 5b. 成立域热图（仅全量）-----------------------------------------------
if (!PREFLIGHT) {
  base_D0 <- agg$D_at_0[agg$r == 0 & agg$c == 0][1]
  agg$r_f <- factor(sprintf("%.2f", agg$r), levels = sprintf("%.2f", sort(unique(agg$r))))
  agg$c_f <- factor(sprintf("%.2f", agg$c), levels = sprintf("%.2f", sort(unique(agg$c))))
  agg$lbl <- sprintf("%+.1f\n%s", agg$D_at_0,
                     ifelse(agg$holds_rate >= 1, "OK", sprintf("%.0f%%", 100 * agg$holds_rate)))
  p <- ggplot(agg, aes(c_f, r_f, fill = D_at_0)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = lbl), size = 3.1, lineheight = 1.05, color = "grey15") +
    scale_fill_gradient2(low = "#4DAF4A", mid = "#FFFFBF", high = "#E41A1C",
                         midpoint = base_D0, name = "mean D(0)\n(pct points)") +
    labs(title = "Step 2 - Support domain of the S2 > S1 ordering",
         subtitle = sprintf("cell = mean D(0) over B = %d resamples; 'OK' = D(t) > 0 for every t in [0, 0.20]\nbaseline (r = 0, c = 0) D(0) = %+.2f; simulation carried out in SCTransform-corrected count space",
                            B, base_D0),
         x = "Ambient fraction c", y = "Doublet fraction r") +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
  pdf(file.path(inv_dir, "07c_domain.pdf"), width = 7.4, height = 6.0)
  print(p); dev.off()
}

# ---- 6. 摘要 ---------------------------------------------------------------
L <- character(0); A <- function(...) L <<- c(L, paste0(...))
A(sprintf("07c 步 2 · 三参数成立域%s —— 结果摘要", if (PREFLIGHT) "（三点预跑）" else ""))
A(sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
A(sprintf("输入 rds md5：%s", unname(tools::md5sum(RDS))))
A(sprintf("输入 all.rds md5：%s", unname(tools::md5sum(ALL))))
A(sprintf("平台口径：UCell 在 SCT counts 层打分（与存量 max|Δ| = 0）"))
A(sprintf("网格：r ∈ {%s} × c ∈ {%s}，B = %d；t 网格 %d 点 [%.3f, %.3f]",
          paste(unique(GRID_GRID$r), collapse = ","), paste(unique(GRID_GRID$c), collapse = ","),
          B, length(GRID_T), min(GRID_T), max(GRID_T)))
A("")
A("【HK】恒等守卫（(r=0,c=0) 逐位复现 07b_invariance.csv）")
A(sprintf("  41 点 × 4 供体 最大绝对偏差 %.3e（tol 1e-9）=> %s", worst, ifelse(HK, "PASS", "FAIL")))
A("")
A("【1】按格点汇总（mean ± 跨 B 次重抽）")
A(paste(sprintf("%5s %5s %10s %10s %22s %9s", "r", "c", "D(0)", "D_max", "支持区间右端[min,max]", "全覆盖率"), collapse = ""))
for (i in seq_len(nrow(agg)))
  A(sprintf("%5.2f %5.2f %+10.3f %+10.3f   %7.3f – %7.3f %8.0f%%",
            agg$r[i], agg$c[i], agg$D_at_0[i], agg$D_max[i],
            agg$support_b_min[i], agg$support_b_max[i], 100 * agg$holds_rate[i]))
A("")
A("【2】决策门 D-C（(r=0.05, c=0.02) 下结论是否仍成立）")
sel <- tab[abs(tab$r - 0.05) < 1e-9 & abs(tab$c - 0.02) < 1e-9, ]
if (nrow(sel)) {
  A(sprintf("  该格点：D(0) = %+.3f（%d/%d 次 D>0）；全覆盖 %d/%d 次",
            mean(sel$D_at_0), sum(sel$D_at_0 > 0), nrow(sel),
            sum(sel$holds_full_grid), nrow(sel)))
} else A("  （本次未含该格点 —— 预跑为对角三点，全量跑后判定）")
A("")
A("【3】交付口径（裁定 R-3f）：只保留 doublet 轴 r（c = 0 列）")
A("  理由：c（ambient）轴实测为阶跃伪影（c = 0.001 即产生全部效应，机制见 docs/76 §4），")
A("        不具剂量反应含义 ⇒ 不计入交付结论，全网格仅作记录（见同目录 07c_domain.csv）。")
del_r <- agg[agg$c == 0, , drop = FALSE]
del_r <- del_r[order(del_r$r), , drop = FALSE]
sel0  <- tab[tab$c == 0, , drop = FALSE]
A(paste(sprintf("%5s %7s %11s %11s %21s %9s",
                "r", "n_rep", "D(0) 均值", "D(0) 最小", "支持右端[min,max]", "全覆盖率"), collapse = ""))
for (i in seq_len(nrow(del_r))) {
  sub <- sel0[abs(sel0$r - del_r$r[i]) < 1e-9, , drop = FALSE]
  A(sprintf("%5.2f %7d %+11.3f %+11.3f   %7.3f – %7.3f %8.0f%%",
            del_r$r[i], nrow(sub), del_r$D_at_0[i], min(sub$D_at_0),
            del_r$support_b_min[i], del_r$support_b_max[i], 100 * del_r$holds_rate[i]))
}
write.csv(del_r[, c("r", "c", "D_at_0", "D_max", "support_b_min", "support_b_max", "holds_rate")],
          file.path(inv_dir, "07c_delivery_r.csv"), row.names = FALSE)
A("")
A(sprintf("用时：%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
A(sprintf("产出：%s", file.path(inv_dir, paste0(tail_name, ".csv"))))
if (!PREFLIGHT) A(sprintf("产出（交付口径 r 轴）：%s", file.path(inv_dir, "07c_delivery_r.csv")))
writeLines(L, file.path(inv_dir, paste0(tail_name, "_summary.txt")))
cat(paste(L, collapse = "\n")); cat("\n===== 07c 完成 =====\n")
