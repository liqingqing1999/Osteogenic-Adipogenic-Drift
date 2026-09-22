#!/usr/bin/env Rscript
# ============================================================================
# fig_scale_diag.R — 开发辅助（helper），**不产出任何图件**
#   用途：打印关键图面板的实际坐标轴范围，供手动核对/统一标尺时参考。
#   不参与论文图件生成链；随仓库保留仅为记录标尺来源。
# ============================================================================
# Diagnostic: print actual axis ranges for the key figure panels.
suppressPackageStartupMessages({ library(Seurat) })
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
proc <- file.path(base_dir, "data/processed")

cat("=== HUMAN drift object ===\n")
seu <- readRDS(file.path(proc, "seurat_gse169396_drift.rds"))
um <- as.data.frame(Embeddings(seu, reduction = "umap.bone"))
cat(sprintf("umap.bone range: x[%.3f, %.3f]  y[%.3f, %.3f]\n",
            min(um[,1]), max(um[,1]), min(um[,2]), max(um[,2])))
cat(sprintf("DriftIndex range: [%.3f, %.3f]\n", min(seu$DriftIndex), max(seu$DriftIndex)))
cat(sprintf("DriftIndex q01/q99: [%.3f, %.3f]\n",
            quantile(seu$DriftIndex, .01), quantile(seu$DriftIndex, .99)))
print(table(seu$sample, seu$group))

cat("\n=== MOUSE bone object ===\n")
m <- readRDS(file.path(proc, "seurat_mouse_metaphysis_bone.rds"))
cat(sprintf("mouse DriftIndex range: [%.3f, %.3f]\n",
            min(m$DriftIndex_mouse), max(m$DriftIndex_mouse)))
cat(sprintf("mouse DriftIndex q01/q99: [%.3f, %.3f]\n",
            quantile(m$DriftIndex_mouse, .01), quantile(m$DriftIndex_mouse, .99)))

cat("\n=== ALL cells (umap.harmony) ===\n")
allc <- readRDS(file.path(proc, "seurat_gse169396_all.rds"))
umh <- as.data.frame(Embeddings(allc, reduction = "umap.harmony"))
cat(sprintf("umap.harmony range: x[%.3f, %.3f]  y[%.3f, %.3f]\n",
            min(umh[,1]), max(umh[,1]), min(umh[,2]), max(umh[,2])))
