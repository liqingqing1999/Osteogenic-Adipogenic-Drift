#!/usr/bin/env Rscript
# ============================================================================
# figS2_stratum_panel.R — 重制 FigS2(A) 中面板：按 densitometric stratum 着色
#   面板改为按密度分层着色（仅 S1/S2 可分层，S3/S4 未分层），
#   不再使用此前的三组 disease-group 着色：
#     S1 = Osteoporotic range (S1) / S2 = Osteopenic range (S2)
#     S3, S4 = Not stratified
#   输出：results/figures/PNG_hires/01_UMAP_integration_compare.png
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
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
out_dir  <- file.path(base_dir, "results/figures/PNG_hires")

LOG <- file.path(base_dir, "logs", "figS2_stratum.log")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, "w"); wr <- function(...) { writeLines(paste0(...), con); flush(con) }
wr("=== FigS2 middle panel re-render ===")

seu <- readRDS(file.path(base_dir, "data/processed/seurat_gse169396_all.rds"))
wr("all cells: ", ncol(seu))

sample_chr <- as.character(seu$sample)
seu$stratum <- ifelse(sample_chr == "S1", "Osteoporotic range (S1)",
                ifelse(sample_chr == "S2", "Osteopenic range (S2)",
                       "Not stratified (S3/S4)"))
seu$stratum <- factor(seu$stratum,
                      levels = c("Osteoporotic range (S1)", "Osteopenic range (S2)",
                                 "Not stratified (S3/S4)"))
writeLines(capture.output(print(table(seu$sample, seu$stratum))), con); flush(con)

col_stratum <- c("Osteoporotic range (S1)" = "#E41A1C",
                 "Osteopenic range (S2)"   = "#FF7F00",
                 "Not stratified (S3/S4)"  = "#BDBDBD")
col_sample  <- c("S1" = "#E41A1C", "S2" = "#FF7F00", "S3" = "#984EA3", "S4" = "#4DAF4A")

p1 <- DimPlot(seu, reduction = "umap.harmony", group.by = "sample", label = FALSE,
              cols = col_sample) +
  ggtitle("By sample") + theme(plot.title = element_text(size = 12, hjust = 0.5))
p2 <- DimPlot(seu, reduction = "umap.harmony", group.by = "stratum", label = FALSE,
              cols = col_stratum) +
  ggtitle("By densitometric stratum") + theme(plot.title = element_text(size = 12, hjust = 0.5))
p3 <- DimPlot(seu, reduction = "umap.harmony", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("By cluster") + theme(plot.title = element_text(size = 12, hjust = 0.5))

ggsave(file.path(out_dir, "01_UMAP_integration_compare.png"), p1 | p2 | p3,
       width = 16, height = 6, dpi = 300)
wr("  01_UMAP_integration_compare.png written")
wr("=== done ===")
close(con)
