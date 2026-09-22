#!/usr/bin/env Rscript
# ============================================================================
# figS1_scale_uniform.R — FigS1（QC 小提琴）pre / post 统一 y 标尺
#   原图 pre 与 post 各自自适应 y 轴，看不出过滤效果；
#   本脚本对三个指标分别取 pre∪post 的最大值作为统一上限。
# 输出：results/figures/PNG_hires/01_QC_violin_prefilter.png
#       results/figures/PNG_hires/01_QC_violin_postfilter.png
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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_sample <- c("S1" = "#E41A1C", "S2" = "#FF7F00", "S3" = "#984EA3", "S4" = "#4DAF4A")

cat("===== FigS1：读取 4 个原始 10x 样本 =====\n")
sids <- c("S1", "S2", "S3", "S4")
obj_list <- lapply(sids, function(sid) {
  d <- file.path(base_dir, "data/raw/10x_standard", sid)
  mtx <- Seurat::Read10X(d)
  obj <- CreateSeuratObject(mtx, project = sid, min.cells = 3, min.features = 200)
  obj[["sample"]] <- sid
  cat(sprintf("  %s: %d cells\n", sid, ncol(obj)))
  obj
})
seu <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
cat(sprintf("合并后: %d cells\n", ncol(seu)))

seu_post <- subset(seu, subset = nFeature_RNA > 500 &
                            nFeature_RNA < 6000 &
                            percent.mt < 20)
cat(sprintf("过滤后: %d cells\n", ncol(seu_post)))

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
ymax <- sapply(feats, function(f)
  max(c(seu[[f]][, 1], seu_post[[f]][, 1]), na.rm = TRUE))
cat("统一 y 上限: ", paste(sprintf("%s=%.0f", feats, ymax), collapse = "  "), "\n")

mk <- function(obj, feat) {
  VlnPlot(obj, features = feat, pt.size = 0, group.by = "sample",
          cols = col_sample) +
    scale_y_continuous(limits = c(0, ymax[[feat]])) +
    labs(title = feat) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"))
}

p_pre  <- wrap_plots(lapply(feats, mk, obj = seu),      ncol = 3)
p_post <- wrap_plots(lapply(feats, mk, obj = seu_post), ncol = 3)

ggsave(file.path(out_dir, "01_QC_violin_prefilter.png"),  p_pre,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(out_dir, "01_QC_violin_postfilter.png"), p_post,
       width = 12, height = 5, dpi = 300)
cat("完成：01_QC_violin_prefilter.png / 01_QC_violin_postfilter.png（统一 y 标尺）\n")
