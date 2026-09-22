#!/usr/bin/env Rscript
# ============================================================================
# fig1_render_uniform.R — 用统一坐标范围 + 统一画布重绘 Fig1 的 3 张 UMAP
#   目的：让 3 张散点图在几何上完全等比（相同 xlim/ylim + coord_equal），
#         先渲染成独立 PNG，再由 compose_fig1_uniform.py 拼成单行 3 面板版
#         （A1 Cell-level DriftIndex / A3 Sample / B Cell state）。
# 输入：data/processed/seurat_gse169396_drift.rds
# 输出：results/figures/fig1_panels/{A1,A3,B}.png
# 用法：Rscript scripts/fig1_render_uniform.R
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
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
out_dir  <- file.path(base_dir, "results/figures/fig1_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(file.path(base_dir, "data/processed/seurat_gse169396_drift.rds"))
cat(sprintf("cells: %d\n", ncol(seu)))

# ---- 1. 提取 UMAP 坐标 -------------------------------------
umap <- as.data.frame(Embeddings(seu, reduction = "umap.bone"))
colnames(umap) <- c("UMAP1", "UMAP2")
umap$sample    <- seu$sample
umap$group     <- seu$group
umap$cell_state <- seu$cell_state
umap$DriftIndex <- seu$DriftIndex

# ---- 2. 统一坐标范围（保留少量边距）--------------------------
# 用全部细胞的坐标范围 + padding，保证 6 张图的散点几何一致
xlim <- range(umap$UMAP1)
ylim <- range(umap$UMAP2)
span <- max(diff(xlim), diff(ylim))
pad  <- 0.06 * span
xlim <- c(min(xlim) - pad, max(xlim) + pad)
ylim <- c(min(ylim) - pad, max(ylim) + pad)
cat(sprintf("Uniform extent: x[%.2f,%.2f] y[%.2f,%.2f]\n", xlim[1], xlim[2], ylim[1], ylim[2]))

# ---- 3. 统一主题 / 画布 / 配色 ---------------------------------------------
base_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title     = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title     = element_text(size = 13),
    axis.text      = element_text(size = 11),
    legend.title   = element_blank(),
    legend.text    = element_text(size = 12),
    legend.position = "right",
    plot.margin    = margin(6, 6, 6, 6)
  )

col_sample <- c("S1" = "#E41A1C", "S2" = "#FF7F00",
                "S3" = "#984EA3", "S4" = "#4DAF4A")
col_state <- c("Uncommitted" = "#BDBDBD", "Osteo-biased" = "#377EB8",
               "Adipo-biased" = "#984EA3", "Hybrid" = "#E41A1C")
col_drift <- c("#377EB8", "white", "#E41A1C")

# Use a wider canvas so the longest legend (B: 4 cell-state names) does not
# squeeze the plot region. With coord_equal, a narrower plot panel becomes
# shorter, making the scatter appear smaller than A1/A3. Pad to 8.6 inch.
DEV_W <- 8.6
DEV_H <- 6.4
DPI   <- 300

# ---- 4. 绘图函数 ------------------------------------------------------------
plot_umap <- function(data, aes_color, title, scale) {
  p <- ggplot(data, aes(x = UMAP1, y = UMAP2))
  if (inherits(scale, "continuous")) {
    p <- p + geom_point(aes(colour = .data[[aes_color]]), size = 0.9, alpha = 0.85) +
         scale_colour_gradientn(colours = col_drift, name = aes_color)
  } else {
    p <- p + geom_point(aes(colour = .data[[aes_color]]), size = 0.9, alpha = 0.85) +
         scale
  }
  p + coord_equal(xlim = xlim, ylim = ylim) +
      labs(x = "UMAP 1", y = "UMAP 2", title = title) +
      base_theme
}

# ---- 5. 渲染 3 张（A1/A3/B —— Fig1 实际使用的面板）-------------------------
panels <- list(
  A1 = list(color = "DriftIndex", title = "Cell-level DriftIndex (AD \u2212 OS)",
            scale = NULL),
  A3 = list(color = "sample", title = "Sample",
            scale = scale_colour_manual(values = col_sample)),
  B  = list(color = "cell_state", title = "Cell state",
            scale = scale_colour_manual(values = col_state))
)

for (nm in names(panels)) {
  spec <- panels[[nm]]
  if (nm == "A1") {
    p <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
         geom_point(aes(colour = DriftIndex), size = 0.9, alpha = 0.85) +
         scale_colour_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                                midpoint = 0) +
         coord_equal(xlim = xlim, ylim = ylim) +
         labs(x = "UMAP 1", y = "UMAP 2", title = spec$title, colour = "DriftIndex") +
         base_theme
  } else {
    p <- plot_umap(umap, spec$color, spec$title, spec$scale)
  }
  f <- file.path(out_dir, paste0(nm, ".png"))
  ggsave(f, p, width = DEV_W, height = DEV_H, dpi = DPI)
  cat(sprintf("  %s -> %s\n", nm, f))
}

cat("Done.\n")
