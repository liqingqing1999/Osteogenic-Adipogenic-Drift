#!/usr/bin/env python
# ============================================================================
# compose_rest_figures.py — 拼装除 Fig1 以外的全部多面板图（Fig2–4, FigS1–S6）
#   来源：results/figures/PNG_hires/*.png
#   说明：Fig1 由 compose_fig1_uniform.py 单独生成（来源 fig1_panels/），
#         本脚本**不含 Fig1**，避免覆盖。
#   布局：target_h=1100, MAX_W=3000, PAD=55（全仓库统一）
# 输出：results/figures/composite/Fig2.png … FigS6.png  (300 dpi)
# ============================================================================
import os
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ---- 自动定位项目根目录（脚本位于 <root>/figure_assembly/；可用环境变量覆盖）----
ROOT = Path(os.environ.get("DRIFT_PROJECT_DIR", Path(__file__).resolve().parent.parent))
SRC = Path(os.environ.get("DRIFT_PANELS_DIR", ROOT / "results" / "figures" / "PNG_hires"))
OUT = Path(os.environ.get("DRIFT_COMPOSITE_DIR", ROOT / "results" / "figures" / "composite"))
os.makedirs(OUT, exist_ok=True)

# ---- 跨平台粗体无衬线字体回退（找不到候选则退回 PIL 内置位图字体）----
_FONT_CANDIDATES = (
    r"C:\Windows\Fonts\arialbd.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/Library/Fonts/Arial Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
)
_FONT_PATH = next((c for c in _FONT_CANDIDATES if os.path.exists(c)), None)


def _bold_font(size):
    """返回粗体字体；找不到 TTF 时退回 PIL 内置位图字体。"""
    return ImageFont.truetype(_FONT_PATH, size) if _FONT_PATH else ImageFont.load_default()

GAP_PANEL = 25
GAP_ROW = 35
PAD = 55
MAX_W = 3000
LETTERS = "ABCDEFGHIJKL"

PLAN = {
    "Fig2": {
        "rows": [
            ["02_DriftIndex_by_group", "02_cell_state_proportion"],
            [["02_DriftIndex_age_dimension", "02_hybrid_pct_age"]],
        ],
        "h": 1100,
    },
    "Fig3": {
        "rows": [
            ["04_trajectory_umap", "04_drift_vs_pseudotime"],
            ["04_genes_along_trajectory"],
        ],
        "h": 1100,
    },
    "Fig4": {
        "rows": [
            ["05_cross_species_driftindex", "05_cross_species_hybrid_pct"],
        ],
        "h": 1100,
    },
    "FigS1": {
        "rows": [
            ["01_QC_violin_prefilter"],
            ["01_QC_violin_postfilter"],
        ],
        "h": 1100,
    },
    "FigS2": {
        "rows": [
            ["01_UMAP_integration_compare"],
            ["01_DotPlot_lineage_markers"],
        ],
        "h": 1100,
    },
    "FigS3": {"rows": [["03_sensitivity_trend"]], "h": 1100},
    "FigS4": {"rows": [["05_mouse_metaphysis_cell_state"]], "h": 1100},
    "FigS5": {
        "rows": [
            ["02b_DriftIndex_by_sample", "02b_hybrid_pct_by_sample"],
        ],
        "h": 1100,
    },
    "FigS6": {
        "rows": [
            ["08_shape_smoother", "08_shape_deciles"],
        ],
        "h": 1100,
    },
}


def open_img(name):
    return Image.open(os.path.join(SRC, name + ".png")).convert("RGB")


def resize_to_h(im, h):
    return im.resize((round(im.width * h / im.height), h), Image.LANCZOS)


def stitch_horizontal(images, gap):
    h = images[0].height
    total_w = sum(im.width for im in images) + gap * (len(images) - 1)
    canvas = Image.new("RGB", (total_w, h), "white")
    x = 0
    for im in images:
        canvas.paste(im, (x, 0))
        x += im.width + gap
    return canvas


def build_cell(spec):
    if isinstance(spec, list):
        ims = [open_img(n) for n in spec]
        h = min(im.height for im in ims)
        ims = [resize_to_h(im, h) for im in ims]
        return stitch_horizontal(ims, GAP_PANEL)
    return open_img(spec)


def draw_label(im, lbl):
    d = ImageDraw.Draw(im)
    r = 26
    font = _bold_font(30)
    d.ellipse([12, 12, 12 + 2 * r, 12 + 2 * r], fill="white", outline="black", width=2)
    bbox = d.textbbox((0, 0), lbl, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    cx, cy = 12 + r, 12 + r
    d.text((cx - tw / 2 - bbox[0], cy - th / 2 - bbox[1]), lbl,
           fill="black", font=font)


def render_row(cell_specs, target_h, letter_iter, draw_labels=True):
    panels = []
    for spec in cell_specs:
        cell = resize_to_h(build_cell(spec), target_h)
        w, h = cell.size
        frame = Image.new("RGB", (w + 2 * PAD, h + 2 * PAD), "white")
        frame.paste(cell, (PAD, PAD))
        panels.append(frame)

    total_w = sum(p.width for p in panels) + GAP_PANEL * (len(panels) - 1)
    row_h = max(p.height for p in panels)
    if total_w > MAX_W:
        scale = MAX_W / total_w
        panels = [p.resize((max(1, round(p.width * scale)),
                            max(1, round(p.height * scale))), Image.LANCZOS)
                  for p in panels]
        total_w = sum(p.width for p in panels) + GAP_PANEL * (len(panels) - 1)
        row_h = max(p.height for p in panels)

    if draw_labels:
        for p in panels:
            draw_label(p, next(letter_iter))

    row_img = Image.new("RGB", (total_w, row_h), "white")
    x = 0
    for p in panels:
        row_img.paste(p, (x, (row_h - p.height) // 2))
        x += p.width + GAP_PANEL
    return row_img


def main():
    for fig, cfg in PLAN.items():
        letter_iter = iter(LETTERS)
        rows = cfg["rows"]
        draw_labels = sum(len(r) for r in rows) > 1
        row_imgs = [render_row(r, cfg["h"], letter_iter, draw_labels) for r in rows]
        canvas_w = max(r.width for r in row_imgs)
        canvas_h = sum(r.height for r in row_imgs) + GAP_ROW * (len(row_imgs) - 1)
        canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
        y = 0
        for r in row_imgs:
            canvas.paste(r, ((canvas_w - r.width) // 2, y))
            y += r.height + GAP_ROW
        canvas.save(os.path.join(OUT, fig + ".png"), dpi=(300, 300))
        print(f"{fig}: {canvas.size}")
    print("Done.")


if __name__ == "__main__":
    main()
