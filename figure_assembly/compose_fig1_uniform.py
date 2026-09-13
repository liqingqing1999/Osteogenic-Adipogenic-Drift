#!/usr/bin/env python
# ============================================================================
# compose_fig1_uniform.py — 把 3 张**等尺寸** UMAP 子图拼成 Fig1（单行 3 面板）
#   Row 1 : A1  A3  B           (Cell-level DriftIndex / Sample / Cell state)
#   面板来源：results/figures/fig1_panels/（由 scripts/fig1_render_uniform.R
#   用统一 xlim/ylim + coord_equal 渲染，三张几何等比）。
#   拼版时统一字母标签 A/B/C，等宽并排。
# 输出：results/figures/composite/Fig1.png  (300 dpi)
# ============================================================================
import os
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ---- 自动定位项目根目录（脚本位于 <root>/figure_assembly/；可用环境变量覆盖）----
ROOT = Path(os.environ.get("DRIFT_PROJECT_DIR", Path(__file__).resolve().parent.parent))
SRC = Path(os.environ.get("DRIFT_PANELS_DIR", ROOT / "results" / "figures" / "fig1_panels"))
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

GAP_PANEL = 30    # panel 间距
GAP_ROW   = 120   # 行间距（保留以兼容旧版，不使用）
TARGET_W  = 3000  # 目标输出宽度（10 inch @ 300 dpi，整页横向 figure 上限）
LABEL_X, LABEL_Y = 18, 18

# 单行布局：3 张 UMAP 等宽并排
ROWS = {
    "Row1": ["A1", "A3", "B"],   # DriftIndex / Sample / Cell state
}
# 每个子图的独立 panel 标签：A / B / C
LABELS = {
    "A1": "A",
    "A3": "B",
    "B":  "C",
}


def open_panel(name):
    return Image.open(os.path.join(SRC, name + ".png")).convert("RGB")


def make_label(im, label):
    """在左上角画独立 panel 字母标签（不覆盖散点）。"""
    im = im.copy()
    d = ImageDraw.Draw(im)
    r = 30
    font = _bold_font(34)
    d.ellipse([LABEL_X, LABEL_Y, LABEL_X + 2 * r, LABEL_Y + 2 * r],
              fill="white", outline="black", width=3)
    bbox = d.textbbox((0, 0), label, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    cx = LABEL_X + r
    cy = LABEL_Y + r
    d.text((cx - tw / 2 - bbox[0], cy - th / 2 - bbox[1]), label,
           fill="black", font=font)
    return im


def build_row(names):
    """把一行内等尺寸 panel 横向拼接，并给每个 panel 打上字母标签。"""
    panels = [make_label(open_panel(n), LABELS[n]) for n in names]
    h = max(p.height for p in panels)
    w = sum(p.width for p in panels) + GAP_PANEL * (len(panels) - 1)
    row = Image.new("RGB", (w, h), "white")
    x = 0
    for p in panels:
        y = (h - p.height) // 2
        row.paste(p, (x, y))
        x += p.width + GAP_PANEL
    return row


def main():
    # 1) 读 panel 原图（1860x1800 @ 300 dpi = 6.2 inch x 6.0 inch）
    panel_imgs = [open_panel(n) for n in ROWS["Row1"]]
    n = len(panel_imgs)

    # 2) 计算每 panel 目标宽度：使总宽 = TARGET_W
    target_panel_w = (TARGET_W - GAP_PANEL * (n - 1)) // n
    panel_h = int(panel_imgs[0].height * target_panel_w / panel_imgs[0].width)

    # 3) 缩放 + 加 label
    labeled = []
    for name, im in zip(ROWS["Row1"], panel_imgs):
        im_s = im.resize((target_panel_w, panel_h), Image.LANCZOS)
        labeled.append(make_label(im_s, LABELS[name]))

    # 4) 横向拼接
    row_w = sum(p.width for p in labeled) + GAP_PANEL * (n - 1)
    row = Image.new("RGB", (row_w, panel_h), "white")
    x = 0
    for p in labeled:
        row.paste(p, (x, 0))
        x += p.width + GAP_PANEL

    out = os.path.join(OUT, "Fig1.png")
    row.save(out, dpi=(300, 300))
    print(f"Fig1: {row.size}")

    # 无标签预览
    old_labels = dict(LABELS)
    for k in LABELS:
        LABELS[k] = ""
    preview_imgs = [im.resize((target_panel_w, panel_h), Image.LANCZOS)
                    for im in panel_imgs]
    pw = sum(p.width for p in preview_imgs) + GAP_PANEL * (n - 1)
    preview = Image.new("RGB", (pw, panel_h), "white")
    x = 0
    for p in preview_imgs:
        preview.paste(p, (x, 0)); x += p.width + GAP_PANEL
    LABELS.clear(); LABELS.update(old_labels)
    out2 = os.path.join(OUT, "Fig1_nolabel_preview.png")
    preview.save(out2, dpi=(300, 300))
    print(f"Fig1_nolabel_preview: {preview.size}")


if __name__ == "__main__":
    main()