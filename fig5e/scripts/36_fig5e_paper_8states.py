#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 5E 论文原图处理链：PDF 提取 → 删除缺失列 → 调色 → 文字重绘

背景
----
论文 Figure 5E 是 8 基因 × 14 状态的 dot plot。本地只覆盖其中 8 个状态
（YS 的 5 个阶段 + 胎儿皮肤 2 组 + TREM2 亚组与之一一对应的 1 组），
缺 AGM / 脑 / 性腺及其 TREM2 亚组共 6 列。为了让"论文原图 vs 本次复现"
在 PPT 里同口径对比，本脚本直接从论文 PDF 的原图上删掉这 6 列，
其余部分（基因名、色标、图例、外框、E 标签）按原样保留。

处理流程
--------
1. 用 PyMuPDF 从论文 PDF 中提取内嵌的 Figure 5 原始位图（1600×2100）
2. 裁出 5E 面板（注意下边界要留到 y≈1680，否则长竖排标签的底部会被切掉）
3. 按点峰线性拟合出的列边界裁剪，保留 8 列后重新拼接（左段 + 列段 + 右段）
4. 调色：文字 → 灰（#555）、点 → 亮度 +26%、边框 → 浅灰
5. 擦除原文字，用 Arial 重绘（字号可精确控制，比缩放位图清晰）
6. 裁掉底部空白

用法
----
    # 论文 PDF 放到 fig5e/data/（或用环境变量指定）
    export PAPER_PDF=/path/to/Science2023_Goh_add7564.pdf
    python fig5e/scripts/36_fig5e_paper_8states.py
    # -> fig5e/figures/fig5e_paper_8states.png

依赖：PyMuPDF（fitz）、Pillow、numpy
输出：brepo/figures/fig5e_paper_8states.png
"""
import os
import numpy as np
import pymupdf
from PIL import Image, ImageDraw, ImageFont

# =============================================================================
# 配置
# =============================================================================
# 路径相对仓库定位：fig5e/scripts -> 仓库根
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
TMP = os.path.join(REPO, ".cache")
OUT_FIG = os.path.join(REPO, "fig5e/figures/fig5e_paper_8states.png")
# 论文全文 PDF：默认放 fig5e/data/，也可用环境变量 PAPER_PDF 指定
PDF = os.environ.get("PAPER_PDF",
                     os.path.join(REPO, "fig5e/data/Science2023_Goh_add7564.pdf"))

# 5E 面板在 Figure 5 原图中的裁剪框（左, 上, 右, 下）
# 注意：下边界 1680 是为保住长竖排标签（"Fetal skin TREM2 Mac."）的底部；
#       裁到 1642 会切掉标签开头的 "Fet"
PANEL_BOX = (10, 1256, 795, 1680)

# 14 列的几何参数（由列方向点峰的线性拟合得到，误差 < 2px）
COL_W = 21.54          # 列宽（像素）
COL1_CX = 125.0        # 第 1 列中心 x
N_COLS = 14

# 论文原图的 14 个状态（从左到右，放大标签区逐字读出）
LABELS_14 = [
    "YS MOP", "YS promonocyte", "YS monocyte", "YS pre Mac.", "YS macrophage",
    "AGM macrophage", "Fetal skin Mac.", "Fetal brain Mac.", "Fetal gonads Mac.",
    "YS TREM2 Mac.", "AGM TREM2 Mac.", "Fetal skin TREM2 Mac.",
    "Fetal brain TREM2 Mac.", "Fetal gonads TREM2 Mac.",
]
# 本地能覆盖的 8 列（0-indexed）：1-5、7、10、12
KEEP = [0, 1, 2, 3, 4, 6, 9, 11]
GENES = ["CD4", "CD14", "C1QA", "CSF1R", "TREM2", "CX3CR1", "P2RY12", "OLFML3"]

# 调色参数
TEXT_GREY = 85         # 文字灰 #555
DOT_BRIGHT = 0.74      # 点：255-(255-v)*0.74 → 亮度约 +26%
BORDER_ADD = 105       # 边框灰底

# 重绘字号
FONT_GENE = 12         # 基因名（斜体）
FONT_STATE = 11        # 状态标签（竖排）
FONT_E = 17            # 面板标签 E（粗体）

# 点图在拼接图中的几何（上/下边框与列中心，实测值）
PLOT_TOP, PLOT_BOTTOM = 74, 255
# 拼接后 8 列的中心（由点峰实测；列间距 21.5±0.5 像素，并非严格等距，
# 故直接写实测值而不是用 COL_W 递推，否则个别列会差 1px）
SPLICED_COL_CENTERS = [126, 147, 169, 190, 212, 234, 255, 277]


# =============================================================================
# 1. 从 PDF 提取 Figure 5 原图（若已存在则跳过）
# =============================================================================
def extract_fig5() -> str:
    """论文 PDF 第 43 页内嵌 1600×2100 的 Figure 5 原始位图。"""
    dst = os.path.join(TMP, "fig5_original.png")
    if os.path.exists(dst):
        print(f"[1] 已存在，跳过提取：{dst}")
        return dst
    doc = pymupdf.open(PDF)
    # 该页内嵌一张 1600×2100 的位图（xref 通过 get_images 定位）
    page = doc[42]
    imgs = page.get_images(full=True)
    target = max(imgs, key=lambda im: im[2] * im[3])   # 取最大的一张
    pix = doc.extract_image(target[0])
    open(dst, "wb").write(pix["image"])
    doc.close()
    print(f"[1] 已提取 Figure 5 原图：{dst} ({pix['width']}x{pix['height']})")
    return dst


# =============================================================================
# 2. 裁出 5E 面板
# =============================================================================
def crop_panel(fig_path: str) -> Image.Image:
    im = Image.open(fig_path).convert("RGB")
    panel = im.crop(PANEL_BOX)
    print(f"[2] 5E 面板：{panel.size[0]}x{panel.size[1]}")
    return panel


# =============================================================================
# 3. 按列裁剪并拼接保留的 8 列
# =============================================================================
def drop_columns(panel: Image.Image) -> np.ndarray:
    a = np.array(panel)
    # 15 个列边界（第 i 列占 [bounds[i], bounds[i+1])）
    bounds = [COL1_CX - COL_W / 2 + COL_W * i for i in range(N_COLS + 1)]

    left = a[:, :int(round(bounds[0]))]                       # 基因名 + E + 左边框 + padding
    cols = [a[:, int(round(bounds[k])):int(round(bounds[k + 1]))] for k in KEEP]
    right = a[:, int(round(bounds[N_COLS])):]                 # padding + 右边框 + 色标 / 图例

    out = np.concatenate([left] + cols + [right], axis=1)
    print(f"[3] 删 {N_COLS - len(KEEP)} 列，保留 {len(KEEP)} 列 -> {out.shape[1]}x{out.shape[0]}")
    print("    保留：" + " / ".join(LABELS_14[k] for k in KEEP))
    print("    删除：" + " / ".join(LABELS_14[i] for i in range(N_COLS) if i not in KEEP))
    return out


# =============================================================================
# 4. 调色：文字变灰、点提亮、边框浅灰
# =============================================================================
def recolor(base: np.ndarray) -> Image.Image:
    H, W, _ = base.shape
    img = base.astype(float).copy()
    lum = img.mean(axis=2)
    sat = img.max(axis=2) - img.min(axis=2)

    # 4a. 文字区（只在这几块区域内改，避免误伤边框与点）
    text_regions = np.zeros((H, W), bool)
    text_regions[55:268, 0:104] = True        # 基因名
    text_regions[272:H, 100:300] = True       # x 轴状态标签
    text_regions[55:270, 430:W] = True        # 色标 / 图例文字
    text_regions[10:70, 0:45] = True          # E 标签
    tm = text_regions & (sat < 40) & (lum < 200)
    img[tm] = TEXT_GREY + img[tm] * 0.30      # 纯黑 -> #555，抗锯齿边缘同步变浅
    print(f"[4] 文字像素调灰：{int(tm.sum())}")

    # 4b. 点图区彩色像素提亮
    dot_region = np.zeros((H, W), bool)
    dot_region[74:272, 104:300] = True
    cm = dot_region & (sat >= 25) & (lum < 250)
    img[cm] = 255 - (255 - img[cm]) * DOT_BRIGHT
    print(f"[4] 点像素提亮：{int(cm.sum())}")

    # 4c. 边框等低饱和暗线轻微提亮（不要比文字更黑）
    bm = (sat < 25) & (lum < 170) & ~text_regions
    img[bm] = img[bm] * 0.35 + BORDER_ADD * 0.65

    return Image.fromarray(np.clip(img, 0, 255).astype(np.uint8))


# =============================================================================
# 5. 擦除原文字并用 Arial 重绘
# =============================================================================
def redraw_text(im: Image.Image) -> Image.Image:
    d = ImageDraw.Draw(im)

    # 5a. 擦除：基因名区、状态标签区（必须从下边框 255 之后开始，否则原标签
    #     顶部会残留成 "IP"/"le"/"c." 之类的碎片）、E 标签区
    d.rectangle([0, 50, 104, 260], fill="white")
    d.rectangle([106, 257, 300, im.size[1]], fill="white")
    d.rectangle([5, 5, 50, 72], fill="white")

    font_gene = ImageFont.truetype("C:/Windows/Fonts/ariali.ttf", FONT_GENE)
    font_state = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", FONT_STATE)
    font_e = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", FONT_E)

    # 5b. 基因名：右对齐到点图左边界（x=102），垂直居中到每一行
    row_h = (PLOT_BOTTOM - PLOT_TOP) / len(GENES)
    for i, g in enumerate(GENES):
        cy = PLOT_TOP + row_h * (i + 0.5)
        bb = d.textbbox((0, 0), g, font=font_gene)
        d.text((102 - (bb[2] - bb[0]), cy - (bb[3] + bb[1]) / 2),
               g, font=font_gene, fill=(TEXT_GREY,) * 3)

    # 5c. 状态标签：竖排（rotate 90 逆时针），底部对齐点图下边框，水平居中到列
    col_centers = SPLICED_COL_CENTERS
    for label, cx in zip([LABELS_14[k] for k in KEEP], col_centers):
        bb = d.textbbox((0, 0), label, font=font_state)
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        tmp = Image.new("RGBA", (w + 6, h + 6), (255, 255, 255, 0))
        ImageDraw.Draw(tmp).text((3 - bb[0], 3 - bb[1]), label,
                                 font=font_state, fill=(TEXT_GREY,) * 3)
        tmp = tmp.rotate(90, expand=True)
        im.paste(tmp, (cx - tmp.width // 2, PLOT_BOTTOM + 4), tmp)

    # 5d. 面板标签 E
    d.text((10, 8), "E", font=font_e, fill=(TEXT_GREY,) * 3)
    print(f"[5] 重绘：基因名 {FONT_GENE}px 斜体 / 状态 {FONT_STATE}px 竖排 / E {FONT_E}px")
    return im


# =============================================================================
# 6. 裁掉底部空白
# =============================================================================
def trim_bottom(im: Image.Image) -> Image.Image:
    a = np.array(im).astype(int)
    nonwhite = np.abs(a - 255).sum(axis=2) > 30
    rows = np.where(nonwhite.any(axis=1))[0]
    bottom = min(rows.max() + 10, im.size[1])
    out = im.crop((0, 0, im.size[0], bottom))
    print(f"[6] 裁底：{im.size[1]} -> {out.size[1]}（比例 {out.size[0]/out.size[1]:.3f}）")
    return out


# =============================================================================
# 主流程
# =============================================================================
def main():
    os.makedirs(TMP, exist_ok=True)
    fig = extract_fig5()
    panel = crop_panel(fig)
    base = drop_columns(panel)
    Image.fromarray(base).save(os.path.join(TMP, "5e_base.png"))   # 中间产物：未调色拼接图
    im = recolor(base)
    im = redraw_text(im)
    im = trim_bottom(im)

    im.save(OUT_FIG)
    print(f"\n完成：{OUT_FIG}  ({im.size[0]}x{im.size[1]})")


if __name__ == "__main__":
    main()
