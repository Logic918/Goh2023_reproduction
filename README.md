# Goh et al. (2023) 人体早期发育细胞图谱 · 复现作图代码

> 复现 **Science 2023《Yolk sac cell atlas reveals multiorgan functions during human early development》**
> 的三张主图：**Figure 2D / Figure 3 / Figure 5E**
> 含数据来源、环境依赖、运行顺序、关键作图口径与已知限制

复现对象：Goh I, Botting RA, Rose A, et al. *Yolk sac cell atlas reveals multiorgan functions during human early development.* **Science** 2023;381:eadd7564.

---

## 1. 复现图一览

| 图 | 内容 | 脚本 | 覆盖度 |
|---|---|---|---|
| **Fig 2D** | 跨组织 marker 点图：34 基因 × 8 个细胞群（YS / 胎肝 / EL / 成人肝 / 胎肾） | `fig2d/scripts/11 → 50 → 40 → 41` | 完整（8/8 列） |
| **Fig 3** | HSPC 分化：A 点图（7 行 × 23 基因）、B early/definitive 比例、E 珠蛋白开关 | `fig3/scripts/03_fig3_reproduce.R` | A/B/E 三个面板（原图 8 个面板） |
| **Fig 5E** | 跨器官巨噬 marker 点图：8 基因 × 8 状态 | `fig5e/scripts/36_fig5e_paper_8states.py` | 8 基因全覆盖；状态 8/14（删除无数据的 6 列） |

---

## 2. 仓库结构

```
Goh2023_reproduction/
├── README.md
├── fig2d/                          # Figure 2D：跨组织 marker 点图
│   ├── scripts/
│   │   ├── 11_adult_liver.R        # ① 成人肝对照列（GSE124395）
│   │   ├── 50_add_kidney.R         # ② 胎肾 h5ad 结构确认
│   │   ├── 40_build_final_db.R     # ③ 核心：8 套数据整合统计
│   │   ├── 41_plot_fig2d_v2.R      # ④ 最终出图 + 原图对比
│   │   ├── 10_fig2d_pipeline.R     # 第一版主流水线（备查）
│   │   ├── 12_plot_fig2d.R         # 第一版点图版式（备查）
│   │   ├── 31_check_alignment.R    # 校验：注释与矩阵对齐
│   │   ├── 32_test_ys_alignment.R  # 校验：YS 内胚层注释
│   │   ├── 34_check_f158.R         # 校验：F158/WE_YS 取数方案
│   │   ├── _verify_geneorder.R     # 校验：基因顺序
│   │   ├── _verify_cellorder.R     # 校验：细胞顺序
│   │   ├── _verify_mapping.R       # 校验：映射关系
│   │   └── _crop_fig2d.R           # 从论文 PDF 裁原图
│   ├── data/                       # 原始数据（git 忽略，见 data/README.md）
│   ├── results/                    # 统计库与缓存
│   └── figures/                    # 出图结果
├── fig3/                           # Figure 3：HSPC 分化
│   ├── scripts/03_fig3_reproduce.R
│   ├── data/  figures/
├── fig5e/                          # Figure 5E：跨器官巨噬点图
│   ├── scripts/36_fig5e_paper_8states.py
│   ├── data/                       # 论文 PDF 放这里
│   └── figures/
├── docs/
│   ├── DIAGNOSTIC_SCRIPTS.md       # 37 个排障脚本清单（未全文收录）
│   └── PIPELINE_SCRIPTS.md         # PPT / 报告生成脚本清单
├── .gitattributes                  # 统一 LF 换行
└── .gitignore
```

---

## 3. 环境依赖

**R 4.6.1**（Fig 2D / Fig 3）＋ **Python 3.11**（Fig 5E）

```r
# ---- R ----
install.packages(c("data.table", "Matrix", "ggplot2", "ggforce",
                   "patchwork", "scales", "grid", "magick", "pdftools"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")   # 读 CELLxGENE 的 h5ad
```

```bash
# ---- Python ----
pip install pymupdf pillow numpy
```

| 包 | 用于 |
|---|---|
| `data.table` / `Matrix` | 大矩阵分块流式读取、稀疏矩阵 |
| `rhdf5` | 读 `.h5ad`（胎肾） |
| `ggplot2` + `ggforce` + `patchwork` + `scales` | 点图、括号、多面板拼图 |
| `magick` / `pdftools` | 裁原图、渲染 PDF 页面 |
| `pymupdf` / `pillow` / `numpy` | 从 PDF 提取原图、删列拼接、调色重绘 |

---

## 4. 数据准备

三张图共用 **E-MTAB-10552**（人卵黄囊 10x + 胎肝 10x）；Fig 2D 另需若干外部数据集。

### 4.1 数据根目录

Fig 2D / Fig 3 的脚本默认读两个数据根目录（**需按本机修改**，见 §7）：

```
C/data/processed/      # 处理后矩阵（mtx + features + metadata）
C/data_out/            # Seurat 对象（.rds）
```

### 4.2 数据集清单

| 数据集 | 来源 | 路径 |
|---|---|---|
| YS 10x（134,765 细胞） | E-MTAB-10552 | `data/processed/E-MTAB-10552_YS_10x/` |
| 胎肝 10x（49,376 细胞） | E-MTAB-10552 | `data/processed/E-MTAB-10552_FL_10x/` |
| 胎肝 CITE-seq | E-MTAB-11613 | `data/processed/E-MTAB-11613_FL_CITEseq/` |
| EL CITE-seq（26,216 细胞） | E-MTAB-11618 | `data/processed/E-MTAB-11618_EL_CITEseq/` |
| YS CITE-seq | E-MTAB-11549 | `data/processed/E-MTAB-11549_YS_CITEseq/` |
| YS Smart-seq2 | E-MTAB-10888 | `data/processed/E-MTAB-10888_YS_SS2/` |
| F158 / WE_YS | 本地 Seurat 对象 | `data_out/01_seurat_f158.rds` |
| **成人肝**（3,417 细胞） | **GSE124395**（Aizarani 2019） | 由 `11_adult_liver.R` **自动下载** |
| **胎肾基质**（8,345 → 1,110 fibroblast） | **Stewart 2019**，经 CELLxGENE | 手动下载 → `fig2d/data/fetal_kidney_stroma.h5ad` |
| **论文全文 PDF** | Science 官网 / 订阅 | 手动放置 → `fig5e/data/Science2023_Goh_add7564.pdf` |

---

## 5. 运行

### 5.1 Figure 2D

先 `setwd()` 到仓库根目录（脚本用相对路径）：

```r
setwd("path/to/Goh2023_reproduction")
```

```bash
Rscript fig2d/scripts/11_adult_liver.R     # ① 成人肝（联网下载 GSE124395）
Rscript fig2d/scripts/50_add_kidney.R      # ② 胎肾 h5ad 结构确认（只读）
Rscript fig2d/scripts/40_build_final_db.R  # ③ 8 套数据整合统计（最久）
Rscript fig2d/scripts/41_plot_fig2d_v2.R   # ④ 出图 + 原图同页对比
```

产出：`results/final_stats.rds`、`results/final_dotplot_stats.csv`、`figures/fig2d_v2_reproduction.pdf|png`、`figures/fig2d_v2_comparison.png`

### 5.2 Figure 3

```bash
Rscript fig3/scripts/03_fig3_reproduce.R
```

产出：`figures/fig3AB.pdf|png`、`figures/fig3A_dotplot_ys.pdf|png`、`figures/fig3B_early_vs_definitive.pdf|png`、`figures/fig3E_globin_switch.pdf|png`、`tables/t15–t17*.csv`

> 首次运行会扫描 YS 矩阵（2.76 亿非零值，约 2.6 分钟），结果缓存到 `data_out/fig3_*_scan.rds`，重跑直接复用。

### 5.3 Figure 5E

> 注意：这个脚本处理的是**论文 PDF 里的原图**（删列 + 调色 + 重绘文字），用于生成与复现图同口径对比的"原图 8 状态版"，不是从表达矩阵画图。

```bash
export PAPER_PDF=/path/to/Science2023_Goh_add7564.pdf
python fig5e/scripts/36_fig5e_paper_8states.py
# -> fig5e/figures/fig5e_paper_8states.png（657 × 385）
```

处理链：PDF 提取 Figure 5 位图 → 裁 5E 面板 → 按列删除 6 个无数据状态 → 调色（文字 #555、点 +26%）→ 擦除重绘文字（12px 斜体基因名 / 11px 竖排状态标签）→ 裁底。

---

## 6. 关键作图口径

### Fig 2D

| 项目 | 取值 |
|---|---|
| 归一化 | CP10K（target = 1e4）+ `log1p` |
| 颜色（组均值） | 基因内 min–max 缩放到 0–1 |
| 点大小（阳性率） | `range = c(0.5, 5.4)` |
| 色阶 | `grey94 → #67000d` |
| 基因 / 列 | 34 个目标基因 × 8 列 |
| 成人肝 / 胎肾细胞数 | 3,417 / 8,345 → 1,110 fibroblast |

### Fig 3

| 项目 | 取值 |
|---|---|
| Panel A | 7 行（6 个 YS 孕龄档 + FL HSC/MPP）× 23 基因，分 4 组（Canonical / Early / Definitive / Patterning） |
| Panel A 配色 | 基因内 min–max，量程只取 6 个 YS 行（FL 测序更深会把 YS 压白） |
| Panel B | 4 个孕龄档（~4 / ~5 / ~7 / ~8 PCW） |
| Panel E | HBZ / HBG1 / BCL11A / ZBTB7A，按孕龄档加权合并后作图 |
| 孕龄换算 | CS → PCW 映射表；胎肝 F## = 受精后 ## 天 |

### Fig 5E

| 项目 | 取值 |
|---|---|
| 色阶 | 论文原图色标实测：浅蓝灰 `#A5B2D9` → 米桃 `#EAC09F` → 暗酒红 `#7C162D` |
| 调色 | 文字 `#555`、点亮度 +26%、边框 `#696969` |
| 文字 | Arial 12px 斜体（基因名）/ 11px 竖排（状态）/ 17px 粗体（面板标签） |

---

## 7. 注意事项

1. **绝对路径需修改**。脚本保留了原始工作机的路径，首次运行前全局替换：
   ```r
   .libPaths("D:/R/library")                                  # -> 你的 R 库路径
   DATA_PROC <- "C:/Users/lenovo/Desktop/C/data/processed"     # -> 你的数据目录（fig2d/40）
   DATA_OUT  <- "C:/Users/lenovo/Desktop/C/data_out"
   ROOT      <- "C:/Users/lenovo/Desktop/C"                    # -> 你的数据根（fig3/03）
   ```
   `_crop_fig2d.R` 里还有论文 PDF 的路径。Python 脚本已改为仓库相对路径，只需指定 `PAPER_PDF`。

2. **大矩阵采用分块流式读取**。YS 10x 的 `mtx` 超过 1 GB，脚本按 500 万行一块读入，16 GB 内存可跑通。

3. **缓存机制**。Fig 2D 的 `40_build_final_db.R` 写 `results/cache_<ds>.rds`；Fig 3 写 `data_out/fig3_*_scan.rds`，并用细胞/基因清单做一致性校验——清单变了会自动重扫。

4. **YS 内胚层不能直接用官方注释**。校验脚本（`31` / `32`）证明官方 `LVL2/LVL3` 与表达矩阵的细胞顺序不一致，Fig 2D 最终改用 `ASGR1` / `SPINK1` 等 marker 打分 + 2-means 分出 Early / Late。

---

## 8. 已知限制

### Fig 2D
- 完整复现（8/8 列），但成人肝列与胎肾列的数据来自外部研究，与原文所用版本可能存在样本差异。

### Fig 3
- **Panel A 只做了 7 行**。原图 13 行中，Gastrulation / AGM / FBM HSC-MPP / iPSC / Def. iPSC 共 5 行来自另外几篇论文的独立数据集（Tyser 2021 等），不属 E-MTAB 数据，本地无法取得。
- **C / F / G / H 未覆盖**：依赖作者提供的 UMAP 坐标与 CellRank 转移矩阵；D 为胚胎照片。

### Fig 5E
- **状态 8/14**。论文 14 列中的 AGM / 脑 / 性腺及其 3 个 TREM2 亚组无本地数据，脚本从原图上删除这 6 列后保留 8 列。
- **本仓库只含"原图处理"脚本**。从表达矩阵重画 Fig 5E 的 R 脚本（点图复现）不在本仓库内。
- 色阶与文字重绘均按论文原图实测校准，属**版式复现**，不改变原图的数值含义。

### 通用
- 37 个排障脚本与 5 个 PPT/报告脚本未全文收录，清单见 `docs/`。
- Fig 3A / 3E 与 Fig 5E 的原始作图由组内其他成员完成。

---

## 9. 引用

```bibtex
@article{goh2023yolk,
  title   = {Yolk sac cell atlas reveals multiorgan functions during human early development},
  author  = {Goh, Issac and Botting, Rachel A. and Rose, Abigail and others},
  journal = {Science},
  volume  = {381},
  pages   = {eadd7564},
  year    = {2023},
  doi     = {10.1126/science.add7564}
}
```

数据来源：

- Aizarani N, et al. *A human liver cell atlas reveals heterogeneity and epithelial progenitors.* Nature 2019. **GSE124395**
- Stewart BJ, et al. *Single-cell transcriptomics of human fetal kidney.* Science 2019.
- 其余数据集见 §4.2 表格中的 ArrayExpress 编号（E-MTAB-*）

---

## 10. 上传到 GitHub / Gitee

本目录已是初始化好的 git 仓库（含提交记录），直接添加远端并推送：

**GitHub**

```bash
# 先在网页端创建一个空仓库（不要勾选 Initialize with README）
git remote add origin https://github.com/<your-name>/Goh2023_reproduction.git
git branch -M main
git push -u origin main
```

**Gitee**

```bash
git remote add origin https://gitee.com/<your-name>/Goh2023_reproduction.git
git branch -M master
git push -u origin master
```

> 数据文件与中间产物已被 `.gitignore` 排除，仓库体积约 300 KB。
