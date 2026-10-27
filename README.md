# Goh et al. (2023) Figure 2D 复现

> 人卵黄囊细胞图谱 · **Figure 2D** 点图的完整复现代码
> 34 个 marker 基因 × 8 个细胞群，含成人肝对照列与胎肾基质列

复现对象：Goh I, Botting RA, Rose A, et al. *Yolk sac cell atlas reveals multiorgan functions during human early development.* **Science** 2023;381:eadd7564. Figure 2D。

---

## 1. 复现目标

论文 Figure 2D 是一张跨组织、跨发育阶段的点图：

- **行**：34 个目标 marker 基因（红系、巨噬、内胚层、基质等谱系）
- **列**：8 个细胞群（YS 早期/晚期内胚层、胎肝 10x / CITE-seq、EL CITE-seq、YS CITE-seq、YS Smart-seq2、成人肝、胎肾基质）
- **颜色**：组均表达量（基因内 min–max 缩放到 0–1）
- **点大小**：阳性率（表达该基因的细胞比例）
- 绿色椭圆：叠加小鼠敲除表型信息

复现图与原图同页对比输出，便于逐像素核对版式。

---

## 2. 仓库结构

```
Goh2023_Fig2D_reproduction/
├── README.md
├── fig2d/
│   ├── scripts/                    # 全部 R 脚本（13 个全文收录）
│   │   ├── 11_adult_liver.R        # ① 成人肝对照列
│   │   ├── 50_add_kidney.R         # ② 胎肾 h5ad 结构确认
│   │   ├── 40_build_final_db.R     # ③ 核心：8 套数据整合统计
│   │   ├── 41_plot_fig2d_v2.R      # ④ 最终出图
│   │   ├── 10_fig2d_pipeline.R     # 第一版主流水线（备查）
│   │   ├── 12_plot_fig2d.R         # 第一版点图版式（备查）
│   │   ├── 31_check_alignment.R    # 校验：注释与矩阵对齐
│   │   ├── 32_test_ys_alignment.R  # 校验：YS 内胚层注释
│   │   ├── 34_check_f158.R         # 校验：F158/WE_YS 取数方案
│   │   ├── _verify_geneorder.R     # 校验：基因顺序
│   │   ├── _verify_cellorder.R     # 校验：细胞顺序
│   │   ├── _verify_mapping.R       # 校验：映射关系
│   │   └── _crop_fig2d.R           # 从论文 PDF 裁原图
│   ├── data/                       # 原始数据（.gitignore 忽略，见 data/README.md）
│   ├── results/                    # 中间产物与统计库
│   └── figures/                    # 出图结果
├── docs/
│   ├── DIAGNOSTIC_SCRIPTS.md       # 37 个排障脚本清单（未全文收录）
│   └── PIPELINE_SCRIPTS.md         # PPT / 报告生成脚本清单
└── .gitignore
```

---

## 3. 环境依赖

**R 4.6.1**（本机实测版本）

```r
# CRAN
install.packages(c("data.table", "Matrix", "ggplot2", "ggforce", "magick", "pdftools"))
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")   # 读 CELLxGENE 的 h5ad
```

| 包 | 用途 |
|---|---|
| `data.table` | 大矩阵分块流式读取、聚合 |
| `Matrix` | 稀疏矩阵 |
| `rhdf5` | 读 `.h5ad`（CELLxGENE 胎肾数据） |
| `ggplot2` + `ggforce` | 点图、椭圆、括号注释 |
| `magick` | 裁切原图、拼对比图 |
| `pdftools` | 渲染论文 PDF 页面 |

---

## 4. 数据准备

脚本默认读取两个数据根目录（**需按本机修改**，见 §7）：

```
C/data/processed/      # 处理后矩阵（mtx + features + metadata）
C/data_out/            # Seurat 对象（.rds）
fig2d/data/            # 本仓库内的小数据文件
```

| 数据集 | 来源 | 本地路径 |
|---|---|---|
| YS 10x（134,765 细胞） | E-MTAB-10552 | `data/processed/E-MTAB-10552_YS_10x/` |
| 胎肝 10x（49,376 细胞） | E-MTAB-10552 | `data/processed/E-MTAB-10552_FL_10x/` |
| 胎肝 CITE-seq | E-MTAB-11613 | `data/processed/E-MTAB-11613_FL_CITEseq/` |
| EL CITE-seq（26,216 细胞） | E-MTAB-11618 | `data/processed/E-MTAB-11618_EL_CITEseq/` |
| YS CITE-seq | E-MTAB-11549 | `data/processed/E-MTAB-11549_YS_CITEseq/` |
| YS Smart-seq2 | E-MTAB-10888 | `data/processed/E-MTAB-10888_YS_SS2/` |
| F158 / WE_YS | 本地 Seurat 对象 | `data_out/01_seurat_f158.rds` |
| **成人肝对照**（3,417 细胞） | **GSE124395**（Aizarani et al., *Nature* 2019） | 由 `11_adult_liver.R` **自动下载**至 `fig2d/data/adult_liver/` |
| **胎肾基质**（8,345 → 1,110 fibroblast） | **Stewart et al., *Science* 2019**，经 CELLxGENE | 需手动下载 → `fig2d/data/fetal_kidney_stroma.h5ad` |

> 成人肝与胎肾两组数据原论文未随稿发布，需单独获取。前者脚本自动下载，后者到 CELLxGENE 下载 `fetal_kidney_stroma.h5ad` 放入 `fig2d/data/`。

---

## 5. 运行

所有脚本用**相对仓库根目录**的路径，先 `setwd()` 到仓库根目录：

```r
setwd("path/to/Goh2023_Fig2D_reproduction")
```

```bash
# ① 成人肝对照列（需联网下载 GSE124395）
Rscript fig2d/scripts/11_adult_liver.R
#   -> results/adult_hepatocyte_stats.rds / _pct.rds

# ② 胎肾 h5ad 结构确认（只读，不落盘）
Rscript fig2d/scripts/50_add_kidney.R
#   -> 打印 obs 分类列分布与 X / raw 层稀疏结构

# ③ 整合 8 套数据，生成统计库（核心，耗时最长）
Rscript fig2d/scripts/40_build_final_db.R
#   -> results/final_stats.rds
#   -> results/final_dotplot_stats.csv
#   -> results/final_group_sizes.csv / database_index.csv
#   -> results/cache_<dataset>.rds（各数据集缓存，可复用加速）

# ④ 出图 + 与原图同页对比
Rscript fig2d/scripts/41_plot_fig2d_v2.R
#   -> figures/fig2d_v2_reproduction.pdf|png
#   -> figures/fig2d_v2_comparison.png
```

**备查脚本**（早期版本，可不跑）：

```bash
Rscript fig2d/scripts/10_fig2d_pipeline.R   # 第一版主流水线
Rscript fig2d/scripts/12_plot_fig2d.R       # 第一版点图版式
```

---

## 6. 关键作图口径（与原文对齐）

| 项目 | 取值 |
|---|---|
| 归一化 | CP10K（target = 1e4）+ `log1p` |
| 颜色（组均值） | 基因内 min–max 缩放到 0–1 |
| 点大小（阳性率） | `range = c(0.5, 5.4)` |
| 色阶 | `grey94 → #67000d` |
| 基因 / 列 | 34 个目标基因 × 8 列 |
| 成人肝细胞数 | 3,417（GSE124395） |
| 胎肾细胞数 | 8,345 → 1,110 fibroblast |
| 基因名 | 斜体 |
| 横轴标签 | 旋转 90° |
| 版式细节 | `coord_cartesian(clip = "off")` 保证括号与竖排标签不被裁 |

---

## 7. 注意事项

1. **绝对路径需修改**。脚本保留了原始工作机的路径，首次运行前请全局替换：
   ```r
   .libPaths("D:/R/library")                    # -> 你的 R 库路径
   DATA_PROC <- "C:/Users/lenovo/Desktop/C/data/processed"   # -> 你的数据目录
   DATA_OUT  <- "C:/Users/lenovo/Desktop/C/data_out"
   ```
   以及 `_crop_fig2d.R` 里论文 PDF 的路径。

2. **大矩阵采用分块流式读取**。YS 10x 的 `mtx` 超过 1 GB，脚本按 500 万行一块读入，内存占用可控（本机 16 GB 可跑通）。

3. **缓存机制**。`40_build_final_db.R` 会把每套数据的统计结果写入 `results/cache_<dataset>.rds`；重跑时若缓存存在会直接读取。改了某套数据的处理逻辑后，删掉对应 cache 再跑。

4. **YS 内胚层不能直接用官方注释**。校验脚本（`31` / `32`）证明官方 `LVL2/LVL3` 与表达矩阵的细胞顺序不一致，最终改用 `ASGR1` / `SPINK1` 等 marker 打分 + 2-means 分出 Early / Late——这是本次复现的关键决定。

---

## 8. 未收录内容

| 内容 | 说明 |
|---|---|
| 37 个排障脚本 | 一次性小工具，不影响出图。清单见 `docs/DIAGNOSTIC_SCRIPTS.md` |
| 5 个 PPT / 报告脚本 | `42` / `43` / `70` / `60` / `80`，用于生成汇报材料。清单见 `docs/PIPELINE_SCRIPTS.md` |
| Fig 3A / 3E / 5E 的代码 | 由组内其他成员完成，未收到脚本 |

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
- 其余数据集见上方表格中的 ArrayExpress 编号（E-MTAB-*）
