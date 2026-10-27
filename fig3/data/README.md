# 数据目录

Figure 3 的脚本读取仓库外的数据根目录（见主 README §4），本目录仅用于放置小文件。

```r
# fig3/scripts/03_fig3_reproduce.R 顶部的路径
ROOT   <- "C:/Users/lenovo/Desktop/C"          # 需改成本机路径
YS_DIR <- file.path(ROOT, "data/processed/E-MTAB-10552_YS_10x")
FL_DIR <- file.path(ROOT, "data/processed/E-MTAB-10552_FL_10x")
```

所需数据：

| 目录 | 来源 | 细胞数 |
|---|---|---|
| `E-MTAB-10552_YS_10x/` | E-MTAB-10552 | 134,765 |
| `E-MTAB-10552_FL_10x/` | E-MTAB-10552 | 49,376 |
