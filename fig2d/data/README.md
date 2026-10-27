# 数据目录

本目录存放体积较小的原始数据；大矩阵请放在仓库外的数据根目录（见主 README §4）。

## 需要手动放入的文件

| 文件 | 来源 | 说明 |
|---|---|---|
| `fetal_kidney_stroma.h5ad` | CELLxGENE（Stewart et al., Science 2019） | 胎肾基质子集，8,345 细胞；`50_add_kidney.R` 读取后取 1,110 个 fibroblast |

## 脚本自动下载

| 目录 | 来源 | 脚本 |
|---|---|---|
| `adult_liver/` | GEO **GSE124395**（Aizarani et al., Nature 2019） | `11_adult_liver.R` 自动下载（genes × cells 文本矩阵，12,621 细胞） |

## 不纳入版本控制

本目录下的数据文件被 `.gitignore` 忽略（`fig2d/data/*` 除本文件外），请自行获取后放入。
