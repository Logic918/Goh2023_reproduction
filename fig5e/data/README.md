# 数据目录

## 需要手动放入

| 文件 | 来源 | 说明 |
|---|---|---|
| `Science2023_Goh_add7564.pdf` | Science 官网（需订阅） | 论文全文 PDF，脚本从中提取内嵌的 Figure 5 原始位图（1600×2100） |

也可用环境变量指定路径，不必放进本目录：

```bash
export PAPER_PDF=/path/to/Science2023_Goh_add7564.pdf
python fig5e/scripts/36_fig5e_paper_8states.py
```

## 输出

脚本产出写到 `fig5e/figures/fig5e_paper_8states.png`（657×385）。
首次运行会在仓库根生成 `.cache/fig5_original.png`（提取出的 Figure 5 位图），已被 `.gitignore` 忽略。
