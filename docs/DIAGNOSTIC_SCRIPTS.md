# 诊断脚本清单（未全文收录）

以下脚本是排障过程中一次性使用的小工具，不影响最终出图结果，
原文档只列清单（说明取自脚本首行注释）。需要时可在原始工作区的 `fig2d/scripts/` 找到全文。

| 脚本 | 行数 | 大小 (KB) | 说明 |
|---|---|---|---|
| `_find_offset.R` | 59 | 2.6 | 定位矩阵列与注释行之 |
| `30_inspect_desktop_data.R` | 55 | 3.0 | ============================== |
| `_verify_alignment_blocks.R` | 53 | 2.6 | 独立验证：注释表的样 |
| `_diag_hep.R` | 48 | 2.3 | FL_10x：官方注释为真值 |
| `_fine_shift.R` | 48 | 2.2 | 去掉粗分块趋势后的精 |
| `_read_team_pptx.R` | 37 | 1.6 | — |
| `_local_shift.R` | 36 | 1.6 | 局部漂移检测：对每个 |
| `_scrublet_test.R` | 34 | 1.5 | 单细胞级决定性检验： |
| `_inspect_data2.R` | 33 | 1.4 | 检查数据（第二版）： |
| `_inspect_data.R` | 31 | 1.5 | 检查各数据文件的元数 |
| `_find_cxg.R` | 30 | 1.5 | — |
| `_check_profile.R` | 28 | 1.6 | — |
| `_check_bc2.R` | 24 | 1.2 | 尝试: metadata V1 -> 去掉 |
| `33_peek_seurat.R` | 24 | 0.9 | — |
| `_check_bc3.R` | 22 | 1.0 | 行块结构：每 5000 行的 |
| `_test_fpar.R` | 22 | 1.1 | — |
| `_assets.R` | 21 | 1.1 | — |
| `_inspect_h5ad.R` | 21 | 0.9 | 分类列的 categories |
| `_render_pages.R` | 21 | 0.7 | 辅助脚本：把 PDF 页面 |
| `_test_order.R` | 21 | 0.8 | — |
| `_diagnose_encoding.R` | 19 | 0.7 | — |
| `_counts.R` | 18 | 0.8 | — |
| `_diag_lib.R` | 18 | 0.7 | — |
| `_extract_pdf.R` | 16 | 0.6 | 首页预览 |
| `_get_refs.R` | 15 | 0.6 | — |
| `_grep_paper.R` | 15 | 0.6 | — |
| `_check_nums.R` | 14 | 0.6 | — |
| `_explore_endoderm.R` | 14 | 0.7 | 分析 YS 内胚层细胞的 |
| `_install_pptx_pkgs.R` | 13 | 0.6 | — |
| `_img_info.R` | 11 | 0.4 | — |
| `_check_bc.R` | 10 | 0.5 | — |
| `_pdf_read.R` | 10 | 0.3 | — |
| `_probe_h5.R` | 10 | 0.5 | — |
| `_dump_nb.R` | 9 | 0.4 | — |
| `_sum.R` | 9 | 0.4 | — |
| `_verify_figures.R` | 9 | 0.4 | — |
| `_summarize_stats.R` | 8 | 0.3 | — |
