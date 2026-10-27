.libPaths("D:/R/library")
suppressMessages({library(rhdf5); library(data.table)})
ROOT <- "C:/Users/lenovo/Desktop/C/data/processed"
h5  <- file.path(ROOT, "E-MTAB-11673_YS_F158/17_04_24_YolkSacRaw_F158_WE_annots.h5ad")

## ---- 1. h5ad obs ------------------------------------------------------------
get_obs <- function(k) {
  p <- paste0("/obs/", k)
  d <- h5read(h5, p)
  if (is.list(d)) {
    if (!is.null(d$categories)) d$categories[d$codes + 1L] else unlist(d)
  } else as.character(d)
}
idx   <- get_obs("_index")
stage <- get_obs("stage"); fid <- get_obs("fetal.ids")
lvl1  <- get_obs("LVL1");  lvl3 <- get_obs("LVL3")
cat("h5ad obs 细胞数:", length(idx), "\n")
cat("fetal.ids:", paste(unique(fid), collapse = ", "), "\n")
cat("stage 分布:\n"); print(table(stage))
cat("LVL1 分布:\n"); print(table(lvl1))
cat("LVL3 top12:\n"); print(head(sort(table(lvl3), decreasing = TRUE), 12))
cat("barcode 示例:", head(idx, 3), "\n")

## ---- 2. YS_10x barcodes 是否包含 h5ad 的 barcode -----------------------------
bc <- fread(file.path(ROOT, "E-MTAB-10552_YS_10x/YS_10x_barcodes.tsv.gz"),
            header = FALSE)[[1]]
cat("\nYS_10x_barcodes 行数:", length(bc), "  示例:", head(bc, 3), "\n")
m_all   <- match(idx, bc)
cat("h5ad barcode 命中 barcodes 文件的比例:",
    sprintf("%.1f%%", 100 * mean(!is.na(m_all))), "\n")
mi <- m_all[!is.na(m_all)]
cat("命中位置: min=", min(mi), " max=", max(mi),
    " 是否单调递增(与 h5ad 行序一致):",
    if (all(diff(mi[order(match(idx[!is.na(m_all)], idx))]) > 0)) "是" else "否", "\n")
cat("命中位置落在 1..134765 的比例:", sprintf("%.1f%%", 100 * mean(mi <= 134765)), "\n")
cat("前 10 个命中位置:", head(mi, 10), "\n")
cat("相邻位置差的前 10 个:", head(diff(mi), 10), "\n")

## ---- 3. metadata 的 orig.dataset / sequencing.type --------------------------
md <- fread(file.path(ROOT, "E-MTAB-10552_YS_10x/YS_10x_core_metadata.csv"))
cat("\norig.dataset 分布:\n"); print(md[, .N, by = orig.dataset][order(-N)])
cat("\nsequencing.type 分布:\n"); print(md[, .N, by = sequencing.type][order(-N)])
cat("\n各 fetal.ids 的 10x 细胞数（前 25）:\n")
print(md[sequencing.type == "10x", .N, by = .(fetal.ids, stage)][order(-N)][1:25])
cat("\n10x 细胞总数:", md[sequencing.type == "10x", .N], "\n")

## ---- 4. 若 h5ad 命中位置集中在某区间，说明 barcodes 文件 = 矩阵列顺序 --------
if (!is.na(m_all[1])) {
  cat("\n=== 结论提示 ===\n")
  cat("若命中位置单调递增且落在 1..134765 内 → barcodes/metadata 的前 N 行即矩阵列顺序\n")
  cat("可据此直接用官方 LVL3/stage 注释复现 YS Endoderm(按 CS 分期) 与 YS Gastrulation\n")
}
