.libPaths("D:/R/library")
suppressMessages({library(Seurat); library(rhdf5); library(data.table)})

## ---- A. Seurat 对象：是否有 raw counts --------------------------------------
o <- readRDS("C:/Users/lenovo/Desktop/C/data_out/01_seurat_f158.rds")
a <- o@assays$RNA
cat("layers:", paste(Layers(a), collapse = ", "), "\n")
for (ly in Layers(a)) {
  m <- LayerData(a, layer = ly)
  cat(sprintf("  layer %-10s dim %s x %s  nnz %s\n", ly, nrow(m), ncol(m),
              format(length(m@x), big.mark = ",")))
}
tg <- c("SERPINA1","APOA2","APOC3","SPINK1","AHSG","AFP","APOA1","ALB","TTR",
        "F2","F5","F10","F3","F7","F12","F11","F8","F9","FGA","FGB","FGG","F13A1",
        "F13B","SERPINC1","PROC","PROS1","SERPINF2","PLAT","PLAU","SERPINE1",
        "SERPINE2","EPO","EGFR","THPO","ASGR1")
cat("目标基因命中:", sum(tg %in% rownames(a)), "/", length(tg), "\n")
cat("缺失:", paste(setdiff(tg, rownames(a)), collapse = ", "), "\n")
cat("meta ENDODERM(LVL3):", sum(o@meta.data$LVL3 == "ENDODERM"),
    " LVL2:", sum(o@meta.data$LVL2 == "ENDODERM"), "\n")
ly_raw <- intersect(c("counts", "Count", "raw"), Layers(a))
if (length(ly_raw)) {
  cat("libsize 中位数(", ly_raw[1], "):",
      median(Matrix::colSums(LayerData(a, layer = ly_raw[1]))), "\n")
}
rm(o, a); invisible(gc())

## ---- B. h5ad：X 的规模（若需要未过滤的原始细胞）-----------------------------
h5 <- "C:/Users/lenovo/Desktop/C/data/processed/E-MTAB-11673_YS_F158/17_04_24_YolkSacRaw_F158_WE_annots.h5ad"
st <- h5ls(h5, recursive = 2)
print(st[grepl("^/X|^/obs/_index|^/var", st$group), ])
cat("\nh5ad X 维度信息:\n")
xgrp <- st[st$group == "/X", ]
print(xgrp)
if ("indptr" %in% xgrp$name) {
  ip <- h5read(h5, "/X/indptr")
  cat("n_obs(indptr-1):", length(ip) - 1, " nnz:", ip[length(ip)], "\n")
}
var_names <- try(h5read(h5, "/var/_index"), silent = TRUE)
if (!inherits(var_names, "try-error")) cat("var 基因数:", length(var_names), "\n")

## ---- C. F158 mtx barcodes 格式 ----------------------------------------------
bc <- fread("C:/Users/lenovo/Desktop/C/data/processed/E-MTAB-11673_YS_F158/F158_YS_barcodes.tsv.gz",
            header = FALSE)[[1]]
cat("\nF158_YS_barcodes:", length(bc), " 示例:", head(bc, 3), "\n")
