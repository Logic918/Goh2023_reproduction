.libPaths("D:/R/library")
suppressMessages({ library(rhdf5); library(Matrix); library(data.table) })

H5 <- "fig2d/data/fetal_kidney_stroma.h5ad"

## CELLxGENE 把分类列存成 categories + codes
read_cat <- function(path) {
  o <- tryCatch(h5read(H5, path), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  if (is.list(o) && !is.null(o$categories))
    as.character(o$categories)[as.integer(o$codes) + 1L]
  else as.character(o)
}

## ---- 1. 注释 ----------------------------------------------------------------
obs <- list()
for (k in c("cell_type", "author_cell_type", "tissue", "development_stage",
            "donor_id", "compartment", "sex", "assay", "cell_state")) {
  v <- read_cat(paste0("obs/", k))
  obs[[k]] <- v
  cat(sprintf("\n== %s (%s 个唯一值) ==\n", k, length(unique(v))))
  print(head(sort(table(v), decreasing = TRUE), 12))
}
cat("\n细胞数:", length(obs$cell_type), "\n")

## ---- 2. 基因名 --------------------------------------------------------------
gn <- read_cat("var/feature_name")
gn_raw <- read_cat("raw/var/feature_name")
cat("\n基因数: X 层", length(gn), " raw 层", length(gn_raw), "\n")
cat("前 6 个:", paste(head(gn_raw, 6), collapse = ", "), "\n")

## ---- 3. 稀疏矩阵结构 --------------------------------------------------------
for (p in c("X", "raw/X")) {
  sh <- tryCatch(h5read(H5, paste0(p, "/shape")), error = function(e) NULL)
  ip <- tryCatch(h5read(H5, paste0(p, "/indptr")), error = function(e) NULL)
  da <- tryCatch(h5read(H5, paste0(p, "/data")), error = function(e) NULL)
  if (is.null(sh)) { cat("\n", p, ": 无数据\n", sep = ""); next }
  cat(sprintf("\n%s: shape=%s  indptr 长度=%s  data 长度=%s  非零占比 %.1f%%\n",
              p, paste(sh, collapse = " x "), length(ip), length(da),
              100 * length(da) / prod(sh)))
  cat("  data 前 8 个:", paste(head(round(as.numeric(da), 3), 8), collapse = ", "), "\n")
}
