.libPaths("D:/R/library")
suppressMessages(library(data.table))
d <- "C:/Users/lenovo/Desktop/C/data/processed/E-MTAB-10552_YS_10x"
g <- fread(file.path(d, "YS_10x_genes.tsv.gz"), header = TRUE)
genes <- g[[1]]
targets <- c("PECAM1","VWF","CDH5","CD34","SPINK2","HBA1","HBB","ALB","TTR","PTPRC","LYZ")
ti <- match(targets, genes)
cat("gene indices:\n"); print(data.table(gene = targets, idx = ti))
m <- fread(file.path(d, "YS_10x_core_metadata.csv"))
cat("row1:", m$V1[1], "|", m$LVL3[1], "\n")
cat("row2:", m$V1[2], "|", m$LVL3[2], "\n")

## 扫描矩阵前 ~8000 行对应的条目（cell 1 与 cell 2）
con <- gzfile(file.path(d, "YS_10x_counts.mtx.gz"), open = "rt")
dummy <- readLines(con, n = 4)   # 4 行头
vals <- list()
repeat {
  lines <- readLines(con, n = 200000, warn = FALSE)
  if (!length(lines)) break
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
              colClasses = list(integer = c("V1","V2"), numeric = "V3"),
              showProgress = FALSE)
  hit <- dt[V1 <= 2 & V2 %in% ti]
  if (nrow(hit)) vals[[length(vals) + 1]] <- hit
  if (max(dt$V1) > 2) break
}
close(con)
res <- rbindlist(vals)
setorder(res, V1, V2)
res[, gene := targets[match(V2, ti)]]
cat("== matrix rows 1-2 expression ==\n")
print(dcast(res, V2 + gene ~ V1, value.var = "V3", fill = 0))
