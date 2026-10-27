.libPaths("D:/R/library")
suppressMessages(library(data.table))
d <- "C:/Users/lenovo/Desktop/C/data/processed/E-MTAB-10552_YS_10x"
g <- fread(file.path(d, "YS_10x_genes.tsv.gz"), header = TRUE)
genes <- g[[1]]
cat("n_cells for housekeeping genes:\n")
hk <- c("MALAT1","EEF1A1","B2M","FTH1","FTL","ACTB","TMSB4X","RPL13A")
print(data.table(gene = hk,
                 n_cells = g$n_cells[match(hk, genes)],
                 total_counts = g$total_counts[match(hk, genes)]))

con <- gzfile(file.path(d, "YS_10x_counts.mtx.gz"), open = "rt")
dummy <- readLines(con, n = 4)
nnz <- integer(30330)
n_lines <- 0
repeat {
  lines <- readLines(con, n = 5000000, warn = FALSE)
  if (!length(lines)) break
  n_lines <- n_lines + length(lines)
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
              colClasses = list(integer = c("V1","V2")), showProgress = FALSE,
              select = c("V1","V2"))
  tt <- tabulate(dt$V2, nbins = 30330)
  nnz <- nnz + tt
  cat("scanned", n_lines, "lines; max cell so far", max(dt$V1), "\n")
  flush.console()
  if (n_lines >= 30000000) break
}
close(con)
## 扫描的细胞数：
n_cells_scanned <- max(dt$V1)
top <- order(-nnz)[1:12]
cat("scanned cells ~", n_cells_scanned, "\n")
cat("top columns by nnz:\n")
print(data.table(col = top, gene_at_col = genes[top], nnz = nnz[top],
                 frac = round(nnz[top] / n_cells_scanned, 3)))
