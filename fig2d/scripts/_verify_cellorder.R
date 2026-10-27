.libPaths("D:/R/library")
suppressMessages(library(data.table))
d <- "C:/Users/lenovo/Desktop/C/data/processed/E-MTAB-10552_YS_10x"
g <- fread(file.path(d, "YS_10x_genes.tsv.gz"), header = TRUE)
genes <- g[[1]]
alb <- match("ALB", genes); hba <- match("HBA1", genes); ptc <- match("PTPRC", genes)
con <- gzfile(file.path(d, "YS_10x_counts.mtx.gz"), open = "rt")
dummy <- readLines(con, n = 4)
alb_flag <- logical(134765); hba_flag <- logical(134765)
n_lines <- 0
repeat {
  lines <- readLines(con, n = 5000000, warn = FALSE)
  if (!length(lines)) break
  n_lines <- n_lines + length(lines)
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
              colClasses = list(integer = c("V1","V2")), showProgress = FALSE,
              select = c("V1","V2"))
  a <- dt[V2 == alb, V1]; if (length(a)) alb_flag[a] <- TRUE
  h <- dt[V2 == hba, V1]; if (length(h)) hba_flag[h] <- TRUE
  if (n_lines >= 40000000) break
}
close(con)
n_cells <- max(which(alb_flag | hba_flag), 1)
cat("scanned up to cell row ~", max(dt$V1), "\n")
blk <- ceiling(seq_along(alb_flag) / 1000)
agg <- data.table(blk = blk[1:max(dt$V1)],
                  alb = alb_flag[1:max(dt$V1)],
                  hba = hba_flag[1:max(dt$V1)])
s <- agg[, .(alb_pct = 100 * mean(alb), hba_pct = 100 * mean(hba)), by = blk]
cat("ALB/HBA1 percentage per 1000-cell block (metadata order: Popescu=liver is rows 10671-23782):\n")
print(s)
