.libPaths("D:/R/library")
suppressMessages({ library(data.table) })

OUT_DIR  <- "fig2d/results"
DATA_DIR <- "fig2d/data/adult_liver"
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)

url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395_Normalhumanlivercellatlasdata.txt.gz"
gz  <- file.path(DATA_DIR, "GSE124395_liver_atlas.txt.gz")

options(timeout = 7200)
if (!file.exists(gz) || file.info(gz)$size < 1e6) {
  cat("[adult] downloading ", url, "\n")
  st <- system2("curl.exe", c("-L", "--ssl-no-revoke", "--retry", "5", "--retry-delay", "5",
                              "-sS", "-o", shQuote(gz), shQuote(url)))
  if (st != 0 || !file.exists(gz) || file.info(gz)$size < 1e6) {
    cat("[adult] download failed (curl status ", st, ")\n"); quit(save = "no", status = 1)
  }
}
cat("[adult] file size: ", file.info(gz)$size, " bytes\n")

fig2d_genes <- c("SERPINA1","APOA2","APOC3","SPINK1","AHSG","AFP","APOA1","ALB","TTR",
                 "F2","F5","F10","F3","F7","F12","F11","F8","F9",
                 "FGA","FGB","FGG","F13A1","F13B",
                 "SERPINC1","PROC","PROS1","SERPINF2","PLAT","PLAU",
                 "SERPINE1","SERPINE2","EPO","EGFR","THPO")
hepato_panel <- c("ALB","AFP","TTR","APOA1","APOA2","APOC3","ASGR1","AHSG",
                  "FGA","FGB","GC","SERPINA1","AMBP","CYP3A7","FETUB")
need <- unique(c(fig2d_genes, hepato_panel))

con <- gzfile(gz, open = "rt")
hdr <- readLines(con, n = 1)
cells <- strsplit(hdr, "\t")[[1]]
cells <- gsub('"', "", cells)
n_c <- length(cells) - 1
cat("[adult] cells: ", n_c, "\n")

keep_vals <- list(); gene_seen <- character(); libsize <- numeric(n_c)
repeat {
  lines <- readLines(con, n = 2000, warn = FALSE)
  if (!length(lines)) break
  for (ln in lines) {
    f <- strsplit(ln, "\t")[[1]]
    gname <- gsub('"', "", f[1])
    v <- as.numeric(gsub('"', "", f[-1]))
    if (length(v) != n_c) v <- as.numeric(gsub('"', "", f[-1])[1:n_c])
    libsize <- libsize + v
    if (gname %in% need) { keep_vals[[gname]] <- v; gene_seen <- c(gene_seen, gname) }
  }
}
close(con)

cat("[adult] target genes found: ", length(intersect(need, gene_seen)), "/",
    length(need), "\n")
mat <- do.call(rbind, keep_vals)
rownames(mat) <- names(keep_vals)
mat[is.na(mat)] <- 0

## CP10K-log
norm <- log1p(sweep(mat, 2, libsize, "/") * 1e4)

## 2-means 判定肝细胞（与 EL 一致的处理）
hp <- intersect(hepato_panel, rownames(norm))
score <- colMeans(norm[hp, , drop = FALSE])
km <- kmeans(scale(score), centers = 2, nstart = 20)
hep_cl <- names(which.max(tapply(score, km$cluster, mean)))
is_hep <- km$cluster == hep_cl
alb <- mat["ALB", ] > 0
cat("[adult] hepatocytes: ", sum(is_hep), "/", n_c,
    " | ALB+ in hep = ", round(100 * mean(alb[is_hep]), 1), "%",
    " | in others = ", round(100 * mean(alb[!is_hep]), 1), "%\n")
cat("[adult] ASGR1-sorted fraction in hep cluster: ",
    round(100 * mean(grepl("^ASGR1_CD45", cells[-1][is_hep])), 1), "%\n")

sel <- norm[, is_hep, drop = FALSE]
adult_means <- rowMeans(sel)
adult_pct   <- rowSums(mat[, is_hep, drop = FALSE] > 0) / sum(is_hep) * 100
idx <- match(fig2d_genes, names(adult_means))
adult_means <- adult_means[idx]; adult_pct <- adult_pct[idx]
names(adult_means) <- fig2d_genes; names(adult_pct) <- fig2d_genes

saveRDS(adult_means, file.path(OUT_DIR, "adult_hepatocyte_stats.rds"))
saveRDS(adult_pct,   file.path(OUT_DIR, "adult_hepatocyte_pct.rds"))
cat("[adult] saved adult_hepatocyte_stats.rds\n")
print(round(adult_means[order(-adult_means)], 2))
