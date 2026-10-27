.libPaths("D:/R/library")
suppressMessages({library(data.table)})
ROOT <- "C:/Users/lenovo/Desktop/C/data/processed"
MTX  <- file.path(ROOT, "E-MTAB-10552_YS_10x/YS_10x_counts.mtx.gz")
GEN  <- file.path(ROOT, "E-MTAB-10552_YS_10x/YS_10x_genes.tsv.gz")
OUT  <- "fig2d/results/ys_colprofile.rds"

## ---- 基因（第 1 列为 symbol，首行为表头）-------------------------------------
gv  <- fread(GEN, header = TRUE, sep = "\t")
sym <- trimws(gv[[1]])
cat("基因数:", length(sym), " 示例:", paste(head(sym, 3), collapse = ", "), "\n")
tg  <- c("HBG1","HBG2","HBB","HBZ","HBE1","HBD","HBM","HBQ1",
         "ALB","AFP","TTR","SPINK1","ASGR1","PTPRC","EPO","APOA1")
tgi <- match(tg, sym); names(tgi) <- tg
cat("命中:", paste(tg[!is.na(tgi)], collapse = ", "),
    " | 缺失:", paste(tg[is.na(tgi)], collapse = ", "), "\n")
tgi <- tgi[!is.na(tgi)]

## ---- 流式扫描（自动识别方向：cells x genes 还是 genes x cells）--------------
con  <- gzfile(MTX, open = "rt")
hd   <- readLines(con, n = 3)
dm   <- as.numeric(trimws(strsplit(hd[3], "\\s+")[[1]]))
nr <- dm[1]; nc <- dm[2]; nnz <- dm[3]
cat(sprintf("矩阵: %d x %d, nnz=%s\n", nr, nc, format(nnz, big.mark = ",")))
## 方向判断：行数 == 细胞数(134765) → cells x genes；行数 == 基因数 → genes x cells
if (nr == 134765L && nc == length(sym)) {
  ccol <- "V1"; gcol <- "V2"; cellfirst <- TRUE
} else if (nc == 134765L && nr == length(sym)) {
  ccol <- "V2"; gcol <- "V1"; cellfirst <- FALSE
} else stop("维度无法判定方向")
cat("方向:", if (cellfirst) "cells x genes（V1=细胞, V2=基因）" else "genes x cells（V1=基因, V2=细胞）", "\n")
NCELL <- 134765L

lib  <- numeric(NCELL)
prof <- matrix(0, nrow = length(tgi), ncol = NCELL, dimnames = list(names(tgi), NULL))
CH <- 5e6; k <- 0; nread <- 0; t0 <- Sys.time()
repeat {
  ln <- readLines(con, n = CH, warn = FALSE)
  if (!length(ln)) break
  nread <- nread + length(ln)
  d <- fread(text = ln, sep = " ", header = FALSE, showProgress = FALSE)
  setnames(d, c("V1", "V2", "V3"))
  cc <- d[[ccol]]; gg <- d[[gcol]]; vv <- d[[3]]
  rm(d, ln)
  a <- data.table(c = cc, v = vv)[, .(s = sum(v)), by = c]
  lib[a$c] <- lib[a$c] + a$s
  rm(a)
  sel <- gg %in% tgi
  if (any(sel)) {
    b <- data.table(gi = gg[sel], c = cc[sel], v = vv[sel])[, .(s = sum(v)), by = .(gi, c)]
    prof[cbind(match(b$gi, tgi), b$c)] <- prof[cbind(match(b$gi, tgi), b$c)] + b$s
    rm(b)
  }
  rm(cc, gg, vv, sel); k <- k + 1
  if (k %% 10 == 0) cat(sprintf("  %d 块, 已读 %s 行 / %.0f 秒\n", k,
      format(nread, big.mark = ","), difftime(Sys.time(), t0, units = "secs")))
}
close(con)
cat("扫描完成:", k, "块, 用时", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "分钟\n")
cat("sanity — 每细胞 UMI 中位数:", median(lib), " 非零细胞:", sum(lib > 0), "/", NCELL, "\n")
stopifnot(median(lib) > 100)   # 防止再次把行列弄反
saveRDS(list(lib = lib, prof = prof, nc = NCELL, nr = nr, nnz = nnz), OUT)

## ---- 检验 2：基因检出数 vs var$n_cells ---------------------------------------
det <- rowSums(prof > 0)
cmp <- data.table(gene = names(tgi), 矩阵检出细胞数 = det,
                  var表n_cells = gv[[2]][tgi])
cmp[, 比值 := round(矩阵检出细胞数 / var表n_cells, 3)]
cat("\n=== 基因检出细胞数 vs 注释表 n_cells ===\n"); print(cmp)
cat("比值中位数:", median(cmp$比值), " (若矩阵是同一批细胞的子集，应≈134765/169494=",
    round(134765 / 169494, 3), ")\n")

## ---- 检验 1：分块指纹比对 -----------------------------------------------------
md <- fread(file.path(ROOT, "E-MTAB-10552_YS_10x/YS_10x_core_metadata.csv"))
hb <- colSums(prof[intersect(c("HBG1","HBG2","HBB","HBZ","HBE1","HBD","HBM","HBQ1"),
                             rownames(prof)), , drop = FALSE])
hbf  <- ifelse(lib > 0, hb / lib, 0)
alb  <- prof["ALB", ] / pmax(lib, 1)
ery  <- hbf > 0.25
cat("\n矩阵红系样细胞比例:", sprintf("%.1f%%", 100 * mean(ery)),
    " | metadata ERYTHROID 比例:", sprintf("%.1f%%", 100 * mean(md$LVL1 == "ERYTHROID")), "\n")

bl <- 2000; nb <- floor(NCELL / bl)
grp <- rep(1:nb, each = bl)
m_ery <- as.vector(tapply(ery[1:(nb * bl)], grp, mean))
r_ery <- as.vector(tapply(md$LVL1[1:(nb * bl)] == "ERYTHROID", grp, mean))
m_alb <- as.vector(tapply(alb[1:(nb * bl)], grp, mean))
cat("\n=== 分块比对（每块 2000 细胞）===\n")
cmp2 <- data.table(块 = 1:nb, 起始行 = (0:(nb - 1)) * bl + 1,
                   矩阵红系比 = round(m_ery, 3), 注释红系比 = round(r_ery, 3),
                   矩阵ALB占比 = round(m_alb, 4))
print(cmp2[seq(1, nb, by = max(1, floor(nb / 34)))])
cat("\n相关性(红系):", sprintf("%.3f", cor(m_ery, r_ery)),
    "  相关性(ALB vs 注释红系比):", sprintf("%.3f", cor(m_alb, r_ery)), "\n")
cat("判定：相关性 > 0.6 支持“第 i 列 = 第 i 行”假设；≈0 说明矩阵与注释行序无关\n")

## 若对齐成立，输出各 stage/LVL3 在矩阵内的分布
cat("\n=== 假设成立时：注释 (前 134,765 行) 的可用构成 ===\n")
md2 <- md[1:NCELL]
print(md2[, .N, by = .(orig.dataset, stage)][order(orig.dataset, stage)][1:20])
cat("\nENDODERM 细胞数(LVL2):", md2[LVL2 == "ENDODERM", .N],
    " | (LVL3):", md2[LVL3 == "ENDODERM", .N], "\n")
