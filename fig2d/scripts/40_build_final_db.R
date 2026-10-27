.libPaths("D:/R/library")
suppressMessages({
  library(data.table); library(Matrix)
})
set.seed(1)

DATA_PROC <- "C:/Users/lenovo/Desktop/C/data/processed"
DATA_OUT  <- "C:/Users/lenovo/Desktop/C/data_out"
RES       <- "fig2d/results"
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## -----------------------------------------------------------------------------
## 1. 基因面板（与原图 Figure 2D 逐字一致）+ marker panel
## -----------------------------------------------------------------------------
fig2d_genes <- list(
  "Metabolic" = c("SERPINA1","APOA2","APOC3","SPINK1","AHSG","AFP",
                  "APOA1","ALB","TTR"),
  "Common pathway"             = c("F2","F5"),
  "Extrinsic pathway"          = c("F10","F3","F7"),
  "Intrinsic pathway"          = c("F12","F11","F8","F9"),
  "Fibrinogen & crosslinkning" = c("FGA","FGB","FGG","F13A1","F13B"),
  "Anticoagulation"            = c("SERPINC1","PROC","PROS1","SERPINF2"),
  "Fibrinolysis"               = c("PLAT","PLAU"),
  "Growth factors"             = c("SERPINE1","SERPINE2","EPO","EGFR","THPO")
)
row_group_big <- c(
  "Metabolic" = "Metabolic", "Common pathway" = "Haemostasis",
  "Extrinsic pathway" = "Haemostasis", "Intrinsic pathway" = "Haemostasis",
  "Fibrinogen & crosslinkning" = "Haemostasis", "Anticoagulation" = "Haemostasis",
  "Fibrinolysis" = "Haemostasis", "Growth factors" = "Growth factors")
green_solid  <- c("F2","F5","F10","F3")
green_hollow <- c("SERPINC1","PROC","PROS1","EPO")

early_markers <- c("PRSS3","GSTA2","LGALS3","APOA1","APOA2","APOC3","TTR","AFP",
                   "ALB","SPINK1","AHSG","FGA","FGB","FGG","F2","F5","F10","F3",
                   "SERPINC1","PROS1","EPO","THPO")
late_markers  <- c("FOS","JUN","JUNB","FOSB","EGR1","ATF3","HSPA1A","HSPA1B",
                   "HSPH1","DNAJB1","BTG2","DUSP1","ZFP36","IER2")
endo_pos <- c("ASGR1","SPINK1","TTR","APOA1","APOA2","APOC3","ALB","AFP",
              "SERPINA1","AHSG","FGA","FGB")
endo_neg <- c("PECAM1","VWF","CDH5","PTPRC","HBA1","HBB","LYZ","C1QA","COL1A1",
              "TAGLN","ACTA2","MSR1")
hepato_panel <- c("ALB","AFP","TTR","APOA1","APOA2","APOC3","ASGR1","AHSG",
                  "FGA","FGB","GC","SERPINA1","AMBP","CYP3A7","FETUB")
all_need <- unique(c(unlist(fig2d_genes), early_markers, late_markers,
                     endo_pos, endo_neg, hepato_panel))
panel34  <- unlist(fig2d_genes, use.names = FALSE)

## -----------------------------------------------------------------------------
## 2. 通用函数
## -----------------------------------------------------------------------------
read_feat <- function(f, header = FALSE, col = 1L) {
  d <- fread(f, header = header, sep = "\t", fill = TRUE,
             colClasses = "character", showProgress = FALSE)
  trimws(d[[col]])
}

## 分块流式提取 mtx.gz（自动识别 cells×genes / genes×cells）
extract_mtx <- function(mtx_gz, genes, keep_genes, n_cells, chunk_lines = 5e6L) {
  g_idx <- match(keep_genes, genes)
  keep_genes <- keep_genes[!is.na(g_idx)]; g_idx <- g_idx[!is.na(g_idx)]
  con  <- gzfile(mtx_gz, open = "rt")
  hd   <- readLines(con, n = 3)
  dline <- trimws(hd[grep("^[0-9]", hd)][1])
  dims <- as.numeric(strsplit(dline, "\\s+")[[1]])
  if (dims[1] == n_cells && dims[2] == length(genes)) {
    ccol <- "V1"; gcol <- "V2"
  } else if (dims[1] == length(genes) && dims[2] == n_cells) {
    ccol <- "V2"; gcol <- "V1"
  } else stop("dims ", dims[1], "x", dims[2], " 与 genes(", length(genes),
              ")/cells(", n_cells, ") 不匹配")

  ## ---- 稠密矩阵分支（如 YS_SS2：只有 "行 列" 两数，无 nnz）-------------------
  if (length(dims) == 2L) {
    log_msg("  稠密矩阵 ", dims[1], " x ", dims[2], "（无 nnz 行，按列主序读入）")
    vals <- scan(con, what = double(), quiet = TRUE)
    D <- matrix(vals, nrow = dims[1], ncol = dims[2])
    rm(vals)
    Dg <- if (ccol == "V1") D else t(D)          # 统一为 cells x genes
    rm(D)
    lib <- rowSums(Dg)
    M <- Matrix(t(Dg[, g_idx, drop = FALSE]), sparse = TRUE)
    dimnames(M) <- list(keep_genes, NULL)
    close(con)
    return(list(counts = M, libsize = lib))
  }
  nnz <- dims[3]
  log_msg("  矩阵 ", dims[1], " x ", dims[2], " nnz=", format(nnz, big.mark = ","),
          if (ccol == "V1") " (cells x genes)" else " (genes x cells)")
  lib <- numeric(n_cells)
  GI <- list(); CI <- list(); XV <- list(); k <- 0L; nread <- 0
  repeat {
    ln <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(ln)) break
    nread <- nread + length(ln)
    d <- fread(text = ln, sep = " ", header = FALSE, showProgress = FALSE)
    setnames(d, c("V1", "V2", "V3"))
    cc <- d[[ccol]]; gg <- d[[gcol]]; vv <- d[[3]]
    rm(d, ln)
    a <- data.table(c = cc, v = vv)[, .(s = sum(v)), by = c]
    lib[a$c] <- lib[a$c] + a$s
    rm(a)
    sel <- gg %in% g_idx
    if (any(sel)) {
      k <- k + 1L
      GI[[k]] <- gg[sel]; CI[[k]] <- cc[sel]; XV[[k]] <- vv[sel]
    }
    rm(cc, gg, vv, sel)
    if (nread >= nnz) break
  }
  close(con)
  M <- if (k) sparseMatrix(i = match(unlist(GI), g_idx), j = unlist(CI),
                           x = unlist(XV), dims = c(length(g_idx), n_cells),
                           dimnames = list(keep_genes, NULL))
       else sparseMatrix(i = integer(), j = integer(), x = numeric(),
                         dims = c(length(g_idx), n_cells),
                         dimnames = list(keep_genes, NULL))
  list(counts = M, libsize = lib)
}

cp10k_log <- function(M, libsize, target = 1e4) {
  sc <- target / pmax(libsize, 1)
  M2 <- t(t(M) * sc); M2@x <- log1p(M2@x); M2
}

cache_or_run <- function(cache_file, expr) {
  if (file.exists(cache_file)) { log_msg("  读缓存 ", basename(cache_file)); readRDS(cache_file) }
  else { v <- force(expr); saveRDS(v, cache_file); v }
}

## =============================================================================
## A. YS_10x（134,765 细胞）
## =============================================================================
log_msg("== A. YS_10x ==")
ys_cache <- file.path(RES, "cache_ys10x.rds")
ys <- cache_or_run(ys_cache, {
  ys_genes <- read_feat(file.path(DATA_PROC, "E-MTAB-10552_YS_10x/YS_10x_genes.tsv.gz"),
                        header = TRUE)
  keep <- intersect(all_need, ys_genes)
  log_msg("  YS 基因表缺失:", paste(setdiff(all_need, ys_genes), collapse = ", "))
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-10552_YS_10x/YS_10x_counts.mtx.gz"),
                   ys_genes, keep, n_cells = 134765)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize,
              genes = ys_genes, n_cells = 134765)
  rm(r); invisible(gc()); out
})
M_ys <- ys$M; n_ys <- ncol(M_ys)
log_msg("  YS_10x: ", n_ys, " 细胞, ", nrow(M_ys), " 基因")

## ---- 关键：官方注释能否按行对齐 ----------------------------------------------
ys_meta <- fread(file.path(DATA_PROC, "E-MTAB-10552_YS_10x/YS_10x_core_metadata.csv"))
prof_file <- file.path(RES, "ys_colprofile.rds")
aligned <- FALSE; align_cor <- NA; align_cell <- NA
if (file.exists(prof_file)) {
  pf <- readRDS(prof_file)
  hb <- colSums(pf$prof[intersect(c("HBG1","HBG2","HBB","HBZ","HBE1","HBD","HBM","HBQ1"),
                                  rownames(pf$prof)), , drop = FALSE])
  hbf <- ifelse(pf$lib > 0, hb / pf$lib, 0)
  ery <- as.numeric(hbf > 0.25)
  ref <- as.numeric(ys_meta$LVL1[seq_len(n_ys)] == "ERYTHROID")
  ## ① 粗分块（2000 细胞）层面：判断是否为同一批细胞、大块顺序是否一致
  bl <- 2000; nb <- floor(n_ys / bl); grp <- rep(1:nb, each = bl)
  m_ery <- as.vector(tapply(ery[1:(nb*bl)], grp, mean))
  r_ery <- as.vector(tapply(ref[1:(nb*bl)], grp, mean))
  align_cor <- cor(m_ery, r_ery)
  ## ② 单细胞层面：先去掉粗分块趋势（滚动均值 2001），再算残差相关
  ##    全局相关会被样本块结构主导，必须去趋势才能反映"块内单细胞顺序"
  rollmean <- function(v, k) {
    n <- length(v); h <- k %/% 2; cs <- cumsum(c(0, v))
    a <- pmax(1, seq_len(n) - h); b <- pmin(n, seq_len(n) + h)
    (cs[b + 1] - cs[a]) / (b - a + 1)
  }
  ery_d <- ery - rollmean(ery, 2001)
  ref_d <- ref - rollmean(ref, 2001)
  align_cell <- cor(ery_d, ref_d)
  ## 判据：粗分块一致 + 去趋势后单细胞一致，两者同时满足才使用官方注释
  aligned <- !is.na(align_cor) && align_cor > 0.6 && !is.na(align_cell) && align_cell > 0.3
  log_msg("  ① 粗分块红系比例相关性 = ", round(align_cor, 3),
          "（同一批细胞、样本块顺序一致）")
  log_msg("  ② 单细胞红系判定相关性 = ", round(align_cell, 3),
          "（<0.3 → 逐细胞注释不可转移）")
  log_msg("  → 官方注释", if (aligned) "可对齐（采用）" else "无法逐细胞对齐（改用 marker 法）")
} else log_msg("  未找到对齐检验结果 → 沿用 marker 法")

## ---- YS 内胚层：官方注释（若可用）或 marker 法 -------------------------------
if (aligned) {
  md_ys <- ys_meta[seq_len(n_ys)]
  is_endo <- md_ys$LVL2 == "ENDODERM" | md_ys$LVL3 == "ENDODERM"
  ys_stage <- md_ys$stage
  log_msg("  官方 ENDODERM 细胞: ", sum(is_endo), "/", n_ys)
} else {
  pos <- intersect(endo_pos, rownames(M_ys)); neg <- intersect(endo_neg, rownames(M_ys))
  pos_score <- colMeans(as.matrix(M_ys[pos, ])); neg_score <- colMeans(as.matrix(M_ys[neg, ]))
  cand <- which(M_ys["ASGR1", ] > 0 & M_ys["SPINK1", ] > 0)
  km_e <- kmeans(cbind(scale(pos_score[cand]), scale(-neg_score[cand])), 2, nstart = 20)
  endo_cl <- names(which.max(tapply(pos_score[cand] - neg_score[cand], km_e$cluster, mean)))
  is_endo <- logical(n_ys); is_endo[cand[km_e$cluster == endo_cl]] <- TRUE
  ys_stage <- NA_character_
  log_msg("  marker 法识别内胚层: ", sum(is_endo))
}

## 早/晚状态（论文：早期代谢/视黄酸活跃，晚期应激基因）
em <- as.matrix(M_ys[intersect(early_markers, rownames(M_ys)), is_endo, drop = FALSE])
lmk <- as.matrix(M_ys[intersect(late_markers, rownames(M_ys)), is_endo, drop = FALSE])
km_s <- kmeans(cbind(scale(colMeans(em)), scale(colMeans(lmk))), 2, nstart = 20)
late_cl <- names(which.max(tapply(colMeans(lmk) - colMeans(em), km_s$cluster, mean)))
ys_late <- logical(sum(is_endo)); ys_late[km_s$cluster == late_cl] <- TRUE

## =============================================================================
## B. F158 / WE_YS（data_out 的 Seurat 对象，官方注释，QC 后 22,462 细胞）
## =============================================================================
log_msg("== B. F158 / WE_YS (data_out) ==")
f158_cache <- file.path(RES, "cache_f158.rds")
f158 <- cache_or_run(f158_cache, {
  suppressMessages(library(Seurat))
  o <- readRDS(file.path(DATA_OUT, "01_seurat_f158.rds"))
  cnt <- LayerData(o@assays$RNA, layer = "counts")
  md  <- o@meta.data
  gn  <- rownames(cnt); keep <- intersect(all_need, gn)
  log_msg("  F158 缺失基因:", paste(setdiff(all_need, gn), collapse = ", "))
  Mc  <- cnt[keep, , drop = FALSE]
  lib <- Matrix::colSums(cnt)
  rm(o, cnt); invisible(gc())
  list(M = cp10k_log(Mc, lib), libsize = lib, meta = md,
       lvl1 = md$LVL1, lvl2 = md$LVL2, lvl3 = md$LVL3, stage = md$stage)
})
log_msg("  F158: ", ncol(f158$M), " 细胞, ENDODERM(LVL3) = ", sum(f158$lvl3 == "ENDODERM"))

## —— 用 F158 官方注释验证 marker 法（ASGR1+SPINK1+）的准确度
v_as <- as.vector(f158$M["ASGR1", ] > 0); v_sp <- as.vector(f158$M["SPINK1", ] > 0)
truth <- f158$lvl3 == "ENDODERM"
marker_call <- v_as & v_sp
tp <- sum(marker_call & truth); fp <- sum(marker_call & !truth); fn <- sum(!marker_call & truth)
marker_eval <- data.table(数据集 = "F158/WE_YS", 官方ENDODERM = sum(truth),
                          marker阳性 = sum(marker_call),
                          精确率 = round(tp / max(sum(marker_call), 1), 3),
                          召回率 = round(tp / max(sum(truth), 1), 3))
log_msg("  marker 法验证: 精确率 ", marker_eval$精确率, " 召回率 ", marker_eval$召回率)
fwrite(marker_eval, file.path(RES, "marker_validation.csv"))

## =============================================================================
## C. FL_10x（官方注释 Hepatocyte）
## =============================================================================
log_msg("== C. FL_10x ==")
fl_cache <- file.path(RES, "cache_fl10x.rds")
fl <- cache_or_run(fl_cache, {
  fl_meta  <- fread(file.path(DATA_PROC, "E-MTAB-10552_FL_10x/FL_10x_barcodes_and_metadata.tsv.gz"))
  fl_genes <- read_feat(file.path(DATA_PROC, "E-MTAB-10552_FL_10x/FL_10x_genes.tsv.gz"), header = TRUE)
  keep <- intersect(all_need, fl_genes)
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-10552_FL_10x/FL_10x_counts.mtx.gz"),
                   fl_genes, keep, n_cells = 49376)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize, meta = fl_meta)
  rm(r); invisible(gc()); out
})
fl_hep <- fl$meta$cell.labels == "Hepatocyte"
log_msg("  FL_10x hepatocytes: ", sum(fl_hep))

## =============================================================================
## D. FL_CITEseq（新增，63,693 细胞）
## =============================================================================
log_msg("== D. FL_CITEseq ==")
flc_cache <- file.path(RES, "cache_flcite.rds")
flc <- cache_or_run(flc_cache, {
  g <- read_feat(file.path(DATA_PROC, "E-MTAB-11613_FL_CITEseq/FL_CITEseq_RNA_features.tsv.gz"))
  keep <- intersect(all_need, g)
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-11613_FL_CITEseq/FL_CITEseq_RNA_matrix.mtx.gz"),
                   g, keep, n_cells = 63693)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize,
              counts = r$counts)
  rm(r); invisible(gc()); out
})
hp <- intersect(hepato_panel, rownames(flc$M))
hpc <- colMeans(as.matrix(flc$M[hp, , drop = FALSE]))
km <- kmeans(scale(hpc), 2, nstart = 20)
flc_hep <- km$cluster == names(which.max(tapply(hpc, km$cluster, mean)))
log_msg("  FL_CITEseq 肝细胞(marker+2-means): ", sum(flc_hep), "/", ncol(flc$M),
        " | ALB+ = ", round(100 * mean(flc$counts["ALB", flc_hep] > 0), 1), "%")

## =============================================================================
## E. EL_CITEseq（26,216 细胞）
## =============================================================================
log_msg("== E. EL_CITEseq ==")
el_cache <- file.path(RES, "cache_el.rds")
el <- cache_or_run(el_cache, {
  g <- read_feat(file.path(DATA_PROC, "E-MTAB-11618_EL_CITEseq/EL_CITEseq_RNA_features.tsv.gz"))
  keep <- intersect(all_need, g)
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-11618_EL_CITEseq/EL_CITEseq_RNA_matrix.mtx.gz"),
                   g, keep, n_cells = 26216)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize, counts = r$counts)
  rm(r); invisible(gc()); out
})
hp <- intersect(hepato_panel, rownames(el$M))
hpc <- colMeans(as.matrix(el$M[hp, , drop = FALSE]))
km <- kmeans(scale(hpc), 2, nstart = 20)
el_hep <- km$cluster == names(which.max(tapply(hpc, km$cluster, mean)))
log_msg("  EL 肝细胞: ", sum(el_hep), " | ALB+ = ",
        round(100 * mean(el$counts["ALB", el_hep] > 0), 1), "%")

## =============================================================================
## F. YS_CITEseq（新增，7,739 细胞；RNA + ADT 验证）
## =============================================================================
log_msg("== F. YS_CITEseq ==")
ysc_cache <- file.path(RES, "cache_yscite.rds")
ysc <- cache_or_run(ysc_cache, {
  g <- read_feat(file.path(DATA_PROC, "E-MTAB-11549_YS_CITEseq/YS_CITEseq_RNA_features.tsv.gz"))
  keep <- intersect(all_need, g)
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-11549_YS_CITEseq/YS_CITEseq_RNA_matrix.mtx.gz"),
                   g, keep, n_cells = 7739)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize, counts = r$counts)
  rm(r); invisible(gc()); out
})
pos <- intersect(endo_pos, rownames(ysc$M)); neg <- intersect(endo_neg, rownames(ysc$M))
ps <- colMeans(as.matrix(ysc$M[pos, , drop = FALSE])); ns <- colMeans(as.matrix(ysc$M[neg, , drop = FALSE]))
cand <- which(as.vector(ysc$M["ASGR1", ] > 0) & as.vector(ysc$M["SPINK1", ] > 0))
if (length(cand) >= 10) {
  km <- kmeans(cbind(scale(ps[cand]), scale(-ns[cand])), 2, nstart = 20)
  cl <- names(which.max(tapply(ps[cand] - ns[cand], km$cluster, mean)))
  ysc_endo <- logical(ncol(ysc$M)); ysc_endo[cand[km$cluster == cl]] <- TRUE
} else ysc_endo <- rep(FALSE, ncol(ysc$M))
log_msg("  YS_CITEseq 内胚层: ", sum(ysc_endo), "/", ncol(ysc$M))

## ADT（CD326/EpCAM）对内胚层的蛋白层面支持
adt_eval <- NULL
try({
  pf <- read_feat(file.path(DATA_PROC, "E-MTAB-11549_YS_CITEseq/YS_CITEseq_protein_features.tsv.gz"))
  ip <- match("CD326", pf)
  if (!is.na(ip)) {
    r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-11549_YS_CITEseq/YS_CITEseq_protein_matrix.mtx.gz"),
                     pf, "CD326", n_cells = 7739)
    cd <- as.vector(r$counts[1, ])
    adt_eval <- data.table(内胚层CD326均值 = round(mean(cd[ysc_endo]), 2),
                           其他细胞CD326均值 = round(mean(cd[!ysc_endo]), 2),
                           内胚层CD326阳性率 = round(100 * mean(cd[ysc_endo] > 0), 1))
    log_msg("  ADT CD326: 内胚层 ", mean(cd[ysc_endo]), " vs 其他 ", mean(cd[!ysc_endo]))
  }
}, silent = TRUE)
if (!is.null(adt_eval)) fwrite(adt_eval, file.path(RES, "adt_validation.csv"))

## =============================================================================
## G. YS_SS2（136 细胞全长，仅作旁证）
## =============================================================================
log_msg("== G. YS_SS2 ==")
ss2_cache <- file.path(RES, "cache_ss2.rds")
ss2 <- cache_or_run(ss2_cache, {
  ## 注意：SS2 的 genes 文件首行是表头（19,735 行 vs 矩阵 19,734 列）
  g <- read_feat(file.path(DATA_PROC, "E-MTAB-10888_YS_SS2/YS_SS2_genes.tsv.gz"), header = TRUE)
  keep <- intersect(all_need, g)
  r <- extract_mtx(file.path(DATA_PROC, "E-MTAB-10888_YS_SS2/YS_SS2_counts.mtx.gz"),
                   g, keep, n_cells = 136)
  out <- list(M = cp10k_log(r$counts, r$libsize), libsize = r$libsize)
  rm(r); invisible(gc()); out
})
log_msg("  YS_SS2: ", ncol(ss2$M), " 细胞（仅旁证，不入主图）")

## =============================================================================
## G2. 胎肾基质 —— 原图 “Fetal kidney fibroblast” 列的数据来源（外部数据）
##     Stewart BJ, et al. Science 2019;365:1461-1466
##     "Spatiotemporal immune zonation of the human kidney"
##     CELLxGENE: "Fetal kidney dataset: stroma"（8,345 细胞，胎儿 9-16 pcw，6 供体）
##     该集合 “Fetal kidney dataset: full” 含 27,197 细胞，
##     与论文 methods 中 kidneys(n=4, k=26,372) 同源（QC 后差 825）
##     下载：https://cellxgene.cziscience.com/e/c52de62a-058d-4d78-a464-bdf552378f43.cxg/
## =============================================================================
log_msg("== G2. 胎肾基质（外部 CELLxGENE）==")
kid_cache <- file.path(RES, "cache_kidney.rds")
kid <- cache_or_run(kid_cache, {
  suppressMessages(library(rhdf5))
  H5 <- "fig2d/data/fetal_kidney_stroma.h5ad"
  ## CELLxGENE 分类列：categories + codes
  rd <- function(p) {
    o <- h5read(H5, p)
    if (is.list(o) && !is.null(o$categories))
      as.character(o$categories)[as.integer(o$codes) + 1L] else as.character(o)
  }
  ct  <- rd("obs/cell_type"); act <- rd("obs/author_cell_type")
  tis <- rd("obs/tissue");    dev <- rd("obs/development_stage")
  don <- rd("obs/donor_id")
  gn  <- rd("raw/var/feature_name")
  ip  <- h5read(H5, "raw/X/indptr")
  idx <- as.integer(h5read(H5, "raw/X/indices")) + 1L
  dat <- as.numeric(h5read(H5, "raw/X/data"))
  n   <- length(ct)
  ci  <- rep(seq_len(n), times = diff(ip))
  ## 全转录组 UMI：把所有非零元素按细胞求和（先算完再丢弃）
  lib <- as.numeric(rowsum(dat, ci, reorder = TRUE))
  rm(ip)
  keep <- which(gn %in% all_need)
  sel  <- which(idx %in% keep)
  cnt  <- sparseMatrix(i = match(idx[sel], keep), j = ci[sel], x = dat[sel],
                       dims = c(length(keep), n))
  rownames(cnt) <- gn[keep]
  rm(idx, dat, ci, sel); invisible(gc())
  out <- list(M = cp10k_log(cnt, lib), libsize = lib, counts = cnt,
              cell_type = ct, author_cell_type = act, tissue = tis,
              dev = dev, donor = don,
              source = "Stewart 2019 Science (CELLxGENE fetal kidney stroma)")
  out
})
log_msg("  胎肾基质: ", ncol(kid$M), " 细胞 | ",
        paste(names(sort(table(kid$cell_type), decreasing = TRUE)),
              sort(table(kid$cell_type), decreasing = TRUE), collapse = " / "))
log_msg("  发育阶段: ", paste(unique(kid$dev), collapse = ", "))
log_msg("  供体: ", paste(sort(unique(kid$donor)), collapse = ", "),
        " | 中位 UMI: ", median(kid$libsize))

## =============================================================================
## H. 成人肝（已有统计）
## =============================================================================
adult_f <- file.path(RES, "adult_hepatocyte_stats.rds")
have_adult <- file.exists(adult_f)
log_msg("== adult: ", if (have_adult) "已就绪" else "缺失", " ==")

## =============================================================================
## 8. 统一数据库 与 分组定义
## =============================================================================
log_msg("== 组装数据库 ==")
db <- list(
  ys10x   = list(M = M_ys,           libsize = ys$libsize,   endo = is_endo,  stage = ys_stage, late = ys_late, aligned = aligned),
  f158    = list(M = f158$M,         libsize = f158$libsize, endo = f158$lvl3 == "ENDODERM", stage = f158$stage, lvl = f158$lvl3),
  fl10x   = list(M = fl$M,           libsize = fl$libsize,   hep = fl_hep),
  flcite  = list(M = flc$M,          libsize = flc$libsize,  hep = flc_hep),
  el      = list(M = el$M,           libsize = el$libsize,   hep = el_hep),
  yscite  = list(M = ysc$M,          libsize = ysc$libsize,  endo = ysc_endo),
  ss2     = list(M = ss2$M,          libsize = ss2$libsize),
  kidney  = list(M = kid$M,          libsize = kid$libsize, cell_type = kid$cell_type,
                 tissue = kid$tissue, dev = kid$dev, donor = kid$donor)
)

## ---- 分组（列）定义 ----------------------------------------------------------
grp <- list()
add_grp <- function(name, M, cells, source, note) {
  if (length(cells) == 0) return(invisible(NULL))
  grp[[name]] <<- list(M = M, cells = cells, source = source, note = note)
}

if (aligned && any(!is.na(ys_stage))) {
  mdy <- ys_meta[seq_len(n_ys)]
  st_map <- list("CS10-11" = c("CS10","CS11"), "CS14-15" = c("CS14","CS15"),
                 "CS17-18" = c("CS17","CS18"), "CS22-23" = c("CS22","CS23"))
  for (nm in names(st_map)) {
    sel <- which(is_endo & ys_stage %in% st_map[[nm]])
    if (length(sel) >= 20)
      add_grp(paste0("YS Endo. ", nm), M_ys, sel, "E-MTAB-10552 YS_10x", "官方注释 LVL2/3=ENDODERM")
  }
  ## 非内胚层对照：CS10-11 期前体/基质群（检验基因面板的特异性；
  ## 原图 “YS Gastrulation” 列用的是 CS7 原肠胚外部数据，本地无，故此处为
  ## 同期非内胚层对照，而非真正的原肠胚细胞）
  gas <- which(ys_stage %in% c("CS10","CS11") &
                 ys_meta$LVL1[seq_len(n_ys)] %in% c("PROGENITOR","STROMA"))
  if (length(gas) >= 50)
    add_grp("YS Prog./Stroma", M_ys, gas, "E-MTAB-10552（含 Wang 2021 外部数据）",
            "CS10-11 前体/基质，非内胚层对照")
} else {
  add_grp("YS Endo. Early", M_ys, which(is_endo)[!ys_late], "E-MTAB-10552 YS_10x",
          "ASGR1+SPINK1+ marker + 2-means")
  add_grp("YS Endo. Late",  M_ys, which(is_endo)[ys_late],  "E-MTAB-10552 YS_10x",
          "ASGR1+SPINK1+ marker + 2-means")
}

add_grp("YS Endo. WE_YS(F158)", f158$M, which(f158$lvl3 == "ENDODERM"),
        "E-MTAB-11673 / data_out Seurat", "官方注释 LVL3=ENDODERM, CS17")

## 胎肾基质：原图 “Fetal kidney fibroblast” 列（外部数据，Stewart 2019 Science）
## 用官方注释的 fibroblast 子集，与原图列名一致
add_grp("Fet. Kidney Fibro.", kid$M, which(kid$cell_type == "fibroblast"),
        "Stewart 2019 Science / CELLxGENE",
        "官方注释 fibroblast（胎肾 9-16 pcw, 6 供体）")

add_grp("Emb. Liver", el$M, which(el_hep), "E-MTAB-11618 EL CITE-seq",
        "15 marker 打分 + 2-means")
add_grp("Fet. Liver", fl$M[, fl_hep, drop = FALSE], seq_len(sum(fl_hep)),
        "E-MTAB-10552 FL_10x", "官方注释 cell.labels=Hepatocyte")
add_grp("Fet. Liver (CITE)", flc$M, which(flc_hep), "E-MTAB-11613 FL CITE-seq",
        "15 marker 打分 + 2-means")
## YS_CITEseq 仅 7,739 细胞、中位 UMI 1,488（很浅），内胚层细胞数过少时
## 不作为主图列，仅保留为 ADT 蛋白层面验证（阈值 150）
if (sum(ysc_endo) >= 150)
  add_grp("YS Endo. (CITE)", ysc$M, which(ysc_endo), "E-MTAB-11549 YS CITE-seq",
          "ASGR1+SPINK1+（ADT CD326 支持）")

log_msg("分组数: ", length(grp), " → ", paste(names(grp), collapse = " | "))

## =============================================================================
## 9. 统计：每组每基因的 表达比例(pct) + CP10K-log 均值
## =============================================================================
stat_list <- list(); gsize <- data.table()
for (g in names(grp)) {
  M <- grp[[g]]$M; sel <- grp[[g]]$cells
  gsize <- rbind(gsize, data.table(group = g, n_cells = length(sel),
                                   source = grp[[g]]$source, note = grp[[g]]$note))
  Ms <- M[, sel, drop = FALSE]
  stat_list[[g]] <- data.table(group = g, gene = rownames(Ms),
                               mean_log = rowMeans(Ms),
                               pct = rowSums(Ms > 0) / ncol(Ms) * 100)
}
stats <- rbindlist(stat_list)

if (have_adult) {
  am <- readRDS(adult_f); ap <- readRDS(file.path(RES, "adult_hepatocyte_pct.rds"))
  stats <- rbind(stats, data.table(group = "adult", gene = names(am),
                                   mean_log = unname(am), pct = unname(ap)))
  gsize <- rbind(gsize, data.table(group = "adult", n_cells = NA_integer_,
                                   source = "GSE124395 Aizarani 2019",
                                   note = "外部数据，marker 圈定肝细胞"))
}

## 基因内 min-max 缩放（对应原图 0-1 legend）
genes_ordered <- unlist(fig2d_genes, use.names = FALSE)
row_group     <- rep(names(fig2d_genes), times = vapply(fig2d_genes, length, 0L))
stats <- stats[!is.na(gene)]
stats[, row_group := row_group[match(gene, genes_ordered)]]
stats[!is.na(mean_log) & !is.na(pct),
      color_scaled := { rng <- range(mean_log, na.rm = TRUE)
                        if (diff(rng) == 0) 0.5 else (mean_log - rng[1]) / diff(rng) },
      by = gene]

fwrite(stats, file.path(RES, "final_dotplot_stats.csv"))
fwrite(gsize, file.path(RES, "final_group_sizes.csv"))

## ---- 数据库索引 --------------------------------------------------------------
idx <- data.table(
  数据集 = c("E-MTAB-10552 YS_10x","E-MTAB-11673 F158/WE_YS","E-MTAB-10552 FL_10x",
             "E-MTAB-11613 FL_CITEseq","E-MTAB-11618 EL_CITEseq",
             "E-MTAB-11549 YS_CITEseq","E-MTAB-10888 YS_SS2",
             "Stewart 2019 胎肾（CELLxGENE）","GSE124395 成人肝"),
  来源文件夹 = c("data/processed","data_out (Seurat) + data/processed","data/processed",
                 "data/processed","data/processed","data/processed","data/processed",
                 "外部 CELLxGENE/HCA","外部 GEO"),
  细胞数 = c(134765, ncol(f158$M), nrow(fl$meta), ncol(flc$M), ncol(el$M),
             ncol(ysc$M), ncol(ss2$M), ncol(kid$M),
             if (have_adult) 3417L else NA),
  在图中的作用 = c("YS 内胚层主数据","缺失样本 WE_YS，官方注释（新增）",
                   "胎肝肝细胞（官方注释）","胎肝肝细胞交叉验证（新增）",
                   "胚胎肝肝细胞","YS 内胚层 ADT 验证（新增）",
                   "全长旁证（不入图）","胎肾基质列（原图缺失项，本次补齐）",
                   "成人肝终点")
)
fwrite(idx, file.path(RES, "database_index.csv"))

saveRDS(list(stats = stats, gsize = gsize, db = db,
             groups = c(names(grp), if (have_adult) "adult"),
             fig2d_genes = fig2d_genes, row_group_big = row_group_big,
             green_solid = green_solid, green_hollow = green_hollow,
             aligned = aligned, align_cor = align_cor, align_cell = align_cell,
             marker_eval = marker_eval, adt_eval = adt_eval, db_index = idx),
        file.path(RES, "final_stats.rds"))

log_msg("分组规模："); print(gsize)
log_msg("数据库构建完成 → final_stats.rds")
