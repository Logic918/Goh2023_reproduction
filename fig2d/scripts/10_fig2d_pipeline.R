.libPaths("D:/R/library")
suppressMessages({
  library(data.table)
  library(Matrix)
})

set.seed(1)

DATA <- "C:/Users/lenovo/Desktop/C/data/processed"
OUT  <- "fig2d"
dir.create(file.path(OUT, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## -----------------------------------------------------------------------------
## 1. 基因面板（逐字取自原图 Figure 2D）
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
  "Metabolic"                  = "Metabolic",
  "Common pathway"             = "Haemostasis",
  "Extrinsic pathway"          = "Haemostasis",
  "Intrinsic pathway"          = "Haemostasis",
  "Fibrinogen & crosslinkning" = "Haemostasis",
  "Anticoagulation"            = "Haemostasis",
  "Fibrinolysis"               = "Haemostasis",
  "Growth factors"             = "Growth factors"
)

## 小鼠纯合敲除出生前致死表型（原图绿色椭圆）
green_solid  <- c("F2","F5","F10","F3")            # 胎肝行使功能之前
green_hollow <- c("SERPINC1","PROC","PROS1","EPO") # 之后

## 早/晚内胚层状态基因（论文正文：早期 = PRSS3/GSTA2/LGALS3 高、视黄酸与
## 脂质代谢活跃；晚期 = cell stress / death 基因）
early_markers <- c("PRSS3","GSTA2","LGALS3","APOA1","APOA2","APOC3","TTR","AFP",
                   "ALB","SPINK1","AHSG","FGA","FGB","FGG","F2","F5","F10","F3",
                   "SERPINC1","PROS1","EPO","THPO")
late_markers  <- c("FOS","JUN","JUNB","FOSB","EGR1","ATF3","HSPA1A","HSPA1B",
                   "HSPH1","DNAJB1","BTG2","DUSP1","ZFP36","IER2")

## 内胚层识别（论文：ASGR1+SPINK1+ endoderm）
endo_pos <- c("ASGR1","SPINK1","TTR","APOA1","APOA2","APOC3","ALB","AFP",
              "SERPINA1","AHSG","FGA","FGB")
endo_neg <- c("PECAM1","VWF","CDH5","PTPRC","HBA1","HBB","LYZ","C1QA","COL1A1",
              "TAGLN","ACTA2","MSR1")

## 胚胎肝肝细胞 marker（EL CITE-seq 无注释）
hepato_panel <- c("ALB","AFP","TTR","APOA1","APOA2","APOC3","ASGR1","AHSG",
                  "FGA","FGB","GC","SERPINA1","AMBP","CYP3A7","FETUB")

all_need <- unique(c(unlist(fig2d_genes), early_markers, late_markers,
                     endo_pos, endo_neg, hepato_panel))

## -----------------------------------------------------------------------------
## 2. 通用函数：分块扫描 .mtx.gz（cells x genes），提取目标基因
## -----------------------------------------------------------------------------
extract_mtx <- function(mtx_gz, genes, keep_genes, n_cells,
                        chunk_lines = 5e6L) {
  n_g <- length(genes)
  g_idx <- match(keep_genes, genes)
  stopifnot(!any(is.na(g_idx)))
  g_logical <- logical(n_g); g_logical[g_idx] <- TRUE

  con <- gzfile(mtx_gz, open = "rt")
  hdr <- readLines(con, n = 3)   # %%MatrixMarket / % / dims
  dims <- as.integer(strsplit(gsub("%.*", "", hdr[3]), "\\s+")[[1]])
  ## 自动判断方向：cells x genes 或 genes x cells
  if (dims[1] == n_cells && dims[2] == n_g) {
    cell_col <- "V1"; gene_col <- "V2"
  } else if (dims[1] == n_g && dims[2] == n_cells) {
    cell_col <- "V2"; gene_col <- "V1"
  } else stop("matrix dims ", dims[1], " x ", dims[2],
              " do not match genes (", n_g, ") / cells (", n_cells, ")")
  nnz <- dims[3]
  log_msg("  matrix ", dims[1], " x ", dims[2], ", nnz = ", nnz,
          if (cell_col == "V1") " (cells x genes)" else " (genes x cells)")

  libsize <- numeric(n_cells)
  gi_all <- list(); ci_all <- list(); x_all <- list(); k <- 0
  total <- 0L
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    total <- total + length(lines)
    dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE,
                colClasses = list(integer = c("V1","V2"), numeric = "V3"),
                showProgress = FALSE)
    setnames(dt, c("gcol", "ccol", "xval"))
    if (gene_col == "V2") setnames(dt, c("ccol", "gcol", "xval"))
    tmp <- dt[, .(s = sum(xval)), by = ccol]
    libsize[tmp$ccol] <- libsize[tmp$ccol] + tmp$s
    dt <- dt[gcol %in% g_idx]
    if (nrow(dt)) {
      k <- k + 1
      gi_all[[k]] <- dt$gcol
      ci_all[[k]] <- dt$ccol
      x_all[[k]]  <- dt$xval
    }
    rm(dt, tmp); invisible(gc(verbose = FALSE))
    log_msg("  scanned ", total, "/", nnz, " entries, kept ",
            sum(vapply(x_all, length, 0L)), " target values")
    if (total >= nnz) break
  }
  close(con)

  gi <- unlist(gi_all); ci <- unlist(ci_all); xx <- unlist(x_all)
  M <- sparseMatrix(i = match(gi, g_idx), j = ci, x = xx,
                    dims = c(length(g_idx), n_cells),
                    dimnames = list(keep_genes, paste0("c", seq_len(n_cells))))
  list(counts = M, libsize = libsize)
}

cp10k_log <- function(M, libsize, target = 1e4) {
  sc <- target / libsize
  M2 <- t(t(M) * sc)
  M2@x <- log1p(M2@x)
  M2
}

## -----------------------------------------------------------------------------
## 3. YS 10x：识别内胚层（ASGR1+/SPINK1+）→ 早/晚状态
## -----------------------------------------------------------------------------
log_msg("== YS 10x ==")
ys_cache <- file.path(OUT, "results", "ys_extract_cache.rds")
if (file.exists(ys_cache)) {
  log_msg("loading YS extract cache ...")
  ys_cache_ls <- readRDS(ys_cache)
  M_ys <- ys_cache_ls$M
} else {
  ys_genes <- fread(file.path(DATA, "E-MTAB-10552_YS_10x/YS_10x_genes.tsv.gz"),
                    header = TRUE)[[1]]
  ys_meta  <- fread(file.path(DATA, "E-MTAB-10552_YS_10x/YS_10x_core_metadata.csv"))
  need_ys  <- intersect(all_need, ys_genes)
  log_msg("genes not found in YS gene list: ",
          paste(setdiff(all_need, ys_genes), collapse = ", "))

  res_ys <- extract_mtx(file.path(DATA, "E-MTAB-10552_YS_10x/YS_10x_counts.mtx.gz"),
                        ys_genes, need_ys, n_cells = 134765)
  M_ys <- cp10k_log(res_ys$counts, res_ys$libsize)
  rm(res_ys); invisible(gc())
  saveRDS(list(M = M_ys), ys_cache)
}

## sanity：metadata 中 ENDODERM 7,656 / 全部 169,494 个细胞（含矩阵外细胞）
ys_meta  <- fread(file.path(DATA, "E-MTAB-10552_YS_10x/YS_10x_core_metadata.csv"))
log_msg("metadata ENDODERM cells = ",
        ys_meta[LVL2 == "ENDODERM", .N], " (of ", nrow(ys_meta),
        " metadata rows; counts matrix holds 134765 cells)")

pos <- intersect(endo_pos, rownames(M_ys))
neg <- intersect(endo_neg, rownames(M_ys))
pos_score <- colMeans(as.matrix(M_ys[pos, ]))
neg_score <- colMeans(as.matrix(M_ys[neg, ]))
## ASGR1 与 SPINK1 共表达（论文定义）
asgr1 <- M_ys["ASGR1", ] > 0
spink1 <- M_ys["SPINK1", ] > 0
cand <- which(asgr1 & spink1)
log_msg("ASGR1+ & SPINK1+ cells: ", length(cand))

## 在候选细胞中 2-means 分离内胚层 vs 其他 ASCR1/SPINK1 散阳性细胞
km_e <- kmeans(cbind(scale(pos_score[cand]), scale(-neg_score[cand])),
               centers = 2, nstart = 20)
endo_cl <- names(which.max(tapply(pos_score[cand] - neg_score[cand],
                                  km_e$cluster, mean)))
is_endo <- cand[km_e$cluster == endo_cl]
log_msg("YS endoderm identified in matrix: ", length(is_endo),
        " (metadata reference: 7656)")

## 早/晚状态
em <- as.matrix(M_ys[intersect(early_markers, rownames(M_ys)), is_endo])
lmk <- as.matrix(M_ys[intersect(late_markers, rownames(M_ys)), is_endo])
early_score <- colMeans(em); late_score <- colMeans(lmk)
km_s <- kmeans(cbind(scale(early_score), scale(late_score)), centers = 2,
               nstart = 20)
late_cl <- names(which.max(tapply(late_score - early_score, km_s$cluster, mean)))
is_late <- km_s$cluster == late_cl
log_msg("early / late endoderm: ", sum(!is_late), " / ", sum(is_late))

M_ys_endo <- M_ys[, is_endo, drop = FALSE]
group_ys <- ifelse(is_late, "YS Endoderm (late)", "YS Endoderm (early)")

## -----------------------------------------------------------------------------
## 4. 胚胎肝肝细胞（EL CITE-seq，marker 打分 + 2-means）
## -----------------------------------------------------------------------------
log_msg("== EL CITE-seq ==")
el_genes <- fread(file.path(DATA, "E-MTAB-11618_EL_CITEseq/EL_CITEseq_RNA_features.tsv.gz"),
                  header = FALSE)
el_genes <- if (ncol(el_genes) >= 2) el_genes[[2]] else el_genes[[1]]
need_el <- intersect(unique(c(all_need, hepato_panel)), el_genes)
res_el <- extract_mtx(file.path(DATA, "E-MTAB-11618_EL_CITEseq/EL_CITEseq_RNA_matrix.mtx.gz"),
                      el_genes, need_el, n_cells = 26216)
M_el_raw <- res_el$counts
M_el <- cp10k_log(M_el_raw, res_el$libsize)
rm(res_el); invisible(gc())

hp <- intersect(hepato_panel, rownames(M_el))
hpc_score <- colMeans(as.matrix(M_el[hp, ]))
km_el <- kmeans(scale(hpc_score), centers = 2, nstart = 20)
hep_cl <- names(which.max(tapply(hpc_score, km_el$cluster, mean)))
is_hep <- km_el$cluster == hep_cl
log_msg("EL hepatocytes (2-means): ", sum(is_hep), "/", ncol(M_el),
        " | ALB+ in hep cluster = ",
        round(100 * mean(M_el_raw["ALB", is_hep] > 0), 1), "%",
        " | in others = ", round(100 * mean(M_el_raw["ALB", !is_hep] > 0), 1), "%")

## -----------------------------------------------------------------------------
## 5. 胎肝肝细胞（FL_10x，自带注释）
## -----------------------------------------------------------------------------
log_msg("== FL 10x ==")
fl_meta <- fread(file.path(DATA, "E-MTAB-10552_FL_10x/FL_10x_barcodes_and_metadata.tsv.gz"))
fl_genes <- fread(file.path(DATA, "E-MTAB-10552_FL_10x/FL_10x_genes.tsv.gz"),
                  header = TRUE)[[1]]
fl_hep_n <- fl_meta[cell.labels == "Hepatocyte", .N]
log_msg("FL hepatocytes (annotated): ", fl_hep_n)
need_fl <- intersect(all_need, fl_genes)
res_fl <- extract_mtx(file.path(DATA, "E-MTAB-10552_FL_10x/FL_10x_counts.mtx.gz"),
                      fl_genes, need_fl, n_cells = 49376)
## 只保留 Hepatocyte 列
keep_idx <- which(fl_meta$cell.labels == "Hepatocyte")
M_fl <- cp10k_log(res_fl$counts[, keep_idx, drop = FALSE],
                  res_fl$libsize[keep_idx])
rm(res_fl); invisible(gc())

## -----------------------------------------------------------------------------
## 6. 成人肝（可选，外部数据）
## -----------------------------------------------------------------------------
adult_file <- file.path(OUT, "results", "adult_hepatocyte_stats.rds")
have_adult <- file.exists(adult_file)
if (have_adult) log_msg("== adult liver stats found ==") else
  log_msg("== adult liver 数据未就绪，该列暂缺 ==")

## -----------------------------------------------------------------------------
## 7. 汇总：每组每基因 表达比例 + CP10K-log 平均表达（基因内 min-max 缩放）
## -----------------------------------------------------------------------------
genes_ordered <- unlist(fig2d_genes, use.names = FALSE)
row_group     <- rep(names(fig2d_genes), times = vapply(fig2d_genes, length, 0L))

group_def <- list(
  "YS Endo. Early"    = list(M = M_ys_endo, cells = which(group_ys == "YS Endoderm (early)")),
  "YS Endo. Late"     = list(M = M_ys_endo, cells = which(group_ys == "YS Endoderm (late)")),
  "Emb. Liver"        = list(M = M_el, cells = which(is_hep)),
  "Fet. Liver"        = list(M = M_fl, cells = seq_len(ncol(M_fl)))
)

stat_list <- list(); group_sizes <- data.table()
for (g in names(group_def)) {
  M <- group_def[[g]]$M; sel <- group_def[[g]]$cells
  group_sizes <- rbind(group_sizes, data.table(group = g, n_cells = length(sel)))
  if (!length(sel)) next
  Ms <- M[, sel, drop = FALSE]
  stat_list[[g]] <- data.table(group = g, gene = rownames(M),
                               mean_log = rowMeans(Ms),
                               pct = rowSums(Ms > 0) / ncol(Ms) * 100)
}
stats <- rbindlist(stat_list)

if (have_adult) {
  adult_means <- readRDS(adult_file)
  adult_pct <- readRDS(file.path(OUT, "results", "adult_hepatocyte_pct.rds"))
  ad <- data.table(group = "adult", gene = names(adult_means),
                   mean_log = unname(adult_means), pct = unname(adult_pct))
  stats <- rbind(stats, ad)
  group_sizes <- rbind(group_sizes, data.table(group = "adult", n_cells = NA))
}

stats[, row_group := row_group[match(gene, genes_ordered)]]
stats[!is.na(mean_log) & !is.na(pct),
      color_scaled := { rng <- range(mean_log, na.rm = TRUE)
                        if (diff(rng) == 0) 0.5 else (mean_log - rng[1]) / diff(rng) },
      by = gene]

fwrite(stats[!is.na(gene)], file.path(OUT, "results", "fig2d_dotplot_stats.csv"))
fwrite(group_sizes, file.path(OUT, "results", "fig2d_group_sizes.csv"))
saveRDS(list(stats = stats, fig2d_genes = fig2d_genes,
             row_group_big = row_group_big, green_solid = green_solid,
             green_hollow = green_hollow, group_sizes = group_sizes,
             have_adult = have_adult, col_order = names(group_def)),
        file.path(OUT, "results", "fig2d_stats.rds"))
log_msg("done.")
print(group_sizes)
