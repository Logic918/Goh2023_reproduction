#!/usr/bin/env Rscript
# =============================================================================
# 03_fig3_reproduce.R   (v2, 2026-10-21)
#
# Science 2023《Yolk sac cell atlas reveals multiorgan functions during human
# early development》(Goh et al., 381: eadd7564) —— Fig. 3 复现
#
# v2 相对 v1 的改动（对齐原图）：
#   Panel A  由「YS 内 2 行 × 12 基因」扩到「7 行 × 23 基因」，基因分 4 组
#            （Canonical / Early / Definitive / Patterning），行按孕龄与
#            early/definitive 拆分，并加 FL HSC/MPP 一行；配色、点样式、
#            分组括号均按原图重做
#   Panel B  由 6 根 CSI 级柱改为原图的 4 根孕龄柱（~4 / ~5 / ~7 / ~8 PCW）
#   Panel E  基因改为原图的 HBZ / HBG1 / BCL11A / ZBTB7A 四个子图（v1 误用
#            了 HBE1/HBG2/HBB/HBD），补 R² 标注
#
# 数据：E-MTAB-10552（ArrayExpress，开放获取）
#   YS_10x  134,765 细胞 × 30,330 基因（矩阵为「细胞 × 基因」方向）
#   FL_10x   49,376 细胞 × 25,756 基因
#
# 原图 Fig. 3 的 A 面板共 13 行，其中 Gastrulation / AGM / FBM HSC-MPP /
# iPSC / Def. iPSC 共 5 行来自另外几篇论文的独立数据集（原肠胚 Tyser 2021、
# AGM、骨髓、iPSC 分化），不属于本研究的 E-MTAB 数据，本地无法取得，故缺失。
# C（UMAP 密度图）、F（FDG）、G/H（CellRank 径向图）依赖作者提供的 UMAP 坐标
# 与 CellRank 转移矩阵，D 为胚胎照片，均不在本脚本范围内。
#
# 运行：Rscript scripts/03_fig3_reproduce.R
# 输出：figures/fig3AB.pdf|png、fig3A_dotplot_ys、fig3B_early_vs_definitive、
#       fig3E_globin_switch.pdf|png（竖排 4 格）；tables/t15-t17
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)      # 稀疏矩阵
  library(data.table)  # 元数据表读写
  library(ggplot2)     # 绘图
  library(patchwork)   # 多面板拼图
  library(scales)      # 百分比坐标轴
  library(grid)        # 箭头/bracket
})

set.seed(2023)
ROOT   <- "C:/Users/lenovo/Desktop/C"
YS_DIR <- file.path(ROOT, "data/processed/E-MTAB-10552_YS_10x")
FL_DIR <- file.path(ROOT, "data/processed/E-MTAB-10552_FL_10x")
FIGDIR <- file.path(ROOT, "figures")
TABDIR <- file.path(ROOT, "tables")
OUTDIR <- file.path(ROOT, "data_out")
for (d in c(FIGDIR, TABDIR, OUTDIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---------------------------------------------------------------- 配色 ----
COL_EARLY <- "#F0C24B"   # 原图 early HSPC：黄
COL_DEF   <- "#5D8C60"   # 原图 definitive HSPC：绿
COL_YS    <- "#E8859B"   # 原图 YS 红系：粉
COL_LIVER <- "#C0392B"   # 原图 EL 红系：红（此处用 FL 代替，见 README）

# =============================================================================
# 0. 工具函数
# =============================================================================

#' 分块扫描 MatrixMarket，只取「目标基因 × 目标细胞」的条目
#' @param mtx_path .mtx.gz；行 = 细胞，列 = 基因
#' @param genes_all 全部基因名（长度 = 列数）
#' @param genes_keep 目标基因
#' @param cell_mask 逻辑向量（长度 = 行数），TRUE 保留；names 为细胞名
scan_mtx_subset <- function(mtx_path, genes_all, genes_keep, cell_mask,
                            chunk = 5e6, label = "") {
  stopifnot(file.exists(mtx_path))
  con <- gzfile(mtx_path, "r"); on.exit(close(con), add = TRUE)
  hdr  <- readLines(con, n = 3L)
  dims <- as.integer(strsplit(hdr[3L], "\\s+")[[1]])
  if (length(cell_mask) != dims[1]) stop("cell_mask 长度 != 矩阵行数")
  if (length(genes_all) != dims[2]) stop("genes_all 长度 != 矩阵列数")

  gid <- match(genes_keep, genes_all); ok <- !is.na(gid)
  if (any(!ok)) warning("矩阵里找不到：", paste(genes_keep[!ok], collapse = ", "))
  gid <- gid[ok]; gname <- genes_keep[ok]

  gmask <- logical(dims[2]); gmask[gid] <- TRUE
  map   <- integer(dims[1]); map[which(cell_mask)] <- seq_len(sum(cell_mask))
  tot <- numeric(dims[1]); buf <- list(); n <- 0; nb <- 0L; t0 <- Sys.time()

  repeat {
    v <- scan(con, what = list(i = integer(), j = integer(), x = numeric()),
              nlines = chunk, quiet = TRUE)
    if (length(v$i) == 0L) break
    n <- n + length(v$i); nb <- nb + 1L
    agg <- rowsum(v$x, v$i, reorder = TRUE); idx <- as.integer(rownames(agg))
    tot[idx] <- tot[idx] + agg[, 1L]
    keep <- cell_mask[v$i] & gmask[v$j]
    if (any(keep)) buf[[length(buf) + 1L]] <- list(i = v$i[keep], j = v$j[keep], x = v$x[keep])
    if (nb %% 10L == 0L)
      cat(sprintf("  [%s] %5.1f%%  %5.0fs\n", label, 100 * n / dims[3],
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  I <- unlist(lapply(buf, `[[`, "i"), use.names = FALSE)
  J <- unlist(lapply(buf, `[[`, "j"), use.names = FALSE)
  X <- unlist(lapply(buf, `[[`, "x"), use.names = FALSE)
  mat <- sparseMatrix(i = map[I], j = match(J, gid), x = X,
                      dims = c(sum(cell_mask), length(gid)),
                      dimnames = list(names(cell_mask)[cell_mask], gname))
  step(label, "扫描完成：", format(n, big.mark = ","), " 个非零值，用时 ",
       round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " 分钟")
  list(mat = mat, tot = tot[cell_mask], cells = names(cell_mask)[cell_mask], genes = gname)
}

read_gene_table <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, colClasses = list(character = 1))
  g <- dt[[1]]; g[!duplicated(g)]
}

#' 带基因清单校验的缓存读取：清单变了就重扫
load_scan <- function(cache, mtx, genes_all, genes_keep, cell_mask, label) {
  if (file.exists(cache)) {
    obj <- readRDS(cache)
    same <- setequal(obj$cells, names(cell_mask)[cell_mask]) &&
      setequal(obj$genes, genes_keep)
    if (same) { step("复用缓存 ", basename(cache), "（", length(obj$cells), " 个细胞）"); return(obj) }
    step("缓存与当前细胞/基因集不一致，重新扫描")
  }
  obj <- scan_mtx_subset(mtx, genes_all, genes_keep, cell_mask, label = label)
  saveRDS(obj, cache); obj
}

# 孕龄：Carnegie stage -> 受精后周数；胎肝 F## = 受精后 ## 天
CS2PCW <- c(CS10 = 4.0, CS11 = 4.3, CS12 = 4.6, CS13 = 4.8, CS14 = 5.0, CS15 = 5.4,
            CS16 = 5.8, CS17 = 6.1, CS18 = 6.5, CS19 = 7.0, CS20 = 7.4, CS21 = 7.6,
            CS22 = 8.0, CS23 = 8.5)
F2PCW  <- c(F32 = 32 / 7, F35 = 35 / 7, F61 = 61 / 7)

# =============================================================================
# 1. 基因面板（照原图 Fig. 3A 的 4 个分组）与全局配色
# =============================================================================
GENE_CANON <- c("CD34", "MLLT3", "SPINK2", "HOPX", "HLF", "RAB27B", "MYB")
GENE_EARLY <- c("DDIT4", "SLC2A3", "RGS16", "LIN28A")
GENE_DEF   <- c("KIT", "ITGA4", "CD74", "PROCR", "EMCN", "GBP4", "ACE")
GENE_PATT  <- c("HOXA7", "HOXA9", "HOXA10", "HOXB7", "HOXB9")
GENE_HSPC  <- c(GENE_CANON, GENE_EARLY, GENE_DEF, GENE_PATT)   # 23 个

GENE_HB    <- c("HBZ", "HBG1", "BCL11A", "ZBTB7A")            # Fig. 3E 四个子图
GENE_ALL   <- unique(c(GENE_HSPC, GENE_HB))

CLS_MAP <- c(EARLY_HSPC_1 = "Early", EARLY_HSPC_2 = "Early",
             DEF_HSPC_1   = "Definitive", DEF_HSPC_2 = "Definitive")

# =============================================================================
# 2. 读 YS 元数据
# =============================================================================
step("读取 YS 元数据 ...")
meta_ys <- fread(file.path(YS_DIR, "YS_10x_core_metadata.csv"), sep = ",",
                 header = TRUE, colClasses = list(character = 1))
setnames(meta_ys, 1, "cell")
bc_ys <- fread(file.path(YS_DIR, "YS_10x_barcodes.tsv.gz"), header = FALSE)[[1]]
n_cell_ys <- 134765L                     # 矩阵行数，第 i 行 = 第 i 个 barcode
meta_ys <- meta_ys[match(bc_ys[seq_len(n_cell_ys)], cell)]
stopifnot(!anyNA(meta_ys$cell))

is_hspc <- !is.na(meta_ys$LVL3) & meta_ys$LVL3 %in% names(CLS_MAP)
is_ery  <- !is.na(meta_ys$LVL2) & meta_ys$LVL2 == "ERYTHROID"
step("YS HSPC ", sum(is_hspc), " 个；红系 ", sum(is_ery), " 个")

cell_mask <- is_hspc | is_ery; names(cell_mask) <- meta_ys$cell

# =============================================================================
# 3. 扫描 YS 矩阵
# =============================================================================
step("扫描 YS 矩阵（2.76 亿非零值，首次约 2.6 分钟）...")
genes_ys <- read_gene_table(file.path(YS_DIR, "YS_10x_genes.tsv.gz"))
ys <- load_scan(file.path(OUTDIR, "fig3_ys_scan.rds"),
                file.path(YS_DIR, "YS_10x_counts.mtx.gz"),
                genes_ys, GENE_ALL, cell_mask, "YS")

# LogNormalize：每细胞总量归一到 1e4 后 log1p
expr_ys <- log1p(sweep(as.matrix(ys$mat), 1, pmax(ys$tot, 1) / 1e4, "/"))
meta_s  <- meta_ys[cell_mask, ]

# =============================================================================
# 4. 扫描 FL 矩阵（红系 + HSC/MPP）
# =============================================================================
step("读取 FL 元数据 ...")
meta_fl <- fread(file.path(FL_DIR, "FL_10x_barcodes_and_metadata.tsv.gz"), sep = "\t",
                 header = TRUE, colClasses = list(character = 1))
setnames(meta_fl, 1, "cell")
fl_lab <- trimws(as.character(meta_fl$cell.labels))
is_ery_fl  <- !is.na(fl_lab) & fl_lab %in% c("Erythroid", "Early Erythroid")
is_hsc_fl  <- !is.na(fl_lab) & fl_lab == "HSC/MPP"
fl_mask <- is_ery_fl | is_hsc_fl; names(fl_mask) <- meta_fl$cell
step("FL 红系 ", sum(is_ery_fl), " 个，HSC/MPP ", sum(is_hsc_fl), " 个")

genes_fl <- read_gene_table(file.path(FL_DIR, "FL_10x_genes.tsv.gz"))
fl <- load_scan(file.path(OUTDIR, "fig3_fl_scan.rds"),
                file.path(FL_DIR, "FL_10x_counts.mtx.gz"),
                genes_fl, GENE_ALL, fl_mask, "FL")
expr_fl <- log1p(sweep(as.matrix(fl$mat), 1, pmax(fl$tot, 1) / 1e4, "/"))
meta_fl_s <- meta_fl[fl_mask, ]

# =============================================================================
# 5. Panel A —— HSPC marker dot plot
#    行 = 孕龄 × early/definitive（6 行）+ FL HSC/MPP（1 行）
# =============================================================================
step("Panel A: dot plot（7 行 × 23 基因）...")

PCW_BIN <- c(CS10 = 4, CS11 = 4, CS12 = 5, CS13 = 5, CS14 = 5, CS15 = 5,
             CS16 = 7, CS17 = 7, CS18 = 7, CS19 = 8, CS20 = 8, CS21 = 8,
             CS22 = 8, CS23 = 8)

grp_ys <- unname(CLS_MAP[meta_s$LVL3])
bin_ys <- unname(PCW_BIN[meta_s$stage])

row_lab <- c("~4 PCW HSPCs", "~5 PCW HSPCs", "~7 PCW HSPCs",
             "~5 PCW HSPCs", "~7 PCW HSPCs", "~8 PCW HSPCs")
row_cls <- c("Early", "Early", "Early", "Definitive", "Definitive", "Definitive")
row_bin <- c(4, 5, 7, 5, 7, 8)

dot_list <- list()
for (i in seq_along(row_lab)) {
  sel <- !is.na(grp_ys) & grp_ys == row_cls[i] & !is.na(bin_ys) & bin_ys == row_bin[i]
  if (sum(sel) < 5) { step("  ", row_lab[i], " (", row_cls[i], ") 细胞过少，跳过：", sum(sel)); next }
  m <- expr_ys[sel, GENE_HSPC, drop = FALSE]
  dot_list[[length(dot_list) + 1L]] <- data.table(
    row = length(dot_list) + 1L, label = row_lab[i], cls = row_cls[i],
    tissue = "YS", n = sum(sel), gene = GENE_HSPC,
    mean = colMeans(m), pct = colMeans(m > 0))
}

# FL HSC/MPP 一行
sel_fl <- is_hsc_fl[fl_mask]
m <- expr_fl[sel_fl, GENE_HSPC, drop = FALSE]
dot_list[[length(dot_list) + 1L]] <- data.table(
  row = length(dot_list) + 1L, label = "FL HSC/MPP", cls = "Definitive",
  tissue = "FL", n = sum(sel_fl), gene = GENE_HSPC,
  mean = colMeans(m), pct = colMeans(m > 0))

dot <- rbindlist(dot_list)
# ---- 颜色标尺 ----
# 每个基因在行间做 min-max 缩放，对应原图 0–1 的 Mean expression 色标。
# 量程只取 6 个 YS 行：FL 是 CD45+ 分选、测序更深，若把 FL 也算进去，
# 每个基因的最大值几乎都被 FL 顶掉，YS 那 6 行会整体压白，HOPX、SPINK2
# 这类本该在 ~4 PCW 最红的基因也会失去对比（原图 HOPX 在 ~4 PCW 就是最深的一格）。
# FL 超出量程的部分截到 1，原始 mean 在 tables/t15 里可查。
ys_rng <- dot[tissue == "YS", .(lo = min(mean), hi = max(mean)), by = gene]
dot <- merge(dot, ys_rng, by = "gene", sort = FALSE)
# 分母设下限 0.05（log 单位）：HOXA9 这类在 YS 里几乎为 0 的基因，量程本来只有
# 0.004，若直接 min-max，0.0046 的小噪声会被拉成最深色，画出假的"高表达"。
dot[, scaled := pmin(pmax((mean - lo) / pmax(hi - lo, 0.05), 0), 1), by = gene]
dot[, gidx := match(gene, GENE_HSPC)]
dot[, gcat := factor(fcase(gene %in% GENE_CANON, "Canonical",
                           gene %in% GENE_EARLY, "Early",
                           gene %in% GENE_DEF,   "Definitive",
                           gene %in% GENE_PATT,  "Patterning"),
                     levels = c("Canonical", "Early", "Definitive", "Patterning"))]
fwrite(dot, file.path(TABDIR, "t15_fig3A_dotplot.csv"))

# ---- 行标签：照原图，行名只写孕龄，early / definitive 用左侧括线标 ----
dot[, row_txt := fcase(tissue == "FL", "FL HSC/MPP",
                       default = label)]                  # label 已是「~4 PCW HSPCs」
dot[, rk := as.integer(row)]                           # 行按创建顺序 1..7
row_levels <- dot[, .SD[1L], by = rk][order(rk), row_txt]   # 不去重：~5 / ~7 各出现两次
nrow_a <- length(row_levels)
dot[, ry := nrow_a - rk + 1]                           # 第 1 行在最上方（y 越大越靠上）
top_y  <- nrow_a + 0.20

cat_bounds <- dot[, .(xmin = min(gidx), xmax = max(gidx), n = .N), by = gcat][order(gcat)]
sep_x <- cat_bounds$xmin[-1] - 0.5          # 分组之间的竖线位置

# ---- 顶部 4 个基因组的括线 ----
mk_bracket <- function(x1, x2, y, depth) {
  data.table(x = c(x1, x1, x2, x2), y = c(y + depth, y, y, y + depth))
}
brk <- rbindlist(lapply(seq_len(nrow(cat_bounds)), function(i)
  mk_bracket(cat_bounds$xmin[i] - 0.26, cat_bounds$xmax[i] + 0.26, top_y + 0.12, 0.26)))
lbl <- data.table(x = (cat_bounds$xmin + cat_bounds$xmax) / 2, y = top_y + 0.54,
                  lab = as.character(cat_bounds$gcat))

# ---- 左侧 early / definitive 括线（原图版式：括号在面板外，文字在旁边）----
side_brk <- function(rows, x = -0.30) {
  y1 <- min(rows); y2 <- max(rows)
  data.table(x = c(x + 0.16, x, x, x + 0.16), y = c(y2, y2, y1, y1))
}
sb <- rbindlist(list(side_brk(c(5, 6, 7)), side_brk(c(2, 3, 4))))
sb_lab <- data.table(x = -0.62, y = c(6, 3), lab = c("Early", "Def."))

pA <- ggplot(dot, aes(gidx, ry)) +
  geom_vline(xintercept = sep_x, colour = "grey80", linewidth = 0.3) +
  # 原图只有一条虚线：YS 的 HSPC 行与其他组织之间
  geom_hline(yintercept = 1.5, linetype = "22", colour = "grey45", linewidth = 0.45) +
  geom_point(aes(size = pct, fill = scaled), shape = 21, colour = "grey25", stroke = 0.28) +
  geom_path(data = brk, aes(x, y), inherit.aes = FALSE, colour = "grey20", linewidth = 0.35) +
  geom_text(data = lbl, aes(x, y, label = lab), inherit.aes = FALSE,
            size = 2.35, colour = "grey10", fontface = "bold") +
  geom_path(data = sb, aes(x, y), inherit.aes = FALSE, colour = "grey25", linewidth = 0.35) +
  geom_text(data = sb_lab, aes(x, y, label = lab), inherit.aes = FALSE, angle = 90,
            size = 2.2, colour = "grey15") +
  scale_fill_gradientn(colours = c("#FFFFFF", "#FCBBA1", "#F67E5A", "#CB181D", "#67000D"),
                       limits = c(0, 1), breaks = c(0, 0.5, 1),
                       name = "Mean expression\nin group",
                       guide = guide_colourbar(order = 2, barwidth = unit(0.55, "cm"),
                                               barheight = unit(2.0, "cm"))) +
  scale_size_area(max_size = 4.6, limits = c(0, 1), breaks = c(.2, .6, 1),
                  labels = c(20, 60, 100),
                  name = "Fraction of cells\nin group (%)",
                  guide = guide_legend(order = 1, direction = "vertical",
                                       override.aes = list(fill = "grey60", colour = "grey25",
                                                           stroke = 0.28))) +
  scale_x_continuous(breaks = seq_along(GENE_HSPC), labels = GENE_HSPC, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_len(nrow_a), labels = rev(row_levels), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.35, length(GENE_HSPC) + 0.65),
                  ylim = c(0.30, top_y + 0.64), clip = "off") +
  labs(tag = "A", x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey35", linewidth = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   face = "italic", size = 7.3),
        axis.text.y = element_text(size = 7.3, margin = margin(r = 20)),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.ticks = element_line(colour = "grey35", linewidth = 0.3),
        legend.position = "right", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 2),
        legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.4),
        legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.12, "cm"),
        legend.box.spacing = unit(0.15, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.005, 0.985),
        plot.margin = margin(10, 6, 4, 50))

# 自检：Early / Definitive 打分能否区分两组（AUC）
auc <- function(x, y) { r <- rank(c(x, y)); n1 <- length(x); n2 <- length(y)
  (sum(r[seq_along(x)]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
sc <- data.table(group = grp_ys[!is.na(grp_ys)],
                 early = rowMeans(expr_ys[!is.na(grp_ys), GENE_EARLY, drop = FALSE]),
                 defin = rowMeans(expr_ys[!is.na(grp_ys), GENE_DEF,   drop = FALSE]))
fwrite(sc, file.path(TABDIR, "t15b_marker_score_by_group.csv"))
step(sprintf("自检 AUC：early 打分 = %.3f；definitive 打分 = %.3f",
             auc(sc[group == "Early", early], sc[group == "Definitive", early]),
             auc(sc[group == "Definitive", defin], sc[group == "Early", defin])))

# =============================================================================
# 6. Panel B —— 4 个孕龄档的 early / definitive 比例
# =============================================================================
step("Panel B: early vs definitive 比例（4 档）...")
hb <- data.table(bin = bin_ys[!is.na(grp_ys)], cls = grp_ys[!is.na(grp_ys)])
hb <- hb[, .(n = .N), by = .(bin, cls)]
hb[, total := sum(n), by = bin][, prop := n / total]
hb[, cls := factor(cls, levels = c("Early", "Definitive"))]
fwrite(hb, file.path(TABDIR, "t16_fig3B_prop.csv"))
print(hb[order(bin, cls)])

pB <- ggplot(hb, aes(factor(bin, levels = sort(unique(bin))), prop * 100, fill = cls)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c(Early = COL_EARLY, Definitive = COL_DEF),
                    name = NULL, breaks = c("Early", "Definitive"),
                    labels = c("Early\nHSPCs", "Definitive\nHSPCs")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_x_discrete(labels = function(x) paste0("~", x, " PCW")) +
  labs(tag = "B", x = NULL, y = "Fraction of cells in population (%)") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey35", linewidth = 0.5),
        axis.text.x = element_text(size = 6.8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.6),
        axis.title.y = element_text(size = 7),
        legend.position = "top", legend.justification = "left",
        legend.key.size = unit(0.28, "cm"), legend.spacing.x = unit(0.08, "cm"),
        legend.text = element_text(size = 6.4, lineheight = 0.9),
        legend.margin = margin(0, 0, 1, 0),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = margin(4, 4, 4, 4))

# =============================================================================
# 7. Panel E —— HBZ / HBG1 / BCL11A / ZBTB7A 随孕龄变化
# =============================================================================
step("Panel E: 珠蛋白及其抑制因子 ...")
ery_in_sub <- is_ery[cell_mask]
cnt_ys <- as.matrix(ys$mat[ery_in_sub, GENE_HB, drop = FALSE])
keep_c <- rowSums(as.matrix(ys$mat[ery_in_sub, c("HBZ", "HBG1"), drop = FALSE])) > 0
g_ys <- data.table(stage = meta_s$stage[ery_in_sub][keep_c],
                   n_total = sum(keep_c))
g_ys[, (GENE_HB) := as.data.table(cnt_ys[keep_c, , drop = FALSE] > 0)]
g_ys <- melt(g_ys, id.vars = c("stage", "n_total"), variable.name = "gene", value.name = "pos")
g_ys <- g_ys[, .(pct = mean(pos), ncell = .N), by = .(stage, gene, n_total)]
g_ys[, `:=`(tissue = "YS", pcw = CS2PCW[stage])]
g_ys <- g_ys[!is.na(pcw)]

ery_fl_in <- is_ery_fl[fl_mask]
cnt_fl <- as.matrix(fl$mat[ery_fl_in, GENE_HB, drop = FALSE])
keep_cf <- rowSums(as.matrix(fl$mat[ery_fl_in, c("HBZ", "HBG1"), drop = FALSE])) > 0
g_fl <- data.table(stage = meta_fl_s$stages[ery_fl_in][keep_cf])
g_fl[, (GENE_HB) := as.data.table(cnt_fl[keep_cf, , drop = FALSE] > 0)]
g_fl <- melt(g_fl, id.vars = "stage", variable.name = "gene", value.name = "pos")
g_fl <- g_fl[, .(pct = mean(pos), ncell = .N), by = .(stage, gene)]
g_fl[, `:=`(tissue = "Liver", pcw = F2PCW[stage])]
g_fl <- g_fl[!is.na(pcw)]

glob <- rbind(g_ys[, .(stage, gene, pct, tissue, pcw, ncell)],
              g_fl[, .(stage, gene, pct, tissue, pcw, ncell)])
GENE_TXT <- c(HBZ = "HBZ (Hb Gower 1)", HBG1 = "HBG1 (HbF)",
              BCL11A = "BCL11A (HBZ repressor)", ZBTB7A = "ZBTB7A (HBZ repressor)")
fwrite(glob, file.path(TABDIR, "t17_globin_switch.csv"))

# 原图横轴是离散孕龄档（~3 / ~4 / ~5 / ~7 / ~8 PCW），每档每组织只有一个点。
# 本数据在 ~5 PCW 有两个胎肝样本（F32 4.57、F35 5.00），直接按连续孕周画会出现
# 两个几乎重叠的方块；所以先归档、再按细胞数加权合并成一个点。
bin_pcw <- function(w) {
  out <- rep(NA_character_, length(w))
  out[w >= 3.5 & w < 4.5] <- "~4"
  out[w >= 4.5 & w < 6.0] <- "~5"
  out[w >= 6.0 & w < 7.5] <- "~7"
  out[w >= 7.5]           <- "~8"
  out
}
glob[, bin := bin_pcw(pcw)]
glob <- glob[!is.na(bin)]
gb <- glob[, .(pct = sum(pct * ncell) / sum(ncell),
               ncell = sum(ncell),
               nsample = uniqueN(stage),
               pcw = mean(pcw)),
           by = .(bin, gene, tissue)]
gb[, `:=`(bin = factor(bin, levels = c("~4", "~5", "~7", "~8")),
          tissue = factor(tissue, levels = c("YS", "Liver")),
          gene = factor(gene, levels = GENE_HB))]
fwrite(gb, file.path(TABDIR, "t17b_globin_switch_binned.csv"))

# 趋势拟合的数字留给正文，图面不标（原图也没有）
r2 <- gb[, .(r2 = summary(lm(pct ~ pcw, data = .SD))$r.squared), by = gene]
print(r2[, .(gene, r2 = round(r2, 3))])

glab <- data.table(gene = factor(GENE_HB, levels = GENE_HB),
                   x = 0.76, y = 46, lab = unname(GENE_TXT[GENE_HB]))

pE <- ggplot(gb, aes(bin, pct * 100, colour = tissue, shape = tissue, fill = tissue)) +
  geom_line(aes(group = tissue), linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(size = ncell), stroke = 0.75) +
  geom_label(data = glab, aes(x, y, label = lab), inherit.aes = FALSE,
             hjust = 0, vjust = 0, size = 2.45, fontface = "italic",
             colour = "grey15", fill = "white", linewidth = 0,
             label.padding = unit(0.04, "lines")) +
  facet_wrap(~ gene, ncol = 1) +
  scale_colour_manual(values = c(YS = COL_YS, Liver = COL_LIVER), guide = "none") +
  scale_fill_manual(values = c(YS = "white", Liver = COL_LIVER), guide = "none") +
  scale_shape_manual(values = c(YS = 21, Liver = 22), name = "Human",
                     labels = c("YS", "Liver (FL)"),
                     guide = guide_legend(
                       direction = "horizontal",
                       override.aes = list(colour = c(COL_YS, COL_LIVER),
                                           fill = c("white", COL_LIVER),
                                           size = 2.6, stroke = 0.75))) +
  scale_size_area(max_size = 4.6, guide = "none") +
  scale_x_discrete(labels = function(x) paste0(x, " PCW"),
                   expand = expansion(add = c(0.40, 0.34))) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0.035, 0)) +
  labs(tag = "E", x = "Approx. post-conception weeks",
       y = "Percentage of Ery\nenriched in Hb (%)") +
  theme_bw(base_size = 8.6) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey45", linewidth = 0.45),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.y = unit(0.11, "cm"),
        axis.text.x = element_text(size = 7.3, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7.1),
        axis.title = element_text(size = 8.0),
        legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 7.8, face = "bold"),
        legend.text = element_text(size = 7.4),
        legend.key.size = unit(0.34, "cm"),
        legend.margin = margin(0, 0, 2, 0),
        legend.box.spacing = unit(0.1, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = margin(4, 6, 4, 4))

# =============================================================================
# 8. 输出（A/B 横向拼成一条，E 单独出竖版，与原图一样是 4 个基因竖排）
# =============================================================================
save_panel <- function(p, name, w, h, widths = NULL) {
  if (!is.null(widths)) p <- p + plot_layout(widths = widths)
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 400,
         device = ragg::agg_png, bg = "white")
}
save_panel(pA | pB, "fig3AB", 9.0, 3.6, widths = c(3.4, 1))
save_panel(pA,  "fig3A_dotplot_ys", 8.2, 3.6)
save_panel(pB,  "fig3B_early_vs_definitive", 2.9, 3.6)
save_panel(pE,  "fig3E_globin_switch", 3.4, 5.6)

step("图已输出到 ", FIGDIR)
cat("\n===== Panel B 数字 =====\n"); print(dcast(hb, bin ~ cls, value.var = "prop"))
