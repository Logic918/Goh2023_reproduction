.libPaths("D:/R/library")
suppressMessages({
  library(ggplot2); library(ggforce); library(patchwork)
  library(magick); library(data.table)
})

OUT <- "fig2d"
S <- readRDS(file.path(OUT, "results", "fig2d_stats.rds"))
stats <- S$stats

col_order <- intersect(c("YS Endo. Early","YS Endo. Late",
                         "Emb. Liver","Fet. Liver","adult"), unique(stats$group))
stats[, group := factor(group, levels = col_order)]
## 基因行顺序与原图一致（SERPINA1 在最上）
stats[, gene := factor(gene, levels = rev(unlist(S$fig2d_genes, use.names = FALSE)))]
stats <- stats[!is.na(group) & !is.na(gene) & !is.na(pct)]
n_g <- length(glevels <- levels(stats$gene))
n_c <- length(col_order)

## ---- 右侧行分组（小括号 + 大括号）----
rg_small <- stats[, .(ymin = min(as.integer(gene)), ymax = max(as.integer(gene))),
                  by = row_group]
rg_small <- rg_small[!is.na(row_group)]
rg_small[, row_group := factor(row_group, levels = names(S$fig2d_genes))]
rg_small[, `:=`(x = n_c + 1.0, xl = n_c + 0.55)]
rg_big <- rg_small[, .(ymin = min(ymin), ymax = max(ymax),
                       big = S$row_group_big[as.character(row_group[1])]),
                   by = row_group][, .(ymin = min(ymin), ymax = max(ymax),
                                       lab = big[1]), by = big]
rg_big[, `:=`(x = n_c + 3.6, xt = n_c + 4.25)]

## ---- 绿色椭圆（小鼠敲除表型）----
gene_levels <- levels(stats$gene)
ellipse_dt <- rbind(
  data.table(gene = intersect(S$green_solid, gene_levels),  lt = "1"),
  data.table(gene = intersect(S$green_hollow, gene_levels), lt = "2")
)
ellipse_dt[, `:=`(y0 = match(gene, gene_levels), x0 = 0.08, a = 0.62, b = 0.30)]

## ---- 底部分组括号 ----
n_ys <- sum(grepl("YS Endo", col_order))
n_hep <- sum(grepl("Liver", col_order))
brackets <- data.table(
  label = c("YS Endoderm", "Hepatocyte"),
  x0 = c(0.6, n_ys + 0.6),
  x1 = c(n_ys + 0.4, n_ys + n_hep + 0.4)
)
brackets[, y := n_g + 1.3]

p <- ggplot(stats, aes(x = group, y = gene)) +
  geom_point(aes(size = pct, color = color_scaled)) +
  geom_ellipse(data = ellipse_dt,
               aes(x0 = x0, y0 = y0, a = a, b = b, angle = 0, linetype = lt),
               inherit.aes = FALSE,
               color = "green3", fill = NA, linewidth = 0.45) +
  scale_size(range = c(0.5, 5.4), breaks = c(20, 60, 100),
             name = "Fraction of cells\nin group (%)") +
  scale_color_gradient(low = "grey94", high = "#67000d",
                       name = "Mean expression\nin group",
                       breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_linetype_manual(values = c("1" = "solid", "2" = "longdash"),
                        labels = c("1" = "before fetal liver function",
                                   "2" = "after fetal liver function"),
                        name = "Mouse knockout phenotype present") +
  guides(size = guide_legend(override.aes = list(shape = 21, color = "grey40",
                                                 fill = "grey70")),
         color = guide_colorbar(barwidth = 0.9, barheight = 4.5),
         linetype = guide_legend(override.aes = list(color = "green3"))) +
  geom_segment(data = rg_small, aes(x = xl, xend = x, y = ymin - 0.35, yend = ymax + 0.35),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = rg_small, aes(x = x + 0.12, y = (ymin + ymax) / 2, label = row_group),
            inherit.aes = FALSE, size = 2.1, hjust = 0, angle = 90, color = "grey15") +
  geom_segment(data = rg_big, aes(x = x, xend = x, y = ymin - 0.4, yend = ymax + 0.4),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = rg_big, aes(x = xt, y = (ymin + ymax) / 2, label = lab),
            inherit.aes = FALSE, size = 3.6, angle = 90, fontface = "bold") +
  geom_segment(data = brackets, aes(x = x0, xend = x1, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = brackets, aes(x = (x0 + x1) / 2, y = y + 1.1, label = label),
            inherit.aes = FALSE, size = 3.6, fontface = "bold.italic") +
  coord_cartesian(clip = "off", xlim = c(0.4, n_c + 5.6), ylim = c(0.4, n_g + 3.2)) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 7.5),
        axis.title = element_blank(),
        plot.margin = margin(12, 35, 16, 5))

ggsave(file.path(OUT, "figures", "fig2d_reproduction.pdf"), p, width = 8.5, height = 8)
ggsave(file.path(OUT, "figures", "fig2d_reproduction.png"), p, width = 8.5, height = 8, dpi = 200)

## ---- 与原图对比 ----
orig <- image_read(file.path(OUT, "results", "fig2d_panelD.png"))
repro <- image_read(file.path(OUT, "figures", "fig2d_reproduction.png"))
h <- 1500
cmp <- image_append(image_join(image_resize(orig, paste0("x", h)),
                               image_resize(repro, paste0("x", h))), stack = TRUE)
image_write(cmp, file.path(OUT, "figures", "fig2d_comparison.png"), format = "png")
cat("figures written.\n")

三、校验与对齐脚本（解释「为什么不用官方注释」）
本次复现的一个关键判断是放弃 E-MTAB-10552 的官方 ENDODERM 注释。下面三个脚本给出了依据：细胞顺序与条形码对不上、相关性检验不达阈，因此改用 mark
