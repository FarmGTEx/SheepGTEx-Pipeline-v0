library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(grid)
library(gridExtra)
library(dendextend)
library(readxl)


# 读取数据
dftis <- read.delim("tissue40.list", header = TRUE)
tislist <- as.list(dftis$Tissue)
tisorder <- read.delim("tissue.order", sep = "\t", header = TRUE)
palette <- setNames(tisorder$Color, tisorder$`Tissue..QTL.`)
lfsr <- readRDS("lfsr_m.s.RDS")
pm <- readRDS("./pm_m.s.RDS")

correlation = "spearman"
metric = "maximum"
method = "complete"
pm_sig <- pm[apply(lfsr < 0.05, 1, any), ]
pm_sig_corr <- cor(pm_sig, method = correlation)

##提取对应关系
#print(palette)
#str(palette)
#tissue_data <- data.frame(
#  Tissue = names(palette),
#  Color = palette,
#  stringsAsFactors = FALSE
#)
#write.csv(tissue_data, "color_tissue.txt", row.names = FALSE)

# 对相关矩阵进行聚类并重新排列

hc <- hclust(dist(pm_sig_corr, method = metric), method = method)
od <- hc$order

pm_sig_corr <- pm_sig_corr[od, od]
nm <- rownames(pm_sig_corr)
row_hc <- hclust(dist(pm_sig_corr, method = metric), method = method)
column_hc <- hclust(dist(t(pm_sig_corr), method = metric), method = method)

# 将聚类结果转换为树状图
row_dend <- as.dendrogram(row_hc)
column_dend <- as.dendrogram(column_hc)
#print(rownames(pm_sig_corr))
#构建新的顺序向量
#nm <- c(1:33,37:42,34:36)
#print(labels(row_dend))
#print(labels(column_dend))

#row_dend <- rotate(row_dend, order = nm)
#column_dend <- rotate(column_dend, order = nm)


#print(labels(row_dend))
sample_names <- colnames(pm_sig_corr)
tissues <- sapply(strsplit(sample_names, "\\."), function(x) x[1])

# 提取在矩阵中匹配到的颜色
matched_colors <- palette[tissues]

# 创建小提琴图注释
ha <- rowAnnotation(
  foo = anno_density(
    pm_sig_corr, 
    type = "violin", 
    gp = gpar(fill = matched_colors),
    border = FALSE,
    width = unit(3, "cm"),
  )
)
hf <- rowAnnotation(
  tissue1 = anno_points(
    x = rep(1, length(rownames(pm_sig_corr))),  # 使用常数来创建点
    axis = FALSE,
    gp = gpar(col = matched_colors[rownames(pm_sig_corr)], fill = matched_colors[rownames(pm_sig_corr)], fontsize = 10),  # 设置颜色和填充色
    pch = 16,  # 设置点的形状为圆形
    size = unit(3, "mm"),  # 设置点的大小
    border = FALSE,
    width = unit(0.25, "cm")
  ))

left_annotation <- c(ha, hf)

bottom_annotation <- HeatmapAnnotation(
  Tissue = anno_simple(
    tissues, 
    col = matched_colors[tissues],
    height = unit(0.3, "cm")
  )
)

ht <- Heatmap(
  pm_sig_corr,
  name = "eQTL",
  cluster_rows = row_dend,
  cluster_columns = column_dend,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = colorRamp2(c(0.6,0.8, 1), c("purple","white","firebrick3")),
  column_dend_side = "bottom",
  bottom_annotation = bottom_annotation,
  left_annotation = left_annotation,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(direction = "vertical"),
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (i > j) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = fill))
    } else if (i == j) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = NA))
      formatted_label <- colnames(pm_sig_corr)[j]
      grid.text(
        formatted_label,
        x, 
        y, 
        gp = gpar(col = "black", fontsize = 9),
        just = "left"
      )
    } else {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = NA))
    }
  },
  width = unit(ncol(pm_sig_corr)* 0.4, "cm"),
  height = unit(nrow(pm_sig_corr)* 0.4, "cm")
)
ht
# 绘制热图
#pdf("./breedspe_heatmap_new_maximum.pdf", width = 12, height = 12)
#draw(ht)
#dev.off()

