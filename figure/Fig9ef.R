library(Mfuzz)
library(readr)
library(tibble)
library(limma)
library(RColorBrewer)
library(tools)

# 设置你的表达文件目录
input_dir <- "../europe.test"
file_list <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# 设置聚类数量
c <- 4

# 设置最小标准差和acore阈值
min_std <- 0.1
min_acore <- 0.5

for (file in file_list) {
  cat("������ 正在处理文件：", file, "\n")
  organ <- file_path_sans_ext(basename(file))  # 组织名

  # 1. 读取数据并处理
  data <- read_tsv(file)
  data <- as.data.frame(data)
  data_without_col1 <- data[, -1]
  rownames(data_without_col1) <- data_without_col1[, 1]
  data_final <- data_without_col1[, -1]

  # 2. 转置并设定时间点标签
  expression_data <- t(data_final)
  time_labels <- c("M6","M6","M6","M6","M6","M6","M6","M6","M6",
                   "M4","M4","M4","M4",
                   "M2","M2","M2","M2",
                   "M1","M1","M1","M1","M1","M1","M1","M1","M1",
                   "M0","M0","M0","M0","M0")
  colnames(expression_data) <- time_labels

  # 3. 表达集构建 + 过滤 + 标准化
  avereps_df <- t(limma::avereps(t(expression_data), ID = colnames(expression_data)))
  gene_tpm <- data.matrix(avereps_df)
  eset <- new("ExpressionSet", exprs = gene_tpm)
  tmp <- filter.std(eset, min.std = min_std)
  gene.s <- standardise(tmp)

  # 4. 聚类
  m <- mestimate(gene.s)
  cl <- mfuzz(gene.s, c = c, m = m)

  # 5. 画图输出 PDF
  color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
  pdf(file.path(input_dir, paste0(organ, "_mfuzz_clusters.pdf")), height = 7, width = 12)
  mfuzz.plot(gene.s, cl,
             mfrow = c(2, 2),
             new.window = FALSE,
             time.labels = colnames(gene.s),
             colo = color.2)
  dev.off()

  # 6. 提取核心基因
  cl.thres <- acore(gene.s, cl = cl, min.acore = min_acore)
  core_dir <- file.path(input_dir, paste0(organ, "_core_genes"))
  dir.create(core_dir, showWarnings = FALSE)

  for (i in seq_along(cl.thres)) {
    core_genes_df <- cl.thres[[i]]
    if (!is.null(core_genes_df) && nrow(core_genes_df) > 0) {
      write.csv(core_genes_df,
                file = file.path(core_dir, paste0("cluster", i, "_core_genes.csv")),
                row.names = TRUE)
    }
  }
}
# 6.1 输出所有聚类中的全部基因信息（含acore值）
all_clusters_dir <- file.path(input_dir, paste0(organ, "_all_clusters"))
dir.create(all_clusters_dir, showWarnings = FALSE)

for (i in 1:c) {
  cluster_genes <- cl$cluster[cl$cluster == i]
  gene_names <- names(cluster_genes)
  acore_values <- cl$membership[gene_names, i]  # 第 i 聚类的 acore 值
  df_cluster <- data.frame(Gene = gene_names, Cluster = i, Acore = acore_values)
  
  write.csv(df_cluster,
            file = file.path(all_clusters_dir, paste0("cluster", i, "_all_genes.csv")),
            row.names = FALSE)
}

