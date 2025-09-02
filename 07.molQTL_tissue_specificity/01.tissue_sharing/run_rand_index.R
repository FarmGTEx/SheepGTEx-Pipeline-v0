library(stats)
library(fossil)


QTL_names <- c("eQTL", "sQTL", "eeQTL", "isoQTL", "stQTL", "aQTL")

# 计算 k=2 的值
k <- 2
rand_results <- matrix(0, nrow = length(QTL_names), ncol = length(QTL_names))  # 存储QTL矩阵对之间的Rand指数

for (i in 1:(length(QTL_names) - 1)) {
  clustering_results_i <- readRDS(paste0("/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/kmeans/", k, "/clustering_results_", QTL_names[i], ".RDS"))
  
  for (j in (i + 1):length(QTL_names)) {
    clustering_results_j <- readRDS(paste0("/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/kmeans/", k, "/clustering_results_", QTL_names[j], ".RDS"))
    
    # 计算 Rand Index
    rand_result <- rand.index(clustering_results_i, clustering_results_j)
    
    # 存入结果矩阵
    rand_results[i, j] <- rand_result
    rand_results[j, i] <- rand_result  # 对称矩阵
    
    # 清除不再使用的对象
    rm(clustering_results_j)
    gc()
  }
  
  # 清除不再使用的对象
  rm(clustering_results_i)
  gc()
}

# 输出 Rand 指数矩阵到文件
colnames(rand_results) <- QTL_names
rownames(rand_results) <- QTL_names
write.csv(rand_results, file = paste0("/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing./kmean/", k, "/rand_index_results_k_", k, ".csv"), row.names = TRUE)

# 清除缓存
rm(rand_results)
gc()
