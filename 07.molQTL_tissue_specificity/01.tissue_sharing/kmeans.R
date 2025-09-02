library(stats)

PCGlnc_pm <- t(read.table("/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/sheep.PCGlnc.gene.merged.tpm.txt", header = TRUE, sep = "\t", row.names = 1))
#sQTL_pm <- t(readRDS("./02.top_pair_sQTL/pm_m.s.RDS"))
#eeQTL_pm <- t(readRDS("./03.top_pair_eeQTL/pm_m.s.RDS"))
#isoQTL_pm <- t(readRDS("./04.top_pair_isoQTL/pm_m.s.RDS"))
#stQTL_pm <- t(readRDS("./05.top_pair_stQTL/pm_m.s.RDS"))
#aQTL_pm <- t(readRDS("./06.top_pair_3aQTL/pm_m.s.RDS"))

QTL_list <- list(PCGlnc_pm)
QTL_names <- c("PCGlnc")

for (i in 1:length(QTL_list)) {
  for (k in 2:20) {
    # 创建保存结果的目录
    output_dir <- paste0("/storage/public/home/2020060185/00.sheep_goatGTEx/01.sheepGTEx/03.QTL/01.eQTL/v6.tissue_sharing/kmeans", k, "/")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # 运行 K-means 聚类
    clustering_result <- kmeans(QTL_list[[i]], centers = k)
    
    # 保存聚类结果到文件
    saveRDS(clustering_result$cluster, paste0(output_dir, "clustering_results_", QTL_names[i], ".RDS"))
  }
}
