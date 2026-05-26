library(fossil)


PCGlnc_pm <- (read.table("./07.expression/median_eqtl_tpm.txt", header = TRUE, sep = "\t", row.names = 1))
PCGlnc_pm_clean <- na.omit(PCGlnc_pm)
nrow(PCGlnc_pm_clean)
PCGlnc_pm <- t(PCGlnc_pm_clean)
PCGlnc_pm <- PCGlnc_pm[order(rownames(PCGlnc_pm)), ]


splicing_pm <- (read.table("./07.expression/median_sqtl_nonorm_counts.txt", header = TRUE, sep = "\t", row.names = 1))
splicing_pm_clean <- na.omit(splicing_pm)
splicing_pm <- t(splicing_pm_clean)
splicing_pm <- splicing_pm[order(rownames(splicing_pm)), ]

exon_pm <- (read.table("./07.expression/median_eeqtl_tpm.txt", header = TRUE, sep = "\t", row.names = 1))
exon_pm_clean <- na.omit(exon_pm)
nrow(exon_pm_clean)
exon_pm <- t(exon_pm_clean)
exon_pm <- exon_pm[order(rownames(exon_pm)), ]

isoform_pm <- (read.table("./07.expression/median_isoqtl_tpm.txt", header = TRUE, sep = "\t", row.names = 1))
isoform_pm_clean <- na.omit(isoform_pm)
isoform_pm <- t(isoform_pm_clean)
isoform_pm <- isoform_pm[order(rownames(isoform_pm)), ]



stability_pm <- (read.table("./07.expression/median_stqtl_tmm.txt", header = TRUE, sep = "\t", row.names = 1))
stability_pm_clean <- na.omit(stability_pm)
nrow(stability_pm_clean)
stability_pm <- t(stability_pm_clean)
stability_pm <- stability_pm[order(rownames(stability_pm)), ]



UTR_pm <- (read.table("./07.expression/median_aqtl_norm.txt", header = TRUE, sep = "\t", row.names = 1))
UTR_pm_clean <- na.omit(UTR_pm)
UTR_pm <- t(UTR_pm_clean)
UTR_pm <- UTR_pm[order(rownames(UTR_pm)), ]

enhancer_pm <- (read.table("./07.expression/median_enqtl_tpm.txt", header = TRUE, sep = "\t", row.names = 1))
enhancer_pm_clean <- na.omit(enhancer_pm)
enhancer_pm <- t(enhancer_pm_clean)
enhancer_pm <- enhancer_pm[order(rownames(enhancer_pm)), ]


# 定义QTL类型
qtl_types <- c("eQTL", "sQTL", "eeQTL", "isoQTL", "stQTL", "3aQTL", "enQTL")
QTL_names <- c("eQTL", "sQTL", "eeQTL", "isoQTL", "stQTL", "3aQTL", "enQTL")
# 创建空列表来存储pm和lfsr数据
QTL_list <- list()
QTL_pm_sig <- list()

# 循环处理每种QTL类型
for (qtl in qtl_types) {
  # 读取pm和lfsr文件
  pm_file <- sprintf("./%02d.top_pair_%s/pm_m.s.RDS", which(qtl_types == qtl), qtl)
  lfsr_file <- sprintf("./%02d.top_pair_%s/lfsr_m.s.RDS", which(qtl_types == qtl), qtl)
  
  pm <- readRDS(pm_file)
  lfsr <- readRDS(lfsr_file)
  
  # 处理数据
  pm_sig <- pm[apply(lfsr < 0.05, 1, any), ]
  pm_sig_t <- t(pm_sig)
  pm_sig_t <- pm_sig_t[order(rownames(pm_sig_t)), ]
  # 存储结果
  QTL_list[[qtl]] <- pm
  QTL_pm_sig[[qtl]] <- pm_sig_t
}

# 如果需要，可以单独命名每个QTL的pm
list2env(QTL_list, envir = .GlobalEnv)



QTL_names <- c("PCGlnc", "splicing", "exon", "isoform", "stability", "enhancer", "3UTR","eQTL", "sQTL", "eeQTL", "isoQTL", "stQTL", "enQTL", "3aQTL")

for (k in 2:20) {
  rand_results <- matrix(0, nrow = length(QTL_names), ncol = length(QTL_names))  
  
  for (i in 1:length(QTL_names)) {
    clustering_results_i <- readRDS(paste0("./kmean/", k, "/clustering_results_new_", QTL_names[i], ".RDS"))
    
    for (j in i:length(QTL_names)) {
      clustering_results_j <- readRDS(paste0("./kmean/", k, "/clustering_results_new_", QTL_names[j], ".RDS"))
      
      # 对 aQTL 进行特殊处理
      if (QTL_names[i] == "3aQTL" || QTL_names[j] == "3aQTL") {
        common_names <- intersect(names(clustering_results_i), names(clustering_results_j))
        clustering_results_i <- clustering_results_i[common_names]
        clustering_results_j <- clustering_results_j[common_names]
      }
      
      tryCatch({
        # 计算 Rand Index
        if (i == j) {
          rand_result <- 1  # 对角线元素为 1
        } else {
          rand_result <- rand.index(clustering_results_i, clustering_results_j)
        }
        
        rand_results[i, j] <- rand_result
        rand_results[j, i] <- rand_result
      }, error = function(e) {
        cat("计算 Rand Index 时出错: ", QTL_names[i], "和", QTL_names[j], "\n")
        cat("错误信息: ", e$message, "\n")
      })
    }
  }
  
  colnames(rand_results) <- QTL_names
  rownames(rand_results) <- QTL_names
  write.csv(rand_results, file = paste0("./kmean/", k, "/rand_index_results_new_k_", k, ".csv"), row.names = TRUE)
}

