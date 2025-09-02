library(data.table)
library(ggplot2)
library(ggrastr)
library(readxl)
library(clusterProfiler)
library(dplyr)
library(IOBR)
library(xlsx)
library(stringr)
library(Rtsne)
library(ggpointdensity)


#Preparing the anno data#
anno<-fread("D:/博士/sheepGTEx/comparative genemic/orthologues_gene_bp.txt",data.table=F)
sheep_orth<-anno$`Gene stable ID`[anno$`Cow homology type`=='ortholog_one2one' & anno$`Human homology type`=='ortholog_one2one' & anno$`Pig homology type`=='ortholog_one2one']
sheep_spe<-anno$`Gene stable ID`[anno$`Cow homology type`=='' & anno$`Human homology type`=='' & anno$`Pig homology type`=='']
sheep_complex <- anno$`Gene stable ID`[
  !(anno$`Gene stable ID` %in% sheep_orth) & 
    !(anno$`Gene stable ID` %in% sheep_spe) 
]
anno$`Gene Classification` <- NA  
anno$`Gene Classification` <- ifelse(
  anno$`Gene stable ID` %in% sheep_orth, "1-1-1-1 gene",
  ifelse(anno$`Gene stable ID` %in% sheep_spe, "sheep specific",
         ifelse(anno$`Gene stable ID` %in% sheep_complex, "complex orthologous", NA)))


######Human#
human_exp<-as.data.frame(fread("G:/gtex/raw read/human/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct"))
human_pc<-as.data.frame(fread('G:/gtex/raw read/human/protein_coding.txt')) #grep PCGs
## sample in common tissues
human_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/human_orh_id_filtered.txt"))
colnames_human_exp <- colnames(human_exp)
cols_to_keep <- c(1, 2, which(colnames_human_exp %in% human_orh_id$Sample))
human_exp_filtered <- human_exp[, cols_to_keep]

## remove replicate data due to verson number
human_exp_processed <- human_exp_filtered
human_exp_processed$Name <- substr(human_exp_processed$Name, 1, 15)
human_exp_processed <- human_exp_processed[, -2]
human_exp_unique <- remove_duplicate_genes(
  human_exp_processed,
  column_of_symbol = "Name",  
  method = "mean"
)
#fwrite(human_exp_unique, file = "D:/博士/sheepGTEx/comparative genemic/e2n/human_exp_unique.txt",sep = "\t",quote = FALSE,row.names = TRUE)
human_exp_unique<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/e2n/human_exp_unique.txt"))
rownames(human_exp_unique) <- human_exp_unique[,1]
human_exp_unique <- human_exp_unique[,-1]

human_matrix_pc <- human_exp_unique[rownames(human_exp_unique) %in% human_pc$`Gene stable ID`, ] #PCG
human_matrix <- human_exp_unique[rownames(human_exp_unique) %in% anno$`Human gene stable ID`, ] #ort
#fwrite(human_matrix_pc, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/human_matrix_pc.txt",sep = "\t",quote = FALSE,row.names = TRUE)
#fwrite(human_matrix, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/human_matrix_ort.txt",sep = "\t",quote = FALSE,row.names = TRUE)

rownames_human_mat <- rownames(human_matrix)
human_ortho_genes <- rownames_human_mat[rownames_human_mat %in% anno$`Human gene stable ID`[anno$`Gene Classification` == "1-1-1-1 gene"]]
writeLines(human_ortho_genes, con = "D:/博士/sheepGTEx/comparative genemic/e2n/test/human_ortho_genes.txt")


human_matrix<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/human_matrix_ort.txt"))
rownames(human_matrix) <- human_matrix[,1]
human_matrix <- human_matrix[,-1]
human_annotation <- anno %>% select(`Human gene stable ID`,`Human chromosome/scaffold start (bp)`,`Human chromosome/scaffold end (bp)`,`Gene Classification`,`Gene stable ID`,)
human_annotation<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/human_annotation.txt"))
#human_meta<-as.data.frame(fread("G:/gtex/raw read/human/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"))
tissuetype <- unique(human_orh_id$Tissue)
sample_names <- colnames(human_matrix)
human_orh_id <- human_orh_id[human_orh_id$Sample %in% sample_names, ]
#write.table(human_orh_id_filtered, file = "D:/博士/sheepGTEx/comparative genemic/human_orh_id_filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE)
human_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/human_orh_id_filtered.txt"))

Results <- data.frame(matrix(nrow = nrow(human_matrix), ncol = length(tissuetype) * 3))
col_names <- c()
for(tissue in tissuetype) {
  col_names <- c(col_names, paste0(tissue, "_median"), paste0(tissue, "_mean"), paste0(tissue, "_sd"))
}
colnames(Results) <- col_names
rownames(Results) <- rownames(human_matrix)
human_annotation$length <- human_annotation$`Human chromosome/scaffold end (bp)` - human_annotation$`Human chromosome/scaffold start (bp)` + 1
for (i in 1:length(tissuetype)) {
  sampleID <- human_orh_id$Sample[human_orh_id$Tissue == tissuetype[i]]
  exp_subset <- human_matrix[, sampleID, drop = FALSE]
  Results[, paste0(tissuetype[i], "_median")] <- apply(exp_subset, 1, median)
  Results[, paste0(tissuetype[i], "_mean")] <- apply(exp_subset, 1, mean)
  Results[, paste0(tissuetype[i], "_sd")] <- apply(exp_subset, 1, sd)
}
Human_ortho_rm <- NULL
for(i in 1:length(tissuetype)){
  tissue <- tissuetype[i]
  summary <- data.frame(
    Gene = rownames(Results),
    median = Results[, paste0(tissue, "_median")],
    mean = Results[, paste0(tissue, "_mean")],
    sd = Results[, paste0(tissue, "_sd")],
    Tissue = tissue,
    Species = 'Human'
  )
  summary$length <- human_annotation$length[match(summary$Gene, human_annotation$`Gene stable ID`)]
  summary$orthology <- human_annotation$`Gene Classification`[match(summary$Gene, human_annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='1-1-1-1 gene')] <- '1-1-1-1 gene'
  summary$orthology[which(summary$orthology=='complex orthologous')] <- 'Complex'
  summary$orthology[is.na(summary$orthology)] <- "Human_spec"
  summary$orthology[summary$orthology != '1-1-1-1 gene' & 
                      summary$orthology != 'Complex'] <- 'Human_spec'
  summary$reads <- summary$median * summary$length
  Human_ortho_rm <- rbind(Human_ortho_rm, summary)
}
fwrite(Human_ortho_rm, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/Human_ortho_sd.txt",sep = "\t",quote = FALSE,row.names = TRUE)


#Pig#
pig_exp<-as.data.frame(fread("G:/gtex/raw read/pig/PigGTEx_v0.Gene.tpm.txt"))
pig_pc<-as.data.frame(fread('G:/gtex/raw read/pig/protein_coding.txt')) #grep PCGs
rownames(pig_exp) <- pig_exp[,1]
pig_exp <- pig_exp[,-1]
## sample in common tissues
#pig_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/pig_id.txt"))
pig_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/n2e/pig_orh_id.txt"))
colnames_pig_exp <- colnames(pig_exp)
cols_to_keep <- c(which(colnames_pig_exp %in% pig_orh_id$Sample))
pig_exp_filtered <- pig_exp[, cols_to_keep]

pig_matrix_pc <- pig_exp_filtered[rownames(pig_exp_filtered) %in% pig_pc$`Gene stable ID`, ] #PCG
pig_matrix <- pig_exp_filtered[rownames(pig_exp_filtered) %in% anno$`Pig gene stable ID`, ] #ort
#fwrite(pig_matrix_pc, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_matrix_pc.txt",sep = "\t",quote = FALSE,row.names = TRUE)
#fwrite(pig_matrix, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_matrix_ort.txt",sep = "\t",quote = FALSE,row.names = TRUE)
rownames_pig_mat <- rownames(pig_matrix)
pig_ortho_genes <- rownames_pig_mat[rownames_pig_mat %in% anno$`Pig gene stable ID`[anno$`Gene Classification` == "1-1-1-1 gene"]]
writeLines(pig_ortho_genes, con = "D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_ortho_genes.txt")

pig_matrix<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_matrix_ort.txt"))
rownames(pig_matrix) <- pig_matrix[,1]
pig_matrix <- pig_matrix[,-1]
pig_annotation <- anno %>% select(`Pig gene stable ID`,`Pig chromosome/scaffold start (bp)`,`Pig chromosome/scaffold end (bp)`,`Gene Classification`,`Gene stable ID`)
#pig_annotation<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/pig_annotation.txt"))
#pig_meta <- read_excel("G:/gtex/raw read/pig/PigGTEx_v0.MetaTable.xlsx")
#pig_meta <- as.data.frame(pig_meta)
tissuetype <- unique(pig_orh_id$Tissue)
#tissuetype <- gsub(" ", ".", tissuetype)
sample_names <- colnames(pig_exp_filtered)
pig_orh_id <- pig_orh_id[pig_orh_id$Sample %in% sample_names, ]
write.table(pig_orh_id, file = "D:/博士/sheepGTEx/comparative genemic/e2n/pig_orh_id.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Results <- data.frame(matrix(nrow = nrow(pig_matrix), ncol = length(tissuetype) * 3))
col_names <- c()
for(tissue in tissuetype) {
  col_names <- c(col_names, paste0(tissue, "_median"), paste0(tissue, "_mean"), paste0(tissue, "_sd"))
}
colnames(Results) <- col_names
rownames(Results) <- rownames(pig_matrix)
pig_annotation$length <- pig_annotation$`Pig chromosome/scaffold end (bp)` - pig_annotation$`Pig chromosome/scaffold start (bp)` + 1
for (i in 1:length(tissuetype)) {
  sampleID <- pig_orh_id$Sample[pig_orh_id$Tissue == tissuetype[i]]
  exp_subset <- pig_matrix[, sampleID, drop = FALSE]
  Results[, paste0(tissuetype[i], "_median")] <- apply(exp_subset, 1, median)
  Results[, paste0(tissuetype[i], "_mean")] <- apply(exp_subset, 1, mean)
  Results[, paste0(tissuetype[i], "_sd")] <- apply(exp_subset, 1, sd)
}
pig_ortho_rm <- NULL
for(i in 1:length(tissuetype)){
  tissue <- tissuetype[i]
  summary <- data.frame(
    Gene = rownames(Results),
    median = Results[, paste0(tissue, "_median")],
    mean = Results[, paste0(tissue, "_mean")],
    sd = Results[, paste0(tissue, "_sd")],
    Tissue = tissue,
    Species = 'Pig'
  )
  summary$length <- pig_annotation$length[match(summary$Gene, pig_annotation$`Gene stable ID`)]
  summary$orthology <- pig_annotation$`Gene Classification`[match(summary$Gene, pig_annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='1-1-1-1 gene')] <- '1-1-1-1 gene'
  summary$orthology[which(summary$orthology=='complex orthologous')] <- 'Complex'
  summary$orthology[is.na(summary$orthology)] <- "Pig_spec"
  summary$orthology[summary$orthology != '1-1-1-1 gene' & 
                      summary$orthology != 'Complex'] <- 'Pig_spec'
  summary$reads <- summary$median * summary$length
  pig_ortho_rm <- rbind(pig_ortho_rm, summary)
}
fwrite(pig_ortho_rm, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_ortho_sd.txt",sep = "\t",quote = FALSE,row.names = TRUE)


#cattle#
cattle_exp <- readRDS("G:/gtex/raw read/cattle/cattleTPMmatrix.rds")
cattle_pc<-as.data.frame(fread('G:/gtex/raw read/cattle/protein_coding.txt')) #grep PCGs
cattle_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/n2e/cattle_orh_id_filtered.txt"))
#rownames(cattle_exp) <- cattle_exp[,1]
#cattle_exp <- cattle_exp[,-1]
colnames_cattle_exp <- colnames(cattle_exp)
cols_to_keep <- c(which(colnames_cattle_exp %in% cattle_orh_id$Sample))
cattle_exp_filtered <- cattle_exp[, cols_to_keep]

cattle_matrix_pc <- cattle_exp_filtered[rownames(cattle_exp_filtered) %in% cattle_pc$`Gene stable ID`, ] #PCG
cattle_matrix <- cattle_exp_filtered[rownames(cattle_exp_filtered) %in% anno$`Cow gene stable ID`, ] #ort
#fwrite(cattle_matrix_pc, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_matrix_pc.txt",sep = "\t",quote = FALSE,row.names = TRUE)
#fwrite(cattle_matrix, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_matrix_ort.txt",sep = "\t",quote = FALSE,row.names = TRUE)

rownames_cattle_mat <- rownames(cattle_matrix)
cattle_ortho_genes <- rownames_cattle_mat[rownames_cattle_mat %in% anno$`Cow gene stable ID`[anno$`Gene Classification` == "1-1-1-1 gene"]]
writeLines(cattle_ortho_genes, con = "D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_ortho_genes.txt")

cattle_matrix<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_matrix_ort.txt"))
rownames(cattle_matrix) <- cattle_matrix[,1]
cattle_matrix <- cattle_matrix[,-1]
cattle_annotation <- anno %>% select(`Cow gene stable ID`,`Cow chromosome/scaffold start (bp)`,`Cow chromosome/scaffold end (bp)`,`Gene Classification`,`Gene stable ID`)
#cattle_annotation<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/cattle_annotation.txt"))
#cattle_meta <- read_excel("G:/gtex/raw read/cattle/Metadata_FarmGTEx_cattle_V0.xlsx")
#cattle_meta <- as.data.frame(cattle_meta)

tissuetype <- unique(cattle_orh_id$Tissue)
sample_names <- colnames(cattle_exp_filtered)
cattle_orh_id <- cattle_orh_id[cattle_orh_id$Sample %in% sample_names, ]
#write.table(cattle_orh_id, file = "D:/博士/sheepGTEx/comparative genemic/e2n/cattle_orh_id_filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Results <- data.frame(matrix(nrow = nrow(cattle_matrix), ncol = length(tissuetype) * 3))
col_names <- c()
for(tissue in tissuetype) {
  col_names <- c(col_names, paste0(tissue, "_median"), paste0(tissue, "_mean"), paste0(tissue, "_sd"))
}
colnames(Results) <- col_names
rownames(Results) <- rownames(cattle_matrix)
cattle_annotation$length <- cattle_annotation$`Cow chromosome/scaffold end (bp)` - cattle_annotation$`Cow chromosome/scaffold start (bp)` + 1
for (i in 1:length(tissuetype)) {
  sampleID <- cattle_orh_id$Sample[cattle_orh_id$Tissue == tissuetype[i]]
  exp_subset <- cattle_matrix[, sampleID, drop = FALSE]
  Results[, paste0(tissuetype[i], "_median")] <- apply(exp_subset, 1, median)
  Results[, paste0(tissuetype[i], "_mean")] <- apply(exp_subset, 1, mean)
  Results[, paste0(tissuetype[i], "_sd")] <- apply(exp_subset, 1, sd)
}
cattle_ortho_rm <- NULL
for(i in 1:length(tissuetype)){
  tissue <- tissuetype[i]
  summary <- data.frame(
    Gene = rownames(Results),
    median = Results[, paste0(tissue, "_median")],
    mean = Results[, paste0(tissue, "_mean")],
    sd = Results[, paste0(tissue, "_sd")],
    Tissue = tissue,
    Species = 'Cattle'
  )
  summary$length <- cattle_annotation$length[match(summary$Gene, cattle_annotation$`Gene stable ID`)]
  summary$orthology <- cattle_annotation$`Gene Classification`[match(summary$Gene, cattle_annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='1-1-1-1 gene')] <- '1-1-1-1 gene'
  summary$orthology[which(summary$orthology=='complex orthologous')] <- 'Complex'
  summary$orthology[is.na(summary$orthology)] <- "Cattle_spec"
  summary$orthology[summary$orthology != '1-1-1-1 gene' & 
                      summary$orthology != 'Complex'] <- 'Cattle_spec'
  summary$reads <- summary$median * summary$length
  cattle_ortho_rm <- rbind(cattle_ortho_rm, summary)
}
fwrite(cattle_ortho_rm, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_ortho_sd.txt",sep = "\t",quote = FALSE,row.names = TRUE)

#sheep#
sheep_exp<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/sheep_merged_e2n_tpm.tsv"))
sheep_exp<-remove_duplicate_genes(sheep_exp, column_of_symbol = "Geneid", method = "mean")
#rownames(sheep_exp) <- sheep_exp[,1]
#sheep_exp <- sheep_exp[,-1]
sheep_orh_id <-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/n2e/sheep_orh_id_filtered.txt"))
colnames_sheep_exp <- colnames(sheep_exp)
cols_to_keep <- c(which(colnames_sheep_exp %in% sheep_orh_id$Sample))
sheep_exp_filtered <- sheep_exp[, cols_to_keep]

sheep_matrix_pc <- sheep_exp_filtered
sheep_matrix <- sheep_exp_filtered[rownames(sheep_exp_filtered) %in% anno$`Gene stable ID`, ] #ort
#fwrite(sheep_matrix_pc, file = "D:/博士/sheepGTEx/comparative genemic/e2n/sheep_matrix_pc.txt",sep = "\t",quote = FALSE,row.names = TRUE)
#fwrite(sheep_matrix, file = "D:/博士/sheepGTEx/comparative genemic/e2n/sheep_matrix_ort.txt",sep = "\t",quote = FALSE,row.names = TRUE)

rownames_sheep_mat <- rownames(sheep_matrix)
sheep_ortho_genes <- rownames_sheep_mat[rownames_sheep_mat %in% anno$`Gene stable ID`[anno$`Gene Classification` == "1-1-1-1 gene"]]
writeLines(sheep_ortho_genes, con = "D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_ortho_genes.txt")

sheep_matrix<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/e2n/sheep_matrix_ort.txt"))
rownames(sheep_matrix) <- sheep_matrix[,1]
sheep_matrix <- sheep_matrix[,-1]

sheep_annotation <- anno %>% select(`Gene stable ID`,`Gene start (bp)`,`Gene end (bp)`,`Gene Classification`,`Gene stable ID`)
#sheep_annotation<-as.data.frame(fread("D:/博士/sheepGTEx/comparative genemic/sheep_annotation.txt"))
tissuetype <- unique(sheep_orh_id$Tissue)
#tissuetype <- gsub(" ", ".", tissuetype)
sample_names <- colnames(sheep_matrix_pc)
sheep_orh_id <- sheep_orh_id[sheep_orh_id$Sample %in% sample_names, ]
#write.table(sheep_orh_id, file = "D:/博士/sheepGTEx/comparative genemic/e2n/sheep_orh_id_filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Results <- data.frame(matrix(nrow = nrow(sheep_matrix), ncol = length(tissuetype) * 3))
col_names <- c()
for(tissue in tissuetype) {
  col_names <- c(col_names, paste0(tissue, "_median"), paste0(tissue, "_mean"), paste0(tissue, "_sd"))
}
colnames(Results) <- col_names
rownames(Results) <- rownames(sheep_matrix)
sheep_annotation$length <- sheep_annotation$`Gene end (bp)` - sheep_annotation$`Gene start (bp)` + 1
for (i in 1:length(tissuetype)) {
  sampleID <- sheep_orh_id$Sample[sheep_orh_id$Tissue == tissuetype[i]]
  exp_subset <- sheep_matrix[, sampleID, drop = FALSE]
  Results[, paste0(tissuetype[i], "_median")] <- apply(exp_subset, 1, median)
  Results[, paste0(tissuetype[i], "_mean")] <- apply(exp_subset, 1, mean)
  Results[, paste0(tissuetype[i], "_sd")] <- apply(exp_subset, 1, sd)
}
sheep_ortho_rm <- NULL
for(i in 1:length(tissuetype)){
  tissue <- tissuetype[i]
  summary <- data.frame(
    Gene = rownames(Results),
    median = Results[, paste0(tissue, "_median")],
    mean = Results[, paste0(tissue, "_mean")],
    sd = Results[, paste0(tissue, "_sd")],
    Tissue = tissue,
    Species = 'Sheep'
  )
  summary$length <- sheep_annotation$length[match(summary$Gene, sheep_annotation$`Gene stable ID`)]
  summary$orthology <- sheep_annotation$`Gene Classification`[match(summary$Gene, sheep_annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='1-1-1-1 gene')] <- '1-1-1-1 gene'
  summary$orthology[which(summary$orthology=='complex orthologous')] <- 'Complex'
  summary$orthology[is.na(summary$orthology)] <- "Sheep_spec"
  summary$orthology[summary$orthology != '1-1-1-1 gene' & 
                      summary$orthology != 'Complex'] <- 'Sheep_spec'
  summary$reads <- summary$median * summary$length
  sheep_ortho_rm <- rbind(sheep_ortho_rm, summary)
}
fwrite(sheep_ortho_rm, file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_ortho_sd.txt",sep = "\t",quote = FALSE,row.names = TRUE)

#####expression heatmap of cv######
sheep_data <- fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_ortho_sd.txt", data.table=F)
cattle_data <- fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_ortho_sd.txt", data.table=F)
pig_data <- fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_ortho_sd.txt", data.table=F)
human_data <- fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/human_ortho_sd.txt", data.table=F)

species_ortho <- fread('D:\\博士\\sheepGTEx\\comparative genemic\\e2n\\test\\filtered_orthologues_gene_bp.txt', sep = '\t')

species_list <- list(
  human = list(data = human_data, col = "Human gene stable ID"),
  cattle = list(data = cattle_data, col = "Cow gene stable ID"), 
  pig = list(data = pig_data, col = "Pig gene stable ID"),
  sheep = list(data = sheep_data, col = "Gene stable ID")
)


for (species_name in names(species_list)) {
  species_list[[species_name]]$data$stable_gene <- species_ortho$`Gene stable ID`[
    match(species_list[[species_name]]$data$Gene, 
          species_ortho[[species_list[[species_name]]$col]])
  ]
}


human_data <- species_list$human$data
cattle_data <- species_list$cattle$data  
pig_data <- species_list$pig$data
sheep_data <- species_list$sheep$data


human_data <- human_data[!is.na(human_data$stable_gene), ]
cattle_data <- cattle_data[!is.na(cattle_data$stable_gene), ]
pig_data <- pig_data[!is.na(pig_data$stable_gene), ]
sheep_data <- sheep_data[!is.na(sheep_data$stable_gene), ]


human_data$cv <- human_data$sd / human_data$mean
cattle_data$cv <- cattle_data$sd / cattle_data$mean
pig_data$cv <- pig_data$sd / pig_data$mean
sheep_data$cv <- sheep_data$sd / sheep_data$mean


plot_and_get_cv_correlation <- function(tissue_name, species1_data, species2_data, 
                                        species1_name, species2_name) {
  
  species1_tissue <- species1_data[species1_data$Tissue == tissue_name, ]
  species2_tissue <- species2_data[species2_data$Tissue == tissue_name, ]
  
  
  merged_data <- merge(species1_tissue[, c("stable_gene", "cv")], 
                       species2_tissue[, c("stable_gene", "cv")], 
                       by = "stable_gene")
  colnames(merged_data)[2:3] <- c(species1_name, species2_name)
  
  valid_data <- merged_data[is.finite(merged_data[[species1_name]]) & 
                              is.finite(merged_data[[species2_name]]), ]
  
  
  if(nrow(valid_data) < 10) {
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste0("Insufficient data"),
               size = 4) +
      theme_void() +
      xlim(0, 1) + ylim(0, 1)
    return(list(plot = p, 
                correlation = NA, 
                p_value = NA, 
                sample_size = nrow(valid_data)))
  }
  
  
  corr_result <- tryCatch({
    cor.test(valid_data[[species1_name]], valid_data[[species2_name]], method = "pearson")
  }, error = function(e) {
    return(list(estimate = NA, p.value = NA))
  })
  
  
  p <- ggplot(valid_data, aes_string(x = species1_name, y = species2_name)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    theme_classic() +
    labs(
      x = paste(species1_name, "CV"),
      y = paste(species2_name, "CV")
    )
  
  
  if(!is.na(corr_result$estimate)) {
    p <- p + annotate("text", x = min(valid_data[[species1_name]]), 
                      y = 0.95 * max(valid_data[[species2_name]]),
                      hjust = 0, vjust = 1, size = 3,
                      label = paste("R =", round(corr_result$estimate, 3), 
                                    "\np =", format.pval(corr_result$p.value, digits = 3)))
  }
  
  return(list(plot = p, 
              correlation = ifelse(is.na(corr_result$estimate), NA, corr_result$estimate), 
              p_value = ifelse(is.na(corr_result$p.value), NA, corr_result$p.value),
              sample_size = nrow(valid_data)))
}


all_tissues <- unique(sheep_data$Tissue)
selected_tissues <- all_tissues[1:16]  


dir.create("D:/博士/sheepGTEx/comparative genemic/e2n/test/plots", 
           showWarnings = FALSE, recursive = TRUE)


correlation_results <- data.frame(
  Tissue = character(),
  Pair = character(),
  R = numeric(),
  P = numeric(),
  SampleSize = integer(),
  stringsAsFactors = FALSE
)


for(tissue in selected_tissues) {
  
  result_sh <- plot_and_get_cv_correlation(tissue, sheep_data, human_data, "Sheep", "Human")
  result_sc <- plot_and_get_cv_correlation(tissue, sheep_data, cattle_data, "Sheep", "Cattle")
  result_sp <- plot_and_get_cv_correlation(tissue, sheep_data, pig_data, "Sheep", "Pig")
  result_hc <- plot_and_get_cv_correlation(tissue, human_data, cattle_data, "Human", "Cattle")
  result_hp <- plot_and_get_cv_correlation(tissue, human_data, pig_data, "Human", "Pig")
  result_cp <- plot_and_get_cv_correlation(tissue, cattle_data, pig_data, "Cattle", "Pig")
  
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Sheep-Human", 
                                          R = result_sh$correlation, 
                                          P = result_sh$p_value,
                                          SampleSize = result_sh$sample_size))
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Sheep-Cattle", 
                                          R = result_sc$correlation, 
                                          P = result_sc$p_value,
                                          SampleSize = result_sc$sample_size))
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Sheep-Pig", 
                                          R = result_sp$correlation, 
                                          P = result_sp$p_value,
                                          SampleSize = result_sp$sample_size))
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Human-Cattle", 
                                          R = result_hc$correlation, 
                                          P = result_hc$p_value,
                                          SampleSize = result_hc$sample_size))
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Human-Pig", 
                                          R = result_hp$correlation, 
                                          P = result_hp$p_value,
                                          SampleSize = result_hp$sample_size))
  
  correlation_results <- rbind(correlation_results,
                               data.frame(Tissue = tissue, 
                                          Pair = "Cattle-Pig", 
                                          R = result_cp$correlation, 
                                          P = result_cp$p_value,
                                          SampleSize = result_cp$sample_size))
  
  
  combined_plot <- plot_grid(
    result_sh$plot, result_sc$plot, result_sp$plot,
    ncol = 3,
    labels = c("Sheep-Human", "Sheep-Cattle", "Sheep-Pig"),
    label_size = 12,
    align = "h"
  )
  
  title <- ggdraw() + 
    draw_label(
      paste("Expression Variation (CV) Correlation in", tissue),
      fontface = "bold",
      size = 14,
      x = 0.5
    )
  
  final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
  ggsave(
    filename = paste0("D:/博士/sheepGTEx/comparative genemic/e2n/test/plots/", tissue, "_cv_correlations.pdf"),
    plot = final_plot,
    width = 15,
    height = 5,
    dpi = 300
  )
}

write.csv(correlation_results, "D:/博士/sheepGTEx/comparative genemic/e2n/test/cv_correlation_results.csv", row.names = FALSE)

# ========== cv heatmap ==========

correlation_matrix <- dcast(correlation_results, Tissue ~ Pair, value.var = "R")
rownames(correlation_matrix) <- correlation_matrix$Tissue
correlation_matrix$Tissue <- NULL
correlation_matrix <- correlation_matrix[complete.cases(correlation_matrix), ]
tissue_means <- rowMeans(correlation_matrix, na.rm = TRUE)
correlation_matrix <- correlation_matrix[order(tissue_means, decreasing = TRUE), ]

sheep_comparisons <- correlation_results %>%
  filter(Pair %in% c("Sheep-Human", "Sheep-Cattle", "Sheep-Pig"))


sheep_matrix <- dcast(sheep_comparisons, Tissue ~ Pair, value.var = "R")
rownames(sheep_matrix) <- sheep_matrix$Tissue
sheep_matrix$Tissue <- NULL

sheep_means <- rowMeans(sheep_matrix, na.rm = TRUE)
sheep_matrix <- sheep_matrix[order(sheep_means, decreasing = TRUE), ]
sheep_matrix_t <- t(as.matrix(sheep_matrix))

tis_col<-read.csv("D:/博士/pic(PDF&TIFF)/tissue_color.txt",sep = "\t")
tis_col$Tissue <- gsub("_", " ", tis_col$Tissue)
tissue_colors <- setNames(tis_col$Color, tis_col$Tissue)
missing_tissues <- setdiff(colnames(sheep_matrix_t), names(tissue_colors))
if(length(missing_tissues) > 0) {
  tissue_colors <- c(tissue_colors, setNames(rainbow(length(missing_tissues)), missing_tissues))
}
col_annotation <- data.frame(Tissue = colnames(sheep_matrix_t))
rownames(col_annotation) <- colnames(sheep_matrix_t)
ann_colors <- list(Tissue = tissue_colors)

ordered_colors <- tissue_colors[colnames(sheep_matrix_t)]
col_anno <- HeatmapAnnotation(
  tissue_dots = anno_points(
    rep(1, ncol(sheep_matrix_t)), 
    pch = 16,  
    size = unit(4, "mm"),
    gp = gpar(col = ordered_colors, fill = ordered_colors)
  ),
  which = "column",
  height = unit(8, "mm"),
  show_legend = FALSE,
  show_annotation_name = FALSE
)
ht <- Heatmap(sheep_matrix_t,
              name = "Correlation",
              col = colorRamp2(c(0.35, 0.525, 0.8), c("purple", "white", "red")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = TRUE,
              column_names_rot = 90,  
              bottom_annotation = col_anno,
              heatmap_legend_param = list(
                title = "Correlation",
                direction = "vertical"
              ))

pdf("D:/博士/sheepGTEx/comparative genemic/e2n/test/final/Sheep_CV_correlation_comparison_point_anno.pdf", 
    width = 7, height = 3)
draw(ht)
dev.off()

#########scatter plot of tau########
library(data.table)
library(ggplot2)
library(ggrastr)
library(readxl)
library(clusterProfiler)
library(dplyr)
library(IOBR)
library(xlsx)
library(stringr)
library(Rtsne)
library(ggpointdensity)

####tau cauculation#####
library(reshape2)

compute_tau <- function(x){
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  if(length(x) <= 1 || max(x) == 0) return(0)
  x_norm <- x / max(x)
  sum(1 - x_norm) / (length(x) - 1)
}

# species files: long per-gene-per-tissue median files and output tau files
base_dir <- "D:/博士/sheepGTEx/comparative genemic/e2n/test"
species_files <- list(
  human  = file.path(base_dir, "Human_ortho_sd.txt"),
  pig    = file.path(base_dir, "pig_ortho_sd.txt"),
  cattle = file.path(base_dir, "cattle_ortho_sd.txt"),
  sheep  = file.path(base_dir, "sheep_ortho_sd.txt")
)

for(sp in names(species_files)){
  long_fp <- species_files[[sp]]
  out_fp  <- file.path(base_dir, paste0(tolower(sp), "_tau.tsv"))
  long_df <- fread(long_fp, data.table = FALSE)
  if(ncol(long_df) == 0) next
  if(!"Gene" %in% colnames(long_df)){
    colnames(long_df)[1] <- "Gene"
  }
  
  wide_df <- tryCatch({
    dcast(long_df, Gene ~ Tissue, value.var = "median")
  }, error = function(e){
    stop(sprintf("Failed to dcast %s: %s", long_fp, e$message))
  })
  
  tau_vec <- apply(wide_df[ , -1, drop = FALSE], 1, compute_tau)
  tau_df <- data.frame(Gene = wide_df$Gene, tau = tau_vec, stringsAsFactors = FALSE)
  fwrite(tau_df, file = out_fp, sep = "\t", quote = FALSE, row.names = FALSE)
}
####

human_tspex<-fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/human_tau.tsv",data.table=F)
sheep_tspex<-fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_tau.tsv",data.table=F)
pig_tspex<-fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/pig_tau.tsv",data.table=F)
cattle_tspex<-fread("D:/博士/sheepGTEx/comparative genemic/e2n/test/cattle_tau.tsv",data.table=F)
human_tspex<-human_tspex[-1,]
sheep_tspex<-sheep_tspex[-1,]
pig_tspex<-pig_tspex[-1,]
cattle_tspex<-cattle_tspex[-1,]
colnames(human_tspex)<-c("Gene","tau")
colnames(sheep_tspex)<-c("Gene","tau")
colnames(pig_tspex)<-c("Gene","tau")
colnames(cattle_tspex)<-c("Gene","tau")

species_ortho<-fread('D:\\博士\\sheepGTEx\\comparative genemic\\e2n\\test\\filtered_orthologues_gene_bp.txt',sep = '\t')

#sheep-cattle
sheep_tspex<-sheep_tspex[match(species_ortho$`Gene stable ID`, sheep_tspex$Gene), ]
human_tspex <- human_tspex[match(species_ortho$`Human gene stable ID`, human_tspex$Gene), ]
human_tspex$Gene <- species_ortho$`Gene stable ID`[match(human_tspex$Gene, species_ortho$`Human gene stable ID`)]
cattle_tspex<-cattle_tspex[match(species_ortho$`Cow gene stable ID`, cattle_tspex$Gene), ]
cattle_tspex$Gene <- species_ortho$`Gene stable ID`[match(cattle_tspex$Gene, species_ortho$`Cow gene stable ID`)]
pig_tspex<-pig_tspex[match(species_ortho$`Pig gene stable ID`, pig_tspex$Gene), ]
pig_tspex$Gene <- species_ortho$`Gene stable ID`[match(pig_tspex$Gene, species_ortho$`Pig gene stable ID`)]

sheep_tspex$tau<-as.numeric(sheep_tspex$tau)
human_tspex$tau<-as.numeric(human_tspex$tau)
cattle_tspex$tau<-as.numeric(cattle_tspex$tau)
pig_tspex$tau<-as.numeric(pig_tspex$tau)


corr <- cor.test(sheep_tspex$tau, cattle_tspex$tau, method='pearson')
sheep_cattle_tau <- cbind(sheep_tspex$tau, cattle_tspex$tau)
colnames(sheep_cattle_tau) <- c("Sheep", "Cattle")
sheep_cattle_tau <- as.data.frame(sheep_cattle_tau)

p1 <- ggplot(sheep_cattle_tau, aes(x=Sheep, y=Cattle)) +
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="turbo", name="Density") +
  stat_smooth(color="blue", method="lm") +
  theme_classic() +
  labs(x="tau in sheep", y="tau in cattle") + 
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),
        axis.text = element_text(color = "black", size = unit(7, "pt")),
        legend.title = element_text(size = unit(7, "pt")), 
        legend.text = element_text(size = unit(6, "pt"))) +
  annotate("text", x = 0.1, y=1, hjust = 0, vjust = 1, size = 2,
           label = paste("R:", round(corr$estimate, 3),
                         "\nP:", round(corr$p.value, 5)),
           color = "black")
p1

ggsave(p1,file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_cattle_tau_new.pdf",dpi=300,width=105,height=80,units = 'mm')


corr<-cor.test(sheep_tspex$tau,human_tspex$tau,method='pearson') #cor=0.5124 pval<2.2e-16
sheep_human_tau<-cbind(sheep_tspex$tau,human_tspex$tau)
colnames(sheep_human_tau)<-c("Sheep","Human")
sheep_human_tau<-as.data.frame(sheep_human_tau)
p2<-ggplot(sheep_human_tau,aes(x=Sheep,y=Human))+
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="turbo", name="Density") +
  stat_smooth(color="blue", method="lm") +
  theme_classic() +
  labs(x="tau in sheep", 
       y="tau in human") + 
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),axis.text = element_text(color = "black", size = unit(7, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))+
  annotate("text",x = 0.1,y=1, hjust = 0, vjust = 1, size =2,
           label = paste("R:", round(corr$estimate, 3),
                         "\nP:", round(corr$p.value, 5)),
           color = "black")
p2
ggsave(p2,file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_human_tau_new.pdf",dpi=300,width=105,height=80,units = 'mm')


corr<-cor.test(sheep_tspex$tau,pig_tspex$tau,method='pearson') #cor=0.5124 pval<2.2e-16
sheep_pig_tau<-cbind(sheep_tspex$tau,pig_tspex$tau)
colnames(sheep_pig_tau)<-c("Sheep","Pig")
sheep_pig_tau<-as.data.frame(sheep_pig_tau)
p3<-ggplot(sheep_pig_tau,aes(x=Sheep,y=Pig))+
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="turbo", name="Density") +
  stat_smooth(color="blue", method="lm") +
  theme_classic() +
  labs(x="tau in sheep", 
       y="tau in Pig") + 
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),axis.text = element_text(color = "black", size = unit(7, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))+
  annotate("text", x = 0.1,y=1, hjust = 0, vjust = 1, size =2,
           label = paste("R:", round(corr$estimate, 3),
                         "\nP:", round(corr$p.value, 5)),
           color = "black")
p3
ggsave(p3,file = "D:/博士/sheepGTEx/comparative genemic/e2n/test/sheep_pig_tau_new.pdf",dpi=300,width=105,height=80,units = 'mm')
