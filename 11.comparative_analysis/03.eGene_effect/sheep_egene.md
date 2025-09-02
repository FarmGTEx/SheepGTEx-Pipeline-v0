````R
cattle_tissues$V1 <- gsub(".nominals.2rd.sig.txt", "", cattle_tissues$V1)
human_tissues$V1 <- gsub(".v10.eGenes.txt.gz", "", human_tissues$V1)
pig_tissues$V1 <- gsub(".cis_qtl_pairs.significant.txt","",pig_tissues$V1)

Adipose blood brain Heart Hypothalamus Ileum (small intestines) Jejunum (small intestines)  Liver Lung Muscle Ovary Pituitary spleen testis uterus

cattle_tissues <- cattle_tissues[c(1,43,3,13,14,15,20,21,27,29,31,34,38,40)]
human_tissues <- human_tissues[c(1,50,13,30,16,43,32,33,35,37,39,44,46,48),]
pig_tissues <- pig_tissues[c(1,5,6,13,14,15,19,20,25,27,29,31,33,34),]
sheep_tissues <- sheep_tissues[c(3,5,6,17,19,20,25,27,30,33,38,44,47,52),]

cattle_tissues <- cattle_tissues[-c(4,10)]
human_tissues <- human_tissues[-c(4,10)]
pig_tissues <- pig_tissues[-c(4,10)]
sheep_tissues <- sheep_tissues[-c(4,10)]

anno <- fread("sheep_esm_ncbi.txt", data.table=F)
ortho_anno <- fread("sheep_ortho_4_species_112.txt",data.table=F)
anno <- anno[!anno$`NCBI gene (formerly Entrezgene) accession` %in% "",]
ortho_anno <- ortho_anno[ortho_anno$`Gene stable ID` %in% anno$`Gene stable ID`,]
ortho_anno <- ortho_anno[ortho_anno$`Human homology type` %in% "ortholog_one2one" & ortho_anno$`Cow homology type` %in% "ortholog_one2one" & ortho_anno$`Pig homology type` %in% "ortholog_one2one",]

human_tissues <- readRDS("human_tissue_list.rds")
cattle_tissues <- readRDS("cattle_tissue_list.rds")
pig_tissues <- readRDS("pig_tissue_list.rds")
sheep_tissues <- readRDS("sheep_tissue_list.rds")

#Number summary
tissues <- rep(sheep_tissues, 2)
tissues[1:14] <- paste0(tissues[1:14], "_ortho")
tissues[15:28] <- paste0(tissues[15:28], "_eGenes")
sum <- array(NA, dim=c(28, 4))
rownames(sum) <- tissues
colnames(sum) <- c("Human", "Cattle", "Pig", "Sheep")
sum <- as.data.frame(sum)

library(ggplot2)
library(reshape2)



for(i in 1:14){
    #load data
    human_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/human-cis-eQTL/",human_tissues[i]), data.table=F)
    cattle_eqtl <- fread(paste0("/faststorage/project/farmgtex/QTL_result/eQTL/",cattle_tissues[i],"/",cattle_tissues[i], ".cis_qtl_pairs.significant.txt"), data.table=F)
    pig_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/pig-cis-eQTL/",pig_tissues[i]), data.table=F)
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/permutation/", sheep_tissues[i], ".cis_qtl_fdr0.05.egenes.txt"), data.table=F)
    
    #human      
    human_eqtl$gene_id <- substr(human_eqtl$gene_id, 1, 15)
    human_eqtl <- human_eqtl[human_eqtl$qval <= 0.05,]
    human_eqtl_ortho <- human_eqtl[human_eqtl$gene_id %in% ortho_anno$`Human gene stable ID`,]
    saveRDS(human_eqtl_ortho, file=paste0("human_ortho_", sheep_tissues[i], "_eGene.rds"))
    sum$Human[i] <- nrow(human_eqtl_ortho)
    sum$Human[i+14] <- nrow(human_eqtl)
    
    #cattle
    cattle_eqtl_ortho <- cattle_eqtl[cattle_eqtl$pheno_id %in% ortho_anno$`Cow gene stable ID`,]
    sum$Cattle[i] <- length(unique(cattle_eqtl_ortho$pheno_id))
    sum$Cattle[i+14] <- length(unique(cattle_eqtl$pheno_id))
    saveRDS(cattle_eqtl_ortho, file=paste0("cattle_ortho_", sheep_tissues[i], "_eGene.rds"))
    
    #pig
    pig_eqtl_ortho <- pig_eqtl[pig_eqtl$phenotype_id %in% ortho_anno$`Pig gene stable ID`,]
    sum$Pig[i] <- length(unique(pig_eqtl_ortho$phenotype_id))
    sum$Pig[i+14] <- length(unique(pig_eqtl$phenotype_id))
    saveRDS(pig_eqtl_ortho, file=paste0("pig_ortho_", sheep_tissues[i], "_eGene.rds"))
    
    #sheep
    sheep_eqtl_ortho <- sheep_eqtl[sheep_eqtl$phenotype_id %in% ortho_anno$`Gene name`,]
    sum$Sheep[i] <- length(unique(sheep_eqtl_ortho$phenotype_id))
    sum$Sheep[i+14] <- length(unique(sheep_eqtl$phenotype_id))
    saveRDS(sheep_eqtl_ortho, file=paste0("sheep_ortho_", sheep_tissues[i], "_eGene.rds"))
}
sum <- sum[!sum$Cattle %in% 0,]
saveRDS(sum, file="eGene_number_sum.rds")
sum_ortho <- sum[1:12,]
sum_egene <- sum[13:24,]
sum_ratio <- array(NA, dim=c(12, 4))
rownames(sum_ratio) <- sheep_tissues
colnames(sum_ratio) <- c("Human", "Cattle", "Pig", "Sheep")
sum_ratio <- as.data.frame(sum_ratio)
sum_ratio$Human <- sum_ortho$Human/sum_egene$Human
sum_ratio$Cattle <- sum_ortho$Cattle/sum_egene$Cattle
sum_ratio$Pig <- sum_ortho$Pig/sum_egene$Pig
sum_ratio$Sheep <- sum_ortho$Sheep/sum_egene$Sheep
sum_ratio$Tissue <- rownames(sum_ratio)
sum_ratio_long <- melt(sum_ratio, id.vars = "Tissue", variable.name = "Species", value.name = "Value")

p <- ggplot(sum_ratio_long, aes(x = Species, y = Tissue, fill = Value)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "red", limits = c(0, max(sum_ratio_long$Value))) +  
  geom_text(aes(label = sprintf("%.2f", Value)), size = 5) +  
  theme_minimal() +
  labs(x = NULL,
       y = NULL,
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave(p, file="eGene_ratio_heatmap.pdf", height=5, width=4)

#effect summary
##human
human_sum <- NULL
for(i in 1:12){
    #load data
    human_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/human-cis-eQTL/",human_tissues[i]), data.table=F)
    human_eqtl$gene_id <- substr(human_eqtl$gene_id, 1, 15)
    human_eqtl <- human_eqtl[human_eqtl$qval <= 0.05,]
    human_eqtl$Tissue <- sheep_tissues[i]
    human_eqtl$z <- human_eqtl$slope/human_eqtl$slope_se
    human_eqtl_ortho <- human_eqtl[human_eqtl$gene_id %in% ortho_anno$`Human gene stable ID`,]
    human_eqtl <- human_eqtl[!human_eqtl$gene_id %in% human_eqtl_ortho$gene_id,]
    human_eqtl$cate <- "spe_eGene"
    human_eqtl_ortho$cate <- "ortho_eGene"
    human_tmp <- rbind(human_eqtl[,c("gene_id", "slope", "slope_se", "z", "Tissue", "cate")], human_eqtl_ortho[,c("gene_id", "slope", "slope_se", "z", "Tissue", "cate")])
    human_sum <- rbind(human_sum, human_tmp)
}
saveRDS(human_sum, file="human_ortho_egene_effect.rds")
##cattle

library(dplyr)
cattle_sum <- NULL
for(i in 1:12){
    #load data
    cattle_eqtl <- fread(paste0("/faststorage/project/farmgtex/QTL_result/eQTL/",cattle_tissues[i],"/",cattle_tissues[i], ".cis_qtl_pairs.significant.txt"), data.table=F)
    cattle_eqtl <- cattle_eqtl %>%
  group_by(pheno_id) %>%
  slice_min(pval_g1, n = 1)
    cattle_eqtl <- as.data.frame(cattle_eqtl)
    cattle_eqtl$Tissue <- sheep_tissues[i]
    cattle_eqtl$z <- cattle_eqtl$beta_g1/cattle_eqtl$beta_se_g1
    cattle_eqtl_ortho <- cattle_eqtl[cattle_eqtl$pheno_id %in% ortho_anno$`Cow gene stable ID`,]
    cattle_eqtl <- cattle_eqtl[!cattle_eqtl$pheno_id %in% cattle_eqtl_ortho$pheno_id,]
    cattle_eqtl$cate <- "spe_eGene"
    cattle_eqtl_ortho$cate <- "ortho_eGene"
    cattle_tmp <- rbind(cattle_eqtl[,c("pheno_id", "beta_g1", "beta_se_g1", "z", "Tissue", "cate")], cattle_eqtl_ortho[,c("pheno_id", "beta_g1", "beta_se_g1", "z", "Tissue", "cate")])
    cattle_tmp <- unique(cattle_tmp)
    cattle_sum <- rbind(cattle_sum, cattle_tmp)
}
saveRDS(cattle_sum, file="cattle_ortho_egene_effect.rds")
##pig
pig_sum <- NULL
for(i in 1:12){
    #load data
    pig_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/pig-cis-eQTL/",pig_tissues[i]), data.table=F)
    pig_eqtl <- pig_eqtl %>%
  group_by(phenotype_id) %>%
  slice_min(pval_nominal, n = 1)
    pig_eqtl <- as.data.frame(pig_eqtl)
    pig_eqtl$Tissue <- sheep_tissues[i]
    pig_eqtl$z <- pig_eqtl$slope/pig_eqtl$slope_se
    pig_eqtl_ortho <- pig_eqtl[pig_eqtl$phenotype_id %in% ortho_anno$`Pig gene stable ID`,]
    pig_eqtl <- pig_eqtl[!pig_eqtl$phenotype_id %in% pig_eqtl_ortho$phenotype_id,]
    pig_eqtl$cate <- "spe_eGene"
    pig_eqtl_ortho$cate <- "ortho_eGene"
    pig_tmp <- rbind(pig_eqtl[,c("phenotype_id", "slope", "slope_se", "z", "Tissue", "cate")], pig_eqtl_ortho[,c("phenotype_id", "slope", "slope_se", "z", "Tissue", "cate")])
    pig_sum <- rbind(pig_sum, pig_tmp)
}
saveRDS(pig_sum, file="pig_ortho_egene_effect.rds")
##sheep
sheep_sum <- NULL
for(i in 1:12){
    #load data
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/permutation/", sheep_tissues[i], ".cis_qtl_fdr0.05.egenes.txt"), data.table=F)
    sheep_eqtl <- as.data.frame(sheep_eqtl)
    sheep_eqtl$Tissue <- sheep_tissues[i]
    sheep_eqtl$z <- sheep_eqtl$slope/sheep_eqtl$slope_se
    sheep_eqtl_ortho <- sheep_eqtl[sheep_eqtl$phenotype_id %in% ortho_anno$`Gene name`,]
    sheep_eqtl <- sheep_eqtl[!sheep_eqtl$phenotype_id %in% sheep_eqtl_ortho$phenotype_id,]
    sheep_eqtl$cate <- "spe_eGene"
    sheep_eqtl_ortho$cate <- "ortho_eGene"
    sheep_tmp <- rbind(sheep_eqtl[,c("phenotype_id", "slope", "slope_se", "z", "Tissue", "cate")], sheep_eqtl_ortho[,c("phenotype_id", "slope", "slope_se", "z", "Tissue", "cate")])
    sheep_sum <- rbind(sheep_sum, sheep_tmp)
}
saveRDS(sheep_sum, file="sheep_ortho_egene_effect.rds")

library(ggplot2)
##sheep
sheep_sum[sheep_sum$cate=="eGene",]
sheep_sum$z <- abs(sheep_sum$z)
sheep_sum$logz <- log2(sheep_sum$z)
p <- ggplot(sheep_sum, aes(x = Tissue, y = logz, fill = cate)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) + 
  scale_fill_manual(values = c("spe_eGene" = "lightgrey", "ortho_eGene" = "red")) +
  theme_classic() + 
  labs(x = NULL,
       y = "log2(|Z-score|)",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, file="sheep_eGene_effect.pdf", height=3, width=5)

##pig
pig_sum$z <- abs(pig_sum$z)
pig_sum$logz <- log2(pig_sum$z)

p <- ggplot(pig_sum, aes(x = Tissue, y = logz, fill = cate)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  
  scale_fill_manual(values = c("spe_eGene" = "lightgrey", "ortho_eGene" = "red")) +  
  theme_classic() +  
  labs(x = NULL,
       y = "log2(|Z-score|)",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave(p, file="pig_eGene_effect.pdf", height=3, width=5)

##cattle
cattle_sum$z <- abs(cattle_sum$z)
cattle_sum$logz <- log2(cattle_sum$z)

p <- ggplot(cattle_sum, aes(x = Tissue, y = logz, fill = cate)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  
  scale_fill_manual(values = c("spe_eGene" = "lightgrey", "ortho_eGene" = "red")) +  
  theme_classic() +  
  labs(x = NULL,
       y = "log2(|Z-score|)",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave(p, file="cattle_eGene_effect.pdf", height=3, width=5)

##human
human_sum$z <- abs(human_sum$z)
human_sum$logz <- log2(human_sum$z)

p <- ggplot(human_sum, aes(x = Tissue, y = logz, fill = cate)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  
  scale_fill_manual(values = c("spe_eGene" = "lightgrey", "ortho_eGene" = "red")) +  
  theme_classic() +  
  labs(x = NULL,
       y = "log2(|Z-score|)",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave(p, file="human_eGene_effect.pdf", height=3, width=5)


colnames(cattle_sum) <- colnames(human_sum)
colnames(sheep_sum) <- colnames(human_sum)
colnames(pig_sum) <- colnames(human_sum)
human_sum$species <- "Human"
cattle_sum$species <- "Cattle"
pig_sum$species <- "Pig"
sheep_sum$species <- "Sheep"
pre <- rbind(human_sum, cattle_sum, pig_sum, sheep_sum)

pre$species <- factor(pre$species, levels = c("Cattle","Sheep","Pig","Human"))
p <- ggplot(pre, aes(x = species, y = logz, fill = cate)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  
  scale_fill_manual(values = c("spe_eGene" = "lightgrey", "ortho_eGene" = "red")) +  
  theme_classic() +  
  labs(x = NULL,
       y = "log2(|Z-score|)",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave(p, file="4species_eGene_effect.pdf", height=3, width=5)
````

````R
#TSS
##human
human_sum <- NULL
for(i in 1:12){
    #load data
    human_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/human-cis-eQTL/",human_tissues[i]), data.table=F)
    human_eqtl$gene_id <- substr(human_eqtl$gene_id, 1, 15)
    human_eqtl <- human_eqtl[human_eqtl$qval <= 0.05,]
    human_eqtl$Tissue <- sheep_tissues[i]
    human_eqtl$z <- human_eqtl$slope/human_eqtl$slope_se
    human_tmp <- human_eqtl[,c("gene_id", "tss_distance", "z", "Tissue")]
    human_sum <- rbind(human_sum, human_tmp)
}
saveRDS(human_sum, file="human_TSS_effect_sum.rds")
##cattle
cattle_sum <- NULL
for(i in 1:12){
    #load data
    cattle_eqtl <- fread(paste0("/faststorage/project/farmgtex/QTL_result/eQTL/",cattle_tissues[i],"/",cattle_tissues[i], ".cis_qtl_pairs.significant.txt"), data.table=F)
    cattle_eqtl <- cattle_eqtl %>%
  group_by(pheno_id) %>%
  slice_min(pval_g1, n = 1)
    cattle_eqtl <- as.data.frame(cattle_eqtl)
    cattle_eqtl$Tissue <- sheep_tissues[i]
    cattle_eqtl$z <- cattle_eqtl$beta_g1/cattle_eqtl$beta_se_g1   
    cattle_tmp <- cattle_eqtl[,c("pheno_id", "start_distance", "z", "Tissue")]
    cattle_tmp <- unique(cattle_tmp)
    cattle_sum <- rbind(cattle_sum, cattle_tmp)
}
saveRDS(cattle_sum, file="cattle_TSS_effect_sum.rds")

##pig
pig_sum <- NULL
for(i in 1:12){
    #load data
    pig_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/pig-cis-eQTL/",pig_tissues[i]), data.table=F)
    pig_eqtl <- pig_eqtl %>%
  group_by(phenotype_id) %>%
  slice_min(pval_nominal, n = 1)
    pig_eqtl <- as.data.frame(pig_eqtl)
    pig_eqtl$Tissue <- sheep_tissues[i]
    pig_eqtl$z <- pig_eqtl$slope/pig_eqtl$slope_se
    pig_tmp <- pig_eqtl[,c("phenotype_id", "tss_distance", "z", "Tissue")]
    pig_sum <- rbind(pig_sum, pig_tmp)
}
saveRDS(pig_sum, file="pig_TSS_effect_sum.rds")

##sheep
sheep_sum <- NULL
for(i in 1:12){
    #load data
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/permutation/", sheep_tissues[i], ".cis_qtl_fdr0.05.egenes.txt"), data.table=F)
    sheep_eqtl <- as.data.frame(sheep_eqtl)
    sheep_eqtl$Tissue <- sheep_tissues[i]
    sheep_eqtl$z <- sheep_eqtl$slope/sheep_eqtl$slope_se
    sheep_tmp <- sheep_eqtl[,c("phenotype_id", "start_distance", "z", "Tissue")]
    sheep_sum <- rbind(sheep_sum, sheep_tmp)
}
saveRDS(sheep_sum, file="sheep_TSS_effect_sum.rds")

colnames(cattle_sum) <- colnames(human_sum)
colnames(sheep_sum) <- colnames(human_sum)
colnames(pig_sum) <- colnames(human_sum)
human_sum$species <- "Human"
cattle_sum$species <- "Cattle"
pig_sum$species <- "Pig"
sheep_sum$species <- "Sheep"
tss <- rbind(human_sum, cattle_sum, pig_sum, sheep_sum)

library(dplyr)
library(ggplot2)
library(tidyr)
tss$abs_tss_kb <- abs(tss$tss_distance) / 1000

get_ecdf_x_at_y <- function(data, target_y) {
  x_vals <- sort(data$abs_tss_kb)
  y_vals <- ecdf(x_vals)(x_vals)
  idx <- which(y_vals >= target_y)[1]
  return(x_vals[idx])
}


x_lines <- tss %>%
  group_by(species) %>%
  summarise(
    x_at_0.5 = get_ecdf_x_at_y(cur_data(), 0.5),
    x_at_0.75 = get_ecdf_x_at_y(cur_data(), 0.75)
  ) %>%
  pivot_longer(
    cols = starts_with("x_at_"),
    names_to = "quantile",
    values_to = "x_value"
  ) %>%
  mutate(y_value = as.numeric(sub("x_at_", "", quantile)))

species_colors <- c(
  "Human" = "red",
  "Cattle" = "green",
  "Sheep" = "blue",
  "Pig" = "brown"
)

p <- ggplot(tss, aes(x = abs_tss_kb, color = species)) + 
  stat_ecdf(geom = "step", size = 1.2) +

  geom_segment(data = x_lines,
               aes(x = x_value, xend = x_value,
                   y = 0, yend = y_value,
                   color = species),
               linetype = "dashed",
              size=0.8) +
  geom_segment(data = x_lines,
               aes(x = 0, xend = x_value,
                   y = y_value, yend = y_value,
                   color = species),
               linetype = "dashed",
              size=0.8) +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "|TSS distance| (kb)",
       y = "Cumulative proportion of lead eQTL") +
  theme_classic()
ggsave(p, file="4species_TSS.pdf", height=3, width=4)
````

````R
sheep_h2 <- fread("omiga_h2.single.txt",data.table=F)
sheep_h2 <- sheep_h2[sheep_h2$pheno_id %in% anno$`NCBI gene (formerly Entrezgene) accession`,]
sheep_h2$ensbid <- anno$`Gene stable ID`[match(sheep_h2$pheno_id, anno$`NCBI gene (formerly Entrezgene) accession`)]
sheep_h2 <- sheep_h2[sheep_h2$ensbid %in% ortho_anno$`Gene stable ID`,]

##sheep
sheep_summary <- sheep_h2 %>%
  group_by(ensbid) %>%
  summarise(
    h2_median = median(h2_g1, na.rm = TRUE),
    h2_sd = sd(h2_g1, na.rm = TRUE)
  ) %>%
  mutate(species = "Sheep")
saveRDS(sheep_summary, file="sheep_h2_summary.rds")

##human
library(dplyr)
library(purrr)
library(readr)
library(stringr)

read_human_hsq <- function(tissue) {
  file_path <- paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/h2hum/hsq/", tissue, ".nofilter.hsq")
  df <- read_tsv(file_path, show_col_types = FALSE) %>%
    mutate(
      Tissue = tissue,
      Gene = str_remove(ID, "\\..*")  
    )
  return(df)
}

human_hsq_all <- map_dfr(human_tissues, read_human_hsq)

human_h2_summary <- human_hsq_all %>%
  group_by(Gene) %>%
  summarise(
    h2_median = median(HSQ, na.rm = TRUE),
    h2_sd = sd(HSQ, na.rm = TRUE),
    n_tissues = n()
  ) %>%
  mutate(species = "Human")
saveRDS(human_h2_summary, file="human_h2_summary.rds")

##pig
library(dplyr)

base_dir <- "/faststorage/project/zhonghaoDK/FarmGTEx/h2pig/output_covar_5or10pc10peer/"

all_data <- do.call(rbind, lapply(sheep_tissues, function(tissue_name) {
  tissue_dir <- file.path(base_dir, tissue_name)
  files <- list.files(path = tissue_dir, pattern = "\\.cis_h2\\.txt$", full.names = TRUE)
  tissue_data <- do.call(rbind, lapply(files, function(file) {
    dat <- fread(file, data.table=F)
    dat$tissue <- tissue_name
    return(dat)
  }))
  
  return(tissue_data)
}))

gene_h2 <- all_data %>%
  group_by(geneid) %>%
  summarise(median_h2 = median(h2, na.rm = TRUE))

saveRDS(gene_h2, file="pig_h2_summary.rds")

##cattle
library(dplyr)

base_dir <- "/faststorage/project/farmgtex/QTL_result/cis_her/"

all_data <- do.call(rbind, lapply(cattle_tissues, function(tissue_name) {
  tissue_dir <- file.path(base_dir, tissue_name)
  file_name <- paste0(tissue_name, ".her_est.txt.gz")
  file_path <- file.path(tissue_dir, file_name)

  dat <- fread(file_path, data.table=F)
  
  dat$tissue <- tissue_name
  
  return(dat)
}))

gene_h2 <- all_data %>%
  group_by(pheno_id) %>%
  summarise(median_h2 = median(h2_g1, na.rm = TRUE))
saveRDS(gene_h2, file="cattle_h2_summary.rds")

####################aFC
##human

human_sum <- NULL
for(i in 1:length(human_tissues)){
  human_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/human-cis-eQTL/", human_tissues[i]), data.table = FALSE)
  human_eqtl$gene_id <- substr(human_eqtl$gene_id, 1, 15)
  human_eqtl <- human_eqtl[human_eqtl$qval <= 0.05, ]
  
  gene_afc_median <- aggregate(afc ~ gene_id, data = human_eqtl, FUN = median)
  gene_afc_median$tissue <- human_tissues[i]
  human_sum <- rbind(human_sum, gene_afc_median)
}

final_human_afc <- aggregate(afc ~ gene_id, data = human_sum, FUN = median)
final_human_afc$afc <- abs(final_human_afc$afc)

##pig
pig_tissues <- gsub(".cis_qtl_pairs.significant.txt","",pig_tissues)

pig_sum <- NULL
for(i in 1:length(pig_tissues)){
  pig_afc <- fread(paste0("/faststorage/project/zhonghaoDK/FarmGTEx/aFC/pig_aFC/", pig_tissues[i], ".log2aFC.txt.gz"), data.table = FALSE)
  pig_afc <- pig_afc[pig_afc$is_eGene == TRUE, ]
  gene_median <- aggregate(log2_aFC ~ phenotype_id, data = pig_afc, FUN = median)
  gene_median$tissue <- pig_tissues[i]
  pig_sum <- rbind(pig_sum, gene_median)
}

final_pig_afc <- aggregate(log2_aFC ~ phenotype_id, data = pig_sum, FUN = median)
final_pig_afc$afc <- abs(final_pig_afc$afc)

##sheep
sheep_sum <- NULL
for(i in 1:length(sheep_tissues)){
  sheep_afc <- fread(paste0("/faststorage/project/zhonghaoDK/FarmGTEx/sheepGTEx/", 
                              sheep_tissues[i], 
                              "/tensorqtl/permutation/", 
                              sheep_tissues[i], ".log2aFC.txt"), 
                      data.table = FALSE)
  sheep_afc <- sheep_afc[sheep_afc$is_eGene == TRUE, ]
  gene_median <- aggregate(log2_aFC ~ phenotype_id, data = sheep_afc, FUN = median)
  gene_median$tissue <- sheep_tissues[i]
  sheep_sum <- rbind(sheep_sum, gene_median)
}
final_sheep_afc <- aggregate(log2_aFC ~ phenotype_id, data = sheep_sum, FUN = median)
final_sheep_afc$afc <- abs(final_sheep_afc$afc)

##cattle
library(data.table)

cattle_sum <- NULL

for(tissue in cattle_tissues) {
  for(chr in c(1:29)) {
    cattle_afc <- paste0("/faststorage/project/farmgtex/QTL_result/eQTL_fold_change/", 
                        tissue, "/", tissue, ".Chr", chr, ".log2aFC.txt")

    cattle_afc <- fread(cattle_afc, data.table = FALSE)
    cattle_afc$tissue <- tissue
    cattle_sum <- rbind(cattle_sum, cattle_afc[, c("pheno_id", "log2_aFC", "qval_g1","tissue")])
  }
}
cattle_sum <- cattle_sum[cattle_sum$qval_g1 <= 0.05, ]
final_cattle_afc <- aggregate(log2_aFC ~ pheno_id, data = cattle_sum, FUN = median)
final_cattle_afc$afc <- abs(final_cattle_afc$afc)


````

````R
##h2 plot
human_h2 <- readRDS("human_h2_summary.rds")
pig_h2 <- readRDS("pig_h2_summary.rds")
cattle_h2 <- readRDS("cattle_h2_summary.rds")
sheep_h2 <- readRDS("sheep_h2_summary.rds")

human_h2 <- human_h2 %>% rename(`Human gene stable ID` = Gene, human_h2 = h2_median)
sheep_h2 <- sheep_h2 %>% rename(`Gene stable ID` = Gene, sheep_h2 = h2_median)
cattle_h2 <- cattle_h2 %>% rename(`Cow gene stable ID` = Gene, cattle_h2 = h2_median)
pig_h2 <- pig_h2 %>% rename(`Pig gene stable ID` = Gene, pig_h2 = h2_median)

combined <- ortho_anno %>%
  inner_join(human_h2, by = "Human gene stable ID") %>%
  inner_join(sheep_h2, by = "Gene stable ID") %>%
  inner_join(cattle_h2, by = "Cow gene stable ID") %>%
  inner_join(pig_h2, by = "Pig gene stable ID")

final_table <- combined %>%
  select(ortho_gene = `Gene stable ID`, human_h2, sheep_h2, cattle_h2, pig_h2)
h2_sum <- final_table

plot_scatter <- function(data, xvar, yvar, xlab = xvar, ylab = yvar, title = NULL) {
  res <- cor.test(data[[xvar]], data[[yvar]])
  r_val <- round(res$estimate, 2)
  p_val <- formatC(res$p.value, format = "e", digits = 2)
  
  if (is.null(title)) {
    title <- paste(xlab, "vs", ylab, ": r =", r_val, ", p =", p_val)
  } else {
    title <- paste(title, "\nr =", r_val, ", p =", p_val)
  }
  
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic()+
    theme(plot.title = element_text(size = 6))
}

h2_sum$human_h2 <- as.numeric(as.character(h2_sum$human_h2))
h2_sum$pig_h2 <- as.numeric(as.character(h2_sum$pig_h2))
h2_sum$sheep_h2 <- as.numeric(as.character(h2_sum$sheep_h2))
h2_sum$cattle_h2 <- as.numeric(as.character(h2_sum$cattle_h2))

ggplot(h2_sum, aes(x = human_h2, y = sheep_h2)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = "Human h2", y = "Sheep h2") +
    theme_classic()

p1 <- plot_scatter(h2_sum, "human_h2", "sheep_h2", xlab = "Human h2", ylab = "Sheep h2")
ggsave(p1, file="human_sheep_h2.pdf", width=3, height=3)
p1 <- plot_scatter(h2_sum, "human_h2", "cattle_h2", xlab = "Human h2", ylab = "Cattle h2")
ggsave(p1, file="human_cattle_h2.pdf", width=3, height=3)
p1 <- plot_scatter(h2_sum, "human_h2", "pig_h2", xlab = "Human h2", ylab = "Pig h2")
ggsave(p1, file="human_pig_h2.pdf", width=3, height=3)
p1 <- plot_scatter(h2_sum, "cattle_h2", "sheep_h2", xlab = "Cattle h2", ylab = "Sheep h2")
ggsave(p1, file="cattle_sheep_h2.pdf", width=3, height=3)
p1 <- plot_scatter(h2_sum, "cattle_h2", "pig_h2", xlab = "Cattle h2", ylab = "Pig h2")
ggsave(p1, file="cattle_pig_h2.pdf", width=3, height=3)
p1 <- plot_scatter(h2_sum, "pig_h2", "sheep_h2", xlab = "Pig h2", ylab = "Sheep h2")
ggsave(p1, file="pig_sheep_h2.pdf", width=3, height=3)

###########aFC plot
human_afc <- readRDS("human_aFC.rds")
pig_afc <- readRDS("pig_aFC.rds")
cattle_afc <- readRDS("cattle_aFC.rds")
sheep_afc <- readRDS("sheep_aFC.rds")
colnames(human_afc) <- colnames(pig_afc) <- colnames(cattle_afc) <- colnames(sheep_afc) <- c("Gene", "afc")
sheep_afc <- sheep_afc %>% 
  left_join(anno, by = c("Gene" = "NCBI gene (formerly Entrezgene) accession")) %>%
  select(`Gene stable ID`, afc) %>%
  rename(sheep_afc = afc)
sheep_afc <- na.omit(sheep_afc)

colnames(human_afc) <- colnames(pig_afc) <- colnames(cattle_afc) <- colnames(sheep_afc) <- c("Gene", "afc")
human_afc <- human_afc %>% rename(`Human gene stable ID` = Gene, human_afc = afc)
sheep_afc <- sheep_afc %>% rename(`Gene stable ID` = Gene, sheep_afc = afc)
cattle_afc <- cattle_afc %>% rename(`Cow gene stable ID` = Gene, cattle_afc = afc)
pig_afc <- pig_afc %>% rename(`Pig gene stable ID` = Gene, pig_afc = afc)
combined_afc <- ortho_anno %>%
  inner_join(human_afc, by = "Human gene stable ID") %>%
  inner_join(sheep_afc, by = "Gene stable ID") %>%
  inner_join(cattle_afc, by = "Cow gene stable ID") %>%
  inner_join(pig_afc, by = "Pig gene stable ID")
afc_sum <- combined_afc %>%
  select(ortho_gene = `Gene stable ID`, human_afc, sheep_afc, cattle_afc, pig_afc)

plot_scatter <- function(data, xvar, yvar, xlab = xvar, ylab = yvar, title = NULL) {
  res <- cor.test(data[[xvar]], data[[yvar]])
  r_val <- round(res$estimate, 2)
  p_val <- formatC(res$p.value, format = "e", digits = 2)
  
  if (is.null(title)) {
    title <- paste(xlab, "vs", ylab, ": r =", r_val, ", p =", p_val)
  } else {
    title <- paste(title, "\nr =", r_val, ",\np =", p_val)
  }
  
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic()+
    theme(plot.title = element_text(size = 6))
}

p1 <- plot_scatter(afc_sum, "human_afc", "sheep_afc", xlab = "Human afc", ylab = "Sheep afc")
ggsave(p1, file="human_sheep_afc.pdf", width=3, height=3)

p1 <- plot_scatter(afc_sum, "human_afc", "cattle_afc", xlab = "Human afc", ylab = "Cattle afc")
ggsave(p1, file="human_cattle_afc.pdf", width=3, height=3)

p1 <- plot_scatter(afc_sum, "human_afc", "pig_afc", xlab = "Human afc", ylab = "Pig afc")
ggsave(p1, file="human_pig_afc.pdf", width=3, height=3)

p1 <- plot_scatter(afc_sum, "cattle_afc", "sheep_afc", xlab = "Cattle afc", ylab = "Sheep afc")
ggsave(p1, file="cattle_sheep_afc.pdf", width=3, height=3)

p1 <- plot_scatter(afc_sum, "cattle_afc", "pig_afc", xlab = "Cattle afc", ylab = "Pig afc")
ggsave(p1, file="cattle_pig_afc.pdf", width=3, height=3)

p1 <- plot_scatter(afc_sum, "pig_afc", "sheep_afc", xlab = "Pig afc", ylab = "Sheep afc")
ggsave(p1, file="pig_sheep_afc.pdf", width=3, height=3)


####afc plot for tissues between human and other species
plot_scatter <- function(data, xvar, yvar, xlab = xvar, ylab = yvar, title = NULL) {
    
  res <- cor.test(data[[xvar]], data[[yvar]])
  r_val <- round(res$estimate, 2)
  p_val <- formatC(res$p.value, format = "e", digits = 2)
  label_text <- paste0(p_tissue, "\nr = ", r_val, "\np = ", p_val)
  
  p <- ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 6))
  
  x_range <- range(data[[xvar]], na.rm = TRUE)
  y_range <- range(data[[yvar]], na.rm = TRUE)
  xpos <- x_range[2]  
  ypos <- y_range[2]  
  
  p + annotate("text", x = xpos, y = ypos, label = label_text,
               hjust = 1, vjust = 1, size = 4)
}

####human vs pig
filtered_ortho <- subset(ortho_anno, 
                         `Human homology type` == "ortholog_one2one" & 
                         `Pig homology type`   == "ortholog_one2one")

for(i in seq_along(human_tissues)) {
  
  h_tissue <- human_tissues[i]
  p_tissue <- pig_tissues[i]
  
  human_sub <- subset(human_sum, tissue == h_tissue)
  pig_sub   <- subset(pig_sum, tissue == p_tissue)
  
  human_ortho <- merge(human_sub, filtered_ortho, 
                       by.x = "gene_id", by.y = "Human gene stable ID")
  
 merge_data <- merge(human_ortho, pig_sub,
                      by.x = "Pig gene stable ID", by.y = "phenotype_id",
                      suffixes = c(".human", ".pig"))
  
  if(nrow(merge_data) > 0) {
    p <- plot_scatter(merge_data, xvar = "afc", yvar = "log2_aFC",
                  xlab = "|log2aFC| in Human",
                  ylab = "|log2aFC| in Pig",
                  title = NULL)
    ggsave(p, file=paste0("aFC_human_", h_tissue, "_pig_", p_tissue, ".pdf"), width = 4, height = 4)
  } else {
    message("tissue pair: human ", h_tissue, " and pig", p_tissue, " cannot find any orthologous genes")
  }
}

###human vs sheep
plot_scatter <- function(data, xvar, yvar, xlab = xvar, ylab = yvar, title = NULL) {
    
  res <- cor.test(data[[xvar]], data[[yvar]])
  r_val <- round(res$estimate, 2)
  p_val <- formatC(res$p.value, format = "e", digits = 2)
  label_text <- paste0(s_tissue, "\nr = ", r_val, "\np = ", p_val)
  
  p <- ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 6))
  
  x_range <- range(data[[xvar]], na.rm = TRUE)
  y_range <- range(data[[yvar]], na.rm = TRUE)
  xpos <- x_range[2]  
  ypos <- y_range[2]  
  
  p + annotate("text", x = xpos, y = ypos, label = label_text,
               hjust = 1, vjust = 1, size = 4)
}

human_ortho <- subset(ortho_anno, `Human homology type` == "ortholog_one2one")
colnames(anno)[colnames(anno) == "NCBI gene (formerly Entrezgene) accession"] <- "sheep_gene_symbol"
sheep_sum2 <- merge(sheep_sum, anno, 
                    by.x = "phenotype_id", by.y = "sheep_gene_symbol", 
                    all.x = TRUE)
names(sheep_sum2)[names(sheep_sum2) == "Gene stable ID"] <- "sheep_gene_id"
human_merge <- merge(human_sum, human_ortho, 
                     by.x = "gene_id", by.y = "Human gene stable ID",
                     all.x = TRUE)
names(human_merge)[names(human_merge) == "Gene stable ID"] <- "sheep_gene_id_ref"
for(i in seq_along(human_tissues)) {
  h_tissue <- human_tissues[i]
  s_tissue <- sheep_tissues[i]
  
  human_sub <- subset(human_merge, tissue == h_tissue)
  sheep_sub <- subset(sheep_sum2, tissue == s_tissue)
  
  merged_data <- merge(human_sub, sheep_sub, 
                       by.x = "sheep_gene_id_ref", by.y = "sheep_gene_id",
                       suffixes = c(".human", ".sheep"))

  if(nrow(merge_data) > 0) {
    p <- plot_scatter(merged_data, xvar = "afc", yvar = "log2_aFC",
                      xlab = "|log2aFC| in Human", 
                      ylab = "|log2aFC| in Sheep",
                      title = paste(h_tissue, "vs", s_tissue))
    ggsave(p, file=paste0("aFC_human_", h_tissue, "_sheep_", s_tissue, ".pdf"), width = 4, height = 4)
  } else {
    message("tissue pair: human ", h_tissue, " and sheep", s_tissue, " cannot find any orthologous genes")
  }
}

####human vs cattle
plot_scatter <- function(data, xvar, yvar, xlab = xvar, ylab = yvar, title = NULL) {
    
  res <- cor.test(data[[xvar]], data[[yvar]])
  r_val <- round(res$estimate, 2)
  p_val <- formatC(res$p.value, format = "e", digits = 2)
  label_text <- paste0(c_tissue, "\nr = ", r_val, "\np = ", p_val)
  
  p <- ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(plot.title = element_text(size = 6))
  
  x_range <- range(data[[xvar]], na.rm = TRUE)
  y_range <- range(data[[yvar]], na.rm = TRUE)
  xpos <- x_range[2]  
  ypos <- y_range[2]  
  
  p + annotate("text", x = xpos, y = ypos, label = label_text,
               hjust = 1, vjust = 1, size = 4)
}


ortho_cow <- subset(ortho_anno, 
                    `Human homology type` == "ortholog_one2one" & 
                    `Cow homology type`   == "ortholog_one2one")

human_merge_cow <- merge(human_sum, ortho_cow, 
                         by.x = "gene_id", by.y = "Human gene stable ID",
                         all.x = TRUE)
names(human_merge_cow)[names(human_merge_cow) == "Cow gene stable ID"] <- "cow_gene_id_ref"

for(i in seq_along(human_tissues)) {
  h_tissue <- human_tissues[i]
  c_tissue <- cattle_tissues[i]

  human_sub <- subset(human_merge_cow, tissue == h_tissue)
  cattle_sub <- subset(cattle_sum, tissue == c_tissue)

  merged_data <- merge(human_sub, cattle_sub, 
                       by.x = "cow_gene_id_ref", by.y = "pheno_id",
                       suffixes = c(".human", ".cow"))
  
  if(nrow(merged_data) > 0) {
    p <- plot_scatter(merged_data, 
                      xvar = "afc", 
                      yvar = "log2_aFC",
                      xlab = "|log2aFC| in Human",
                      ylab = "|log2aFC| in Cattle",
                      title = paste(h_tissue, "vs", c_tissue))
     ggsave(p, file=paste0("aFC_human_", h_tissue, "_cattle_", c_tissue, ".pdf"), width = 4, height = 4)
  } else {
    message("tissue pair: human ", h_tissue, " and sheep", s_tissue, " cannot find any orthologous genes")
  }
}
````

