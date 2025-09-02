````R
library(data.table)
library(dplyr)
anno <- fread("sheep_esm_ncbi.txt", data.table=F)
ortho_anno <- fread("sheep_ortho_4_species_112.txt",data.table=F)
anno <- anno[!anno$`NCBI gene (formerly Entrezgene) accession` %in% "",]
ortho_anno <- ortho_anno[ortho_anno$`Gene stable ID` %in% anno$`Gene stable ID`,]
ortho_anno <- ortho_anno[ortho_anno$`Human homology type` %in% "ortholog_one2one" & ortho_anno$`Cow homology type` %in% "ortholog_one2one" & ortho_anno$`Pig homology type` %in% "ortholog_one2one",]

human_tissues <- readRDS("human_tissue_list.rds")
cattle_tissues <- readRDS("cattle_tissue_list.rds")
pig_tissues <- readRDS("pig_tissue_list.rds")
sheep_tissues <- readRDS("sheep_tissue_list.rds")

ances_anno <- fread("four_species_ancestral_new.tsv", data.table=F)
sheep_cattle_ances_anno <- ances_anno[!ances_anno$cattle %in% "N",]
sheep_human_ances_anno <- ances_anno[!ances_anno$human %in% "N",]
sheep_pig_ances_anno <- ances_anno[!ances_anno$pig %in% "N",]

sheep_human_sum <- NULL
for(i in 1:length(sheep_tissues)){
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/nominal/", sheep_tissues[i], ".cis_qtl_pairs.all_chr.txt.gz"), data.table=F)
    human_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/human-cis-eQTL/",human_tissues[i]), data.table=F)
    human_eqtl <- human_eqtl %>%
  mutate(human_key = paste0(chr, "_", variant_pos))

human_eqtl <- human_eqtl %>%
  mutate(human_key = paste0(chr, "_", variant_pos))

sheep_human_ances_anno <- sheep_human_ances_anno %>%
  mutate(human_key = paste0("chr", human_chrom, "_", human_pos))

sheep_eqtl_top <- sheep_eqtl %>%
  group_by(variant_id) %>%
  filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(variant_id, .keep_all = TRUE)

human_eqtl_top <- human_eqtl %>%
  group_by(human_key) %>%
  filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(human_key, .keep_all = TRUE)

shared_eqtls <- inner_join(
  sheep_eqtl_top,
  sheep_human_ances_anno,
  by = c("variant_id" = "Chrom_Pos")
) %>%
inner_join(
  human_eqtl_top,
  by = "human_key",
  suffix = c("_sheep", "_human")
)

final_results <- shared_eqtls %>%
  mutate(
    zscore_sheep = slope_sheep / slope_se_sheep,
    zscore_human = slope_human / slope_se_human
  ) %>%
  select(
    sheep_variant_id = variant_id_sheep,
    human_variant_id = variant_id_human,
    sheep_phenotype_id = phenotype_id,
    human_gene_name = gene_name,
    zscore_sheep,
    zscore_human
  )
 final_results$tissue <- sheep_tissues[i]
 final_results$species <- "Human"
 sheep_human_sum <- rbind(sheep_human_sum, final_results)
 }   

sheep_pig_sum <- NULL
for(i in 1:length(sheep_tissues)){
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/nominal/", sheep_tissues[i], ".cis_qtl_pairs.all_chr.txt.gz"), data.table=F)
    pig_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/pig-cis-eQTL/",pig_tissues[i]), data.table=F)
    
  sheep_eqtl_top <- sheep_eqtl %>%
  group_by(variant_id) %>%
  filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(variant_id, .keep_all = TRUE)

pig_eqtl_top <- pig_eqtl %>%
  mutate(pig_key = str_extract(variant_id, "^[0-9]+_[0-9]+")) %>%
  group_by(pig_key) %>%
  filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(pig_key, .keep_all = TRUE)

ancestral_map_with_pig_key <- sheep_pig_ances_anno %>%
  mutate(pig_key = paste(pig_chrom, pig_pos, sep = "_"))

shared_eqtls_pig <- inner_join(
  sheep_eqtl_top,
  ancestral_map_with_pig_key,
  by = c("variant_id" = "Chrom_Pos")
) %>%
inner_join(
  pig_eqtl_top,
  by = "pig_key",
  suffix = c("_sheep", "_pig")
)

final_results_pig <- shared_eqtls_pig %>%
  mutate(
    zscore_sheep = slope_sheep / slope_se_sheep,
    zscore_pig = slope_pig / slope_se_pig
  ) %>%
  select(
    sheep_variant_id = variant_id_sheep,
    pig_variant_id = variant_id_pig,
    sheep_phenotype_id = phenotype_id_sheep,
    pig_phenotype_id = phenotype_id_pig,
    zscore_sheep,
    zscore_pig
  )
final_results_pig$tissue <- sheep_tissues[i]
final_results_pig$species <- "Pig"
sheep_pig_sum <- rbind(sheep_pig_sum, final_results_pig)
}

sheep_cattle_sum <- NULL
for(i in 1:length(sheep_tissues)){
    sheep_eqtl <- fread(paste0("/home/zhonghaobai/zhonghaoDK/FarmGTEx/sheepGTEx/", sheep_tissues[i], "/tensorqtl/nominal/", sheep_tissues[i], ".cis_qtl_pairs.all_chr.txt.gz"), data.table=F)
    cattle_eqtl <- try( data.table::fread(paste0("/faststorage/project/farmgtex/QTL_result/eQTL/",cattle_tissues[i],"/",cattle_tissues[i], ".cis_qtl_pairs.significant.txt"), data.table=F),
  silent = TRUE
)

if (inherits(cattle_eqtl, "try-error") || nrow(cattle_eqtl) == 0){ 
  message(paste("Skipping tissue:", cattle_tissues[i], " - File not found, unreadable, or empty."))
} else {

  sheep_eqtl_top <- sheep_eqtl %>%
    group_by(variant_id) %>%
    filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(variant_id, .keep_all = TRUE)

  cattle_eqtl_top <- cattle_eqtl %>%
    rename(
      slope = beta_g1,
      slope_se = beta_se_g1,
      pval_nominal = pval_g1,
      phenotype_id = pheno_id
    ) %>%
    group_by(variant_id) %>%
    filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(variant_id, .keep_all = TRUE)

  ancestral_map_with_cattle_key <- sheep_cattle_ances_anno %>%
    mutate(cattle_key = paste(cattle_chrom, cattle_pos, sep = "_"))
  
  shared_eqtls_cattle <- inner_join(
    sheep_eqtl_top,
    ancestral_map_with_cattle_key,
    by = c("variant_id" = "Chrom_Pos")
  ) %>%
  inner_join(
    cattle_eqtl_top,
    by = c("cattle_key" = "variant_id"),
    suffix = c("_sheep", "_cattle")
  )
  
  if (nrow(shared_eqtls_cattle) > 0) {
    final_results_cattle <- shared_eqtls_cattle %>%
      mutate(
        zscore_sheep = slope_sheep / slope_se_sheep,
        zscore_cattle = slope_cattle / slope_se_cattle
      ) %>%
      select(
        sheep_variant_id = variant_id,
        cattle_variant_id = cattle_key,
        sheep_phenotype_id = phenotype_id_sheep,
        cattle_phenotype_id = phenotype_id_cattle,
        zscore_sheep,
        zscore_cattle
      )
  } else {
    final_results_cattle <- data.frame()
  }
}
final_results_cattle$tissue <- sheep_tissues[i]
final_results_cattle$species <- "Cattle"
sheep_cattle_sum <- rbind(sheep_cattle_sum, final_results_cattle)
}

###summary
all_tissue_correlations <- data.frame(
  tissue = character(),
  cor_human = numeric(),
  pval_human = numeric(),
  cor_pig = numeric(),
  pval_pig = numeric(),
  cor_cattle = numeric(),
  pval_cattle = numeric()
)

for(i in 1:length(sheep_tissues)){
    current_tissue <- sheep_tissues[i]
final_results_human <- sheep_human_sum[sheep_human_sum$tissue %in% sheep_tissues[i],]
final_results_cattle <- sheep_cattle_sum[sheep_cattle_sum$tissue %in% sheep_tissues[i],]
final_results_pig <- sheep_pig_sum[sheep_pig_sum$tissue %in% sheep_tissues[i],]

cor_h <- NA; pval_h <- NA
cor_p <- NA; pval_p <- NA
cor_c <- NA; pval_c <- NA

  corr_test_human <- cor.test(
    abs(final_results_human$zscore_sheep), 
    abs(final_results_human$zscore_human)
  )
  cor_h <- corr_test_human$estimate
  pval_h <- corr_test_human$p.value

  corr_test_pig <- cor.test(
    abs(final_results_pig$zscore_sheep), 
    abs(final_results_pig$zscore_pig)
  )
  cor_p <- corr_test_pig$estimate
  pval_p <- corr_test_pig$p.value

  corr_test_cattle <- cor.test(
    abs(final_results_cattle$zscore_sheep), 
    abs(final_results_cattle$zscore_cattle)
  )
  cor_c <- corr_test_cattle$estimate
  pval_c <- corr_test_cattle$p.value

new_row <- data.frame(
  tissue = current_tissue,
  cor_human = cor_h,
  pval_human = pval_h,
  cor_pig = cor_p,
  pval_pig = pval_p,
  cor_cattle = cor_c,
  pval_cattle = pval_c
)

all_tissue_correlations <- rbind(all_tissue_correlations, new_row)
}

sorted_correlations <- all_tissue_correlations %>%
  mutate(mean_cor = rowMeans(select(., cor_human, cor_pig, cor_cattle), na.rm = TRUE)) %>%
  arrange(desc(mean_cor)) %>%
  select(-mean_cor)

heatmap_data <- sorted_correlations %>%
  select(cor_human, cor_pig, cor_cattle)

correlation_matrix <- t(as.matrix(heatmap_data))

rownames(correlation_matrix) <- c("Human", "Pig", "Cattle")
colnames(correlation_matrix) <- sorted_correlations$tissue

tissue_color_map <- tis_col %>%
  select(tissue, color) %>%
  distinct() %>%
  { setNames(.$color, .$tissue) }

ordered_colors <- tissue_color_map[colnames(correlation_matrix)]

col_annotation <- HeatmapAnnotation(
  Tissue = anno_points(
    rep(1, ncol(correlation_matrix)), 
    pch = 16,
    size = unit(4, "mm"),
    gp = gpar(col = ordered_colors)
  ),
  show_annotation_name = FALSE,
  height = unit(8, "mm")
)

color_fun <- colorRamp2(c(-0.6, 0, 0.6), c("purple", "white", "#d90429"))

# Generate the heatmap
ht <- Heatmap(
  correlation_matrix,
  name = "Correlation",
  col = color_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  bottom_annotation = col_annotation,
  na_col = "grey",
  column_names_rot = 90,
  row_names_side = "left",
  heatmap_legend_param = list(title = "Correlation")
)

pdf("sheep_4sp_ances_correlation_heatmap.pdf", width = 10, height = 4)
draw(ht)
dev.off()
````

