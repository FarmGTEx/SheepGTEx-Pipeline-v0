###1.counts to TMM
library(data.table)
library(dplyr)
library(stringr)
library(edgeR)
library(limma)

sample_info = fread("sample_info.txt", header = TRUE, sep = "\t", data.table = FALSE)
tissues = sample_info %>%
  filter(Stage %in% c("Adult", "Lamb", "Prenatal")) %>%
  distinct(Tissue_s)
counts = fread("F_PCGlnc.gene.count.txt", data.table = FALSE, header = TRUE)
gene = fread("all_gene.list", data.table = FALSE, header = TRUE)
counts = counts %>% filter(Geneid %in% gene$Gene_ID)
names(counts)[1] = "gene_id"

for (t in unique(tissues$Tissue_s)) {
  cat("Analyzing ", t, "\n")
  tissue_samples = subset(sample_info, Tissue_s == t & Stage %in% c("Adult", "Lamb", "Prenatal"))
  samples_list = tissue_samples$Biosample
  tissue_counts = counts %>% select(gene_id, all_of(samples_list))
  rownames(tissue_counts) = tissue_counts$gene_id
  tissue_counts$gene_id = NULL
  if (!all(tissue_samples$Biosample %in% colnames(tissue_counts))) {
    stop("Sample names in tissue_samples do not match column names in tissue_counts")
  }
  tissue_counts = tissue_counts[apply(tissue_counts, 1, function(x) sum(x < 6) / length(x) <= 0.8), ]
  y = DGEList(counts = tissue_counts, group = tissue_samples$Stage)
  keep = filterByExpr(y)
  y = y[keep, , keep.lib.sizes = FALSE]
  y = calcNormFactors(y, method = "TMM")
  count_norm = cpm(y) %>% as.data.frame()
  colnames(count_norm) = colnames(tissue_counts)
  write.table(count_norm, file = paste0("TMM/TMM_", t, ".txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

###2.wilcoxon rank-sum test
libpath = "/storage/public/apps/software/R/R-4.3.0/lib64/R/library"
library(tidyr, lib.loc=libpath)
library(dplyr, lib.loc=libpath)
library(data.table, lib.loc=libpath)
library(SmartSVA)
library(limma)
library(edgeR)
output_folder_path <- "./results/"
gene_info <- fread("../gene_tss.txt")
stage_info <- fread("../sample_info.txt", header = TRUE, sep = "\t", data.table = FALSE)
tissues <- fread("../remain_tissue.txt", header = FALSE)[[1]]
for (t in tissues) {
  cat("Analyzing ", t, "\n")
  tmm_matrix <- fread(paste0("../TMM/TMM_", t, ".txt"))
  tmm_matrix <- as.data.frame(tmm_matrix)
  rownames(tmm_matrix) <- tmm_matrix$Geneid
  tmm_matrix$Geneid <- NULL
  stage_info_tissue <- stage_info %>% filter(Tissue_s == t)
  stage_counts <- stage_info_tissue %>%
    group_by(Stage_n) %>%
    summarise(SampleCount = n(), .groups = 'drop') %>%
    filter(SampleCount > 10)
  if (nrow(stage_counts) < 2) {
    cat("Skipping tissue:", t, "- Less than 2 valid stages\n")
    next
  }
  valid_stages <- stage_counts$Stage_n
  cat("Tissue:", t, "- Valid stages:", paste(valid_stages, collapse = ", "), "\n")

  stage_info_tissue <- stage_info_tissue %>% filter(Stage_n %in% valid_stages)
  stage_conditions <- factor(stage_info_tissue$Stage_n)
  tmm_matrix <- tmm_matrix[, colnames(tmm_matrix) %in% stage_info_tissue$Biosample, drop = FALSE]
  stage_conditions <- stage_conditions[match(colnames(tmm_matrix), stage_info_tissue$Biosample)]
  valid_genes <- gene_info$gene[gene_info$CHROM %in% c(paste0("chr", 1:26), "chrX")]
  filtered_tmm <- tmm_matrix[rownames(tmm_matrix) %in% valid_genes, ]
  mod <- model.matrix(~ stage_conditions)
  rownames(mod) <- colnames(filtered_tmm)
  n.sv <- max(num.sv(filtered_tmm, mod), 1)
  print(n.sv)
  sv <- smartsva.cpp(as.matrix(filtered_tmm), mod = mod, n.sv = n.sv)
  surrogateVariables <- sv$sv
  rownames(surrogateVariables) <- colnames(filtered_tmm)
  tmm_adjusted <- removeBatchEffect(as.matrix(filtered_tmm), covariates = surrogateVariables, design = mod)
  adjusted_tmm_file <- paste0("./TMM/", t, "_TMM_Adjusted.txt")
  write.table(tmm_adjusted, file = adjusted_tmm_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  pairwise_results <- combn(levels(stage_conditions), 2, function(stage_pair) {
    stage1 <- stage_pair[1]
    stage2 <- stage_pair[2]
    group1_samples <- which(stage_conditions == stage1)
    group2_samples <- which(stage_conditions == stage2)
    log2_fold_changes <- log2(rowMeans(tmm_adjusted[, group2_samples, drop = FALSE]) /
                              rowMeans(tmm_adjusted[, group1_samples, drop = FALSE]))
    pvalues <- apply(tmm_adjusted, 1, function(gene_expression) {
      data <- cbind.data.frame(
        gene = gene_expression[c(group1_samples, group2_samples)],
        group = factor(c(rep(stage1, length(group1_samples)), rep(stage2, length(group2_samples))))
      )
      wilcox.test(gene ~ group, data)$p.value
    })
    data.frame(
      GeneID = rownames(tmm_adjusted),
      Stage1 = stage1,
      Stage2 = stage2,
      log2FoldChange = log2_fold_changes,
      PValue = pvalues,
      FDR = p.adjust(pvalues, method = "fdr"),
      Bonferroni = p.adjust(pvalues, method = "bonferroni")
    )
  }, simplify = FALSE)
  all_pairwise_results <- do.call(rbind, pairwise_results)
  stage_means <- sapply(levels(stage_conditions), function(stage) {
    rowMeans(tmm_adjusted[, which(stage_conditions == stage), drop = FALSE])
  })
  preferred_stage <- apply(stage_means, 1, function(means) {
    names(which.max(means))
  })
  output_file_base <- paste0(output_folder_path, t)
  write.table(all_pairwise_results, file = paste0(output_file_base, "_all_results.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  preferred_stage_results <- data.frame(
    GeneID = rownames(tmm_adjusted),
    PreferredStage = preferred_stage
  )
  write.table(preferred_stage_results, file = paste0(output_file_base, "_prefer_stage.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

###3.GSEA
###Adult and Prenatal
library(clusterProfiler)
library(org.Hs.eg.db)
library(metap)
files <- list.files(path = "./01input/", pattern = ".*_bonfer005\\.txt", full.names = TRUE)
tissue_names <- gsub("_bonfer005\\.txt", "", basename(files))
gsea_results <- list()
for (i in seq_along(files)) {
  message("Processing file: ", basename(files[i]))
  data <- read.table(files[i], header = TRUE, sep = "\t")
  data <- data[!is.na(data$log2FoldChange), ]
  #data <- subset(data, Chromosome != "chrX")
  data <- data[order(-data$log2FoldChange),]
  geneSymbols <- data$GeneID
  geneList <- as.numeric(data$log2FoldChange)
  names(geneList) <- geneSymbols
  gseaResults <- gseGO(geneList = geneList,
                       ont = "BP",
                       keyType = "SYMBOL",
                       OrgDb = org.Hs.eg.db,
                       nPerm = 1000,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pAdjustMethod = "BH",
                       verbose = FALSE,
                       pvalueCutoff = 1)
  gsea_df <- as.data.frame(gseaResults@result)
  gsea_df$Tissue <- tissue_names[i]
  gsea_results[[tissue_names[i]]] <- gsea_df
}
combined_results <- do.call(rbind, gsea_results)
colnames(combined_results)
term_counts <- combined_results %>%
  group_by(ID) %>%
  summarise(tissue_count = length(unique(Tissue))) %>% 
  filter(tissue_count >= 5)
valid_terms <- term_counts$ID
filtered_combined_results <- combined_results %>%
  filter(ID %in% valid_terms)
write.csv(combined_results, "./output/Combined_GSEAwithTissue.csv")
write.csv(filtered_combined_results, "./output/5TissueCombined_GSEA.csv")
##Meta_anlaysis
library(dplyr)
data <- read.csv("./output/5TissueCombined_GSEA.csv")
filtered_data <- data %>%
  mutate(direction = ifelse(enrichmentScore > 0, "Adult_biased", "Prenatal_biased")) %>%
  group_by(ID) %>%
  mutate(
    dominant_direction = ifelse(sum(enrichmentScore > 0) > sum(enrichmentScore < 0), "Adult_biased", "Prenatal_biased"),
    pvalue_adjusted = ifelse(direction != dominant_direction, 1/pvalue, pvalue)
  ) %>%
  ungroup()
fisher_results_with_direction <- filtered_data %>%
  group_by(ID) %>%
  summarise(
    fisher_statistic = -2 * sum(log(pvalue)),
    n = n(),
    dominant_direction = dplyr::first(dominant_direction),
    genes = dplyr::first(core_enrichment)
  ) %>%
  mutate(
    fisher_pvalue = pchisq(fisher_statistic, df = 2 * n, lower.tail = FALSE)
  ) %>%
  ungroup()
fisher_results_corrected_with_direction <- fisher_results_with_direction %>%
  mutate(
    adjust_pvalue = p.adjust(fisher_pvalue, method = "fdr")
  )
significant_terms_with_direction <- fisher_results_corrected_with_direction %>%
  filter(adjust_pvalue <= 0.05)

print(significant_terms_with_direction)
write.csv(significant_terms_with_direction, "./output/Meta_analysis_sigresult.fdr005.csv")

###Adult and Lamb
library(clusterProfiler)
library(org.Hs.eg.db)
library(metap)

files <- list.files(path = "./input/", pattern = ".*_stage\\.txt", full.names = TRUE)
tissue_names <- gsub("_stage\\.txt", "", basename(files))
gsea_results <- list()
for (i in seq_along(files)) {
  message("Processing file: ", basename(files[i]))
  data <- read.table(files[i], header = TRUE, sep = "\t")
  data <- data[!is.na(data$log2FoldChange), ]
  #data <- subset(data, Chromosome != "chrX")
  data <- data[order(-data$log2FoldChange),]
  geneSymbols <- data$GeneID
  geneList <- as.numeric(data$log2FoldChange)
  names(geneList) <- geneSymbols
  gseaResults <- gseGO(geneList = geneList,
                       ont = "BP",
                       keyType = "SYMBOL",
                       OrgDb = org.Hs.eg.db,
                       nPerm = 1000,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pAdjustMethod = "BH",
                       verbose = FALSE,
                       pvalueCutoff = 1)

  gsea_df <- as.data.frame(gseaResults@result)
  gsea_df$Tissue <- tissue_names[i]
  gsea_results[[tissue_names[i]]] <- gsea_df
}
combined_results <- do.call(rbind, gsea_results)
colnames(combined_results)
term_counts <- combined_results %>%
  group_by(ID) %>%
  summarise(tissue_count = length(unique(Tissue))) %>% 
  filter(tissue_count >= 9)
valid_terms <- term_counts$ID
filtered_combined_results <- combined_results %>%
  filter(ID %in% valid_terms)

write.csv(combined_results, "./12output/Combined_GSEAwithTissue.csv")
write.csv(filtered_combined_results, "./12output/9TissueCombined_GSEA.csv")

###Meta_anlaysis
library(dplyr)
data <- read.csv("./12output/9TissueCombined_GSEA.csv")

filtered_data <- data %>%
  mutate(direction = ifelse(enrichmentScore > 0, "Lamb_biased", "Adult_biased")) %>%
  group_by(ID) %>%
  mutate(
    dominant_direction = ifelse(sum(enrichmentScore > 0) > sum(enrichmentScore < 0), "Adult_biased", "Lamb_biased"),
    pvalue_adjusted = ifelse(direction != dominant_direction, 1/pvalue, pvalue)
  ) %>%
  ungroup()
fisher_results_with_direction <- filtered_data %>%
  group_by(ID) %>%
  summarise(
    fisher_statistic = -2 * sum(log(pvalue)),
    n = n(),
    dominant_direction = dplyr::first(dominant_direction),
    genes = dplyr::first(core_enrichment)
  ) %>%
  mutate(
    fisher_pvalue = pchisq(fisher_statistic, df = 2 * n, lower.tail = FALSE)
  ) %>%
  ungroup()
fisher_results_corrected_with_direction <- fisher_results_with_direction %>%
  mutate(
    adjust_pvalue = p.adjust(fisher_pvalue, method = "fdr")
  )
significant_terms_with_direction <- fisher_results_corrected_with_direction %>%
  filter(adjust_pvalue <= 0.05)

print(significant_terms_with_direction)

write.csv(significant_terms_with_direction, "./12output/Meta_analysis_sigresult.fdr005.csv")
