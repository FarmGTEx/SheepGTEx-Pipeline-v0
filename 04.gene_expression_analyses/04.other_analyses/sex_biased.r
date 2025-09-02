###1.counts to TMM
library(data.table)
library(dplyr)
library(stringr)
library(edgeR)
library(limma)

sex_samples = fread("../sex_deg.samples.list", header = TRUE, sep = "\t", data.table = FALSE)
sex_samples = sex_samples[, 1:3]
tissues = sex_samples %>%
  filter(Sex == "0" | Sex == "1") %>%
  group_by(Tissue_s, Sex) %>%
  tally() %>%
  filter(n >= 10) %>%
  as.data.frame() %>%
  subset(Tissue_s %in% names(subset(table(.$Tissue_s), table(.$Tissue_s) == 2)))
counts = fread("../F_PCGlnc.gene.count.txt", data.table = FALSE, header = TRUE)
gene = fread("../all_gene.list", data.table = FALSE, header = TRUE)
counts = counts %>% filter(Geneid %in% gene$Gene_ID)
names(counts)[1] = "gene_id"

for (t in unique(tissues$Tissue_s)) {
  cat("Analyzing ", t, "\n")
  samples = subset(sex_samples, Tissue_s == t & (Sex == "0" | Sex == "1"))[,"Biosample"]
  tissue_counts = counts %>% select(gene_id, all_of(samples))
  rownames(tissue_counts) = tissue_counts$gene_id
  tissue_counts$gene_id = NULL
  tissue_samples = sex_samples %>% filter(Biosample %in% samples) %>% select(Biosample, Sex) %>% as.data.frame()
  rownames(tissue_samples) = tissue_samples$Biosample
  tissue_samples$Biosample = NULL
  tissue_samples$Sex = factor(tissue_samples$Sex)
  if (!all(rownames(tissue_samples) %in% colnames(tissue_counts))) {
    stop("Sample names in tissue_samples do not match column names in tissue_counts")
  }
tissue_counts = tissue_counts[apply(tissue_counts, 1, function(x) sum(x < 6) / length(x) <= 0.8), ]
  y = DGEList(counts = tissue_counts, group = tissue_samples$Sex)
  keep = filterByExpr(y)
  y = y[keep, , keep.lib.sizes = FALSE]
  y = calcNormFactors(y, method = "TMM")
  count_norm = cpm(y) %>% as.data.frame()
  colnames(count_norm) = colnames(tissue_counts)
  write.table(count_norm, file = paste0("TMM_new/TMM_", t, ".txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

###2.wilcoxon rank-sum test
libpath = "/storage/public/apps/software/R/R-4.3.0/lib64/R/library"
library(tidyr, lib.loc=libpath)
library(dplyr, lib.loc=libpath)
library(data.table, lib.loc=libpath)
require(SmartSVA)
library(limma)
library(edgeR)
#t <- as.character(commandArgs(trailingOnly = TRUE)[1])

output_folder_path <- paste0("./results_new/")
gene_info <- fread("../gene_tss.txt")
conditions <- fread("../sexinfo_updated.txt", header = TRUE, sep = "\t", data.table = FALSE)
conditions[, 4:6] <- NULL
tissues <- fread("sex_tissue.txt", header =F)[[1]]
for (t in tissues) {
cat("Analyzing ", t, "\n")
tmm_matrix <- fread(paste0( "./TMM_new/TMM_", t, ".txt"))
tmm_matrix <- as.data.frame(tmm_matrix)
rownames(tmm_matrix) = tmm_matrix$Geneid
tmm_matrix$Geneid = NULL

sexinfo <- conditions %>% filter(Tissue_s == t)
sexinfo_conditions <- factor(sexinfo$Sex)
if (length(unique(sexinfo_conditions)) < 2) next
valid_genes <- gene_info$gene[gene_info$CHROM %in% c(paste0("chr", 1:26), "chrX")]
filtered_tmm <- tmm_matrix[rownames(tmm_matrix) %in% valid_genes, ]
mod <- model.matrix(~ sexinfo_conditions)
rownames(mod) <- colnames(filtered_tmm)
n.sv <- num.sv(filtered_tmm, mod)
sv <- smartsva.cpp(as.matrix(filtered_tmm), mod = mod, n.sv = n.sv)
surrogateVariables <- sv$sv
rownames(surrogateVariables) <- colnames(filtered_tmm)
tmm_adjusted <- removeBatchEffect(as.matrix(filtered_tmm), covariates = surrogateVariables, design = mod)
pvalues <- sapply(1:nrow(tmm_adjusted), function(i) {
data <- cbind.data.frame(gene = as.numeric(t(tmm_adjusted[i,])), sexinfo_conditions)
    wilcox.test(gene ~ sexinfo_conditions, data)$p.value
})
fdr <- p.adjust(pvalues, method = "fdr")
bonferroni <- p.adjust(pvalues, method = "bonferroni")
conditionsLevel <- levels(sexinfo_conditions)
print(conditionsLevel)
dataCon1 <- tmm_adjusted[, which(sexinfo_conditions == conditionsLevel[1])]
dataCon2 <- tmm_adjusted[, which(sexinfo_conditions == conditionsLevel[2])]
foldChanges <- log2(rowMeans(dataCon2) / rowMeans(dataCon1))
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr, Bonferroni = bonferroni)
rownames(outRst) <- rownames(tmm_adjusted)
outRst <- na.omit(outRst)
outRst$Chromosome <- gene_info$CHROM[match(rownames(outRst), gene_info$gene)]
outRst$Direction <- ifelse(outRst$log2foldChange > 0, conditionsLevel[2], conditionsLevel[1])
significant_results_fdr <- outRst[outRst$FDR < 0.05, ]
significant_results_bonferroni <- outRst[outRst$Bonferroni < 0.05, ]
large_effect_genes_fdr <- significant_results_fdr[abs(significant_results_fdr$log2foldChange) > 1, ]
large_effect_genes_bonferroni <- significant_results_bonferroni[abs(significant_results_bonferroni$log2foldChange) > 1, ]
output_file_base <- paste0(output_folder_path)
write.table(outRst, file = paste0(output_file_base, t, ".txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(significant_results_fdr, file = paste0(output_file_base, t, "_FDR.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(significant_results_bonferroni, file = paste0(output_file_base, t,  "_Bonferroni.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(large_effect_genes_fdr, file = paste0(output_file_base, t, "_LargeEffect_FDR.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(large_effect_genes_bonferroni, file = paste0(output_file_base, t, "_LargeEffect_Bonferroni.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

###3.GSEA
##gsea
library(clusterProfiler)
library(org.Hs.eg.db)
library(metap)
files <- list.files(path = "./wilcox_results/", pattern = ".*_Bonferroni\\.tsv", full.names = TRUE)
tissue_names <- gsub("_Bonferroni\\.tsv", "", basename(files))
gsea_results <- list()
for (i in seq_along(files)) {
  message("Processing file: ", basename(files[i]))
  data <- read.table(files[i], header = TRUE, sep = "\t")
  data <- subset(data, Chromosome != "chrX")
  data <- data[order(-data$log2foldChange),]
  geneSymbols <- data$Geneid
  geneList <- as.numeric(data$log2foldChange)
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
  filter(tissue_count >= 10)
valid_terms <- term_counts$ID
filtered_combined_results <- combined_results %>%
  filter(ID %in% valid_terms)
write.csv(combined_results, "./output/Combined_GSEAwithTissue.csv")
write.csv(filtered_combined_results, "./output/10TissueCombined_GSEA.csv")
##Meta_anlaysis
library(dplyr)
data <- read.csv("./output/1206/10TissueCombined_GSEA.csv")
filtered_data <- data %>%
  mutate(direction = ifelse(enrichmentScore > 0, "male_biased", "female_biased")) %>%
  group_by(ID) %>%
  mutate(
    dominant_direction = ifelse(sum(enrichmentScore > 0) > sum(enrichmentScore < 0), "male_biased", "female_biased"),
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
write.csv(significant_terms_with_direction, "./output/1206/Meta_analysis_sigresult.fdr005.csv")





