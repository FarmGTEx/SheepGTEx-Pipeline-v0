library(qvalue)
library(data.table)
library(dplyr)

ARGS <- commandArgs(trailingOnly = TRUE)
file_trans = ARGS[1] # trans-QTL result file
out_prefix = ARGS[2] # output file prifx
fdr_thresholds = as.numeric(ARGS[3])
# fdr_thresholds = 0.05

#Read data
trans = fread(file_trans)
trans$pval_multiple = trans$pval_g1*1000000

# run p.adjust (BH) on pvalues for each gene
trans_gene <- trans %>%  group_by(pheno_id) %>% mutate(pval_adj_BH_gene = p.adjust(pval_multiple, method = "BH"))

# run qvalue on pvalues for best signals
eqtl <- trans_gene %>%  group_by(pheno_id) %>% slice_min(pval_multiple) %>% sample_n(1)
pval_adj_BH = p.adjust(eqtl$pval_multiple, method = "BH")

#Determine significance threshold
set0 = eqtl[which(pval_adj_BH <= fdr_thresholds),]
set1 = eqtl[which(pval_adj_BH > fdr_thresholds),]

#eqtl$qval = Q$qvalues
eqtl$qval = pval_adj_BH
#eqtl$pval_nominal_threshold = nthresholds
eqtl$is_eGene = (eqtl$qval < fdr_thresholds) & (eqtl$pval_adj_BH_gene < fdr_thresholds)
eqtl_subset <- eqtl %>% select(tissue, pheno_id, qval, is_eGene)
trans_sig <- left_join(trans_gene, eqtl_subset, by = c("tissue", "pheno_id")) %>% filter((is_eGene == TRUE) & (pval_adj_BH_gene < fdr_thresholds))

#Output
fwrite(eqtl, paste0(out_prefix, ".txt.gz"),sep="\t")
fwrite(trans_sig, paste0(out_prefix, ".sig.txt.gz"),sep="\t")
cat("Done\n")
