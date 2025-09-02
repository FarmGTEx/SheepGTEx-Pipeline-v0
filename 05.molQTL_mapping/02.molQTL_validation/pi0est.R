library(qvalue)
library(data.table)
library(dplyr)


ARGS <- commandArgs(trailingOnly = TRUE)
discover_path=ARGS[1] # path of tensorqtl results of discovery data (permutation)
validat_path=ARGS[2] # path of tensorqtl results of validation data (nominal)
discover_name=ARGS[3] # name of discovery data
validat_name=ARGS[4] # name of validation data
tissue=ARGS[5] # tissue name
outprefix=ARGS[6] # output the pi1 value and correlation of effect size (slope or z-score)

permfile = paste0(discover_path, "/permutation/", tissue, ".cis_qtl_fdr0.05.txt.gz")
dfp = fread(permfile) %>% as.data.frame
dfp = dfp[dfp$is_eGene ==TRUE,]
dfp['z'] = dfp['slope']/dfp['slope_se']

# extract p value and slope value
dfp_pval <- c()
dfp_slope <- c()
dfp_zscore <- c()
dfn_pval <- c()
dfn_slope <- c()
dfn_zscore <- c()
chrom0 = NaN
for(i in 1:nrow(dfp)){
  rowp <- dfp[i, ]
  chrom = strsplit(rowp$variant_id, "_")[[1]][1]
  
  if (chrom != chrom0) {
    nominalfile = paste0(validat_path, "/nominal/", tissue, ".cis_qtl_pairs.chr", chrom, ".txt.gz")
    dfn = fread(nominalfile) %>% as.data.frame
    chrom0 = chrom
  }
  rown = dfn[(dfn$phenotype_id==rowp$phenotype_id)&(dfn$variant_id==rowp$variant_id),]
  
  if (dim(rown)[1]) {
    dfp_pval <- c(dfp_pval, rowp$pval_nominal)
    dfp_slope <- c(dfp_slope, rowp$slope)
    dfp_zscore <- c(dfp_zscore, rowp$slope/rowp$slope_se)
    dfn_pval <- c(dfn_pval, rown$pval_nominal)
    dfn_slope <- c(dfn_slope, rown$slope)
    dfn_zscore <- c(dfn_zscore, rown$slope/rown$slope_se)
  } else {
    cat(rowp$phenotype_id, rowp$variant_id, "is not exist in validation data!\n")
  }
}

# calculate pi1 value and correlation of effect size
q_values <- try(qvalue(p=dfn_pval)$qvalues, silent = TRUE)
if ("try-error" %in% class(q_values)) {
    # if qvalue fails this is probably because there are no p-values greater than
    # 0.95 (the highest lambda value
    # if so add a single p-value of 1 to try to combat the problem
    qobj <- qvalue(p=c(dfn_pval, 1))
    print(qobj$pvalues[length(qobj$pvalues)])
    pval = data.frame(qobj$pvalues[-length(qobj$pvalues)])
} else {
    qobj <- qvalue(p=dfn_pval)
    pval = data.frame(qobj$pvalues)
}
pi1 = 1 - qobj$pi0
rs = cor(dfp_slope, dfn_slope, use = "complete.obs")
rz = cor(dfp_zscore, dfn_zscore, use = "complete.obs")
df <- data.frame(slope1=unlist(dfp_slope), slope2=unlist(dfn_zscore),
                 zscore1=unlist(dfp_zscore), zscore2=unlist(dfn_zscore))

print(qobj$pvalues[length(qobj$pvalues)])
fwrite(pval, paste0(outprefix, ".pval.txt"), col.names=FALSE, sep="\t")
fwrite(list(discover_name, validat_name, tissue, pi1, rs, rz), paste0(outprefix, ".pi1.txt"), col.names=FALSE, sep="\t")
fwrite(df, paste0(outprefix, ".effect.txt"), sep="\t")
