library(data.table)
library(dplyr)
library(purrr)
library(coloc)

ARGS <- commandArgs(trailingOnly = TRUE)
file1 = ARGS[1] # nominal file of eQTL of tissue1
file2 = ARGS[2] # nominal file of eQTL of tissue2
genefile = ARGS[3] # list of significant genes in at least one QTL
#tissue = ARGS[4] # tissue name
chrom = ARGS[4] # chromosome name
outfile = ARGS[5] # output file of coloc
if (!file.exists(file1)) { stop("Can not find the nominal file1") }
if (!file.exists(file2)) { stop("Can not find the nominal file2") }

# input qtl files
nominal1 = fread(file1, header = T)
nominal2 = fread(file2, header = T)
genelist = fread(genefile, header = F)$V1

df = inner_join(nominal1, nominal2, by = c("phenotype_id", "variant_id"))
df$phenotype_id.y = df$phenotype_id

# extract genes
df = df[phenotype_id %in% genelist]
sample_size = round(max(df$ma_count.x / df$af.x)/2)

# calculate the pph4
unique_phenotypes <- unique(df$phenotype_id.y)
for (phenotype in unique_phenotypes) {
  df_target <- subset(df, phenotype_id.y == phenotype)
  # ignore genes with missing data
  df_test <- as.data.frame(df_target)
  if (any(is.na(df_test[c("slope.x", "pval_nominal.x", "slope.y", "pval_nominal.y")]))) { next }
  
  gene = unique(df_target$phenotype_id)
  df1_coloc = list(beta=df_target$slope.x, pvalues=df_target$pval_nominal.x,
                   snp=df_target$variant_id, type="quant", N=rep(sample_size, nrow(df_target)),
                   MAF=as.numeric(ifelse(df_target$af.x <= 0.5, df_target$af.x, 1-df_target$af.x)))
  df2_coloc = list(beta=df_target$slope.y, pvalues=df_target$pval_nominal.y,
                   snp=df_target$variant_id, type="quant", N=rep(sample_size, nrow(df_target)),
                   MAF=as.numeric(ifelse(df_target$af.y <= 0.5, df_target$af.y, 1-df_target$af.y)))
  print(phenotype)
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = gene
  results$phenotype = phenotype
  results$chrom = chrom
  #results$tissue = tissue
  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf phenotype1 phenotype2 chrom tissue
  fwrite(results, outfile, append=TRUE, col.names=FALSE, sep="\t")
}
