#!/opt/app/languages/R-3.6.3/bin/Rscript
# 2022-01-08
# Edited by https://github.com/3UTR/3aQTL-pipe/blob/main/src/curate_pheno_geno_covariates.R

library(optparse)
library(tidyr)
library(dplyr)
library(impute)

# -- global variable
option_list <- list(
  make_option(c("-p", "--pheno_data"),type = "character", default = "Dapars2_res.all_chromosomes.txt", action = "store", help = "Proivde the merged ouput from DaPars2, Dapars2_res.all_chromosomes.txt in default"),
  make_option(c("-t","--tss_pos"),type = "character", action = "store", help = "the TSS position file of genes with header of gene, CHROM, POS. 1-based coordinate."),
  make_option(c("-s","--sample_participant"),type = "character", action = "store", help = "File listing sample IDs to participant IDs to include (no header), formated as: sample\tparticipant"),
  make_option(c("-o","--outprefix"),type = "character", action = "store", help = "Output file prefix")
  #make_option(c("-g","--geno_pca"),type = "character", default = "./Matrix_eQTL/genotype_pca.eigenvec", action = "store", help = "the eigenvector of genotpye pca analysis, ./Matrix_eQTL/genotype_pca.eigenvec in default"),
  #make_option(c("-c","--known_covs"),type = "character", default = "NA", action = "store", help = "input a text file contains known covariates if available, NA in defaut"),
  #make_option(c("-n","--top_N_pca"),type = "integer", default = "5", action = "store", help = "specify the top N PCA on gt would be used, defaut value = 5")
)
opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))
apa_res_file <- opt$pheno_data
tss_file <- opt$tss_pos
sample_file <- opt$sample_participant
outprefix <- opt$outprefix


cat("Start phenotype matrix\n","Open APA results file:",apa_res_file,"\n")
# --------------- prepare phenotype matrix ------------------
pdui_mat <- read.table(apa_res_file, stringsAsFactors=FALSE, header=TRUE,sep="\t",check.names=FALSE)
pdui_mat %>% separate(., col="Gene", into=c("transcript_id", "gene_id", "chr", "strand"), sep = "\\|") -> pdui_mat
pdui_mat <- pdui_mat[,-c(5,6,7)]

# extract samples 
sample_to_participant_s <- read.table(sample_file, sep='\t', header=FALSE, stringsAsFactors=FALSE)
pdui_mat_ann = pdui_mat[1:4]
pdui_mat_pdui = pdui_mat[5:ncol(pdui_mat)]
pdui_mat_pdui = pdui_mat_pdui[, sample_to_participant_s$V1]
# change sample IDs to participant IDs
colnames(pdui_mat_pdui) = sample_to_participant_s$V2
pdui_mat.sel = cbind(pdui_mat_ann,pdui_mat_pdui)

# extract genes
tss_df <- read.table(tss_file, stringsAsFactors=FALSE, header=TRUE,sep="\t",check.names=FALSE)[, c("gene", "CHROM", "POS")]
genes <- as.vector(tss_df$gene)
pdui_mat.sel <- pdui_mat.sel[pdui_mat.sel$gene_id %in% genes, ]
pdui_mat.sel %>% unite("ID",transcript_id:strand, sep = "|") -> pdui_mat.sel

rownames(pdui_mat.sel) <- pdui_mat.sel$ID
pdui_mat.sel$ID <- NULL
pdui_mat.sel <- as.matrix(pdui_mat.sel)

cat('remove genes with more than 50% entries missing and individuals with more than 80% missing data.\n')
pdui_mat.sel <- pdui_mat.sel[, colMeans(is.na(pdui_mat.sel)) <= 0.8]
pdui_mat.sel <- pdui_mat.sel[rowMeans(is.na(pdui_mat.sel)) <= 0.5,]
dim(pdui_mat.sel)
class(pdui_mat.sel) <- 'numeric'

#impute missing values in PDUI matrix
mat.ds <- pdui_mat.sel
mat_impute <- impute.knn(mat.ds)
df_w <- as.data.frame(mat_impute$data)

#remove constant phenotypes
constant <- apply(df_w, 1, function(i) length(unique(i)) == 1)
constant_names <- rownames(df_w)[constant]
cat(length(constant_names),"constant phenotypes were removed:\n",constant_names,"\n")
df_w <- df_w[!constant,]
dim(df_w)

#quantile normalization
for(gene in 1:nrow(df_w)){
  mat = df_w[gene,]
  mat = apply(mat,1,rank,ties.method = "average")
  mat = qnorm(mat / (ncol(df_w)+1))
  df_w[gene,] = mat
}

# output BED file for tensorQTL
df_w = as.data.frame(df_w)
norm_pdui_mat <- cbind(rownames(df_w),df_w)
y <- colnames(norm_pdui_mat)[-1]
colnames(norm_pdui_mat) <- c("ID",y)
norm_pdui_mat %>% separate(., col="ID", into=c("transcript_id", "gene", "chr", "strand"), sep = "\\|") -> tmp_df
merged_df <- tss_df %>% inner_join(tmp_df, by = "gene")
merged_df %>% unite("ID",transcript_id,gene,chr,strand, sep = "|") -> norm_pdui_mat.sel
norm_pdui_mat.sel$chrom = as.numeric(sub("X", "99", sub("^chr", "", norm_pdui_mat.sel$CHROM)))
sorted_indices <- order(norm_pdui_mat.sel$chrom, norm_pdui_mat.sel$POS)
norm_pdui_mat.sel <- norm_pdui_mat.sel[sorted_indices, ]
norm_pdui_mat.sel$chrom = NULL
norm_pdui_mat.sel %>%
  mutate(start = POS-1 ) %>%
  rename("#chr"=CHROM, "end" = POS) %>%
  select("#chr", start, end, ID, everything()) -> norm_pdui_mat.sel
cat("Output normalized phenotype matrix for tensorQTL:\n")
write.table(norm_pdui_mat.sel,file=paste0(outprefix,".bed"),row.names=F,col.names=T,quote=F,sep="\t")
