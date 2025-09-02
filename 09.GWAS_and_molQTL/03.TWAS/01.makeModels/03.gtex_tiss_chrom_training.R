"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
work_path <- argv[2]

source("./gtex_v7_nested_cv_elnet.R")

setwd(work_path)
snp_annot_file <- "./output/snp_annot.chr" %&% chrom %&% ".txt"
gene_annot_file <- "./../gene_annot.txt"
genotype_file <- "./output/genotype.chr" %&% chrom %&% ".txt"
expression_file <- "./output/transformed_expression.txt"
covariates_file <- "./output/covariates.txt"
prefix <- "Model_training"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


