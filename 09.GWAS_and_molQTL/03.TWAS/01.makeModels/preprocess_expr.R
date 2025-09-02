# Load packages
library(tidyverse)
library(dplyr)
library(RSQLite)

"%&%" <- function(a,b) paste(a,b, sep='')
argv <- commandArgs(trailingOnly = TRUE)
type <- argv[1]
tissue <- argv[2]

gene_exp_transpose_file = type %&% "/" %&% tissue %&% "/output/transformed_expression.txt"
covariates_file = type %&% "/" %&% tissue %&% "/data/covariates.txt"
gene_exp_transpose = read.table(file = gene_exp_transpose_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
covariates = read.table(file = covariates_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)

#Set the column names for the PEER factors (covariates) as the subject IDs
colnames(covariates) = rownames(gene_exp_transpose)

# write out a covariates matrix
write.table(covariates, file = type %&% "/" %&% tissue %&% "/output/covariates.txt", sep = "\t", row.names = TRUE)

## Make a copy of the transposed gene expression dataframe so that we can replace the values with the residuals of the multiple linear regressions.
expression = gene_exp_transpose

# This loops through all the columns of the transposed gene expression which correspond to each gene,
# for each gene it runs linear regression on the PEER factor covariates. Then it sets the residuals to the new expression for that gene.

for (i in 1:length(colnames(gene_exp_transpose))) {
    fit = lm(gene_exp_transpose[,i] ~ t(as.matrix(covariates)))
    expression[,i] <- fit$residuals
  }

# Write out the residual expression file
write.table(expression, file = type %&% "/" %&% tissue %&% "/output/residuals_expression.txt", sep = "\t", row.names = TRUE)
