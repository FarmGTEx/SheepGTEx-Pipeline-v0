library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)


ARGS <- commandArgs(trailingOnly = TRUE)
phenopcfile = ARGS[1] # phenotype PCs
genopcfile = ARGS[2] # genotype PCs
covarfile = ARGS[3] # known covariates
prefix = ARGS[4]

# input top 40 phenotype PCs
pheno_PCs <- fread(phenopcfile, header =T, data.table = F)
rownames(pheno_PCs) <- pheno_PCs$V1
pheno_PCs$V1 = NULL
pheno_PCsTop <- pheno_PCs[, 1:40]

# input top 40 genotype PCs
geno_PCs <- fread(genopcfile, header =T, data.table = F)
rownames(geno_PCs) <- geno_PCs$V1
geno_PCs$V1 = NULL
geno_PCsTop <- geno_PCs[, 1:40]
geno_PCsTop <- geno_PCsTop[rownames(pheno_PCsTop),]

# input other known covariates
otherCovariates<-fread(covarfile, header = T, data.table = F)
rownames(otherCovariates) <- otherCovariates$Sample
otherCovariates <- otherCovariates[rownames(pheno_PCsTop),]
otherCovariates$Sample = NULL
otherCovariates <- otherCovariates %>%
  select(where(~ n_distinct(.) > 1))
dim(otherCovariates)

# pairwise comparasion of R2 using lm function
df1 <- cbind(pheno_PCsTop, geno_PCsTop)
df2 <- cbind(pheno_PCsTop, geno_PCsTop, otherCovariates)
names(df2) = str_c(" ", names(df2))
df2 <- df2 %>%
  mutate_if(~ !is.numeric(.), as.factor)
cols_df1 <- names(df1)
cols_df2 <- names(df2)
## adjust R2 using lm function
column_combinations <- crossing(col1 = names(df1), col2 = names(df2))
result_df <- pmap_dfr(
  column_combinations, 
  ~ {
    lm_fit <- lm(df1[[..1]] ~ df2[[..2]])
    adj_r_squared <- summary(lm_fit)$adj.r.squared
    tibble(name = paste0(..1, ..2), adj_r_squared = adj_r_squared)
  }
)

result_df <- result_df %>% separate_wider_delim(name, delim=' ', names=c("y", "x"))
result_df$adj_r_squared[result_df$adj_r_squared < 0] <- 0
fwrite(result_df, file = paste0(prefix, "_adjr2.tsv"), sep = "\t", quote = F, col.names = T, row.names=F)
