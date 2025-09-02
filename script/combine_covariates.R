library(data.table, lib.loc = "/storage/public/home/2020060185/anaconda3/envs/rnaseq/lib/R/library")
library(dplyr, lib.loc = "/storage/public/home/2020060185/anaconda3/envs/rnaseq/lib/R/library")
library(PCAForQTL)

# input data
ARGS <- commandArgs(trailingOnly = TRUE)
phenopcfile = ARGS[1] # phenotype top PCs
genopcfile = ARGS[2] # genotype top PCs
covarfile = ARGS[3] # other known covariates such as sex
outfile = ARGS[4] # output combined covariates file
if (!file.exists(phenopcfile)) { stop("Can not find the phenotype top PCs file") }
if (!file.exists(genopcfile)) { stop("Can not find the genotype top PCs file") }
if (!file.exists(covarfile)) { stop("Can not find the known covariates file") }

# input phenotype PCs
pheno_PCsTop<-fread(phenopcfile, header =T, data.table = F)
rownames(pheno_PCsTop) <- pheno_PCsTop$V1
pheno_PCsTop$V1 = NULL
dim(pheno_PCsTop)

# input genotype PCs
geno_PCsTop<-fread(genopcfile, header =T, data.table = F)
rownames(geno_PCsTop) <- geno_PCsTop$V1
geno_PCsTop <- geno_PCsTop[rownames(pheno_PCsTop),]
geno_PCsTop$V1 = NULL
dim(geno_PCsTop)

# input other known covariates
otherCovariates<-fread(covarfile, header = T, data.table = F)

# compare and combine covariates
if (dim(otherCovariates)[2]>1) {
  rownames(otherCovariates) <- otherCovariates$Sample
  otherCovariates <- otherCovariates[rownames(pheno_PCsTop),]
  otherCovariates$Sample = NULL
  otherCovariatesFiltered<-PCAForQTL::filterKnownCovariates(otherCovariates,geno_PCsTop, unadjustedR2_cutoff=0.9, verbose = T)
  knownCovariates<-cbind(otherCovariatesFiltered,geno_PCsTop)
  dim(knownCovariates)
} else {
  knownCovariates<-geno_PCsTop
}

## compare combine other covariates and genotype PCs with expression PCs
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,pheno_PCsTop, unadjustedR2_cutoff=0.9, verbose = T)
dim(knownCovariatesFiltered)
covariatesToUse<-cbind(knownCovariatesFiltered,pheno_PCsTop)
dim(covariatesToUse)
fwrite(as.data.frame(t(covariatesToUse)), file = outfile, sep = "\t", quote = F, row.names=T)

