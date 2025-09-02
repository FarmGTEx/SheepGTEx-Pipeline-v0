libpath = "/storage/public/apps/software/R/R-4.3.0/lib64/R/library"
library(data.table, lib.loc = libpath)
library(limma, lib.loc = libpath)
require(PCAForQTL)

# input data
ARGS <- commandArgs(trailingOnly = TRUE)
pheno_bed = ARGS[1] # phenotype bed file. The first four columns are: #chr, start, end, and gene_id.
sample_file = ARGS[2] # File listing sample IDs to participant IDs to include (no header), formated as: sample\tparticipant\ttissue")
out_file = ARGS[3] # output file name


pheno_mat <- as.data.frame(fread(pheno_bed, header=T, check.names=F, stringsAsFactors=F)) #The first four columns are: #chr, start, end, and gene_id.
# extract samples 
sample_to_participant_s <- read.table(sample_file, sep='\t', header=FALSE, stringsAsFactors=FALSE)
pheno_ann = pheno_mat[1:6]
pheno_expr = pheno_mat[7:ncol(pheno_mat)]
pheno_expr = pheno_expr[, sample_to_participant_s$V1]
# change sample IDs to participant IDs
colnames(pheno_expr) = sample_to_participant_s$V2
pheno_mat = cbind(pheno_ann,pheno_expr)

#remove constant phenotypes
#constant <- apply(pheno_expr, 1, function(i) length(unique(i)) == 1)
#constant_names <- rownames(pheno_expr)[constant]
#cat(length(constant_names),"constant phenotypes were removed:\n",constant_names,"\n")
#pheno_mat <- pheno_mat[!constant,]
dim(pheno_mat)
pheno_ann = pheno_mat[1:6]
pheno_expr = pheno_mat[7:ncol(pheno_mat)]

#estimate phenotype PCs
expr <- t(pheno_expr)
dim(expr)
#pheno_prcomp<-prcomp(expr,center=TRUE, scale.=TRUE) #This should take less than a minute.
pheno_prcomp<-prcomp(expr,center=TRUE, scale.=FALSE) #This should take less than a minute.
pheno_PCs<-pheno_prcomp$x
dim(pheno_PCs)
pheno_importance<-summary(pheno_prcomp)$importance
pheno_PVEs<-pheno_importance[2,]
sum(pheno_PVEs) #Theoretically, this should be 1.
##choose K for PCs
pheno_RunElbow<-PCAForQTL::runElbow(prcompResult=pheno_prcomp)
print(paste0("phenotype Elbow: ", pheno_RunElbow))
pheno_K_elbow<-pheno_RunElbow
pheno_PCs<-as.data.frame(pheno_PCs)
pheno_PCsTop<-pheno_PCs[,1:pheno_K_elbow]
##remove batch effact
pheno_rmbatch <- as.data.frame(removeBatchEffect(pheno_expr, covariates=pheno_PCsTop))

# output BED file after remove batch effect
pheno_mat = cbind(pheno_ann,pheno_rmbatch)
write.table(pheno_mat,file=out_file,row.names=F,col.names=T,quote=F,sep="\t")
