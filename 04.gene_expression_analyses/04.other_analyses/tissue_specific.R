library(tidyverse)
library(data.table)
library(ggsci)
suppressWarnings(library(edgeR, quietly = T))
args = commandArgs(trailingOnly=TRUE)

counts = fread("./F_PCGlnc.gene.count.txt", data.table=F)
row.names(counts) = counts$Geneid
counts[,1] = NULL
DF=as.data.frame(t(counts))
sheep_meta=as.data.frame(fread("./samplelist_final.txt", header=T, sep="\t"))
tissues=unique(sheep_meta$Tissue_s)
sheep_exp=DF
tissues= args[1]
tissue_cata<-unique(sheep_meta$Tissue_categories[sheep_meta$Tissue_s==tissues])
tissue_cata_sample <- sheep_meta$samples[sheep_meta$Tissue_categories==tissue_cata]
Tissue_s_sample<-sheep_meta$samples[sheep_meta$Tissue_s==tissues]
remove_sample<-tissue_cata_sample[!tissue_cata_sample%in%Tissue_s_sample]
sheep_meta_data<-sheep_meta[!sheep_meta$samples%in%remove_sample,]
Meta_data_sheep_sort<-sheep_meta_data[order(sheep_meta_data$Tissue_s),]
sheep_expression<-t(sheep_exp)
sheep_expression[1:10,1:10]
sheep_expression_sort<-sheep_expression[, match(Meta_data_sheep_sort$samples, colnames(sheep_expression))]
sheep_expression_sort[1:10,1:10]
tissue_index<-which(Meta_data_sheep_sort$Tissue_s==tissues)
df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_sheep_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
table(df)
design <- model.matrix(~-1+factor(df))
dim(design)
colnames(design) <- c("Others","Tissue")
table(design)
rownames(design)=colnames(sheep_expression_sort)
conditions = factor((as.data.frame(design))$Tissue) 

y <- DGEList(counts=sheep_expression_sort,group=conditions)
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)
write.table(count_norm, file = "./count_norm.txt", row.names = TRUE)

pvalues <- sapply(1:nrow(count_norm),function(i){
data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
p=wilcox.test(gene~conditions, data)$p.value
return(p)
})

fdr=p.adjust(pvalues,method = "fdr")
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.table(outRst, file=paste("tissue_specificity_Wilcoxon_rank_sum/", args[1], ".txt", sep=""), sep="\t", quote=F, row.names = T)
