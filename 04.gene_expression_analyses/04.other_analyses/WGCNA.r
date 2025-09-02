#!/usr/bin/Rscript
###1.WGCNA
library(WGCNA)
library(ggplot2)
library(data.table)
library(dplyr)
options(stringsAsFactors = FALSE)
library(GO.db)
library(ggplot2)
library(org.Hs.eg.db)
library("clusterProfiler")
args = commandArgs(trailingOnly=TRUE)
##data
datExpr0 = fread(file=(paste0("../",args[1],"/results/",args[1],"_adjust_TMM.txt")), header = T, data.table = F)
row.names(datExpr0) = datExpr0[,1]
datExpr0[,1] = NULL
std = apply(datExpr0, 1, sd)
datExpr0 = datExpr0[std!=0,]
datExpr0 = as.data.frame(t(datExpr0))
num_datExpr0 = datExpr0 %>% mutate(across(1:ncol(datExpr0), as.numeric))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
        if (sum(!gsg$goodGenes)>0)
                printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0)
                printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
                datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

png(file.path(args[1],paste0("01.softpower_",args[1],".png")),width=10,height=6,unit="in",res=600)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
        main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red");
        abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
        main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
while (!is.null(dev.list())) dev.off()
save(datExpr0, file = file.path(args[1],paste0("01.datExpr0_",args[1],".RData")))

cor <- WGCNA::cor
net = blockwiseModules(datExpr0,
                       power = 6, TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       saveTOMs = F,corType = "pearson",
                       saveTOMFileBase = "TOM-blockwise",
                       verbose = 3)
cor<-stats::cor
## Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms;
save(net, MEs, moduleLabels, moduleColors, geneTree,file=file.path(args[1], paste0("02.net_",args[1],".RData")))

module_gene<-data.frame(geneID=colnames(datExpr0),module=moduleColors,MEs=paste0("ME",net$colors))
write.csv(module_gene,file.path(args[1], paste0("02.gene_module_",args[1],".csv")),row.names=F)

df<-as.data.frame(table(module_gene$module))
colnames(df)<-c("Module","Gene_counts")
write.csv(df,file.path(args[1],paste0("02.module_gene_counts_",args[1],".csv")),row.names=F)
pdf(file.path(args[1],paste0("02.module_gene_counts_",args[1],".pdf")),width=20,height=15)

p<-ggplot(data = df, mapping = aes(x = Module,y=Gene_counts)) + geom_bar(stat = 'identity')+
        geom_text(mapping = aes(label = Gene_counts),vjust = -1)+
        theme(axis.text.x = element_text(angle = 45,vjust = 0.5,size = 12),plot.margin = margin(t = 1,r = 1,b = 1,l = 1,unit = "cm"))
print(p)
while (!is.null(dev.list()))  dev.off()


pdf(file.path(args[1],file=paste0("02.Eigengene_heatmap_",args[1],".pdf")),width=20,height=20)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
while (!is.null(dev.list()))  dev.off()


res_mat = data.frame(MEs)
sample_order = colnames(t(datExpr0))
row.names(res_mat) = sample_order
res_mat = res_mat[,colnames(res_mat)!="MEgrey"]
write.csv(res_mat, file.path(args[1],paste0("02_separate_eigen_",args[1], ".csv")), quote = F)

for (j in unique(module_gene$module)) {
        row.names(module_gene) <- module_gene[,1]
        module <- module_gene[which(row.names(module_gene) %in% module_gene$geneID[module_gene$module == j]),]
        module1<- module[,1]
        module1 <- as.character(module1)
        module_trance <- bitr(module1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T) #ID转换
        module_result_go <- enrichGO(module_trance$ENTREZID,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",ont = "ALL",readable = T) #GO分析
        module_result <- as.data.frame(module_result_go)
        a = paste0("03.gene_annotation_", j, ".xls")
        write.table(module_result,file.path(args[1],a),quote = F,sep = "\t")
###prepare files of cytoscape
    probes = colnames(datExpr0)
    inModule = (moduleColors==j);
    modProbes = probes[inModule];
    ### Recalculate topological overlap
    TOM = TOMsimilarityFromExpr(datExpr0, power = sft$power);
    #Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = file.path(args[1],paste("04_CytoscapeInput-edges-", paste(j, collapse="-"), ".txt", sep="")),
      nodeFile = file.path(args[1],paste("04_CytoscapeInput-nodes-", paste(j, collapse="-"), ".txt", sep="")),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes,
      nodeAttr = moduleColors[inModule])
}

###2.unanno
library(org.Hs.eg.db)
library("clusterProfiler")
result_df <- data.frame()
args = commandArgs(trailingOnly=TRUE)

read.table (paste0(args[1],"/02.gene_module_",args[1],".csv"), header = T, sep = ",", comment.char = "") -> gene_module
gene_module=as.data.frame(gene_module)
row.names(gene_module)=gene_module[,1]
for (t in unique(gene_module$module)) {
  gene_list <- gene_module %>% filter((geneID %in% gene_module$geneID[gene_module$module == t]))
  module_trance <- bitr(gene_list$geneID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = T) ##
  annot <- AnnotationDbi::select(org.Hs.eg.db, keys = module_trance$ENTREZID, columns = "GO") ##
  annot_gene_count <- length(unique(annot[,1]))
  total_gene_count <- nrow(gene_list)
  annot_ratio <- annot_gene_count / total_gene_count
  unannot_ratio <- (total_gene_count - annot_gene_count) / total_gene_count
  module_result <- data.frame(tissue = args[1],module = t, annot_count = annot_gene_count, total_count = total_gene_count, annot_ratio = annot_ratio, unannot_ratio = unannot_ratio)
  result_df <- rbind(result_df, module_result)
  unanno_gene <- as.data.frame(setdiff(gene_module$geneID, module_trance$SYMBOL))
  write.csv(unanno_gene, file.path(args[1],paste0("05.unannot_",t,".csv")), quote = F, row.names = F, sep = "\t")
#  write.csv(result_df, file.path(args[1],paste0("00.annot_",args[1],".csv")), quote = F, row.names = F, sep = "\t")
}

###3.module count
library(ggplot2)
color <- read.table("tissue_color.txt", header = T, sep = "\t", comment.char = "")
Tissues <- unique(color$Tissue)
color$Tissue <- factor(color$Tissue, levels=Tissues)
Tissuecolors <- color[match(Tissues,color$Tissue),"Color"]
read.table ("03.module_gene_counts.txt", header = F, sep = "\t", comment.char = "") -> df
colnames(df)=c("count","tissue","aa")
df$tissue <- factor(df$tissue, levels = Tissues)
df$aa <- as.factor(df$aa)
p<-ggplot(data = df, mapping = aes(x = tissue,y=count,fill=tissue)) + geom_bar(stat = 'identity')+
  geom_text(mapping = aes(label = count),vjust = 1.5)+
  ylab("Number of modules")+
  xlab("Tissue")+
  theme_classic() +
  scale_fill_manual(values = Tissuecolors)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 10))+
  theme(axis.text.y = element_text(size = 10,hjust = 1)) +
  theme(axis.title = element_text(size = 20))
pdf("03.module_gene_counts.pdf",width=20,height=15)
print(p)
while (!is.null(dev.list()))  dev.off()

###4.gene count
library(ggplot2)
library(ggpubr)
color <- read.table("tissue_color.txt", header = T, sep = "\t", comment.char = "")
Tissues <- unique(color$Tissue)
color$Tissue <- factor(color$Tissue, levels=Tissues)
Tissuecolors <- color[match(Tissues,color$Tissue),"Color"]
data = read.table("04.all_tissue_module_gene.txt", sep = "\t")
colnames(data)= c("modul", "N_number", "tissue")
data$tissue = factor(data$tissue, levels=Tissues)

p<-ggplot(data, aes(x=tissue, y=N_number,color=tissue)) +
        geom_point() +
        geom_boxplot(aes(group=tissue, x=tissue), size=0.8, width=0.8, alpha=0) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 10,hjust = 1,angle = 90)) +
        theme(axis.text.y = element_text(size = 12,hjust = 1)) +
        theme(axis.title = element_text(size = 20))+
        ylab("Number of genes")+
        xlab("Tissue")+
        scale_y_continuous(breaks=seq(0,20000,1000))+
        scale_color_manual(values=c(Tissuecolors))
pdf("04.gene_counts.pdf",width=20,height=15)
print(p)
while (!is.null(dev.list()))  dev.off()
