# load library
library(data.table)
library("variancePartition")

# https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html
# load data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
#data(varPartData)

ARGS <- commandArgs(trailingOnly = TRUE)
geneExpr = ARGS[1]
info = ARGS[2]
out_rdata_prefix = ARGS[3] 


geneExpr = as.data.frame(fread(geneExpr, sep = "\t", header = T))
rowname_geneExpr <- geneExpr[,1]
geneExpr = geneExpr[,-c(1:6)]
geneExpr[geneExpr < 0] <- 0

info = as.data.frame(fread(info, sep = "\t", header = T))
rownames(info) <- info$ID
info[c("ID", "Individual", "Tissue", "Breed", "Group", "Sex", "Stage", "BioProject")] <- lapply(info[c("ID", "Individual", "Tissue", "Breed", "Group", "Sex", "Stage", "BioProject")], factor)

# calculate Coefficient of Variation (CV), log scale and remove rows with variance = 0
rownames(geneExpr) <- rowname_geneExpr
geneExpr <- geneExpr[, rownames(info)]
geneExpr_mean = as.data.frame(apply(geneExpr, MARGIN=1, mean))
geneExpr_sd = as.data.frame(apply(geneExpr, MARGIN=1, sd))
geneExpr_cv = geneExpr_sd / geneExpr_mean
names(geneExpr_mean)[1] <- "mean"
names(geneExpr_sd)[1] <- "sd"
names(geneExpr_cv)[1] <- "cv"
geneExpr = log(geneExpr+0.25)
geneExpr <- geneExpr[geneExpr_sd != 0, ]


# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
#form <- ~ Age + (1 | Individual) + (1 | Tissue) + (1 | Batch)
form <- ~ (1 | Tissue) + (1 | BioProject) + (1 | Individual) + (1 | Breed) + (1 | Sex) + (1 | Stage)
form0 <- ~ Tissue + BioProject + Individual + Breed + Sex + Stage
# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form0, info)

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,
#     a linear mixed model is used
# If all variables are modeled as fixed effects,
#       a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable
#       to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix,
# and EList output by voom() in the limma package,
# or an ExpressionSet
varPart <- fitExtractVarPartModel(geneExpr, form, info)

# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
vp <- merge(vp, geneExpr_mean, by = "row.names", all = FALSE)
row.names(vp) <- vp$Row.names
vp <- vp[, -1]
vp <- merge(vp, geneExpr_sd, by = "row.names", all = FALSE)
row.names(vp) <- vp$Row.names
vp <- vp[, -1]
vp <- merge(vp, geneExpr_cv, by = "row.names", all = FALSE)
row.names(vp) <- vp$Row.names
vp <- vp[, -1]

rm(geneExpr)
out_rdata = paste0(out_rdata_prefix, ".RData")
save.image(out_rdata)

#load(out_rdata)

# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1:10, 1:7])

# Plot correlation matrix
# between all pairs of variables
#plotCorrMatrix(C)

# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart(vp[, 1:7]) + theme_classic()
#fig <- plotVarPart(vp)
#ggsave(file, fig)

# variance partition
annot = as.data.frame(fread("gene.annot", sep = "\t", header = F, col.names=c("gene","annotation")))
rownames(annot) <- annot[,1]
#annot <- annot[, -1]
vp <- merge(vp, annot, by = "row.names", all = FALSE)
row.names(vp) <- vp$Row.names
vp <- vp[, -1]

ggplot(vp) +
  geom_point(aes(x = Tissue, y = Breed, color = annotation), size = 0.2) +
  geom_point(data = subset(vp, annotation == 'lncRNA'), aes(x = Tissue, y = Sex, color = annotation), size = 0.2) +
  theme_classic() +
  xlab("Tissue variance") + xlim(0, 1) +
  ylab("Breed variance") + ylim(0,1) +
  scale_color_manual(values = c( "protein_coding" = "#d99433", "lncRNA" = "#64a846"))

ggplot(vp) +
  geom_point(aes(x = Tissue, y = Individual, color = annotation), size = 0.2) +
  geom_point(data = subset(vp, annotation == 'lncRNA'), aes(x = Tissue, y = Individual, color = annotation), size = 0.2) +
  theme_classic() +
  xlab("Tissue variance") + xlim(0, 1) +
  ylab("Individual variance") + ylim(0,1) +
  scale_color_manual(values = c( "protein_coding" = "#d99433", "lncRNA" = "#64a846"))

# extract top and bottom 1000 genes for GO enrichment
vp <- vp[order(vp$cv,decreasing=T),]
top = rownames(vp[1:2000, ])
write.table(as.data.frame(top), file = "top_variable_gene.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(top), file = "top_variable_gene.rmbatch.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
bottom <- rownames(vp[(nrow(vp) - 2000 + 1):nrow(vp), ])
write.table(as.data.frame(bottom), file = "bottom_variable_gene.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(bottom), file = "bottom_variable_gene.rmbatch.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)

go_top <- enrichGO(
  gene          = top,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable = T)
go_bottom <- enrichGO(
  gene          = bottom,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable = T)

go_top = as.data.frame(go_top)
go_bottom = as.data.frame(go_bottom)
write.csv(as.data.frame(go_top), "go_top.csv", row.names = FALSE)
write.csv(as.data.frame(go_bottom), "go_bottom.csv", row.names = FALSE)

# plot
#go_top <- fread("go_top.csv", header = TRUE)
go_top <- go_top[order(go_top$FoldEnrichment,decreasing=T),]
go_top_top10 = go_top[1:10, ]
#go_bottom <- fread("go_bottom.csv", header = TRUE)
go_bottom <- go_bottom[order(go_bottom$FoldEnrichment,decreasing=T),]
go_bottom_top10 = go_bottom[1:10, ]

p1 <- #ggplot(go_top_top10,aes(GeneRatio,fct_reorder(factor(Description),GeneRatio))) + 
  ggplot(go_top_top10,aes(FoldEnrichment,fct_reorder(factor(Description),FoldEnrichment))) + 
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_color_gradient(low="blue",high ="red") + 
  labs(color=expression(-log[10](p.adjust)),size="Count", shape="Ontology",
       x="Fold Enrichment",y="GO term",title="GO enrichment of high variance genes") + 
  theme_bw()
p2 <- #ggplot(go_bottom_top10,aes(GeneRatio,fct_reorder(factor(Description),GeneRatio))) + 
  ggplot(go_bottom_top10,aes(FoldEnrichment,fct_reorder(factor(Description),FoldEnrichment))) + 
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_color_gradient(low="blue",high ="red") + 
  labs(color=expression(-log[10](p.adjust)),size="Count", shape="Ontology",
       x="Fold Enrichment",y="GO term",title="GO enrichment of low variance genes") + 
  theme_bw()

library("forcats")
p1 / p2