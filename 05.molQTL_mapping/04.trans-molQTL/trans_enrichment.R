library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

setwd('C:/Users/au783873/OneDrive - 西北农林科技大学/备份/博后/GTEx/results/03.molQTL/eQTL')

genelist <- fread("trans_hotsplot.genelist", header = FALSE)

go_results <- enrichGO(
  gene          = genelist$V1,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable = T)

go_results = as.data.frame(go_results)
write.csv(as.data.frame(go_results), "trans_hotsplot.go.csv", row.names = FALSE)
