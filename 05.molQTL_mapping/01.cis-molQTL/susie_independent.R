library(data.table)
library(ggplot2)
library(ggpubr)

#setwd('E:\\OneDrive - 西北农林科技大学\\备份\\博后\\GTEx\\results\\03.molQTL\\eQTL')
# https://github.com/gandallab/devBrain_xQTL/blob/master/code/04-QTL/cis-eQTL/susie.ipynb
ARGS <- commandArgs(trailingOnly = TRUE)
infile = ARGS[1] # columns: tissue	phenotype_id	cs_count	num_rank
outfile = ARGS[2]
tis = ARGS[3]

type = strsplit(outfile, split = "\\.")[[1]][1]
df <- fread(infile)
if (!is.na(tis)) { df <- subset(df, tissue == tis) }
ggplot(df, aes(x = cs_count, y = num_rank)) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  geom_count(color = "#648FFF") +
  #   geom_count(aes(color = ..n..), size=5) +
  scale_size(range = c(1, 20), "Gene count") + 
  #     scale_size_area() +
  #   geom_hex() +
  labs(x = "Number of SuSiE credible sets", y = paste0("Number of independent cis-",type, "s")) +
  #   ggtitle("7491 genes in common") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_x_continuous(limits = c(1, 10), breaks = seq(1,10,by=1)) +
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,by=1)) +
  stat_cor(method = "spearman", label.x = 1.25, label.y = 10, cor.coef.name = "rho", size = 6)
ggsave(outfile, width = 8, height = 8)
