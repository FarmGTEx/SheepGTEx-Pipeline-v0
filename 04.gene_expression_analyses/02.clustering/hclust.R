library("data.table")
library(dplyr)
library(ggplot2)
library("cowplot")
library(ggtree)
library("ggtreeExtra")
library(ggnewscale)

##input data
data <- fread("../sheep.PCGlnc.gene.merged.tpm.txt", sep = "\t", header = T)
data=as.data.frame(data)
##Keep valid data
indivs=as.vector(as.matrix(fread("../select.samplelist", header=F)))
data=data[,-c(1:6)]
data<- data[rowSums(data)>0,]
data = data[,indivs]

##Log transformation of data frame
DF=apply(log(data+0.25), MARGIN = 2, FUN = scale)
DF_sd = apply(DF, MARGIN=1, sd)
DF <- DF[order(DF_sd,decreasing=T),]
res <- cor(DF[1:5000,],method="pearson")
dist <- as.dist(1 - res)

##hclust
hc=hclust(dist, method="complete")
rm(data, indivs, DF, DF_sd, res, dist)
save.image("hclust.RData")

load("hclust.RData")

##group and color
indivs=as.vector(hc$labels)
color <- read.table("sample_color.txt", header = T, sep = "\t", comment.char = "")
color = color[color$sample %in% indivs,]
Systems <- unique(color$system)
Tissues <- unique(color$type)
color$system <- factor(color$system, levels=Systems)
color$type <- factor(color$type, levels=Tissues)
Systemcolors <- color[match(Systems,color$system),"color1"]
Tissuecolors <- color[match(Tissues,color$type),"color2"]
color$newlabel = paste0(color$sample, " ", color$type, " ", color$system)
color$Abundance = abs(rnorm(n=length(hc$labels), mean=10, sd=0.00000001))

#https://github.com/YuLab-SMU/ggtreeExtra/issues/2
##circular tree
###plot tree
p <- open_tree(ggtree(hc, layout="circular", color="#6C6E72") , 10)
### add label annotation
#p <- p +
#  geom_fruit(
#    data=color[, c("sample", "newlabel")],
#    geom=geom_text,
#    mapping=aes(y=sample, x=newlabel, label=newlabel), 
#    size=0.5, pwidth=0, offset=0.02, hjust=0
#  )
### add system information
p <- p +
  geom_fruit(data=color[, c("sample", "system", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=system),
             pwidth=0.1, offset=0.01,
             stat="identity", width=1,
             orientation="y"
  ) +
  scale_fill_manual(values=Systemcolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+
  labs(fill="System")
### add tissue information
p <- p + 
  new_scale_fill() +
  geom_fruit(data=color[, c("sample", "type", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=type),
             pwidth=0.1, offset=0,
             stat="identity", width=1,
             orientation="y"
  ) +
  scale_fill_manual(values=Tissuecolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
  labs(fill="Tissue/cell type") 

pdf("hclust.cirtree.pdf", width=50, height=50)
print(p)
dev.off()

##normal tree
###plot tree
p = ggtree(hc, hang=0, color="#6C6E72")
###rotate branches in ggtree
p <- p + geom_text(aes(label=node))
p <- p %>% rotate(8201) %>% rotate(8244) %>% rotate(8209) %>% rotate(8263) %>%
  rotate(8181) %>% rotate(8178) %>% rotate(8180) %>% rotate(8212) %>% rotate(8227) %>%
  rotate(8205) %>% rotate(8266) %>% rotate(8204)
### add label only
#p <- p %<+% color +
#  geom_tippoint(aes(colour=type),
#                alpha=0) +
#  geom_tiplab(aes(colour=type),
#              align=TRUE,
#              linetype=3,
#              size=0.5,
#              linesize=0.2,
#              show.legend=TRUE
#)
### add label annotation
p <- p +
  geom_fruit(
    data=color[, c("sample", "newlabel")],
    geom=geom_text,
    mapping=aes(y=sample, x=newlabel, label=newlabel), 
    size=0.5, pwidth=0, offset=0.001, hjust=0
  )
### add system information
p <- p + 
  geom_fruit(data=color[, c("sample", "system", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=system),
             pwidth=0.05, offset=0.01,
             stat="identity",
             orientation="y"
  ) +
  scale_fill_manual(values=Systemcolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+
  labs(fill="System")
### add tissue information
p <- p + 
  new_scale_fill() +
  geom_fruit(data=color[, c("sample", "type", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=type),
             pwidth=0.05, offset=0,
             stat="identity",
             orientation="y"
  ) +
  scale_fill_manual(values=Tissuecolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
  labs(fill="Tissue/cell type") 

pdf("hclust.tree.pdf", width=50, height=100)
print(p)
dev.off()

##dendrogram tree
###plot tree
p = ggtree(hc, hang=0, root.position=8, color="#6C6E72")
###rotate branches in ggtree
p <- p %>% rotate(8201) %>% rotate(8244) %>% rotate(8209) %>% rotate(8263) %>%
  rotate(8181) %>% rotate(8178) %>% rotate(8180) %>% rotate(8212) %>% rotate(8227) %>%
  rotate(8205) %>% rotate(8266) %>% rotate(8204)
#### add label annotation
#p <- p +
#  geom_fruit(
#    data=color[, c("sample", "newlabel")],
#    geom=geom_text,
#    mapping=aes(y=sample, x=newlabel, label=newlabel), 
#    size=0.5, pwidth=0, offset=0.001, hjust=0
#  )
### add system information
p <- p +
  geom_fruit(data=color[, c("sample", "system", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=system),
             pwidth=0.05, offset=0.01,
             stat="identity",
             orientation="y"
  ) +
  scale_fill_manual(values=Systemcolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+
  labs(fill="System")
### add tissue information
p <- p + 
  new_scale_fill() +
  geom_fruit(data=color[, c("sample", "type", "Abundance")],
             geom=geom_bar,
             mapping=aes(y=sample, x=Abundance, fill=type),
             pwidth=0.05, offset=0,
             stat="identity",
             orientation="y"
  ) +
  scale_fill_manual(values=Tissuecolors,
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
  labs(fill="Tissue/cell type")

p <- ggdraw()+
  draw_plot(ggplotify::as.ggplot(p, angle = -90),
            width = 0.7,height = 0.7,
            hjust = -0.2,
            vjust = -0.1)

pdf("hclust.dentree.pdf", width=15, height=15)
print(p)
dev.off()

