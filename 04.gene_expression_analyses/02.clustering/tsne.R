library("data.table")
library(Rtsne)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)

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
DF <- t(DF[1:5000,])

##tsne
tsne <- Rtsne(DF, theta=0, partial_pca=T, verbose=T, check_duplicates=F, num_threads=4)
tSNE_df <- tsne$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2")
tSNE_df$sample <- rownames(DF)
rm(data, indivs, DF, DF_sd)
save.image("tSNE.RData")

load("tSNE.RData")

## color
color <- read.table("sample_color.txt", header = T, sep = "\t", comment.char = "")
Systems <- unique(color$system)
Tissues <- unique(color$type)
color$system <- factor(color$system, levels=Systems)
color$type <- factor(color$type, levels=Tissues)
## combine
tSNE_reuslt <- inner_join(tSNE_df, color, by = "sample")
tSNE_reuslt$system = factor(tSNE_reuslt$system, levels=Systems)
tSNE_reuslt$type = factor(tSNE_reuslt$type, levels=Tissues)

## plot
p <- ggplot(tSNE_reuslt, aes(label=paste0(sample," ",system," ",type)))+
  geom_point(aes(x=tSNE1,y=tSNE2,color=type),size=1)+
  scale_color_manual(values=c(
    "Heart"="#bc58e3",
    "Artery"="#a82ed7",
    "Rumen"="#fc9891",
    "Reticulum"="#d97c68",
    "Omasum"="#f8ae81",
    "Abomasum"="#f18264",
    "Small intestine"="#cb8119",
    "Duodenum"="#eb9d63",
    "Jejunum"="#eb951c",
    "Ileum"="#ce9639",
    "Large intestine"="#c2723d",
    "Cecum"="#daaa6c",
    "Colon"="#f2c063",
    "Rectum"="#efe0a8",
    "Liver"="#ad8c8b",
    "Embryo"="#9fe2d8",
    "Adipose"="#ffc4f1",
    "Subcutaneous adipose"="#fb9df1",
    "Intermuscular adipose"="#d083c8",
    "Perirenal adipose"="#d974cf",
    "Caul adipose"="#B849A0",
    "Tail adipose"="#db8fca",
    "Pericardial adipose"="#e01fae",
    "Adrenal"="#a15866",
    "Thyroid"="#f359d1",
    "Mammary gland"="#fed9d0",
    "Ovary"="#69d28c",
    "Ovarian follicle"="#c1ecd0",
    "Corpus luteum"="#88eca9",
    "Granulosa cell"="#badac5",
    "Oocyte"="#58a75d",
    "Uterus"="#8bc774",
    "Cervix"="#b6d7a9",
    "Oviduct"="#79ffaa",
    "Placentome"="#54ab89",
    "Blood"="#e11b1b",
    "PBMC"="#ec1360",
    "Bone marrow"="#d84b4b",
    "Lymph node"="#c33a11",
    "Spleen"="#962932",
    "Thymus"="#ff3a32",
    "Milk"="#ddc3bd",
    "Skin"="#d09dc5",
    "Hair follicle"="#ef9cd8",
    "Horn"="#ba4567",
    "Soft horn"="#a25d73",
    "Testis"="#7ef351",
    "Epididymis"="#70d24b",
    "Muscle"="#a180ca",
    "Longissimus muscle"="#754fa4",
    "Infraspinatus muscle"="#ba99e2",
    "Biceps muscle"="#8232cd",
    "Skeletal muscle"="#9800ff",
    "Brain"="#fcf61e",
    "Cerebrum"="#d4dc23",
    "Cerebral cortex"="#dcd71a",
    "Cerebral medulla"="#f2c905",
    "Hippocampus"="#f1cb05",
    "Pineal"="#807120",
    "Cerebellum"="#efd80b",
    "Hypothalamus"="#825e19",
    "Pituitary"="#f1d95d",
    "Adenohypophysis"="#eddf61",
    "Brainstem"="#f4d578",
    "Medulla oblongata"="#d0b35b",
    "Optic chiasm"="#fece01",
    "Optic nerve"="#fee062",
    "Trachea"="#afe0f7",
    "Lung"="#36b5f1",
    "Periosteum"="#E2EFDA",
    "Kidney"="#4f3136",
    "Salivary gland"="#fdbbb7",
    "Esophagus"="#fb857d",
    "Pancreatic islet"="#bf0092",
    "Amnion"="#bce718",
    "Lymph"="#ef6a43",
    "Tonsil"="#c85b37",
    "Retina"="#a98156",
    "Bone"="#EDE4D4"))+
  theme(legend.title =element_blank())+labs(x="tSNE1",y="tSNE2") +
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"),
        axis.title.x=element_text(family="Times",size = 15,face="plain"))+#设置y轴标题的字体属性
  theme(legend.text=element_text(family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16))+
  theme(legend.title=element_text(family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=16))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- ggplotly(p)
saveWidget(
  widget = p, #the plotly object
  file = "tsne.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)
