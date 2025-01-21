data <- read.table("celltypedata.txt",sep = "\t",header = T,row.names = 1)
data <- data[rowSums(data) > 1,]
data<-as.matrix(t(data))
##################################mosaic+bar#############################
colnames(data) <- paste0(colnames(data)," cell")
library("pheatmap")
pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         color = c("#f7e9ed","#f6b6c6"),
         border = "white",
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = c("45")
         )

df <- data.frame(celltype = colnames(data),nmatrix = colSums(data))
library(ggplot2)
ggplot(df, aes(celltype, nmatrix)) +
  geom_col(position= position_dodge(0),fill = "#DDF1D0") +
  theme_classic() + 
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        # axis.line = element_blank(),
        # axis.ticks =  element_blank(),
        legend.position = "none")+
  labs(x="",y="")


df <- data.frame(cancer = rownames(data),nmatrix = rowSums(data))
ggplot(df, aes(cancer, nmatrix)) +
  geom_col(fill = "#DDF1D0") +
  theme_classic() + 
  theme(
        # axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.line = element_blank(),
        # axis.ticks =  element_blank(),
        legend.position = "none")+
  labs(x="",y="")+
  coord_flip()

##################################Deconvolved cell proportion#############################
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
list <- list.files()
df1=data.frame()
for(cancer in CANCER[5:23]){
  data = read.table(paste0("TAPE_",cancer,".txt"),header=F,sep = "\t")
  tumor<-read.table(paste0("TCGA_",cancer,"_CS.txt"), 
                    header = TRUE,row.names = 1,sep="\t")
  df <-data.frame(Composition=colMeans(data))
  df$cancer = cancer
  df$CellType=names(tumor)[order(names(tumor))]
  df1 = rbind(df1,df)
}
# write.table(df1,"all.cell.frection.txt",quote = F,sep = "\t")
library(tidyverse)
df1 = read.table("all.cell.frection.txt",header = T,sep = "\t")
df1 <- df1[order(df1$cancer,decreasing = T), ]
color_list <- c("#b2cbe6", "#9FA8DA",
  "#FCE4EC", "#FFB3C8","#e6a0af","#EF9A9A", "#E57373","#EF5350","#B71C1C",
  "#EC407A", "#C2185B","#880E4F","#CE93D8",  "#cb86b5","#9C27B0", "#7B1FA2", 
  "#D4F69F", "#A1FFA1","#91FF91","#76d4a4", "#81C784", "#2E7D32",
   "#FFFDE7","#FFECB3","#DCE775","#C0CA33", "#FFEE58","#FDD835","#F9A825",
  "#B2EBF2", "#4DD0E1","#00BCD4","#0097A7", "#006064",
   "#D7CCC8","#A1887F"
  )

df1$CellType <- factor(df1$CellType, levels = 
 c("AC_Like_Malignant", "Malignant",
   "B","Mast", "MonoMacro","Neutrophils", "NK", "pDC","Plasma",
   "CD4Tconv", "CD8T", "CD8Tex", "DC",  "Promonocyte","Tprolif", "Treg",
   "Pericytes", "Astrocyte","Neuron","Oligodendrocyte", "SMC","Vascular",
   "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts","Myocyte","Acinar", "Ductal",
   "Hepatic progenitor","HSC", "Progenitor","OPC", "Endometrial stromal cells",
    "EryPro", "GMP"
            ))
p1 = ggbarplot(df1, x = "cancer", y = "Composition",
          size = 0, fill = "CellType", color = "transparent") +
  theme(legend.position = "bottom") +
  theme(legend.position = "right") +
  labs(x = "", y = "", title = "") +
  coord_flip() +
  scale_fill_manual(values = setNames(color_list, levels(df1$CellType)))
p1
####### cor verify ########
library(tidyverse)
library(ggpubr)
library(ggplot2)
cor <- read.table("cor.txt", sep="\t", header=T)
ggplot(cor) + 
  geom_col(aes(x = TCGA, y = -log(p), fill = TCGA),
           position = "dodge2",
           show.legend = TRUE,
           alpha = 0.9) +
  coord_polar() +
  scale_fill_manual(values = c("#CCCCFF", "#FF6666", "#FEBB81", "#CCFF00",
                               "#FFCCCC", "#FF99CC", "#D44842", "#CCCC99",
                               "#99CC66", "#66CC99", "#CC99CC")) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1)
    )
  ) +labs(x="")+
  theme_minimal() +
  theme(
    # axis.title = element_blank(),
    # axis.ticks = element_blank(),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(color = "gray12", size = 12),
    legend.position = "none"
  )
########marker verify#########
library("readxl")
library("data.table")
marker = read_excel("LUSC_marker.xlsx",col_names  =T)
markergene =data.frame(marker[!duplicated(marker$gene),])

setwd("F:/Desktop/毕业设计/inputdata/TAPE-tumor")

list <- list.files()
list=list[grep("LUSC",list)]
data1 = fread(list[1],header = T,sep = ",")
markergene= intersect(colnames(data1),markergene$gene)#65marker
marker=marker[marker$gene %in% markergene,]
marker=marker[!duplicated(marker$gene),]


list=list[c(1,5,6,8)]
data=data.frame()
for (i in list) {
  data1 = fread(i, header = T, sep = ",")
  data1 = data.frame(data1)[,marker$gene]
  sample_rows <- sample(1:nrow(data1), 100)
  data1_sampled <- data1[sample_rows, ]
  data = rbind(data, data1_sampled)
}

data1 = data.frame(t(data))
data1 <- data.frame(lapply(data1, as.numeric))

rownames(data1)=marker$gene
colnames(data1) = paste0("Sample", 1:ncol(data1))

library(vegan) 
data2 <- decostand(data1, "standardize", MARGIN = 2)
annotation_row = data.frame(CellType = factor(marker$cell))
rownames(annotation_row) = marker$gene
CellType_color <- c("#FCE4EC","#D4F69F", "#A1FFA1","#e6a0af") 
names(CellType_color) <- c("B","Endothelial","Epithelial","MonoMacro") 
celltype=rep(c("B","Endothelial","Epithelial","MonoMacro"),each=100)
annotation_col = data.frame(CellType = factor(celltype))
rownames(annotation_col) = paste0("Sample", 1:ncol(data2))        

ann_colors <- list(CellType=CellType_color)
library(pheatmap)
pheatmap(data2,
         gaps_col = c(100,200,300),
         gaps_row = c(7,8,20),
         cluster_rows = F, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         border = F,
         legend = T,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#FFFBFC",  "#b14e53"))(10)
)
