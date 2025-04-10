

##################################bar#############################
library("pheatmap")

data <- read.table("celltypedata.txt",sep = "\t",header = T,row.names = 1)
data <- data[rowSums(data) > 1,]
data<-as.matrix(t(data))

pheatmap(data,
         cluster_cols = F,
         cluster_rows = F,
         legend = F,
         color = c("#f7e9ed","#f6b6c6"),
         border = "white",
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = c("90")
         )

###############################################################
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(tidyverse)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
setwd("E://InPut//Pred")
list <- list.files()
df1=data.frame()
for(cancer in CANCER){
  data = read.csv(paste0(cancer,"_Pred.csv"),header=T)
  df <-data.frame(Composition=colMeans(data))
  df$cancer = cancer
  df$CellType=colnames(data)
  df1 = rbind(df1,df)
}
color_list <- c("#ffe0e9", "#4c8dae", "#8ca1c4","#2e4e7e","#f47983", "#ff4777","#d07b76","#ff9b60",
                "#e8ac9d","#ffb1c9","#48c0a3", "#a4e2c6", 
                "#f2ecde",
                "#FFfcc6", "#ffeda8","#fff897","#ecfc00")
df1$CellType <- factor(df1$CellType, levels = 
             c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono.Macro","NK",
               "Neutrophils","Plasma","Tprolif", "Treg",
               "Malignant",
               "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts"
                ))
p1 = ggbarplot(df1, x = "cancer", y = "Composition",
          size = 0, fill = "CellType", color = "transparent") +
  theme(legend.position = "bottom") +
  theme(legend.position = "right") +
  labs(x = "", y = "", title = "") +
  coord_flip() +
  scale_fill_manual(values = setNames(color_list, levels(df1$CellType)))


####### ######
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



######## marker#########
library("readxl")
library("data.table")
library(dplyr)
library(stringr)
library(ggplot2)
library(rstatix)
library(ggpubr)

genes_to_check = c(
    'CD19',"MS4A1","CD79A", #B
    'MZB1',"IGHG1","JCHAIN", #plasma
    'CD3D','CD3E','CD3G','CD4',"IG7R","FOXP3", #T
    "CD1C", "ITGAE","ITGAM",#DC
    'CPA3' ,'KIT' ,"TPSB2",#mast
    "ITGAX" ,"CD14","CD68", #MONO/MACRO
    "ACTA2","TNS1","CDH11",#Myofibroblasts
    'PECAM1', 'CLDN5','RAMP2', #end
    'PGC', 'PGA3','MUC5AC', #epi
    'COL1A1',"PGC","PDGFRA",#fib
    "PRF1","KLRD1","NKG2D", #NK
    "FCGR3B","FUT3","CEACAM8"#nei
 )

input<-"E:\\InPut/fjj"
fjj<-list.files(input,pattern = "\\.csv$",full.names = T)
df1<-data.frame()
for (file in fjj) {
  cancer = strsplit(basename(file), "_")[[1]][1]
  cell = str_extract(file, "(?<=_).+(?=\\.)")
  data<- fread(file,header = T,sep = ",")
  rownames(data) <- data[[1]] 
  data <- data[, -1]          
  subdata <- data %>%slice_sample(n = 40) 
  subdata <-  as.data.frame(subdata)
  missing_genes <- setdiff(genes_to_check, colnames(subdata))
  if (length(missing_genes) > 0) { for (gene in missing_genes) { subdata[[gene]] <- 0  } }
  subdata <-  subdata[,genes_to_check]
  subdata$cell = cell
  df1 = rbind(df1,subdata)
}
rio::export(df1,"E:\\InPut/fjjsub_marker_exp.txt")
data_long <- melt(df1, id.vars = "cell")   
colnames(data_long) <- c("cell_type", "gene", "expression")

gene_keep = c(
  'CD19',"MS4A1","CD79A", #B
  'MZB1', #plasma
  'CD3D','CD3E','CD3G',"FOXP3", #T
  "KLRD1", #NK
  "CD1C", #DC
  'CPA3' ,'KIT' ,#mast
  "ITGAX" ,"CD68", #MONO/MACRO
  "FCGR3B",#NEI
  'PECAM1', 'CLDN5','RAMP2', #end
  'PGC', #epi
  'COL1A1',"PDGFRA",#fib
  "ACTA2"#Myofibroblasts
)
mydata <-data_long[data_long$gene%in%gene_keep,]
mydata$cell_type <- factor(mydata$cell_type, levels = c("B","Plasma","CD4Tconv","CD8T", "CD8Tex", "Tprolif", "Treg",
                                                        "NK","DC","Mast", "Mono_Macro","Neutrophils" ,
                                                        "Endothelial", "Epithelial","Malignant",
                                                        "Fibroblasts", "Myofibroblasts"))
mydata$gene <- factor(mydata$gene, levels = rev(gene_keep))

result <- mydata %>%group_by(cell_type, gene) %>%
  summarise(mean_exp = mean(expression, na.rm = TRUE))%>%
  mutate(sig=case_when(mean_exp >=quantile(mean_exp, 0.9, na.rm = TRUE)  ~ "up", T~ "down"))


p=ggplot(result, aes(x = cell_type, y = gene, size = mean_exp, color=mean_exp)) + 
  geom_point()+
  geom_vline(aes(xintercept = 2.5), linetype = "dashed", color = "grey40") + 
  geom_vline(aes(xintercept = 8.5), linetype = "dashed", color = "grey40") +
  geom_vline(aes(xintercept = 10.5), linetype = "dashed", color = "grey40") +  
  geom_vline(aes(xintercept = 12.5), linetype = "dashed", color = "grey40") + 
  geom_vline(aes(xintercept = 15.5), linetype = "dashed", color = "grey40") + 
  geom_hline(aes(yintercept = 3.5), linetype = "dashed", color = "grey40") +  
  geom_hline(aes(yintercept = 7.5), linetype = "dashed", color = "grey40") + 
  geom_hline(aes(yintercept = 10.5), linetype = "dashed", color = "grey40") +
  geom_hline(aes(yintercept = 13.5), linetype = "dashed", color = "grey40") + 
  geom_hline(aes(yintercept = 14.5), linetype = "dashed", color = "grey40") +  
  geom_hline(aes(yintercept = 18.5), linetype = "dashed", color = "grey40") + 
  theme_classic() +
  scale_color_gradient(low = "#d6ecf0",high="#f20c00")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

