

##################################反卷积马赛克+bar#############################
library("pheatmap")

data <- read.table("F:/Desktop/HMU/2024文章/毕业设计/inputdata/celltypedata.txt",
                   sep = "\t",header = T,row.names = 1)
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

# df <- data.frame(celltype = colnames(data),nmatrix = colSums(data))
# library(ggplot2)
# ggplot(df, aes(celltype, nmatrix)) +
#   geom_col(position= position_dodge(0),fill = "grey80") +
#   theme_classic() + 
#   theme(axis.text.x =element_text(angle = 45,hjust = 1),
#         # axis.text.x = element_blank(),
#         # axis.text.y = element_blank(),
#         # axis.line = element_blank(),
#         # axis.ticks =  element_blank(),
#         legend.position = "none")+
#   labs(x="",y="")
##################################反卷积细胞比例#############################
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
                "#e8ac9d","#ffb1c9","#48c0a3", "#a4e2c6",  #免疫细胞12
                "#f2ecde",
                "#FFfcc6", "#ffeda8","#fff897","#ecfc00")
df1$CellType <- factor(df1$CellType, levels = 
             c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono.Macro","NK",
               "Neutrophils","Plasma","Tprolif", "Treg",#免疫细胞10
               "Malignant",
               "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts"#功能性细胞4
                ))
p1 = ggbarplot(df1, x = "cancer", y = "Composition",
          size = 0, fill = "CellType", color = "transparent") +
  theme(legend.position = "bottom") +
  theme(legend.position = "right") +
  labs(x = "", y = "", title = "") +
  coord_flip() +
  scale_fill_manual(values = setNames(color_list, levels(df1$CellType)))
device_size <- dev.size()
width <- device_size[1]  # 获取宽度
height <- device_size[2]  # 获取高度
pdf("F:/Desktop/毕业设计/plot/文章图/R1细胞比例图.pdf",
    width = width, height = height)
p1
dev.off()
###################TAPE免疫细胞比例
# immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono.Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
# immdf1<-df1 
# immdf1$CellType <- as.character(immdf1$CellType)
# immdf1$Composition[is.na(immdf1$CellType) | !immdf1$CellType %in% immune_cell] <- 0
# color_list <- c("#FCE4EC", "#E545C7", "#9C27B0","#7B1FA2","#FFA1A1", "#c3a3b4","#DA0050",
#                 "#a66a53","#FF9500","#FF0000","#8E81ED", "#8D3EF9")
# immdf1$CellType <- factor(immdf1$CellType, levels = 
#                          c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono.Macro","NK","Neutrophils","Plasma","Tprolif", "Treg"))
# ggbarplot(immdf1, x = "cancer", y = "Composition",
#           size = 0, fill = "CellType", color = "transparent",
#           position = position_fill()) +  # 使用 position_fill() 来显示百分比
#   theme(legend.position = "right") +  # 只需要设置一次图例位置
#   labs(x = "", y = "", title = "") +
#   coord_flip() +
#   scale_fill_manual(values = setNames(color_list, levels(immdf1$CellType))) +
#   scale_y_continuous(labels = scales::percent)  # 将 y 轴标签格式化为百分比
# ############################cibersort#############################
# #不能有重复的基因名和缺失/负值数据
# library(CIBERSORT)
# library(data.table)
# library(rio)
# library(dplyr)
# library(limma)
# library(tibble)
# library(patchwork)
# source("F:\\Desktop\\HMU\\2024文章\\毕业设计\\code\\CIBERSORT.R")
# sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
# TCGA<-list.files("E:/InPut/TCGAbulk/",pattern = "\\.RData$",full.names = T)
# for (file in TCGA[11:23]) {
#   cancer=sub(".*-(.*)\\..*", "\\1", basename(file))
#   data<- get(load(paste0("E:/InPut/TCGAbulk/TCGA-",cancer ,".RData")))
#   data <- data[,substr(colnames(data),14,15) < 10]   
#   data <- data[rowSums(data) > 0,]
#   data <- data.frame(avereps(data))
#   results <-cibersort(sig_matrix ,data, perm = 20, QN = F)
#   #perm为置换次数。用于估算单个样本免疫浸润的p值。
#   #QN为分位数归一化。一般芯片数据设置为T，测序数据设置为F。
#   df=data.frame(B=results[,1]+results[,2],
#                 Plasma=results[,3],
#                 CD4Tconv=results[,5]+results[,6]+results[,7]+results[,8],
#                 CD8T=results[,4]+results[,10],
#                 DC=results[,17]+results[,18],
#                 Treg=results[,9],
#                 NK=results[,11]+results[,12],
#                 Mono_Macro=results[,13]+results[,14]+results[,15]+results[,16],
#                 Mast=results[,19]+results[,20],
#                 Neutrophils=results[,22],
#                 Eosinophils=results[,21])
#   write.table(df, file = paste0("F:\\Desktop\\HMU\\2024文章\\毕业设计\\inputdata\\ciber\\",cancer,"_ciber.txt"), sep = "\t", quote = FALSE)
#   message("Saved: ", basename(file))
# }
# #可视化
# setwd("F:\\Desktop\\HMU\\2024文章\\毕业设计\\inputdata\\ciber\\")
# list <- list.files()
# df1=data.frame()
# for(cancer in CANCER[1:23]){
#   data = read.table(paste0(cancer,"_ciber.txt"),header=T)
#   df <-data.frame(Composition=colMeans(data))
#   df$cancer = cancer
#   df$CellType=colnames(data)
#   df1 = rbind(df1,df)
# }
# color_list <- c("#FCE4EC", "#E545C7", "#9C27B0","#7B1FA2","#FFA1A1", "#c3a3b4","#DA0050",
#                 "#783333","#FF9500","#FD002B","#8E81ED", "#8D3EF9",  #免疫细胞12
#                 "pink")
# df1$CellType <- factor(df1$CellType, levels = 
#                          c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono_Macro","NK","Neutrophils","Plasma","Tprolif", "Treg",#免疫细胞10
#                            "Eosinophils"))
# p1 = ggbarplot(df1, x = "cancer", y = "Composition",
#                size = 0, fill = "CellType", color = "transparent") +
#   theme(legend.position = "bottom") +
#   theme(legend.position = "right") +
#   labs(x = "", y = "", title = "") +
#   coord_flip() +
#   scale_fill_manual(values = setNames(color_list, levels(df1$CellType)))
# p1

####### 验证每个癌反卷积的谱之间cor相关性是不是显著的 ######@###
library(tidyverse)
library(ggpubr)
library(ggplot2)
cor <- read.table("F:/Desktop/毕业设计/inputdata/cor.txt", sep="\t", header=T)
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
######## 验证单个癌症免疫细胞marker是否在对应反卷积表达谱高表达#########
# library("readxl")
# library("data.table")
# marker = read_excel("F:/Desktop/毕业设计/inputdata/LUSC_marker.xlsx",col_names  =T)
# markergene =data.frame(marker[!duplicated(marker$gene),])
# 
# setwd("F:/Desktop/毕业设计/inputdata/TAPE-tumor")
# #随便读一个LUSC的反卷积，将自查maker在表达谱中有表达的maker筛选出来
# list <- list.files()
# list=list[grep("LUSC",list)]
# data1 = fread(list[1],header = T,sep = ",")
# markergene= intersect(colnames(data1),markergene$gene)#65marker
# marker=marker[marker$gene %in% markergene,]
# marker=marker[!duplicated(marker$gene),]
# 
# #刚自查的是B endothelial epithelial mono 表达谱取出来整合到一起
# list=list[c(1,5,6,8)]
# data=data.frame()
# for (i in list) {
#   data1 = fread(i, header = T, sep = ",")
#   data1 = data.frame(data1)[,marker$gene]
#   sample_rows <- sample(1:nrow(data1), 100)
#   data1_sampled <- data1[sample_rows, ]
#   data = rbind(data, data1_sampled)
# }
# 
# data1 = data.frame(t(data))
# data1 <- data.frame(lapply(data1, as.numeric))
# 
# rownames(data1)=marker$gene
# colnames(data1) = paste0("Sample", 1:ncol(data1))
# 
# library(vegan) 
# data2 <- decostand(data1, "standardize", MARGIN = 2)#于对数据进行标准化处理
# annotation_row = data.frame(CellType = factor(marker$cell))
# rownames(annotation_row) = marker$gene
# CellType_color <- c("#FCE4EC","#D4F69F", "#A1FFA1","#e6a0af") 
# names(CellType_color) <- c("B","Endothelial","Epithelial","MonoMacro") #类型颜色
# celltype=rep(c("B","Endothelial","Epithelial","MonoMacro"),each=100)
# annotation_col = data.frame(CellType = factor(celltype))
# rownames(annotation_col) = paste0("Sample", 1:ncol(data2))        
# 
# ann_colors <- list(CellType=CellType_color)
# library(pheatmap)
# pheatmap(data2,
#          gaps_col = c(100,200,300),
#          gaps_row = c(7,8,20),
#          cluster_rows = F, 
#          cluster_cols = F,
#          show_rownames = T,
#          show_colnames = F,
#          annotation_row = annotation_row,
#          annotation_col = annotation_col,
#          border = F,
#          legend = T,
#          annotation_colors = ann_colors,
#          color = colorRampPalette(c("#FFFBFC",  "#b14e53"))(10)
# )

######## 验证癌症免疫细胞marker是否在对应反卷积表达谱高表达#########
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
    "FCGR3B","FUT3","CEACAM8"#中性粒
 )

input<-"E:\\InPut/fjj"
fjj<-list.files(input,pattern = "\\.csv$",full.names = T)
df1<-data.frame()
for (file in fjj) {
  cancer = strsplit(basename(file), "_")[[1]][1]
  cell = str_extract(file, "(?<=_).+(?=\\.)")
  data<- fread(file,header = T,sep = ",")
  rownames(data) <- data[[1]]  # 将第一列设置为行名
  data <- data[, -1]           # 删除第一列
  subdata <- data %>%slice_sample(n = 40) 
  subdata <-  as.data.frame(subdata)
  missing_genes <- setdiff(genes_to_check, colnames(subdata))
  if (length(missing_genes) > 0) { for (gene in missing_genes) { subdata[[gene]] <- 0  } }
  subdata <-  subdata[,genes_to_check]
  subdata$cell = cell
  df1 = rbind(df1,subdata)
}
rio::export(df1,"E:\\InPut/fjjsub_marker_exp.txt")
data_long <- melt(df1, id.vars = "cell")    # 将数据转换为长格式
colnames(data_long) <- c("cell_type", "gene", "expression")

gene_keep = c(
  'CD19',"MS4A1","CD79A", #B
  'MZB1', #plasma
  'CD3D','CD3E','CD3G',"FOXP3", #T
  "KLRD1", #NK
  "CD1C", #DC
  'CPA3' ,'KIT' ,#mast
  "ITGAX" ,"CD68", #MONO/MACRO
  "FCGR3B",#中性粒
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
  geom_vline(aes(xintercept = 2.5), linetype = "dashed", color = "grey40") +  # 添加垂直线
  geom_vline(aes(xintercept = 8.5), linetype = "dashed", color = "grey40") +  # 添加垂直线
  geom_vline(aes(xintercept = 10.5), linetype = "dashed", color = "grey40") +  # 添加垂直线
  geom_vline(aes(xintercept = 12.5), linetype = "dashed", color = "grey40") +  # 添加垂直线
  geom_vline(aes(xintercept = 15.5), linetype = "dashed", color = "grey40") +  # 添加垂直线
  geom_hline(aes(yintercept = 3.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  geom_hline(aes(yintercept = 7.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  geom_hline(aes(yintercept = 10.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  geom_hline(aes(yintercept = 13.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  geom_hline(aes(yintercept = 14.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  geom_hline(aes(yintercept = 18.5), linetype = "dashed", color = "grey40") +  # 添加水平线
  theme_classic() +
  scale_color_gradient(low = "#d6ecf0",high="#f20c00")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
device_size <- dev.size()
width <- device_size[1]  # 获取宽度
height <- device_size[2]  # 获取高度
pdf(paste0("F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/","R1C.pdf"),
    width = width, height = height)
print(p)
dev.off()
# ggplot(mydata, aes(x = cell_type, y = expression, fill = cell_type)) +
#   geom_violin(trim = T, scale = "area") +
#   geom_boxplot(width = 0.1, alpha = 0.5,outlier.shape = NA) +  # 添加箱线图
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
#   labs(title = "Marker Gene Expression in Different Cell Subtypes",
#        x = "Cell Subtype",
#        y = "Expression Level") +
#   facet_wrap(~ gene, scales = "free_y",ncol = 2)+  # 按基因分面显示
#   scale_fill_manual(values = color_list)
  # scale_color_manual(values = color_list)
