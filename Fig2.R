##########################Data were pooled and each cancer was run separately#############################################
library(data.table)
library(dplyr)
library(tidyverse)
BLCA<-get(load("D:/TCGA/TCGA-BLCA_RNASeq_Counts.RData"))
BLCA_tumor <- BLCA[,substr(colnames(BLCA),14,15) < 10] 
BLCA_adj <- BLCA[,substr(colnames(BLCA),14,15) >= 10] 
GTEX=fread("gtex_RSEM_gene_fpkm")
GTEX=as.data.frame(GTEX) #60698
GTEX_phe=fread("GTEX_phenotype")
GTEX_phe=as.data.frame(GTEX_phe)

GTEX_pro=fread("probeMap_gencode.v23.annotation.gene.probemap")
GTEX_pro=as.data.frame(GTEX_pro)
GTEX_pro=GTEX_pro[,c(1,2)]
GTEX=merge(GTEX_pro,GTEX,by.y = "sample",by.x = "id")
rownames(GTEX_phe)=GTEX_phe$Sample
GTEX_phe=GTEX_phe[,-1]
colnames(GTEX_phe)=c("body_site_detail (SMTSD)","primary_site","gender","patient","cohort")

GTEX_phe1=filter(GTEX_phe,primary_site=="Cervix Uteri")
GTEX1=GTEX[,c("id","gene",intersect(rownames(GTEX_phe1),colnames(GTEX)))]#
GTEX1=distinct(GTEX1,gene,.keep_all = T)# 
rownames(GTEX1)=GTEX1$gene
GTEX1=GTEX1[,-c(1,2)]
GTEX1=log2(2^GTEX1-0.001+1)

same_gene=intersect(row.names(GTEX1),row.names(BLCA))
data=cbind(BLCA_tumor[same_gene,],BLCA_adj[same_gene,],GTEX1[same_gene,])


#library(limma)
data2=normalizeBetweenArrays(data)
data2=as.data.frame(data2)
data2=data2[apply(data2, 1, var, na.rm=TRUE) != 0,]
samp_info=c(rep(c("tumor","adj","GTEX"),c(ncol(BLCA_tumor),ncol(BLCA_adj),length(GTEX1))))
samp_info=as.data.frame(samp_info)
colnames(samp_info)="batch"
samp_info$batch <- as.factor(samp_info$batch)
#library(sva)
design_matrix <- model.matrix(~ batch, data = samp_info)

adjusted_matrix <- ComBat(dat = data2, batch = samp_info$batch,par.prior = TRUE, prior.plots = FALSE)

#####################################
data=cbind(BLCA_tumor,BLCA_adj)
adjusted_matrix<-normalizeBetweenArrays(data)
library(GSVA)
BLCA_CS<- gsva(expr=data,
               gset.idx.list=geneset,
               method="ssgsea",
               kcdf="Gaussian",
               verbose=T, 
               mx.diff=T)
BLCA_CS<-as.data.frame( BLCA_CS[2,]-BLCA_CS[1,])
colnames(BLCA_CS)<-"CS"
BLCA_CS$info<-c(rep(c("tumor","adj"),c(ncol(BLCA_tumor),ncol(BLCA_adj))))
BLCA_CS$cancer<-rep("BLCA",nrow(BLCA_CS))

#####################################
geneset<-get(load("geneset.RData"))
#library("GSVA")
BLCA_CS<- gsva(expr=adjusted_matrix,
               gset.idx.list=geneset,
               method="ssgsea",
               kcdf="Gaussian",
               verbose=T, 
               mx.diff=T)
BLCA_CS<-as.data.frame( BLCA_CS[2,]-BLCA_CS[1,])
colnames(BLCA_CS)<-"CS"
BLCA_CS$info<-samp_info$batch
BLCA_CS$cancer<-rep("BLCA",nrow(samp_info))


#df<-data.frame()
df<-rbind(df,BLCA_CS)
##################################gtex adj tumor#############################
library(ggplot2)
library(ggprism)
library(ggpubr) 
library(ggbreak)
library(dplyr)
library(ggsignif) 

df$info <- factor(df$info, levels = c("tumor", "adj", "GTEX"))
mean_data <- df %>% group_by(info, cancer) %>%summarise(mean_len = mean(CSS), .groups = "drop")  
p=ggplot(df,aes(x = info, y = CSS,  fill = info,color = info)) + 
  theme(legend.position = "top")+
  geom_violin(width = 0.8)+
  geom_point(data = mean_data, aes(x = info, y = mean_len), shape = 21, size = 1.8, fill = "black") +
  geom_line(data = mean_data, aes(x = info, y = mean_len, group = cancer), size = 0.7,color="black") +
  labs(x = "",  y = "CSS")+
  facet_wrap(~cancer,nrow=1,strip.position = "bottom")+
  theme_classic()+
  scale_fill_manual(values = c('#dc3023','#ffc080','#5861ac'))+
  scale_color_manual(values = c('#dc3023','#ffc080','#5861ac'))+
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.ticks=element_blank())+
  stat_compare_means(aes(label=..p.signif..),
                     method = "aov",
                     label.y = 0.95)


#######################CSS  circus###################################################
library('ComplexHeatmap')
library('circlize')
library(gridBase)
library(dplyr)

df=data.frame(cancer=c(rep("BLCA",414), rep("BRCA",1109),rep("CESC",306),rep("CHOL",36),rep("COAD",480),rep("DLBC",48),
                       rep("ESCA",162),rep("GBM",169),rep("HNSC",502),rep("KICH",65),rep("KIRC",539),rep("LAML",151),
                       rep("LIHC",374),rep("LUAD",535),rep("LUSC",502),rep("OV",379),rep("PAAD",178),rep("PRAD",499),rep("SKCM",471),
                       rep("STAD",375),rep("THCA",510),rep("UCEC",552),rep("UVM",80)))
allCSS <- read.table("fjjCSS.txt",header=T,sep="\t",row.names = 1)
matrixCSS <- reshape(allCSS[,1:3], idvar = "sample", timevar = "cell", direction = "wide")
matrixCSS <- matrixCSS[,c(1:13,15:17)]
colnames(matrixCSS) <- sub("CSS\\.", "", colnames(matrixCSS))
rownames(matrixCSS) = matrixCSS$sample
matrixCSS = matrixCSS[,-1]

matrixCSS = matrixCSS[,c("B","CD4Tconv","CD8T","CD8Tex","NK","Mast","MonoMacro",
                         "Plasma","Tprolif","Treg","Endothelial","Epithelial","Fibroblasts","Myofibroblasts","Malignant")]
# matrixCSS[is.na(matrixCSS)] <- 0  
mycol= colorRamp2(c(-1, 0, 1),c("#238b49","white","#e74a32"))
ann_row = df
row.names(ann_row) = rownames(matrixCSS)
ann_row <- as.matrix(ann_row)
circos.par(gap.after=c(rep(1,22),30),
           cell.padding = c(0.05, 0, 0, 0)) 
circos.heatmap(matrixCSS,
               col=mycol, 
               split =ann_row,
               # show.sector.labels = T,
               # rownames.cex=0.9,
               cluster=F,
               rownames.side = c("none"),
               na.col = "grey98", 
               track.height = 0.75)
lg=Legend(title="CSS",col_fun=mycol,direction = c("horizontal"))
grid.draw(lg)
circos.clear()
colnames(matrixCSS)
########################density map########################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
allCSS <- read.table("fjjCSS.txt",header=T,sep="\t",row.names = 1)
# unique(allCSS$cell[allCSS$cancer=="ESCA"])
# ESCA=CSS[allCSS$cancer=="ESCA",c(2,3)]%>%filter(cell%in%c("B","CD4Tconv","CD8T","MonoMacro","Plasma","Tprolif"))
# ggdensity(ESCA, x = "CSS", 
#           rug = TRUE, xlab = "ESCA_CSS",
#           color = "cell",
#           fill = "cell",bw=0.1) +
#   guides(color="none")+ 
#   labs(fill=" ")+
#   scale_color_manual(values = c("#F06292","#8BC34A","#4DD0E1","#FFEE58","#9575CD","#8D6E63"))+
#   scale_fill_manual(values = c("#F06292","#8BC34A","#4DD0E1","#FFEE58","#9575CD","#8D6E63")) +
#   theme(legend.position = "right")

p1=ggdensity(allCSS[allCSS$cell=="B",], x = "CSS", 
          rug = TRUE, xlab = "B cells",ylab = "",color = "#FCE4EC",fill = "#FCE4EC",bw=0.1) +
          guides(color="none")+ labs(fill=" ")
p2=ggdensity(allCSS[allCSS$cell=="CD4Tconv",], x = "CSS", 
          rug = TRUE, xlab = "CD4Tconv cells",ylab = "",color = "#EC407A",fill = "#EC407A",bw=0.1) +
          guides(color="none")+ labs(fill=" ")+theme(legend.position = "right")
p3=ggdensity(allCSS[allCSS$cell=="CD8T",], x = "CSS", 
          rug = TRUE, xlab = "CD8T cells",ylab = "",color = "#C2185B",fill = "#C2185B",bw=0.1) +
          guides(color="none")+ labs(fill=" ")
p4=ggdensity(allCSS[allCSS$cell=="CD8Tex",], x = "CSS", 
          rug = TRUE, xlab = "CD8Tex cells",ylab = "",color = "#880E4F",fill = "#880E4F",bw=0.1) +
          guides(color="none")+ labs(fill=" ")
p5=ggdensity(allCSS[allCSS$cell=="MonoMacro",], x = "CSS", 
          rug = TRUE, xlab = "MonoMacro cells",ylab = "",color = "#FFB3C8",fill = "#FFB3C8",bw=0.1) +
          guides(color="none")+ labs(fill=" ")
p6=ggdensity(allCSS[allCSS$cell=="Plasma",], x = "CSS", 
          rug = TRUE, xlab = "Plasma cells",ylab = "",color = "#B71C1C",fill = "#B71C1C",bw=0.1) +
          guides(color="none")+ labs(fill=" ")
p7=ggdensity(allCSS[allCSS$cell=="Treg",], x = "CSS", 
          rug = TRUE, xlab = "Treg cells",ylab = "",color = "#7B1FA2",fill = "#7B1FA2",bw=0.1) +
          guides(color="none")+ labs(fill=" ")


p1 / p2 / p3 / p4 / p5 / p6 / p7

#####################################bulk  verify############################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
##################### OV
OV_verify_GSE128700<-read.csv("OV_verify_GSE128700.csv")
rownames(OV_verify_GSE128700)<-OV_verify_GSE128700$name
Ensembl_ID <- OV_verify_GSE128700$name
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
gene_symbol<-gene_symbol[which(!duplicated(gene_symbol$SYMBOL)),]
data<-OV_verify_GSE128700[gene_symbol$ENSEMBL,]
data$symbol<-gene_symbol$SYMBOL
df<-data.frame(data$G400_1_1,data$G400_1_2,data$G400_1_3,
               data$G400_2_1,data$G400_2_2,data$G400_2_3)
rownames(df)<-data$symbol
data1=data.frame(sample=c("G400_1_1","G400_1_2","G400_1_3",
                          "G400_2_1","G400_2_2","G400_2_3"),
                 group=c(rep("con",3),rep("sen",3)))

CSS <- gsva(gsvaParam(as.matrix(df),geneset,kcdf="Gaussian"))
data1$CSS=CSS=CSS[2,]-CSS[1,]

##################### COAD
COAD<-read.csv("COAD_verify_GSE226998.csv")[,1:5]
colnames(COAD)<-COAD[1,]
COAD<-COAD[-1,]
COAD<-COAD[which(!duplicated(COAD$`Gene Symbol`)),]
df<-COAD[,2:5]
df=as.data.frame(lapply(df,as.numeric))
rownames(df)<-COAD$`Gene Symbol`
data2=data.frame(sample=colnames(df),
                 group=c(rep("con",2),rep("sen",2)))
CSS <- gsva(gsvaParam(as.matrix(df),geneset,kcdf="Gaussian"))
data2$CSS=CSS[2,]-CSS[1,]
##################### LAML
LAML<-read.table("LAML_verify_GSE186200.txt",header = T)
rownames(LAML)<-LAML$GeneID
gene_symbol <- bitr(LAML$GeneID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
gene_symbol<-gene_symbol[which(!duplicated(gene_symbol$SYMBOL)),]
data<-LAML[gene_symbol$ENSEMBL,]


df<-data[,2:10]
rownames(df)<-gene_symbol$SYMBOL
data3=data.frame(sample=colnames(df),
                 group=c(rep("con",3),rep("sen",3),rep("con",3)))
CSS <- gsva(gsvaParam(as.matrix(df),geneset,kcdf="Gaussian"))
data3$CSS=CSS[1,]
##################### BRCA
GSE222984<-read.table("BRCA_verify_GSE222984.txt",header = T,sep = "\t")
df<-GSE222984[,3:11]
rownames(df)<-GSE222984$Gene.symbol
data4=data.frame(sample=colnames(df),
                 group=c(rep("con",3),rep("sen",6)))
CSS <- gsva(gsvaParam(as.matrix(df),geneset,kcdf="Gaussian"))
data4$CSS=CSS[2,]-CSS[1,]
##################### GBM
GSE74304<-read.table("GBM_verify_GSE74304.txt",header = T,sep = "\t")
rownames(GSE74304)<-GSE74304$ID_REF
GPL<-rio::import("GPL570.txt")
GPL$`Gene Symbol`<-gsub("[/].*","",GPL$`Gene Symbol`)
GPL<-GPL[which(!duplicated(GPL$`Gene Symbol`)),]
GSE74304<-GSE74304[GPL$ID,]
df<-GSE74304[,2:11]
rownames(df)<-GPL$`Gene Symbol`

data5=data.frame(sample=colnames(df),
                 group=c(rep("con",2),rep("sen",3),rep("con",2),rep("sen",3)))
CSS <- gsva(gsvaParam(as.matrix(df),geneset,kcdf="Gaussian"))
data5$CSS=CSS[2,]-CSS[1,]


library("ggplot2")
library("ggprism")
library(ggpubr) 
library(ggbreak)
library(ggsignif)
library("ggbeeswarm")
p1 <- ggplot(data1,aes(x=group,y=CSS,color=CSS))+
  geom_quasirandom(width =0.2,size=5)+
  labs(x="GSE128700",y="CSS")+
  scale_color_gradient(low="#5390b5",high ="#d56e5e")+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),
                     method = "t.test")
p2 <- ggplot(data2,aes(x=group,y=CSS,color=CSS))+
  geom_quasirandom(width =0.2,size=5)+
  labs(x="GSE226998",y="")+
  scale_color_gradient(low="#5390b5",high ="#d56e5e")+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),
                     method = "t.test")


p3 <- ggplot(data3,aes(x=group,y=CSS,color=CSS))+
  geom_quasirandom(width =0.2,size=5)+
  labs(x="GSE186200",y="")+
  scale_color_gradient(low="#5390b5",high ="#d56e5e")+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),
                     method = "t.test")


p4 <- ggplot(data4,aes(x=group,y=CSS,color=CSS))+
  geom_quasirandom(width =0.2,size=5)+
  labs(x="GSE222984",y="")+
  scale_color_gradient(low="#5390b5",high ="#d56e5e")+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),
                     method = "t.test")

p5 <- ggplot(data5,aes(x=group,y=CSS,color=CSS))+
  geom_quasirandom(width =0.2,size=5)+
  labs(x="GSE74304",y="")+
  scale_color_gradient(low="#5390b5",high ="#d56e5e")+
  theme_classic()+
  stat_compare_means(aes(label=..p.signif..),
                     method = "t.test")



ggarrange(p1,p2,p3,p4,p5,ncol=5,common.legend = T, legend="right")



##########co -sen###########
CSSgene = get(load("geneset.RData"))
load( "F:/Desktop/毕业设计/inputdata/IMMfjjcss.RData")
IMMfjjcss$group = IMMfjjcss$quantile_CSS_Group
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LIHC","LUAD","LUSC",
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
fjjcss <- read.table("fjjCSS.txt",sep = "\t",header=T,row.names=1)
fjjcss <- fjjcss %>%
  group_by(cell, cancer) %>% mutate(
    q1 = quantile(CSS, 0.25, na.rm = TRUE),
    q3 = quantile(CSS, 0.75, na.rm = TRUE),
    dif = max(CSS) - min(CSS),
    file = paste(cancer,cell),
    quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
  )

cancer_dfs <- list()
unique_cancers <- unique(fjjcss$cancer)
for (cancer_type in unique_cancers) {
  df <- fjjcss %>%ungroup() %>% 
    filter(cancer == cancer_type) %>%
    dplyr::select(sample, cell, CSS) %>%
    pivot_wider( names_from = cell,values_from = CSS,
                 values_fill = list(group = NA) )
  cancer_dfs[[as.character(cancer_type)]] <- df
}
# cancer_dfs <- cancer_dfs[-8]  
result <- data.frame(Column1 = character(),
                     Column2 = character(),
                     p = numeric(),
                     cor = numeric(),
                     stringsAsFactors = FALSE,
                     cancer = character())
for (cancer in names(cancer_dfs)) {
  data <- cancer_dfs[[cancer]]
  columns_to_check <-colnames(data)[-c(1,2)]

  for (i in 1:(length(columns_to_check) - 1)) {
    for (j in (i + 1):length(columns_to_check)) {
      col1 <- columns_to_check[i]
      col2 <- columns_to_check[j]
      cor = cor.test(data[[col1]], data[[col2]])$estimate
      p_value <- cor.test(data[[col1]], data[[col2]])$p.value
      result <- rbind(result, data.frame(Column1 = col1,
                                         Column2 = col2,
                                         p_cor = p_value,
                                         cor = cor,
                                         cancer = cancer))
    }
  }
}
cell_both_high_cor = result
# cell_counts <- result %>%filter(abs(cor)>=0.75&p_value<=0.05)%>%
#   pivot_longer(cols = c(Column1, Column2), names_to = "Column", values_to = "CellType") %>%
#   count(CellType, name = "Frequency")
# rio::export(res,"outputdata/cell_both_high_cor.tsv")

cancer_dfs <- list()
unique_cancers <- unique(fjjcss$cancer)
for (cancer_type in unique_cancers) {
  df <- fjjcss %>%
    filter(cancer == cancer_type) %>%
    select(sample, cell, quantile_CSS_Group) %>%
    pivot_wider(
      names_from = cell,
      values_from = quantile_CSS_Group,
      values_fill = list(quantile_CSS_Group = NA)
    )

  cancer_dfs[[as.character(cancer_type)]] <- df
}
# cancer_dfs <- cancer_dfs[-8]  
result <- data.frame(Column1 = character(),
                     Column2 = character(),
                     cramer_v = numeric(),
                     p = numeric(),
                     stringsAsFactors = FALSE,
                     cancer = character())
for (cancer in names(cancer_dfs)) {
  data <- cancer_dfs[[cancer]]
  columns_to_check <-colnames(data)[-c(1,2)]
  for (i in 1:(length(columns_to_check) - 1)) {
    for (j in (i + 1):length(columns_to_check)) {
      col1 <- columns_to_check[i]
      col2 <- columns_to_check[j]
      contingency_table <- table(data[[col1]], data[[col2]])
      stats <- vcd::assocstats(contingency_table)
      cramev <- stats$cramer  # Cramér's V
      p_value <- stats$chisq_tests[2, 3] 
      result <- rbind(result, data.frame(Column1 = col1,
                                         Column2 = col2,
                                         p_cramev = p_value,
                                         cramer_v = cramev,
                                         cancer = cancer))
    }
  }
}
cell_both_high_cramerv = result
res = merge(cell_both_high_cor,cell_both_high_cramerv)
res1 = res[res$p_cor<=0.05 & res$p_cramev<=0.05 ,]
cell_del <- res1 %>%
  pivot_longer(cols = c(Column1, Column2), names_to = "Column", values_to = "CellType") %>%
  count(CellType, name = "Frequency")%>%
  filter(Frequency<=5)%>%select(CellType)
res2 = res1[!(res1$Column1 %in% cell_del$CellType | res1$Column2 %in% cell_del$CellType), ]

##########network###

library(igraph)
res3<-res2
edgesum <- res2 %>%
  gather(key = "point", value = "node", Column1, Column2) %>%
  # filter(abs(cor) > 0.5) %>% 
  group_by(node) %>%
  summarise(sumcor = sum(abs(cor))) %>%
  as.data.frame()
res3$cor <-abs(res3$cor)
cancer_col <- data.frame(
  cancer = unique(res2$cancer)[order(unique(res2$cancer))] ,
  cancer_col = c(
    "#CCCCFF",#BLCA
    "#6699CC",#BRCA
    "#663366",#CESC
    "#FF6666",#CHOL
    "#FFFF00",#COAD
    "#003399",#DLBC
    "#FEBB81",#ESCA
    "#CCFF00",#GBM
    "#CC6699",#HNSC
    "#FFFF99",#KICH
    "#DB5C68",#KIRC
    "#CC3399",#LAML
    "#FFCCCC",#LIHC
    "#FF99CC",#LUAD
    "#D44842",#LUSC
    "#CCCC99",#OV
    "#666699",#PAAD
    "#99CC66",#PRAD
    "#66CC99",#SKCM
    "#339999",#STAD
    "#CC99CC",#THCA
    "#FAC228",#UCEC
    "#CCDDDD"#UVM
  )
)
cell_col = data.frame(
  cell =  c("AC_Like_Malignant", "Malignant", 
            "B","Mast", "MonoMacro","Neutrophils", "NK", "pDC","Plasma",
            "CD4Tconv", "CD8T", "CD8Tex", "DC",  "Promonocyte","Tprolif", "Treg",
            "Pericytes", "Astrocyte","Neuron","Oligodendrocyte", "SMC","Vascular",
            "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts","Myocyte","Acinar", "Ductal",
            "Hepatic progenitor","HSC", "Progenitor","OPC", "Endometrial stromal cells", 
            "EryPro", "GMP"
  ),
  cell_col = c("#b2cbe6", "#9FA8DA",
               "#FCE4EC", "#FFB3C8","#e6a0af","#EF9A9A", "#E57373","#EF5350","#B71C1C",
               "#EC407A", "#C2185B","#880E4F","#CE93D8",  "#cb86b5","#9C27B0", "#7B1FA2", 
               "#D4F69F", "#A1FFA1","#91FF91","#76d4a4", "#81C784", "#2E7D32",
               "#FFFDE7","#FFECB3","#DCE775","#C0CA33", "#FFEE58","#FDD835","#F9A825", 
               "#B2EBF2", "#4DD0E1","#00BCD4","#0097A7", "#006064",
               "#D7CCC8","#A1887F"
  )
)
res3 = merge(res3,cancer_col)
res3 <- res3 %>%mutate(cancer_col = ifelse(cor <= 0.5, "grey35", cancer_col))
res3 <- res3[order(res3$cramer_v),]
network <- graph_from_data_frame(d = cbind(res3$Column1,res3$Column2,res3$cor),
                                 directed = FALSE)


node_size <- numeric(length(V(network)$name))
names(node_size) <-V(network)$name

for (i in seq_along(V(network)$name)) {
  node_name <- V(network)$name[i]
  match_index <- match(node_name, edgesum$node)
  if (!is.na(match_index)) {
    node_size[i] <- edgesum$sumcor[match_index]
  } else {
    node_size[i] <- 1 
  }
}

cellcol <- NULL
for (i in seq_along(V(network)$name)) {
  node_name <- V(network)$name[i]
  match_index <- match(node_name, cell_col$cell)
  if (!is.na(match_index)) {
    cellcol[i] <- cell_col$cell_col[match_index]
  } else {
    cellcol[i] <- "white"  
  }
}

E(network)$width <-(exp(res3$cor) - 1)^3
plot(network, 
     vertex.size = node_size, 
     layout = layout.circle,  
     vertex.color = cellcol,  
     vertex.label.cex = 1.5,  
     vertex.label.color = "black", 
     vertex.frame.color = "transparent",  
     edge.lty=1,  
     edge.curved = 0.3, 
     edge.color = res3$cancer_col)

library(ggplot2)
cancer_col$y.infor<-1:23
tl <- ggplot()+ geom_col(cancer_col,mapping=aes(cancer,y.infor,fill = cancer))+ 
  scale_fill_manual(values = cancer_col$cancer_col)+
  guides(fill = guide_legend(ncol = 1))
