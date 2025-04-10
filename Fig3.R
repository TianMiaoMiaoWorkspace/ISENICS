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
GTEX1=distinct(GTEX1,gene,.keep_all = T)
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
write.table(df,"threegroup.txt")
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


#################CSS##################
library(GSVA)
library(purrr)
library(dplyr)
library(tidyverse)
cd8gene<-c("IFNG","TNF","GZMB")
Mp<-c( "LGALS1", "PLP2", "S100A6", "HLA-DPA1", "HLA-DPB1",  "TSPO", "VIM", "HNRNPDL", "S100A4","HLA-F", "NPC2", "PFDN5", "ZBTB38", "HLA-DRB1","CNBP","AHNAK","ANKRD12","CIRBP")
Mn<-c( "H3F3A", "PLAC8", "GYPC", "CD27", "NOP10", "FKBP1A", "GIMAP4", "SLC25A5", "SERBP1", "SELL", "OSTC", "SOX4", "ABRACL", "UCP2", "RGS10", "PFN1", "HMGN1", "LIMD2", "CD7", "PGLS", "FYB1", "EVL", "TMIGD2", "IL27RA", "ADSL", "CD52", "MARCKSL1", "DENND2D", "ILF2", "LIMS2", "ARPC2", "SMC4", "MYC", "C9orf16", "ATP5T1C", "ZNF22", "TSPAN14", "CCR7")
siAge<-list(Mn=Mn,Mp=Mp)
SASP=read.table("F:\\Desktop\\major_revision\\all_SASP.csv",header=T,sep=",")
a=unlist(geneset)
SenMayo<-rio::import("F:\\Desktop\\major_revision/SenMayo.xlsx",sep = "\t",header = T)
SenMayo<-list(SenMayo$`Gene(human)`)
senCID<-get(load("F:\\Desktop\\major_revision\\senCID.RData"))
a = intersect(unlist(siAge),unlist(geneset))
b = intersect(SenMayo$`Gene(human)`,unlist(geneset))
c = intersect(senCID,unlist(geneset))
input<-"E:\\InPut/fjjtest/"
fjj<-list.files(input,pattern = "\\.csv$",full.names = T)
df1<-data.frame()
for (file in fjj) {
  cancer = strsplit(basename(file), "_")[[1]][1]
  cell = str_extract(file, "(?<=_).+(?=\\.)")
  data<- fread(file,header = T,sep = ",")
  data = data.frame(data)
  
  SASPexp <- data.frame(matrix(0, nrow = nrow(data), ncol = length(setdiff(SASP$x, colnames(data))), dimnames = list(NULL, setdiff(SASP$x, colnames(data)))))
  if (length(intersect(SASP$x, colnames(data))) > 0) {
    SASPexp <- cbind(data[, intersect(SASP$SASP_factor, colnames(data))], SASPexp)}
  SASPexp <- SASPexp[, match(SASP$x, colnames(SASPexp))]
  
  cd8geneexp <- data %>% dplyr::select(intersect(cd8gene,colnames(data))) 
  for (gene in setdiff(cd8gene, colnames(cd8geneexp))) {cd8geneexp[[gene]] <- 0 }
  cd8geneexp <- cd8geneexp[, cd8gene, drop=FALSE]
  
  sample <- data[[1]]  
  rownames(data)<-sample
  data <- t(data[, -1])
  CSS <- gsva(gsvaParam(data,geneset,kcdf="Gaussian"))
  SENMAYO <- gsva(gsvaParam(data,SenMayo,kcdf="Gaussian"))
  # SIAGE <-gsva(gsvaParam(data,siAge,kcdf="Gaussian"))
  tryCatch({
    SIAGE <- gsva(gsvaParam(data, siAge, kcdf = "Gaussian"))
  }, error = function(e) {
    message("GSVA failed: ", e$message)
    SIAGE <<- NULL
  })
  if (is.null(SIAGE)) {
    SIAGE <- matrix(rep(0, ncol(data)*2), nrow = 2)
    rownames(SIAGE) <- c("Mn","Mp")
  } else if (!"Mp" %in% rownames(SIAGE)) {
    SIAGE <- rbind(SIAGE, Mp = rep(0, ncol(SIAGE)))
  } 
  
  df=data.frame(sample=sample,
                cell=cell,
                cancer=cancer,
                CSS=CSS[2,]-CSS[1,],
                SENMAYO=SENMAYO[1,],
                SIAGE=SIAGE[2,]-SIAGE[1,],
                p21=data["CDKN1A",],
                p16=data["CDKN2A",],
                SA_β_Gal=data["GLB1",])
  if (is.null(SIAGE)) {
    SIAGE <- matrix(rep(0, ncol(data)*2), nrow = 2)
    rownames(SIAGE) <- c("Mn","Mp")
  }
  df = cbind(df,SASPexp,cd8geneexp)
  df1 = rbind(df1,df)
}
df1$cell <- gsub("Mono.Macro", "Mono_Macro", df1$cell)
table(df1$cell,df1$cancer)

###########################################################################
library('ComplexHeatmap')
library('circlize')
library(gridBase)
library(dplyr)

df=data.frame(cancer=c(rep("BLCA",414), rep("BRCA",1109),rep("CESC",306),rep("CHOL",36),rep("COAD",480),rep("DLBC",48),
                       rep("ESCA",162),rep("GBM",169),rep("HNSC",502),rep("KICH",65),rep("KIRC",539),rep("LAML",151),
                       rep("LIHC",374),rep("LUAD",535),rep("LUSC",502),rep("OV",379),rep("PAAD",178),rep("PRAD",499),rep("SKCM",471),
                       rep("STAD",375),rep("THCA",510),rep("UCEC",552),rep("UVM",80)))
allCSS <- df1
matrixCSS <- reshape(allCSS[,c("sample","cell","SENMAYO")], idvar = "sample", timevar = "cell", direction = "wide")
colnames(matrixCSS) <- sub("SENMAYO\\.", "", colnames(matrixCSS))
rownames(matrixCSS) = matrixCSS$sample
matrixCSS = matrixCSS[,-1]
matrixCSS = matrixCSS[,c("B","CD4Tconv","CD8T","CD8Tex",
                         "DC","Mast","Mono_Macro","Neutrophils",
                         "NK","Plasma","Tprolif","Treg",
                         "Endothelial","Epithelial","Fibroblasts","Myofibroblasts",
                         "Malignant")]
matrixCSS_transformed <-  sign(matrixCSS) * abs(matrixCSS)^(1/2) 

mycol= colorRamp2(c(-1, -0.2, 0, 0.2,1), c("#00796B", "#B2DFDB","white", "#FFCDD2" ,"#C62828"))

# mycol= colorRamp2(c(-1,0,1),c("green","white","red"))
ann_row = df
row.names(ann_row) = rownames(matrixCSS_transformed)
ann_row <- as.matrix(ann_row)
circos.par(gap.after=c(rep(1,22),30),
           cell.padding = c(0.05, 0, 0, 0))
circos.heatmap(matrixCSS_transformed,
               col=mycol ,
               split =ann_row,
               # show.sector.labels = T,
               # rownames.cex=0.9,
               cluster=F,
               rownames.side = c("none"),
               na.col = "grey98", 
               track.height = 0.75)
lg=Legend(title="SENMAYO",col_fun=mycol,direction = c("horizontal"))
grid.draw(lg)
circos.clear()
colnames(matrixCSS)#自己去AI添加


##################################SASP##################################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
library(pheatmap)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
immcell=c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono_Macro","NK",
       "Neutrophils","Plasma","Tprolif", "Treg")

immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono_Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
SASP=read.table("F:\\Desktop\\major_revision\\all_SASP.csv",header=T,sep=",")
IMMfjjcss <- allCSS %>% filter(cell %in% immune_cell)

CSS_SASPcor <- IMMfjjcss %>%
  group_by(cell, cancer) %>%
  filter(if_any(11:100, ~ . > quantile(., 0.25) + 1.5 * IQR(.))) %>%
  summarise(across(11:100,  ~ if (sum(!is.na(.)) > 1) cor(CSS, ., use = "pairwise.complete.obs")
                              else NA_real_,.names = "{.col}"), .groups = "drop")


color.pals = c(
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
color_list <- c("#ffe0e9", "#4c8dae", "#8ca1c4","#2e4e7e","#f47983", "#ff4777","#d07b76","#ff9b60",
                "#e8ac9d","#ffb1c9","#48c0a3", "#a4e2c6")
##########CSS_SASPcor
cor_results <- CSS_SASPcor[, -c(1, 2)] %>% 
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
cor_results <-abs(cor_results)
rownames(cor_results)<-paste(CSS_SASPcor$cancer,CSS_SASPcor$cell)

row_annotation <- data.frame(
  cell = factor(CSS_SASPcor$cell),
  cancer = factor(CSS_SASPcor$cancer)
)
rownames(row_annotation) <- rownames(cor_results)


annotation_colors <- list(
  cell = cellcolor,
  cancer = setNames(color.pals, unique(row_annotation$cancer))
)
pheatmap(cor_results,
         main = "CSS_SASPcor",  
         border_color=NA,
         color = colorRampPalette(c("white", "red"))(50),  
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         cluster_rows = F,  
         cluster_cols = F, 
         show_rownames = F,
         show_colnames = TRUE, 
         fontcolor = "black",
         annotation_row = row_annotation,
         annotation_colors = annotation_colors)  
positive_count = rowSums(cor_results> 0, na.rm = TRUE)
non_na_count = rowSums(!is.na(cor_results))
positive_ratio = data.frame(positive_count / non_na_count )
ggplot(positive_ratio, aes(x = positive_count.non_na_count)) +
  geom_density(fill="grey80",alpha = 0.5) + 
  theme_classic() 


SenMayo<-rio::import("F:\\Desktop\\major_revision/SenMayo.xlsx",sep = "\t",header = T)
senCID<-rio::import("F:\\Desktop\\major_revision\\senCIDseneset.txt",header = F)
con<-intersect(SenMayo$`Gene(human)`,senCID$V1)
geneset<-get(load("F:/Desktop/HMU/2024文章/毕业设计/inputdata/geneset.RData"))
con<-intersect(unlist(geneset),senCID$V1)


#######################cd8gene##################
library(dplyr)
library(ggplot2)
library(data.table)
library("ggprism")
library(ggpubr)
library(ggbreak)
library(ggsignif)
library(ggbeeswarm)
library(tidyr)
library(dplyr)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
cd8gene<-c("IFNG","TNF","GZMB")
allCSS <- read.table("E:/InPut/反卷积CSS.txt",header=T,sep="\t")
CD8T_gene <- allCSS %>%filter(cell=="CD8T")%>%select(sample,cell,cancer,CSS)%>%
  group_by(cancer) %>% mutate(
    q1 = quantile(CSS, 0.25, na.rm = TRUE),
    q3 = quantile(CSS, 0.75, na.rm = TRUE),
    dif = max(CSS) - min(CSS),
    quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
  )

setwd("E:\\InPut\\TCGAbulk")
df=data.frame()
for (cancer in CANCER) {
  data <- fread(paste0(cancer,".txt"), header = TRUE, sep = "\t")
  data = data.frame(data)[,c("V1","IFNG","TNF","GZMB")]
  sample = CD8T_gene[CD8T_gene$cancer==cancer,"sample"]
  rownames(data) = data$V1
  if(length(sample$sample)>0){
    data = data[sample$sample,]
    data$cancer = cancer
    CSS=CD8T_gene[CD8T_gene$cancer==cancer,"CSS"]
    data$CSS=CSS$CSS
    quantile_CSS_Group=CD8T_gene[CD8T_gene$cancer==cancer,"quantile_CSS_Group"]
    data$quantile_CSS_Group=quantile_CSS_Group$quantile_CSS_Group
    
    df=rbind(df,data)
    print(paste0("OK ",cancer,"!"))
  }
}
rio::export(df,"F:\\Desktop\\major_revision\\CD8T.txt")
cd8gene<-c("IFNG","TNF","GZMB")
#,"CESC","UVM","CHOL"
data=df[df$cancer%in%c("BRCA"),]
data=data[data$IFNG>0,]
data=data[data$TNF>0,]
data=data[data$GZMB>0,]
data=data[data$quantile_CSS_Group!="Medium",]
long_data <- data %>%
  pivot_longer(
    cols = c(IFNG, TNF, GZMB), 
    names_to = "gene",
    values_to = "value"  
  ) %>%
  select(cancer, quantile_CSS_Group, gene, value,cancer)  

ggplot(long_data, aes(x = quantile_CSS_Group, y = value, fill = quantile_CSS_Group)) +
  geom_violin(trim = TRUE,position = "dodge") +   
  theme_minimal() +  
  stat_compare_means(aes(group = quantile_CSS_Group), 
                     method = "wilcox.test",label = "p.signif",
                     label.y = c(2, 100))+
  facet_wrap(~ cancer+gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("#d56e5e","#5390b5")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SenCID")

################################################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)

color_list <- c("#ffe0e9", "#4c8dae", "#8ca1c4","#2e4e7e","#f47983", "#ff4777","#d07b76",
                "#e8ac9d","#ff9b60","#ffb1c9","#48c0a3", "#a4e2c6",  
                "#FFfcc6", "#ffeda8","#fff897","#ecfc00","#f2ecde")
cellorder=c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono_Macro","Neutrophils","NK",
                    "Plasma","Tprolif", "Treg","Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts","Malignant")
allCSS$cell <- factor(allCSS$cell, levels =rev(cellorder) )  

ggplot(allCSS, aes(x = CSS, y = cell,fill = cell)) +
  geom_density_ridges(rel_min_height = 0.005)+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  theme_ridges()+ 
  scale_fill_manual(values = rev(color_list))

#####################################CSSbulk4############################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(limma)
library(purrr)
library(dplyr)
library(tidyverse)
library(data.table)
input<-"E:\\InPut\\verifyTPM\\免疫fjj"
fjj<-list.files(input,pattern = "\\.csv$",full.names = T)
df1<-data.frame()
for (file in fjj) {
  cancer = strsplit(basename(file), "_")[[1]][1]
  cell = str_extract(file, "(?<=_).+(?=\\.)")
  data<- fread(file,header = T,sep = ",")
  data = data.frame(data)
  sample <- data[[1]] 
  rownames(data)<-sample
  data <- t(data[, -1])
  tryCatch({
    SIAGE <- gsva(gsvaParam(data, siAge, kcdf = "Gaussian"))
  }, error = function(e) {
    message("GSVA failed: ", e$message)
    SIAGE <<- NULL
  })
  if (is.null(SIAGE)) {
    SIAGE <- matrix(rep(0, ncol(data)*2), nrow = 2)
    rownames(SIAGE) <- c("Mn","Mp")
  } else if (!"Mp" %in% rownames(SIAGE)) {
    SIAGE <- rbind(SIAGE, Mp = rep(0, ncol(SIAGE)))
  } 
  CSS <- gsva(gsvaParam(data,geneset,kcdf="Gaussian"))
  SENMAYO <- gsva(gsvaParam(data,SenMayo,kcdf="Gaussian"))
  SenCID <- gsva(gsvaParam(data,list(senCID),kcdf="Gaussian"))
  df=data.frame(sample=sample,
                cell=cell,
                cancer=cancer,
                SENMAYO=SENMAYO[1,],
                SenCID = SenCID[1,],
                SIAGE = SIAGE[2,]-SIAGE[1,],
                CSS=CSS[2,]-CSS[1,]
                )
  df1 = rbind(df1,df)
}
df1$cell <- gsub("Mono.Macro", "Mono_Macro", df1$cell)
# table(df1$cell,df1$cancer)
con=c("data.G400_1_1","data.G400_1_2","data.G400_1_3",
      "TPM_Con_1","TPM_Con_2","TPM_Con_3",
      "SW620.con1","SW620.con2",
      "GSM1916644","GSM1916645","GSM1916649","GSM1916650")

sen=c("data.G400_2_1","data.G400_2_2","data.G400_2_3",
      "TPM_Doxo_1","TPM_Doxo_2","TPM_Doxo_3","TPM_Abe_1","TPM_Abe_2","TPM_Abe_3",
      "SW620.sen1","SW620.sen2",
      "GSM1916646","GSM1916647","GSM1916648","GSM1916651","GSM1916652","GSM1916653")
df1$group <- case_when(
  df1$sample %in% con ~ "con",
  df1$sample %in% sen ~ "sen"
)
# rio::export(df1,"E:\\InPut\\verifyTPM\\verify.txt")
###################
library("ggplot2")
library("ggprism")
library(ggpubr)
library(ggbreak)
library(ggsignif)
library("ggbeeswarm")
df1<-rio::import("E:\\InPut\\verifyTPM\\verify.txt")

p1 = ggplot(df1, aes(x = group, y = CSS, fill = group)) +
  geom_violin(trim = TRUE,position = "dodge",scale = "area") +
  theme_classic() +  
  stat_compare_means(aes(group = group), method = "wilcox.test",label = "p.signif")+
  facet_wrap(~ cancer, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(title = "CSS")

p2 = ggplot(df1, aes(x = group, y = SenCID, fill = group)) +
  geom_violin(trim = TRUE,position = "dodge",scale = "area") +  
  theme_classic() +  
  stat_compare_means(aes(group = group), method = "wilcox.test",label = "p.signif")+
  facet_wrap(~ cancer, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SenCID")
p3 = ggplot(df1, aes(x = group, y = SENMAYO, fill = group)) +
  geom_violin(trim = TRUE,position = "dodge",scale = "area") + 
  theme_classic() + 
  stat_compare_means(aes(group = group), method = "wilcox.test",label = "p.signif")+
  facet_wrap(~ cancer, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SENMAYO")
p4 = ggplot(df1, aes(x = group, y = SIAGE, fill = group)) +
  geom_violin(trim = TRUE,position = "dodge",scale = "area") +  
  theme_classic() +  
  stat_compare_means(aes(group = group), method = "wilcox.test",label = "p.signif")+
  facet_wrap(~ cancer, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SIAGE")
p1/p2/p3/p4
#####################

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LIHC","LUAD","LUSC",####去掉LAML在IOBR找不到
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
fjjcss <- read.table("E:/InPut/反卷积CSS.txt",sep = "\t",header=T)
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

cancer_dfs <- list()
unique_cancers <- unique(fjjcss$cancer)
for (cancer_type in unique_cancers) {
  df <- fjjcss %>%
    filter(cancer == cancer_type) %>%
    dplyr::select(sample, cell, quantile_CSS_Group) %>%
    pivot_wider(
      names_from = cell,
      values_from = quantile_CSS_Group,
      values_fill = list(quantile_CSS_Group = NA)
    )
 
  cancer_dfs[[as.character(cancer_type)]] <- df
}

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
  filter(Frequency<=5)%>%dplyr::select(CellType)#
res2 = res1[!(res1$Column1 %in% cell_del$CellType | res1$Column2 %in% cell_del$CellType), ]
rio::export(res2,"outputdata/cell_both_high.tsv")

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
  cell =   c("B","CD4Tconv","CD8T", "CD8Tex", "DC","Mast", "Mono_Macro","NK","Neutrophils","Plasma","Tprolif", "Treg",
             "Malignant",
             "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts"
  ),
  cell_col =c("#ffe0e9", "#4c8dae", "#8ca1c4","#2e4e7e","#f47983", "#ff4777","#d07b76","#ff9b60",
                            "#e8ac9d","#ffb1c9","#48c0a3", "#a4e2c6", 
                            "#f2ecde",
                            "#FFfcc6", "#ffeda8","#fff897","#ecfc00")
  
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
