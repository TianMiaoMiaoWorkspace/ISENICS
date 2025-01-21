library(R.utils)
library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(harmony)
library(cowplot)
library(clustree)
library(patchwork)
library(GSVA)
library(pROC)

folder_path <- "GSE203115"
subfolders <- list.dirs(folder_path, full.names = TRUE, recursive = TRUE)
for (subfolder in subfolders) {
  gz_files <- list.files(subfolder, pattern = "\\.gz$", full.names = TRUE)
  if (length(gz_files) > 0) {
    for (gz_file in gz_files) {
      gunzip(gz_file, overwrite = TRUE)  
    }
  }
}
setwd("GSE203115")
dir="GSE203115/raw/"
samples=list.files( dir )
sceList = lapply(samples,function(pro){ 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  sce =CreateSeuratObject(counts =  tmp ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 200 )
  return(sce)
}) 
View(sceList)
#####merge
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  ) 
names(sce.all@assays$RNA@layers)
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
########
phe = sce.all@meta.data
table(phe$orig.ident)
# View(phe)
phe$group = case_when(phe$orig.ident=="GSE203115ESCC1"~"responder",
                      phe$orig.ident=="GSE203115ESCC2"~"Nonresponder",
                      phe$orig.ident=="GSE203115ESCC3"~"responder") 
sce.all@meta.data = phe
saveRDS(sce.all,"GSE203115.rds")
###############CSS###############
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
plot1<-FeatureScatter(sce.all,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(sce.all,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1+plot2
sce.all <- subset(sce.all, subset = nFeature_RNA > 200 & 
                    nFeature_RNA < 7000 &nCount_RNA < 60000 & percent.mt < 60)
sce.all <- NormalizeData(sce.all)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 3000)
sce.all <- ScaleData(sce.all, vars.to.regress = c('nCount_RNA'))
sce.all <- RunPCA(sce.all)

sce.all<-JackStraw(sce.all,num.replicate = 100)
sce.all<-ScoreJackStraw(sce.all)
JackStrawPlot(sce.all)
ElbowPlot(sce.all)#10

seuratObj <- RunHarmony(sce.all, "orig.ident")
names(seuratObj@reductions)

seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F ) 
seuratObj <- RunTSNE(seuratObj, dims = 1:15, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "tsne",label=F ) 

sce.all.filt=seuratObj

sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:15) 

sce.all.filt.all=sce.all.filt

for (res in c(0.01, 0.05, 0.1, 0.2, 0.3,0.4, 0.5,0.8,1)) {
  sce.all.filt.all=FindClusters(sce.all.filt.all, #graph.name = "CCA_snn", 
                                resolution = res, algorithm = 1)
}
colnames(sce.all.filt.all@meta.data)
apply(sce.all.filt.all@meta.data[,grep("RNA_snn",colnames(sce.all.filt.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 5, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.05") + 
                   ggtitle("louvain_0.05"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1") , DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
p1_dim

p2_dim=plot_grid(ncol = 4, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.4") + 
                   ggtitle("louvain_0.4"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
                   ggtitle("louvain_0.5"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"))
p2_dim

clustree(sce.all.filt.all@meta.data, prefix = "RNA_snn_res.")
table(sce.all.filt.all@active.ident) 
saveRDS(sce.all.filt.all,"F:\\Desktop\\毕业设计\\testdata\\GSE203115.rds")

sce.all.int = readRDS("F:\\Desktop\\毕业设计\\testdata\\GSE203115.rds") 
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
colnames(sce.all.int@meta.data) 
sce.all.int@meta.data <- subset(sce.all.int@meta.data, select = c(orig.ident, nCount_RNA, nFeature_RNA, group, RNA_snn_res.0.5))

scRNA=sce.all.int
genes_to_check = c('EPCAM','KRT19','KRT18',  #上皮Epithelial
                   'COL1A1' , 'DCN', 'COL1A2','COL3A1','BGN','POSTN','C1R','CFD',  #成纤维Fibroblasts
                   'CD3D', 'CD3E', 'CD3G','CD8A', 'CD4','CD2', #T
                   'MS4A1','CD19', 'CD79A',  #B
                   'KLRD1', 'NCR1',  #NK
                   'MZB1',#浆细胞
                   'CPA3' ,'KIT', 'TPSAB1','TPSB2',#肥大
                   'FLT1', 'PECAM1', 'VWF','RAMP2',  #内皮 Endothelial
                   'CD68','FCGR3A', 'TNF','IL1B','CLEC7A','IRF8','CD163', 'IL10RA','TGFB1'#巨噬
) 
DotPlot(scRNA, features = unique(genes_to_check),
        assay='RNA'  )  + coord_flip()

mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',"pink","yellow",
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)


tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.5",label = T,label.box = T)
tsne
celltype=data.frame(ClusterID=0:16,
                    celltype= 0:16) 
celltype[celltype$ClusterID %in% c( 4,7),2]='MonoMacro'
celltype[celltype$ClusterID %in% c(5,9),2]='Endothelial'
celltype[celltype$ClusterID %in% c(8),2]='Plasma'
celltype[celltype$ClusterID %in% c(10),2]='B'
celltype[celltype$ClusterID %in% c(0,1,12),2]='T'
celltype[celltype$ClusterID %in% c(3,6,14,16),2]='Fibroblasts'
celltype[celltype$ClusterID %in% c(2,13,15),2]='Epithelial'
# cluster11.markers <- FindMarkers(scRNA, ident.1 = 11, min.pct = 0.5)
# head(cluster11.markers, 5)
celltype[celltype$ClusterID %in% c(11),2]='Mast'

scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)

th=theme(axis.text.x = element_text(angle = 45,  vjust = 0.5, hjust=0.5)) 
DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.8,alpha = 0.5,
        group.by = "celltype",label = T)
colnames(scRNA@meta.data) 
DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.2,alpha = 0.1,
        group.by = "group") 


CSSgene = get(load("F:/Desktop/毕业设计/inputdata/geneset.RData"))
sc=as.matrix(scRNA@assays$RNA$counts)
CSS = gsva(gsvaParam(sc,CSSgene, kcdf = "Gaussian"))
scRNA@meta.data$CSS<- CSS[2,]-CSS[1,]
phe <-scRNA@meta.data
phe <-phe%>%mutate(
  q1 = quantile(CSS, 0.25, na.rm = TRUE),
  q3 = quantile(CSS, 0.75, na.rm = TRUE),
  dif = max(CSS) - min(CSS),
  quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
)
scRNA@meta.data$quantile_CSS_Group <- phe$quantile_CSS_Group
saveRDS(scRNA, "F:\\Desktop\\毕业设计\\testdata\\GSE203115/GSE203115.rds")

scRNA_filtered <- subset(scRNA, subset = quantile_CSS_Group != "Medium" &celltype %in%c("T","B","MonoMacro","Plasma"))
scRNA_filtered$Rgroup <- ifelse(scRNA_filtered$group == "responder", 1, 0)
saveRDS(scRNA_filtered, "F:\\Desktop\\毕业设计\\testdata\\GSE203115/GSE203115_scRNA_filtereds.rds")



##############NRCSS##############
scRNA_filtered<-readRDS( "GSE203115_scRNA_filtereds.rds")

DimPlot(scRNA_filtered, reduction = "tsne",pt.size = 0.4,alpha = 0.7,
        group.by = "celltype")
DimPlot(scRNA_filtered, reduction = "tsne",cols = c("#d56e5e", "#5390b5"),pt.size = 0.4,alpha = 0.9,
        group.by = "quantile_CSS_Group") 
DimPlot(scRNA_filtered, reduction = "tsne",cols = c("#bccf90","#f6a09a"),pt.size = 0.4,alpha = 0.9,
        group.by = "group")

plotdata<-scRNA_filtered@meta.data
plotdata[1:2,]

plotdata_summary <- plotdata %>%
  group_by(celltype, quantile_CSS_Group, group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(celltype, quantile_CSS_Group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
ggplot(plotdata_summary, aes(x = quantile_CSS_Group, y = percentage, fill =group,colour =quantile_CSS_Group )) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Type", y = "Percentage", fill = "Group") +
  theme_classic() +
  scale_fill_manual(values = c("#bccf90","#f6a09a")) +
  scale_color_manual(values = c("#d56e5e", "#5390b5")) +
  facet_wrap(~celltype, nrow = 1, strip.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(title = "Group")) 
###################IMMdeg cox  lasso#######################
library(tidyverse)
library(dplyr)
library(limma)

sc=data.frame(scRNA_filtered@assays$RNA$counts)
sampleinfo=data.frame(sampleinfo=scRNA_filtered$quantile_CSS_Group,
                      sample=colnames(sc))
design <- model.matrix(~ sampleinfo$sampleinfo)  
fit <- lmFit(sc, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, p.value=0.05, lfc=1, number=Inf)  # FC卡2 LFC为1
gene = rownames(results)#605immdegs

### COX
pfilter <- 0.05   
uniresult <- data.frame(gene = character(),
                        Estimate = numeric(),
                        StdErr = numeric(),
                        zValue = numeric(),
                        pvalue = numeric(),
                        stringsAsFactors = FALSE)
for(i in gene[1:length(gene)]){   

  unicox <- glm(scRNA_filtered$Rgroup ~ t(sc[i,]), family = binomial(link = "logit"), data = sc)  
  unisum <- summary(unicox)  
  pvalue <- unisum$coefficients[2, 4] 
  if(pvalue < pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(gene = i,
                             Estimate = unisum$coefficients[2, 1],  
                             StdErr = unisum$coefficients[2, 2], 
                             zValue = unisum$coefficients[2, 3], 
                             pvalue = unisum$coefficients[2, 4]  ))
  }} 
write.csv(uniresult,file = "COX.csv",row.names = F)

############lasso
library(glmnet)
x <- t(sc[uniresult$gene,])  
x <- scale(x)
y <- scRNA_filtered$Rgroup  

lasso_model <- cv.glmnet(x, y, alpha = 1, family = "binomial")  

best_lambda <- lasso_model$lambda.min 

coef_matrix <- coef(lasso_model, s = "lambda.min") 
coef_matrix <- as.matrix(coef_matrix)

selected_features <- rownames(coef_matrix)[coef_matrix != 0][-1]
selected_features#####197genes
write.csv(selected_features,file = "lasso.csv",row.names = F)

selected_features<-rio::import("lasso.csv")
ModeGenelist<-intersect(unlist(CSSgene),selected_features$x)
ModeGenelist<-ModeGenelist[-13]#

#########ModeGenelist########
library(caret)
library(DALEX)
library(ingredients)
Dataset1=cbind(Dataset1$Var,Dataset1[,ModeGenelist])
colnames(Dataset1)[1]<-c("Var")

Dataset2=cbind(Dataset2$Var,Dataset2[,ModeGenelist])
colnames(Dataset2)[1]<-c("Var")

fit<-glm(Var ~ .,family = binomial("logit"),data = Dataset1)
explain_titanic_glm <- explain(fit,Dataset1[,-1],Dataset1[,1])
plot(feature_importance(explain_titanic_glm, B = 1))

##########################Mime1 AUC+ROC #########
library(caret)
library(Mime1)
# scRNA_filtered@meta.data[1:3,]
scRNA_filtered$Rgroup <- ifelse(scRNA_filtered$group == "responder", "Y", "N")
trainIndex <- createDataPartition(scRNA_filtered$Rgroup, p = 0.7, list = FALSE)
trainData <- scRNA_filtered[, trainIndex]#1924cells  18218 features
testData <- scRNA_filtered[, -trainIndex]#824cells   18218 features
Dataset1 <- data.frame(ID = colnames(trainData), Var = trainData$Rgroup, t(trainData@assays$RNA$data), stringsAsFactors = FALSE)
Dataset2 <- data.frame(ID = colnames(testData), Var = testData$Rgroup, t(testData@assays$RNA$data), stringsAsFactors = FALSE)

list_train_vali_Data <- list(Dataset1=Dataset1,Dataset2=Dataset2)
list_train_vali_Data$Dataset1 [1:2,1:4]#ID 
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$Dataset1,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = ModeGenelist,
                                    methods = methods,
                                    seed = 5201314,
                                    cores_for_parallel = 16)
#nb':Naive Bayes algorithm. 'svmRadialWeights': Support Vector Machine (SVM). 'rf': Random Forest. 'kknn': K-nearest Neighbors.
#'adaboost': AdaBoost Classification Trees. 'LogitBoost':Boosted Logistic Regressions. 'cancerclass': Cancerclass.
save(res.ici,file="Mime1res.RDa")

#####
auc_vis_category_all(res.ici,dataset = c("Dataset2"),order= c("Dataset2"),
                     c("#bddd22","#ffc773","#70f3ff","#ffb3a7","#cca4e3","#a4e2c6","#f9906f"))
library(pROC)
true_labels <- factor(list_train_vali_Data$Dataset2$Var, levels = c("N", "Y"))

roc_data_list <- list()

for (model_name in names(res.ici$model)) {
  model_prob <- predict(res.ici$model[[model_name]], list_train_vali_Data$Dataset2, type = "prob")[,2]
  model_roc <- roc(true_labels, model_prob)
  roc_data_list[[model_name]] <- data.frame(
    tpr = model_roc$sensitivities,
    fpr = 1 - model_roc$specificities,
    model = model_name
  )}
roc_data <- do.call(rbind, roc_data_list)
head(roc_data)
library(ggplot2)
ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 0.8) +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  theme_classic() +
  scale_color_manual(values = c("#70f3ff","#ffb3a7","#a4e2c6","#f9906f","#ffc773","#bddd22"))+
  theme(legend.title = element_blank())

##########################other marker ROC #########
##########CD28_2746,CD274(PDL1)_8327,CD152(CTLA4)_2747,IFNG_11107
TrainRawData <- data.frame(trainData@assays$RNA$data, stringsAsFactors = FALSE)

genes <- c("CD274", "HAVCR2", "STAT1","CD28","CTLA4","IFNG")
roc_results <- data.frame(tpr = numeric(0), fpr = numeric(0), model = character(0))
for (singlegene in genes) {
  TrainData <- TrainRawData[singlegene, ] %>% t() %>% as.data.frame()
  TrainData$group <- trainData$group
  roc <- roc(response = TrainData$group, 
             predictor = TrainData[, 1],
             levels = c('responder', 'Nonresponder'))
  temp_roc <- data.frame(
    tpr = roc$sensitivities,
    fpr = 1 - roc$specificities,
    model = singlegene
  ) 
  roc_results <- rbind(roc_results, temp_roc)}


ggplot(roc_results, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 0.8) +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  theme_classic() +
  scale_color_manual(values = c("#9b59b6", "#e74c3c", "#2ecc71", "#3498db", "#f39c12", "#1abc9c"))+
  theme(legend.title = element_blank())+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") 


############scRNA——NR R############
scRNA<-get(load("scRNA.RData"))
imm <- c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono/Macro","NK","Plasma","Tprolif","Treg")
scRNA_immsubset <- subset(scRNA, subset = Celltype..major.lineage. %in% imm)
table(scRNA_immsubset$Celltype..major.lineage.)

library(Matrix)
data_list <- data_list <- list(scRNA_immsubset@assays$RNA$data.1, scRNA_immsubset@assays$RNA$data.2, scRNA_immsubset@assays$RNA$data.3,
                               scRNA_immsubset@assays$RNA$data.4, scRNA_immsubset@assays$RNA$data.5, scRNA_immsubset@assays$RNA$data.6,
                               scRNA_immsubset@assays$RNA$data.7, scRNA_immsubset@assays$RNA$data.8, scRNA_immsubset@assays$RNA$data.9,
                               scRNA_immsubset@assays$RNA$data.10, scRNA_immsubset@assays$RNA$data.11, scRNA_immsubset@assays$RNA$data.12,
                               scRNA_immsubset@assays$RNA$data.13, scRNA_immsubset@assays$RNA$data.14, scRNA_immsubset@assays$RNA$data.15,
                               scRNA_immsubset@assays$RNA$data.16, scRNA_immsubset@assays$RNA$data.17, scRNA_immsubset@assays$RNA$data.18,
                               scRNA_immsubset@assays$RNA$data.19, scRNA_immsubset@assays$RNA$data.20, scRNA_immsubset@assays$RNA$data.21,
                               scRNA_immsubset@assays$RNA$data.22, scRNA_immsubset@assays$RNA$data.23)
all_genes <- unique(unlist(lapply(data_list, rownames)))
all_cells <- unique(unlist(lapply(data_list, colnames)))
combined_data <- Matrix(0, nrow = length(all_genes), ncol = length(all_cells), sparse = TRUE)
rownames(combined_data) <- all_genes
colnames(combined_data) <- all_cells
for (data in data_list) {
  genes <- rownames(data)
  cells <- colnames(data)
  gene_indices <- match(genes, all_genes)
  cell_indices <- match(cells, all_cells)
  combined_data[gene_indices, cell_indices] <- data
}
scRNA_immsubset@assays$RNA$combined_data<-combined_data


testdata<-data.frame(t(scRNA_immsubset@assays$RNA$combined_data[ModeGenelist,]))
model_prob <- predict(res.ici$model$adaboost, testdata, type = "raw")
scRNA_immsubset$respond<-model_prob


plotdata<-scRNA_immsubset@meta.data[,c(6,7,10)]
colnames(plotdata)<-c("quantile_CSS_Group","celltype","group")
plotdata_summary <- plotdata %>%
  group_by(celltype, quantile_CSS_Group, group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(celltype, quantile_CSS_Group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
ggplot(plotdata_summary, aes(x = quantile_CSS_Group, y = percentage, fill =group,colour =quantile_CSS_Group )) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Type", y = "Percentage", fill = "Group") +
  theme_classic() +
  scale_fill_manual(values = c("#bccf90","#f6a09a")) +
  scale_color_manual(values = c("#d56e5e", "#5390b5")) +
  facet_wrap(~celltype, nrow = 1, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(title = "Group")) 
