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
library(rio)
library(DoubletFinder)
library(scCDC)

###############seurat GSE203115###############
folder_path <- ".\\GSE203115"
subfolders <- list.dirs(folder_path, full.names = TRUE, recursive = TRUE)
for (subfolder in subfolders) {
  gz_files <- list.files(subfolder, pattern = "\\.gz$", full.names = TRUE)
  if (length(gz_files) > 0) {
    for (gz_file in gz_files) {
      gunzip(gz_file, overwrite = TRUE) 
    }
  }
}

dir="./GSE203115/raw/"
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
#####merge sample
do.call(rbind,lapply(sceList, dim))
GSE203115=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  )

########添加meta信息
phe = GSE203115@meta.data
table(phe$orig.ident)
# View(phe)
phe$response = case_when(phe$orig.ident=="GSE203115ESCC1"~"R",
                      phe$orig.ident=="GSE203115ESCC2"~"NR",
                      phe$orig.ident=="GSE203115ESCC3"~"R")
GSE203115@meta.data = phe

saveRDS(GSE203115,".\\GSE203115.rds")

###################质量控制#####################
GSE203115<-import("GSE203115.txt")
GSE203115<-GSE203115%>%column_to_rownames("V1")
GSE203115<- CreateSeuratObject(count, min.cells = 3, min.features = 200)  # 创建Seurat对象

GSE203115<-readRDS("GSE203115.rds")


count=GSE203115
res<-add_percent_features(count,"human")
res<-QC_mad_filter(res)
VlnPlot(res, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb"),
        ncol = 5, pt.size = 0.1, combine = TRUE) +theme(plot.title = element_text(size = 10))
res<-QC_SC_FilterContamination(res,nCount_RNA_filter = 500,
                               nFeature_RNA_filter = 250,
                               percent.mt_filter = 12, ####difined by vlnplot
                               percent.rb_filter=50,
                               percent.hb_filter=2)
best.pc<-QC_best.pc(res,pc.contribution=90)
#dim.usage difined by best.pc
dim.usage=21

res<-QC_RM_doublet(res,dim.usage,doublet_rate=8)
res<-QC_RNA_corrected(res,dim.usage)
res <- subset(res, Doublet == "Singlet")
res@meta.data <- select(res@meta.data, -Doublet)
saveRDS(res,"QCGSE203115.RData")



#### same work flow for GSE145281 and GSE207422 were omitted

# GSE145281:14474->12242
# GSE207422:83230->79204
# GSE203115:12762->9380


###################细胞注释#####################

GSE207422<-readRDS("GSE207422.RData")
GSE207422$Doublet<-NULL
GSE207422$orig.ident<-GSE207422$sample
GSE207422$sample<-NULL
GSE207422$source<-"GSE207422"

GSE145281<-readRDS("GSE145281.RData")
GSE145281$response <- case_when(
  GSE145281$orig.ident %in% c("R1" ,"R2" ,"R3", "R4" ,"R5") ~ "R",
  GSE145281$orig.ident %in% c("NR1" ,"NR2" ,"NR3", "NR4" ,"NR5") ~ "NR"
)
GSE145281$source<-"GSE145281"

GSE203115<-readRDS("GSE203115.rds")
GSE203115$source<-"GSE203115"
GSE203115$response = case_when(GSE203115$orig.ident=="GSE203115ESCC1"~"R",
                               GSE203115$orig.ident=="GSE203115ESCC2"~"NR",
                               GSE203115$orig.ident=="GSE203115ESCC3"~"R")
# seurat_list<-list(GSE203115=GSE203115,GSE145281=GSE145281,GSE207422=GSE207422)
# sce.all <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
sce.all<-GSE207422
scRNA<-NormalizeData(sce.all)%>%FindVariableFeatures(nfeatures = 3000)%>% ScaleData()
scRNA<-RunPCA(scRNA,verbose=F)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP (Before Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))

scRNA<-RunHarmony(scRNA,group.by.vars="orig.ident")
scRNA<-RunUMAP(scRNA,reduction = "harmony",dims=pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", pt.size = 0.1,raster=FALSE) +
  ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))


###############222222细胞注释计算CSS###############
# library(dplyr)
# library(Seurat)
# library(SingleR)
# library(celldex)
scRNA<- FindNeighbors(scRNA,dims = 1:30) 
scRNA<- FindClusters(scRNA, resolution = 0.5, algorithm = 1)
p1=DimPlot(scRNA,reduction = "umap", group.by = "RNA_snn_res.0.5" ,label =T) 
###########single R
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <-BlueprintEncodeData()
anno <- SingleR(scRNA@assays$RNA@counts,
                ref = list(BP=bpe.se,HPCA=hpca.se),
                labels = list(bpe.se$label.fine,hpca.se$label.main),
                clusters = scRNA@meta.data$RNA_snn_res.0.5)
plotScoreHeatmap(anno,scores.use = c(0),clusters = anno@rownames,show_colnames = T)
celltype = data.frame(ClusterID=rownames(anno),
                      celltype=anno$labels,
                      stringsAsFactors = F)
scRNA@meta.data$singleR = celltype[match(scRNA@meta.data$RNA_snn_res.0.5,celltype$ClusterID),'celltype']
p2=DimPlot(scRNA,reduction = "umap", group.by = "singleR" ) 
p1|p2
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


p3=DotPlot(scRNA, features = unique(genes_to_check),
        assay='RNA'  )  + coord_flip()

p1|p2|p3
#####203115
if(T){
  celltype=data.frame(ClusterID=0:19,
                      celltype= 0:19) 
  celltype[celltype$ClusterID %in% c(11),2]='B'
  celltype[celltype$ClusterID %in% c(9),2]='Plasma'
  celltype[celltype$ClusterID %in% c(0,1,2,16),2]='T'
  celltype[celltype$ClusterID %in% c(13),2]='Mast'
  celltype[celltype$ClusterID %in% c(3,14),2]='Mono/Macro'
  celltype[celltype$ClusterID %in% c(7),2]='Myofibroblasts'
  celltype[celltype$ClusterID %in% c(4,18),2]='Endothelial'
  celltype[celltype$ClusterID %in% c(6,8,10,19),2]='Epithelial'
  celltype[celltype$ClusterID %in% c(5,12,17),2]='Fibroblasts'
  # celltype[celltype$ClusterID %in% c(0,2),2]='NK'
  celltype[celltype$ClusterID %in% c(15),2]='Neutrophils'
}
if(T){#145281
  celltype=data.frame(ClusterID=0:18,
                      celltype= 0:18) 
  celltype[celltype$ClusterID %in% c(10),2]='B'
  celltype[celltype$ClusterID %in% c(2,4,5,9,11,12,14,17),2]='T'
  celltype[celltype$ClusterID %in% c(13),2]='Mast'
  celltype[celltype$ClusterID %in% c(0,1,3,6,7,8,15,16,18),2]='Mono/Macro'
  celltype[celltype$ClusterID %in% c(13),2]='NK'
}
if(T){#207422
  celltype=data.frame(ClusterID=0:26,
                      celltype= 0:26) 
  celltype[celltype$ClusterID %in% c(2,26),2]='B'
  celltype[celltype$ClusterID %in% c(10),2]='Plasma'
  celltype[celltype$ClusterID %in% c(0,1,4,5,6,13,20,24),2]='T'
  celltype[celltype$ClusterID %in% c(19),2]='Mast'
  celltype[celltype$ClusterID %in% c(8,11,12,15,23),2]='Mono/Macro'
  celltype[celltype$ClusterID %in% c(25),2]='Endothelial'
  celltype[celltype$ClusterID %in% c(7,16,17,18,22),2]='Epithelial'
  celltype[celltype$ClusterID %in% c(21),2]='Fibroblasts'
  celltype[celltype$ClusterID %in% c(3,9,14),2]='Neutrophils'
}
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)


p4=DimPlot(scRNA, reduction = "umap",group.by = "celltype",raster=FALSE,label =T,
        cols=c("#ffe0e9", "#FFfcc6", "#ffeda8","#fff897",
               "#ff4777","#d07b76","#ecfc00",
               "#e8ac9d","#ff9b60","#ffb1c9","#8ca1c4")) 
p2|p4
#######################计算全体CSS
library(irGSEA)
CSSgene = get(load("F:/Desktop/HMU/2024文章/毕业设计/inputdata/geneset.RData"))
immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono/Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
csp=list(geneset$csp)
csn=list(geneset$csn)
CSSp <- irGSEA.score(object = scRNA,
                     seeds = 123, ncores = 15,  # 设置随机种子和并行核心数
                     custom = T, 
                     geneset = csp, 
                     msigdb = FALSE,
                     species = "Homo sapiens", 
                     geneid = "symbol", 
                     method = "AUCell", kcdf = "Gaussian", minGSSize = 1, maxGSSize = 2000)
CSSn <- irGSEA.score(object = scRNA,
                     seeds = 123, ncores = 15,  # 设置随机种子和并行核心数
                     custom = T, 
                     geneset = csn, 
                     msigdb = FALSE,
                     species = "Homo sapiens", 
                     geneid = "symbol", 
                     method = "AUCell", kcdf = "Gaussian", minGSSize = 1, maxGSSize = 2000)
CSS=CSSp@assays$AUCell@counts[1,]-CSSn@assays$AUCell@counts[1,]
scRNA$CSS<- CSS

CSS=scRNA$CSS
q1 = quantile(CSS, 0.25, na.rm = TRUE)
q3 = quantile(CSS, 0.75, na.rm = TRUE)
scRNA$quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))

saveRDS(scRNA, "E:\\InPut\\GSE207422_cell.rds")
scRNA<-readRDS("E:\\InPut\\GSE203115_cell.rds")

################整合3个效果最好###################
GSE207422<-readRDS("E:\\InPut\\GSE207422_cell.rds")
GSE145281<-readRDS("E:\\InPut\\GSE145281_cell.rds")
GSE203115<-readRDS("E:\\InPut\\GSE203115_cell.rds")

seurat_list<-list(GSE203115=GSE203115,GSE145281=GSE145281,GSE207422=GSE207422)
sce.all <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
sce.all<-  saveRDS(sce.all,"E:\\InPut\\免疫应答\\all3_cell.rds")
#######筛选细胞
all<-  readRDS("E:\\InPut\\免疫应答\\all3_cell.rds")
table(all@meta.data$celltype)
immune_cell=c("B","T","Plasma",'Mast',"Mono/Macro",'NK','Neutrophils')
all <- subset(all, subset = celltype != "NK")##264太少/104208
all <- subset(all, subset = quantile_CSS_Group != "Medium" &celltype %in%immune_cell)
all<-NormalizeData(all)%>%FindVariableFeatures(nfeatures = 3000)%>% ScaleData()
all<-RunPCA(all,verbose=F)
ElbowPlot(all, ndims = 50)
pc.num=1:30
all <- RunUMAP(all, reduction = "pca", dims = pc.num)
DimPlot(all, reduction = "umap", group.by = "source",cols=c("#f47983","#4b5cc4","#ffc773")) +
  ggtitle("UMAP (Before Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))
all<-RunHarmony(all,group.by.vars="source")
all<-RunUMAP(all,reduction = "harmony",dims=pc.num)
DimPlot(all, reduction = "umap", group.by = "source",cols=c("#f47983","#4b5cc4","#ffc773")) +
  ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))
#####这是对每个source分别划分亚型后计算CSS分组的
p1=DimPlot(all, reduction = "umap",cols = c("#d56e5e", "#5390b5"),pt.size = 0.1,
        group.by = "quantile_CSS_Group") 
p2=DimPlot(all, reduction = "umap",cols = c("#bccf90","#f6a09a"),pt.size = 0.1,
        group.by = "response")
p3=DimPlot(all, reduction = "umap", cols=c("#ffe0e9", "#ff4777","#d07b76","#e8ac9d","#ffb1c9","#8ca1c4"),pt.size = 0.1,
        group.by = "celltype")

plotdata<-all@meta.data
plotdata_summary <- plotdata %>%
  group_by(celltype, quantile_CSS_Group, response) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(celltype, quantile_CSS_Group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
p4=ggplot(plotdata_summary, aes(x = quantile_CSS_Group, y = percentage, fill =response,colour =quantile_CSS_Group )) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Type", y = "Percentage", fill = "Group") +
  theme_classic() +coord_flip()+
  scale_fill_manual(values = c("#bccf90","#f6a09a")) +
  scale_color_manual(values = c("#d56e5e", "#5390b5")) +
  facet_wrap(~celltype, ncol = 1, strip.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 如果需要，旋转x轴标签以便阅读
  guides(fill = guide_legend(title = "Group")) # 自定义图例标题
p4
(p3|p4)/(p1|p2)
##################卡方检验##################
all<-  readRDS("E:\\InPut\\免疫应答\\all3_cell.rds")
table(all@meta.data$celltype)
immune_cell=c("B","T","Plasma",'Mast',"Mono/Macro",'NK','Neutrophils')
all <- subset(all, subset = celltype != "NK")##264太少/104208
all <- subset(all, subset = quantile_CSS_Group != "Medium" &celltype %in%immune_cell)
table(all@meta.data$celltype)

chi_result <- chisq.test(all$quantile_CSS_Group,all$response)
col_palette <- colorRampPalette(c("#a3d900", "white", "#c32136"))(20)#标准化残差热图
corrplot(chi_result$residuals, is.cor = FALSE, method = "circle",col = col_palette)
mosaic(table(all$quantile_CSS_Group,all$response), shade = TRUE, legend = TRUE)
###################333333CSSIMMdeg与113diff交集CSSgene交集的单因素cox和lasso筛选模型基因#######################
scRNA_filtered=all#####43143->43055
library(tidyverse)
library(dplyr)
library(limma)
library(Matrix) 
sc=scRNA_filtered@assays$RNA@counts
CSSgene = get(load("F:/Desktop/HMU/2024文章/毕业设计/inputdata/geneset.RData"))
CSSgene =unlist(CSSgene)
diff<-read.table("F:/Desktop/HMU/2024文章/毕业设计/outputdata/sig_113files_DEG.txt",sep = "\t")
# #无法分配大小为8.9 Gb的向量
# 按数据集拆分子集

groups <- unique(scRNA_filtered$celltype)
fit_list <- list()
for (grp in groups) {
  sub_cells <- which(scRNA_filtered$celltype == grp)
  sc_sub <- sc[, sub_cells]
  sampleinfo_sub <-   data.frame(
    sampleinfo = scRNA_filtered$quantile_CSS_Group[sub_cells],  
    sample = colnames(sc_sub)
  )
  sampleinfo_sub$sampleinfo <- factor(sampleinfo_sub$sampleinfo)
  design_sub <- model.matrix(~ sampleinfo_sub$sampleinfo)
  fit_sub <- lmFit(sc_sub, design_sub)
  fit_list[[grp]] <- fit_sub
}


res=data.frame()
for (i in 1:6) {
  fit=fit_list[[1]]
  fit <- eBayes(fit)
  results <- topTable(fit, coef=2, p.value=0.05, lfc=1, number=Inf)  # FC卡2 LFC为1
  res=rbind(res,results)
}
gene = rownames(res)#1350immdegs
mygene=intersect(diff$symbol,gene)#67immdegs
scRNA_filtered$Rgroup <- ifelse(scRNA_filtered$response == "R", 1, 0)
### 单因素COX回归分析
pfilter <- 0.05   
uniresult <- data.frame(gene = character(),
                        Estimate = numeric(),
                        StdErr = numeric(),
                        zValue = numeric(),
                        pvalue = numeric(),
                        stringsAsFactors = FALSE)
for(i in mygene[2:length(mygene)]){   
  # 使用逻辑回归分析每个基因对 Rgroup 的影响
  unicox <- glm(scRNA_filtered$Rgroup ~ sc[i,], family = binomial(link = "logit"), data = sc)  
  unisum <- summary(unicox)  
  pvalue <- unisum$coefficients[2, 4] 
  if(pvalue < pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(gene = i,
                             Estimate = unisum$coefficients[2, 1],  # 回归系数
                             StdErr = unisum$coefficients[2, 2],  # 标准误
                             zValue = unisum$coefficients[2, 3],  # z值
                             pvalue = unisum$coefficients[2, 4]  ))
  }print(i)
  } 
write.csv(uniresult,file = "单因素COX分析结果.csv",row.names = F)
uniresult<-rio::import("单因素COX分析结果.csv")
############lasso
library(glmnet)
x <- t(sc[uniresult$gene,])  # 自变量：基因表达数据
x <- scale(x)
y <- scRNA_filtered$Rgroup  # 目标变量：分类标签（例如：0 和 1）
# 使用cv.glmnet进行LASSO回归并选择最佳lambda
lasso_model <- cv.glmnet(x, y, alpha = 1, family = "binomial")  # alpha = 1表示LASSO回归,而alpha = 0则表示Ridge回归
#查看最佳lambda值
best_lambda <- lasso_model$lambda.min  # 最佳的lambda值（最小的交叉验证误差）
#查看选中的特征
coef_matrix <- coef(lasso_model, s = "lambda.min") # 查看在最佳lambda下被选中的变量系数
coef_matrix <- as.matrix(coef_matrix)
# 查看非零系数的变量名称
selected_features <- rownames(coef_matrix)[coef_matrix != 0][-1]
selected_features#####62genes
write.csv(selected_features,file = "lasso结果.csv",row.names = F)
selected_features<-rio::import("lasso结果.csv")

ModeGenelist<-intersect(selected_features$x,CSSgene )
#########ModeGenelist重要度计算########
load("F:/Desktop/HMU/2024文章/毕业设计/testdata/GSE203115/Mime1res.RDa")
library(caret)
library(DALEX)
library(ingredients)
# Dataset1=cbind(Dataset1$Var,Dataset1[,ModeGenelist])
# colnames(Dataset1)[1]<-c("Var")
trainData <- scRNA_filtered[ModeGenelist, trainIndex]#30140cells  27839 features
Dataset1 <- data.frame(Var = trainData$Var,
                       t(trainData@assays$RNA@data), 
                       stringsAsFactors = FALSE)
Dataset1$Var=ifelse(Dataset1$Var=="Y",1,0)

#单一的线性逻辑回归模型
fit<-glm(Var ~ .,family = binomial("logit"),data = Dataset1)
explain_titanic_glm <- explain(fit,Dataset1[,-1],Dataset1[,1])
plot(feature_importance(explain_titanic_glm, B = 1)) 

##########################Mime1免疫应答分类器 AUC+ROC #########
library(caret)
library(Mime1)
library(tidyverse)
library(dplyr)
library(Matrix)
setwd("E:\\InPut\\免疫应答")
load("E:/InPut/免疫应答/.RData")
selected_features<-rio::import("lasso结果.csv")
ModeGenelist <-c("ETS1","GABARAP","HMGB1","TUBA1B","TUBB","RPS16","PTMA")
ModeGenelist <-selected_features$x
scRNA_filtered$Var = case_when(scRNA_filtered$response=="R"~ "Y",scRNA_filtered$response=="NR"~ "N")

# scRNA<-subset(scRNA_filtered, source == "GSE203115")#####数据太大了，仅取一个子集
trainIndex <- createDataPartition(scRNA_filtered$Var, p = 0.6, list = FALSE)
trainData <- scRNA_filtered[, trainIndex]#30140cells  27839 features
testData <- scRNA_filtered[, -trainIndex]#12915cells   27839 features
Dataset1 <- data.frame(ID = colnames(trainData), Var = factor(trainData$Var),
                       t(trainData@assays$RNA@data), 
                       stringsAsFactors = FALSE)

# Dataset1 <- Dataset1[c(1:30,30135:30140),]
Dataset2 <- data.frame(ID = colnames(testData), Var = factor(testData$Var),
                       t(testData@assays$RNA@data), 
                       stringsAsFactors = FALSE)

# Dataset2<- Dataset2[c(1:10,12910:12915),]
list_train_vali_Data <- list(Dataset1=Dataset1,Dataset2=Dataset2)

# load("F:/Desktop/HMU/2024文章/毕业设计/testdata/Mime1/Example.ici.Rdata")
# load("F:/Desktop/HMU/2024文章/毕业设计/testdata/Mime1/genelist.Rdata")
# install.packages("F:/Download/Mime-fe3309253bacb468b6df0184e8365514803a295b.zip", repos = NULL, type = "win.binary")
# trace("ML.Dev.Pred.Category.Sig", edit = TRUE)###手动更改参数method
# methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost')
res.ici <- mymime(train_data = list_train_vali_Data$Dataset1,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = ModeGenelist,
                                    methods = c('nb','rf','kknn','LogitBoost','adaboost'),
                                    seed = 5201314,
                                    cores_for_parallel = 16)

#nb':Naive Bayes algorithm. 'svmRadialWeights': Support Vector Machine (SVM). 'rf': Random Forest. 'kknn': K-nearest Neighbors.
#'adaboost': AdaBoost Classification Trees. 'LogitBoost':Boosted Logistic Regressions. 'cancerclass': Cancerclass.
save(res.ici,file="Mime1res.RDa")

#####可视化
auc_vis_category_all(res.ici,dataset = c("Dataset2"),order= c("Dataset2"),
                     c("#bddd22","#ffc773","#70f3ff","#ffb3a7","#cca4e3"))
library(pROC)
res.ici[["auc"]][["Dataset2"]]
true_labels <- factor(list_train_vali_Data$Dataset2$Var, levels = c("N", "Y"))
# 创建一个空的数据框用于存储所有模型的 ROC 数据
roc_data_list <- list()
# 循环遍历每个模型进行 ROC 曲线计算
methods = c('nb','rf','kknn','LogitBoost','adaboost')
for (model_name in methods) {
  model_prob <- predict(res.ici$model[[model_name]], list_train_vali_Data$Dataset2, type = "prob")[,2]
  model_roc <- roc(true_labels, model_prob)
  roc_data_list[[model_name]] <- data.frame(
    tpr = model_roc$sensitivities,
    fpr = 1 - model_roc$specificities,
    model = model_name
  )}##cancerclass画不出来
roc_data <- do.call(rbind, roc_data_list)
head(roc_data)
library(ggplot2)
ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 0.8) +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  theme_classic() +
  scale_color_manual(values = c("#70f3ff","#ffb3a7","#a4e2c6","#f9906f","#ffc773"))+
  theme(legend.title = element_blank())

##########################其他免疫应答maker单基因ROC #########
library(pROC)
trainIndex <- createDataPartition(scRNA_filtered$Var, p = 0.7, list = FALSE)
trainData <- scRNA_filtered[, trainIndex]#30140cells  27839 features
Dataset1 <- data.frame(Var =factor(trainData$Var),
                       t(trainData@assays$RNA@data), 
                       stringsAsFactors = FALSE)
# 基因列表
genes <- c("CD274", "HAVCR2", "STAT1","CD28","CTLA4","IFNG")
roc_results <- data.frame(tpr = numeric(0), fpr = numeric(0), model = character(0))
for (singlegene in genes) {
  TrainData <- Dataset1[,singlegene] %>% as.data.frame()
  TrainData$Var <- Dataset1$Var
  roc <- roc(response = TrainData$Var, 
             predictor = TrainData[, 1],
             levels = c('Y', 'N'))
  print(paste(singlegene, auc(roc))) 
  temp_roc <- data.frame(
    tpr = roc$sensitivities,
    fpr = 1 - roc$specificities,
    model = singlegene
  ) 
  roc_results <- rbind(roc_results, temp_roc)}

# 可视化
ggplot(roc_results, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 0.8) +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  theme_classic() +
  scale_color_manual(values = c("#9b59b6", "#e74c3c", "#2ecc71", "#3498db", "#f39c12", "#1abc9c"))+
  theme(legend.title = element_blank())+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") 


