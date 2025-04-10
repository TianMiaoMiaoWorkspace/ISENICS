################################CSS#############################
library(GSVA)
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
library(parallel)
library(rio)
library(irGSEA)
load("inputdata/geneset.RData")
csp=list(geneset$csp)
csn=list(geneset$csn)
immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono/Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
# CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")

input<-"E:\\InPut/singlecell/"
h5<-list.files(input,pattern = "\\.h5$",full.names = T)
scRNAlist=list()
cellchat_list=list()
for (file in h5[1:22]) {
  cancer = strsplit(basename(file), "_")[[1]][1]
  data<- Read10X_h5(file)
  data<-(2^data - 1) * 10
  sc<-as(as.matrix(data), "dgCMatrix")
  scRNAlist[[cancer]] <- CreateSeuratObject(counts = sc,project = cancer)
  # CSS <- gsva(gsvaParam(sc,geneset,kcdf="Gaussian"))
  CSSp <- irGSEA.score(object = scRNAlist[[cancer]],
                       seeds = 123, ncores = 15, 
                       custom = T, 
                       geneset = csp, 
                       msigdb = FALSE,
                       species = "Homo sapiens", 
                       geneid = "symbol", 
                       method = "AUCell", kcdf = "Gaussian", minGSSize = 1, maxGSSize = 2000)
  CSSn <- irGSEA.score(object = scRNAlist[[cancer]],
                       seeds = 123, ncores = 15,  
                       custom = T, 
                       geneset = csn, 
                       msigdb = FALSE,
                       species = "Homo sapiens", 
                       geneid = "symbol", 
                       method = "AUCell", kcdf = "Gaussian", minGSSize = 1, maxGSSize = 2000)
  CSS=CSSp@assays$AUCell@counts[1,]-CSSn@assays$AUCell@counts[1,]
  meta <-import(paste0("E:\\InPut/celltype/",cancer,"_celltype.txt"))
  colnames(meta)<-"cell"
  meta$CSS <- CSS
  meta = meta %>% group_by(cell) %>% 
    mutate(
      q1 = quantile(CSS, 0.25, na.rm = TRUE),
      q3 = quantile(CSS, 0.75, na.rm = TRUE),
      quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
    )
  sc = sc[,which(meta$quantile_CSS_Group!="Medium" & meta$cell %in% immune_cell)]
  scmeta<-meta%>%filter(cell %in% immune_cell &quantile_CSS_Group!="Medium")
  scRNAlist[[cancer]] <- CreateSeuratObject(counts = sc,project = cancer)
  scRNAlist[[cancer]] $cell <-scmeta$cell
  scRNAlist[[cancer]] $CSS  <-scmeta$CSS
  # scRNAlist[[cancer]] $cancer  <-scmeta$cancer
  scRNAlist[[cancer]] $quantile_CSS_Group  <-scmeta$quantile_CSS_Group
  ###################
  high<-subset(scRNAlist[[cancer]],quantile_CSS_Group=="High")
  low<-subset( scRNAlist[[cancer]],quantile_CSS_Group=="Low")
  ##################high
  if (length(unique(high$cell)) > 1 & min(table(high$cell)) > 2)  {
  
    cellchat <- createCellChat(object = high, group.by = "cell")
    cellchat@DB <- CellChatDB.human 
    cellchat <- subsetData(cellchat)  
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
    cellchat <- computeCommunProbPathway(cellchat)#pathway level
    cellchat <- aggregateNet(cellchat)
    high = cellchat
  } else {
 
    high  <- createCellChat(object = high,meta = high@meta.data,group.by = "cell")
    
  }
  #######################low
  if  (length(unique(low$cell)) > 1 & min(table(low$cell)) > 2) {

    cellchat <- createCellChat(object = low, group.by = "cell")
    cellchat@DB <- CellChatDB.human 
    cellchat <- subsetData(cellchat)  
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
    cellchat <- computeCommunProbPathway(cellchat)#pathway level
    cellchat <- aggregateNet(cellchat)
    low = cellchat
  } else {
   
    low  <- createCellChat(object = low,meta = low@meta.data,group.by = "cell")
    
  }
  ###########################
  cellchat_list[[cancer]]<-mergeCellChat(list(Low = low,High = high), add.names = names(list(Low = low,High = high)))
}
for (i in CANCER) {
  scRNAlist[[i]]$cancer =  i
}
save(scRNAlist,file="E:/InPut/scRNAlist.RData")
save(cellchat_list,file="E:/InPut/cellchat_list.RData")
###############UMAP#################
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
  "#FF99CC",#LUAD   NSCLC
  # "#D44842",#LUSC
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
library(Seurat)
library(dplyr)
library(CellChat)
library(data.table)
library(harmony)
load("E:/InPut/scRNAlist.RData")
load("E:/InPut/cellchat_list.RData")
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","NSCLC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
###########
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:length(scRNAlist)], add.cell.ids = NULL, project = "MergedData")
a=scRNA@meta.data
scRNA<-NormalizeData(scRNA)%>%FindVariableFeatures(nfeatures = 3000)%>% ScaleData()

scRNA<-RunPCA(scRNA,verbose=F)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "cancer", cols = color.pals) +
  ggtitle("UMAP (Before Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))

scRNA<-RunHarmony(scRNA,group.by.vars="cancer")
scRNA<-RunUMAP(scRNA,reduction = "harmony",dims=pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "cancer", cols = color.pals, pt.size = 0.1) +
  ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))

scRNA<-RunUMAP(scRNA,reduction = "harmony",dims=pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "cancer", cols = color.pals, pt.size = 0.1) +
  ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))

p1=DimPlot(scRNA,group.by = "cell",cols =color_list)
p2=DimPlot(scRNA,group.by = "quantile_CSS_Group",cols =c("#B04F48", "#104F80"))
p1|p2

save(scRNA,file = "E:/InPut/scRNA.RData")
##############################cellchat############################################
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
scRNA<-get(load("E:/InPut/scRNA.RData"))
high <- scRNA[,which(scRNA$quantile_CSS_Group=="High")]
low  <- scRNA[,which(scRNA$quantile_CSS_Group=="Low")]
##################high
expr_data <- GetAssayData(high, slot = "data") 
cellchat <- createCellChat(object = high,
                           meta = high@meta.data,
                           group.by = "cell")
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)#
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
high = cellchat
##################################low
expr_data <- GetAssayData(low, slot = "data")
cellchat <- createCellChat(object = low,
                           meta = low@meta.data,
                           group.by = "cell")
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
low = cellchat
save(low,file = "E:/InPut/scRNA_cellchat_low.RData")
save(high,file = "E:/InPut/scRNA_cellchat_high.RData")

cellchat <- mergeCellChat(list(CSS_low = low,CSS_high = high),  add.names = names(list(CSS_low = low,CSS_high = high)))
save(cellchat,file = "E:/InPut/scRNA_cellchat.RData")
load("E:/InPut/scRNA_cellchat.RData")
netVisual_heatmap(cellchat,
                  comparison = c(1,2),
                  color.heatmap = c("#104F80", "#B04F48"),
                  font.size = 15,
                  measure = "weight",
                  color.use = color_list[1:12],
                  cluster.cols = T)

#############
rankNet(cellchat, mode = "comparison",
        stacked = T, 
        comparison = c(1, 2),
        color.use = c("#104F80","#B04F48"), 
        do.stat = TRUE)


gg1=netVisual_bubble(cellchat, sources.use = 7, 
                     comparison = c(1, 2), 
                     targets.use = c(7),
                     thresh = 0.01,
                     max.dataset = 2,
                     grid.on = F,
                     title.name = "Increased in CSS_high", 
                     angle.x = 45, remove.isolate = T)
gg2=netVisual_bubble(cellchat, sources.use = 3, 
                     comparison = c(1, 2),
                     targets.use = c(3),
                     thresh = 0.01,
                     max.dataset = 1,
                     grid.on = F,
                     title.name = "Decreased in CSS_high",
                     angle.x = 45, remove.isolate = T)
ggarrange(gg1, gg2, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right")


##########################cancer MonoMacro#################################
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","NSCLC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")

cellchat_list<-get(load("E:/InPut/cellchat_list.RData"))
cancer=CANCER[22]
cancer
cellchat=cellchat_list[[cancer]]
gg1 <- compareInteractions(cellchat, show.legend = F,
                           group = c(1,2),
                           size.text = 8,
                           color.use = c("#104F80","#B04F48"))
gg2 <- compareInteractions(cellchat, show.legend = F,
                           group = c(1,2),
                           size.text = 8,
                           color.use = c("#104F80","#B04F48"), 
                           measure = "weight")
ggarrange(gg1,gg2)
#####################
cellchatdata<-rio::import("cellchat.txt")
#weight
cellchatdata$diff  = cellchatdata$`High weight` -cellchatdata$`Low weight`
cellchatdata$height = 0.5
p1 = ggplot(cellchatdata, aes(x = V1, y = height, color = diff, size = height)) +
  geom_point(alpha = 1) +  
  scale_color_gradientn(
    colors = c("#104F80", "white", "#B04F48"), 
    values = scales::rescale(c(0,
                               -min(cellchatdata$diff) / (max(cellchatdata$diff) - min(cellchatdata$diff)),
                               1))
  ) +  theme_classic() +theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),   
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.ticks = element_blank(),  
    axis.line = element_blank()    
  ) 
#count
cellchatdata$diff  = cellchatdata$`High count` -cellchatdata$`Low count`
cellchatdata$height = 0.5
# summary(interaction$diff)
p2 = ggplot(cellchatdata, aes(x = V1, y = height, color = diff, size = height)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colors = c("#388E3C", "white", "#B04F48"),
    values = scales::rescale(c(0,
                               -min(cellchatdata$diff)/(max(cellchatdata$diff)-min(cellchatdata$diff)),
                               1)
    ) ) +
  theme_classic()+theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),   
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),    
    axis.line = element_blank()   
  ) 




cellchat_list<-get(load("E:/InPut/cellchat_list.RData"))
df = data.frame()
for(cancer in CANCER){
  cellchat=cellchat_list[[cancer]]
  a=data.frame(cellchat@net[["High"]][["count"]]-
                 cellchat@net[["Low"]][["count"]] )
  if ("Mono.Macro" %in% colnames(a)) {
    a = select(a, "Mono.Macro") 
    a$cancer = cancer 
    a$cell = rownames(a)
    df = rbind(df, a) 
  } else {
    next  
  }
}
unique_cells <- unique(df$cell)
all_combinations <- expand.grid(cancer = CANCER, cell = unique_cells)
all_combinations$Mono.Macro = 0
combined_df <- all_combinations %>%
  left_join(df, by = c("cancer", "cell")) %>%
  mutate(Mono.Macro = coalesce(Mono.Macro.y, Mono.Macro.x, 0)) %>%
  select(cancer, cell, Mono.Macro)

 p3=ggplot(combined_df, aes(x = cancer, y = cell, fill = Mono.Macro)) +
  geom_tile(na.rm = TRUE,color = "grey30") +                            
  scale_fill_gradientn(
    colors = c("#104F80","white","#B04F48"),
    values = scales::rescale(c(0,
                               -min(filtered_df$Mono.Macro)/(max(filtered_df$Mono.Macro)-min(filtered_df$Mono.Macro)),
                               1)
    ) )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(
    axis.title.x = element_blank(), 
    axis.ticks = element_blank(),    
    axis.line = element_blank()    
  )
p3

 all <- combined_df %>%
   group_by(cancer) %>%         
   summarise(cell = "all", Mono.Macro = sum(Mono.Macro))%>%
   mutate(height = 0.5)
p4=ggplot(all, aes(x = cancer, y = height,fill = Mono.Macro)) +
   geom_col(color = "grey50") +
   scale_fill_gradientn(
     colors = c("#104F80","white","#B04F48"),
     values = scales::rescale(c(0,
                                -min(all$Mono.Macro)/(max(all$Mono.Macro)-min(all$Mono.Macro)),
                                1)) )+
   theme_classic()+theme(
     axis.text.x = element_blank(), 
     axis.text.y = element_blank(),   
     axis.title.x = element_blank(),  
     axis.title.y = element_blank(),  
     axis.ticks = element_blank(),   
     axis.line = element_blank()    
   )

###############################UMAP####################################
 library(Seurat)
 library(dplyr)
 library(CellChat)
 library(data.table)
 color_list <- c("#ffe0e9", "#4c8dae", "#8ca1c4","#2e4e7e","#f47983", "#ff4777","#d07b76","#ff9b60",
                 "#e8ac9d","#ffb1c9","#48c0a3", "#a4e2c6")
 load("E:/InPut/scRNAlist.RData")
 CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","NSCLC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
 plots_list <- list()
 for (cancer in CANCER) {
   scRNA <- scRNAlist[[cancer]]
   scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
   scRNA <- RunPCA(scRNA, verbose = FALSE)
   scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:20)
   p1 <- DimPlot(scRNA, group.by = "cell", cols = color_list) + 
     theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + 
     ggtitle(cancer)
   p2 <- DimPlot(scRNA, group.by = "quantile_CSS_Group", 
                 cols = c("#B04F48", "#104F80")) + 
     theme(legend.position = "none", plot.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
   
   plots_list[[cancer]] <- list(p1 = p1, p2 = p2)
 }
 


###################################################################
library(rio)
library(ggplot2)
library(monocle)
library(tidyverse)
library(Biobase)
library(ggpubr)
library(dplyr)
library(Seurat)
load("E:/InPut/scRNAlist.RData")
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","NSCLC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
plots_list <- list()
 for (cancer in CANCER[7:22]) {
   scRNA <- scRNAlist[[cancer]]
   pd <- new("AnnotatedDataFrame", as.data.frame(scRNA@meta.data))
   sampleNames(pd) <-colnames(scRNA)
   fd <- data.frame(gene_short_name = rownames(scRNA), row.names = rownames(scRNA))
   fd <- new("AnnotatedDataFrame", fd)
   cds <- newCellDataSet(GetAssayData(object = scRNA, slot = "counts", assay = "RNA"),
                         phenoData = pd, featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
   pheno_data <- pData(cds)
   sample_names <- rownames(pheno_data)
   if(nrow(pheno_data) >0){
     cds_filtered <- cds[, sample_names]
     cds_filtered <- estimateSizeFactors(cds_filtered)
     cds_filtered <- estimateDispersions(cds_filtered)
     cds_filtered <- detectGenes(cds_filtered, min_expr = 0.01)
     markers <- differentialGeneTest(cds_filtered, 
                                     fullModelFormulaStr = "~quantile_CSS_Group",
                                     cores = 4 ,relative_expr = T)
     ordering_genes <-rownames(subset(markers , qval< 0.01 ))
     cds_filtered <- setOrderingFilter(cds_filtered, ordering_genes = ordering_genes)
     cds_filtered <- reduceDimension(cds_filtered, 
                                     method = 'PCA', 
                                     # auto_param_selection = F,ndims=10,
                                     max_components = 2, 
                                     verbose = TRUE)
     cds_filtered <-orderCells(cds_filtered)

   p1 <- plot_cell_trajectory(cds_filtered, color_by = "quantile_CSS_Group",cell_size = 1) + 
     theme_bw(base_rect_size = 1.5)+
     # labs(title=cancer)+
     scale_color_manual(values = c("#B04F48","#104F80"))+
     theme(
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank(),
       axis.text = element_blank(),
       axis.ticks = element_blank(),
       legend.position = "none",
       axis.title = element_blank(),
       plot.title = element_text(hjust = 0.5))
   p2 <- plot_cell_trajectory(cds_filtered, color_by = "Pseudotime",cell_size = 1) + 
     theme_bw(base_rect_size = 1.5)+
     # labs(title=cancer)+
     theme(
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank(),
       axis.text = element_blank(),
       axis.ticks = element_blank(),
       # legend.position = "none",
       axis.title = element_blank(),
       plot.title = element_text(hjust = 0.5))+
     scale_color_gradientn(colors = c("grey", "orange"))+
     guides(color = guide_colorbar(barwidth = 1, barheight = 3))
   
   plots_list[[cancer]] <- list(p1 = p1, p2 = p2)
   }
   print(paste("ok!",cancer))
 }  
