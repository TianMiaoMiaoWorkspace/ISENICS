################################计算CSS#############################
library(GSVA)
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
library(parallel)
library(rio)
library(irGSEA)
load("F:/Desktop/HMU/2024文章/毕业设计/inputdata/geneset.RData")
csp=list(geneset$csp)
csn=list(geneset$csn)
immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono/Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
# CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
#"LUSC和LUAD是同一个"共22cancer，删去LUSC 14需特殊处理
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
                       seeds = 123, ncores = 15,  # 设置随机种子和并行核心数
                       custom = T, 
                       geneset = csp, 
                       msigdb = FALSE,
                       species = "Homo sapiens", 
                       geneid = "symbol", 
                       method = "AUCell", kcdf = "Gaussian", minGSSize = 1, maxGSSize = 2000)
  CSSn <- irGSEA.score(object = scRNAlist[[cancer]],
                       seeds = 123, ncores = 15,  # 设置随机种子和并行核心数
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
    )#保留免疫细胞非medium
  sc = sc[,which(meta$quantile_CSS_Group!="Medium" & meta$cell %in% immune_cell)]
  scmeta<-meta%>%filter(cell %in% immune_cell &quantile_CSS_Group!="Medium")
  scRNAlist[[cancer]] <- CreateSeuratObject(counts = sc,project = cancer)
  scRNAlist[[cancer]] $cell <-scmeta$cell
  scRNAlist[[cancer]] $CSS  <-scmeta$CSS
  # scRNAlist[[cancer]] $cancer  <-scmeta$cancer
  scRNAlist[[cancer]] $quantile_CSS_Group  <-scmeta$quantile_CSS_Group
  ###################细胞互作
  high<-subset(scRNAlist[[cancer]],quantile_CSS_Group=="High")
  low<-subset( scRNAlist[[cancer]],quantile_CSS_Group=="Low")
  ##################high
  if (length(unique(high$cell)) > 1 & min(table(high$cell)) > 2)  {
    # 只有多个细胞类型时才运行 CellChat
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
    print(paste(cancer,"high  CellChat 计算。"))
    high  <- createCellChat(object = high,meta = high@meta.data,group.by = "cell")#以细胞亚型分组
    
  }
  #######################low
  if  (length(unique(low$cell)) > 1 & min(table(low$cell)) > 2) {
    # 只有多个细胞类型时才运行 CellChat
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
    print(paste(cancer,"low 跳过 CellChat 计算。"))
    low  <- createCellChat(object = low,meta = low@meta.data,group.by = "cell")#以细胞亚型分组
    
  }
  ###########################合并两个chellchat########https://www.jianshu.com/p/49a0a0b50987
  cellchat_list[[cancer]]<-mergeCellChat(list(Low = low,High = high), add.names = names(list(Low = low,High = high)))
}
for (i in CANCER) {
  scRNAlist[[i]]$cancer =  i
}
save(scRNAlist,file="E:/InPut/scRNAlist.RData")
save(cellchat_list,file="E:/InPut/cellchat_list.RData")
###############衰老和非衰老细胞的UMAP图#################
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
###########数据整合
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:length(scRNAlist)], add.cell.ids = NULL, project = "MergedData")
a=scRNA@meta.data
scRNA<-NormalizeData(scRNA)%>%FindVariableFeatures(nfeatures = 3000)%>% ScaleData()
#标准化函数 SCTransform 函数可替代 NormalizeData、ScaleData 和 FindVariableFeatures 三个函数
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
#####尝试
# scRNA<-RunHarmony(scRNA,group.by.vars="cell")
# scRNA<-RunUMAP(scRNA,reduction = "harmony",dims=pc.num)
# DimPlot(scRNA, reduction = "umap", group.by = "quantile_CSS_Group", cols = color.pals, pt.size = 0.1) +
#   ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))
# DimPlot(scRNA, reduction = "umap", group.by = "cell", cols = color.pals, pt.size = 0.1) +
#   ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))
#####尝试结束

scRNA<-RunUMAP(scRNA,reduction = "harmony",dims=pc.num)
DimPlot(scRNA, reduction = "umap", group.by = "cancer", cols = color.pals, pt.size = 0.1) +
  ggtitle("UMAP (After Batch Correction)") + theme(plot.title = element_text(hjust = 0.5))

p1=DimPlot(scRNA,group.by = "cell",cols =color_list)
p2=DimPlot(scRNA,group.by = "quantile_CSS_Group",cols =c("#B04F48", "#104F80"))
p1|p2

save(scRNA,file = "E:/InPut/scRNA.RData")
##############################单细胞整合cellchat############################################
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
scRNA<-get(load("E:/InPut/scRNA.RData"))
high <- scRNA[,which(scRNA$quantile_CSS_Group=="High")]
low  <- scRNA[,which(scRNA$quantile_CSS_Group=="Low")]
##################high
expr_data <- GetAssayData(high, slot = "data")  # 使用标准化数据
cellchat <- createCellChat(object = high,
                           meta = high@meta.data,
                           group.by = "cell")#以细胞亚型分组
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) #将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- identifyOverExpressedGenes(cellchat)#识别过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#识别过表达配体受体对
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) #推断细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
high = cellchat
##################################low
expr_data <- GetAssayData(low, slot = "data")  # 使用标准化数据
cellchat <- createCellChat(object = low,
                           meta = low@meta.data,
                           group.by = "cell")#以细胞亚型分组
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) #将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- identifyOverExpressedGenes(cellchat)#识别过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#识别过表达配体受体对
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) #推断细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
low = cellchat
save(low,file = "E:/InPut/scRNA_cellchat_low.RData")
save(high,file = "E:/InPut/scRNA_cellchat_high.RData")
###########################合并两个chellchat
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
###############################互作网络图
# object.list=list(CSS_low = cellchat@net$CSS_low,CSS_high = cellchat@net$CSS_high)
# weight.max <- getMaxWeight(list(CSS_low = low,CSS_high = high),attribute = c("idents","count"))
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(list(CSS_low = low,CSS_high = high))) {
#   netVisual_circle(list(CSS_low = low,CSS_high = high)[[i]]@net$count,
#                    weight.scale = T, label.edge= F,
#                    edge.weight.max = weight.max[2],
#                    edge.width.max = 5, edge.color =color_list,
#                    title.name = paste0("Number of interactions - ",
#                                        names(list(CSS_low = low,CSS_high = high))[i]))
# }
#############比较每个信号通路的整体信息流
rankNet(cellchat, mode = "comparison",
        stacked = T, 
        comparison = c(1, 2),
        color.use = c("#104F80","#B04F48"), 
        do.stat = TRUE)
###############比较与每个细胞群相关的传出（或传入）信号
# library(ComplexHeatmap)
# library(ggpubr)
# object.list=list(CSS_low = low,CSS_high = high)
# for(i in 1:2){
#   object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
#   pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
#   ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
#   ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# }
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#########识别上调和下调的信号配体对
#sources.use = 7 表示巨噬细胞
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

device_size <- dev.size()
width <- device_size[1]  # 获取宽度
height <- device_size[2]  # 获取高度
pdf("F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/R5细胞互作mono受体配体2.pdf",
    width = width, height = height)
gg2
dev.off()
##########################单个cancer MonoMacro细胞互作#################################
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
#####################逐个手动记录细胞互作数目和强度的变化
cellchatdata<-rio::import("F:\\Desktop\\HMU\\2024文章\\毕业设计\\inputdata\\cellchat.txt")
#weight
cellchatdata$diff  = cellchatdata$`High weight` -cellchatdata$`Low weight`
cellchatdata$height = 0.5
p1 = ggplot(cellchatdata, aes(x = V1, y = height, color = diff, size = height)) +
  geom_point(alpha = 1) +  # 添加散点图
  scale_color_gradientn(
    colors = c("#104F80", "white", "#B04F48"),  # 颜色列表
    values = scales::rescale(c(0,
                               -min(cellchatdata$diff) / (max(cellchatdata$diff) - min(cellchatdata$diff)),
                               1))
  ) +  theme_classic() +theme(
    axis.text.x = element_blank(),    # 去掉横坐标文字
    axis.text.y = element_blank(),    # 去掉纵坐标文字
    axis.title.x = element_blank(),   # 去掉横坐标标题
    axis.title.y = element_blank(),   # 去掉纵坐标标题
    axis.ticks = element_blank(),     # 去掉坐标轴刻度
    axis.line = element_blank()       # 去掉坐标轴线
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
    axis.text.x = element_blank(),    # 去掉横坐标文字
    axis.text.y = element_blank(),    # 去掉纵坐标文字
    axis.title.x = element_blank(),   # 去掉横坐标标题
    axis.title.y = element_blank(),   # 去掉纵坐标标题
    axis.ticks = element_blank(),     # 去掉坐标轴刻度
    axis.line = element_blank()       # 去掉坐标轴线
  ) 

device_size <- dev.size()
width <- device_size[1]  # 获取宽度
height <- device_size[2]  # 获取高度
pdf("F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/R5单细胞细胞互作子图1.pdf",
    width = width, height = height)
p1 / p2 
dev.off()


cellchat_list<-get(load("E:/InPut/cellchat_list.RData"))
df = data.frame()
for(cancer in CANCER){
  cellchat=cellchat_list[[cancer]]
  a=data.frame(cellchat@net[["High"]][["count"]]-
                 cellchat@net[["Low"]][["count"]] )
  if ("Mono.Macro" %in% colnames(a)) {
    a = select(a, "Mono.Macro") 
    a$cancer = cancer  # 添加癌症类型列
    a$cell = rownames(a)
    df = rbind(df, a)  # 绑定到 df
  } else {
    next  # 如果不存在，则跳过当前循环
  }
}
unique_cells <- unique(df$cell)# 生成所有 cancer 和 cell 的组合
all_combinations <- expand.grid(cancer = CANCER, cell = unique_cells)
all_combinations$Mono.Macro = 0
combined_df <- all_combinations %>%
  left_join(df, by = c("cancer", "cell")) %>%
  mutate(Mono.Macro = coalesce(Mono.Macro.y, Mono.Macro.x, 0)) %>%
  select(cancer, cell, Mono.Macro)

# filtered_df <- combined_df %>%
#   group_by(cell) %>%  # 按 cell 分组
#   filter(sum(Mono.Macro == 0) < 21) %>%
#   ungroup()   # 解除分组

 p3=ggplot(combined_df, aes(x = cancer, y = cell, fill = Mono.Macro)) +
  geom_tile(na.rm = TRUE,color = "grey30") +                               # 创建热图
  scale_fill_gradientn(
    colors = c("#104F80","white","#B04F48"),
    values = scales::rescale(c(0,
                               -min(filtered_df$Mono.Macro)/(max(filtered_df$Mono.Macro)-min(filtered_df$Mono.Macro)),
                               1)
    ) )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(
    axis.title.x = element_blank(),   # 去掉横坐标标题
    axis.title.y = element_blank(),   # 去掉纵坐标标题
    axis.ticks = element_blank(),     # 去掉坐标轴刻度
    axis.line = element_blank()       # 去掉坐标轴线
  )
p3

 all <- combined_df %>%
   group_by(cancer) %>%             # 按 cell 列分组
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
     axis.text.x = element_blank(),    # 去掉横坐标文字
     axis.text.y = element_blank(),    # 去掉纵坐标文字
     axis.title.x = element_blank(),   # 去掉横坐标标题
     axis.title.y = element_blank(),   # 去掉纵坐标标题
     axis.ticks = element_blank(),     # 去掉坐标轴刻度
     axis.line = element_blank()       # 去掉坐标轴线
   )
 device_size <- dev.size()
 width <- device_size[1]  # 获取宽度
 height <- device_size[2]  # 获取高度
 pdf("F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/R5热图p4.pdf",
     width = width, height = height)
 p4
 dev.off()
###############################泛癌UMAP####################################
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
 


 pdf(file="F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/Sup1.pdf", width=width, height=height)
 ggarrange(plotlist = unlist(plots_list, recursive = FALSE), 
           ncol = 4, nrow = 44 / 4)
 dev.off()
###############################泛癌拟时序分析####################################
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
     cds_filtered <- estimateSizeFactors(cds_filtered)#估计因子大小
     cds_filtered <- estimateDispersions(cds_filtered)#估计离散度
     cds_filtered <- detectGenes(cds_filtered, min_expr = 0.01)#过滤低质量细胞
     markers <- differentialGeneTest(cds_filtered, 
                                     fullModelFormulaStr = "~quantile_CSS_Group",
                                     cores = 4 ,relative_expr = T)
     ordering_genes <-rownames(subset(markers , qval< 0.01 ))
     cds_filtered <- setOrderingFilter(cds_filtered, ordering_genes = ordering_genes)
     cds_filtered <- reduceDimension(cds_filtered, #降维
                                     method = 'PCA', 
                                     # auto_param_selection = F,ndims=10,
                                     max_components = 2, 
                                     verbose = TRUE)
     cds_filtered <-orderCells(cds_filtered)
     saveRDS(cds_filtered,file = paste0("F:/Desktop/HMU/2024文章/毕业设计/outputdata/monocle3/",cancer,"_cds.rds"))
  
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
 
 
library(patchwork)
device_size <- dev.size()
width <- device_size[1]  # 获取宽度
height <- device_size[2]  # 获取高度
pdf(file="F:/Desktop/HMU/2024文章/毕业设计/plot/文章图/Sup拟时序.pdf", width=width, height=height)
ggarrange(plotlist = unlist(plots_list, recursive = FALSE), 
          ncol = 4, nrow = 44 / 4, widths = c(1, 1.7, 1, 1.7))
dev.off()

#######################meta信息免疫治疗等分组#################################
library(ggplot2)
library(reshape2)
library(dplyr)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC",
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
cancer=CANCER[1]
meta <- read.table(paste0("F:/Desktop/毕业设计/inputdata/单细胞/", cancer, "_meta.tsv"), header = TRUE, sep = "\t")
meta <- meta[which(meta$quantile_CSS_Group != "Medium"&
                     meta$Celltype..malignancy.=="Immune cells"), ]  # 过滤数据
table(meta$Treatment)
treat<-meta[which(meta$Treatment %in%c("Anti-PD-L1")),#c("aCTLA4","aPD1","aPD1+aCTLA4")
            c("Celltype..major.lineage.","quantile_CSS_Group")]
untreat<-meta[which(meta$Treatment=="Untreated"),#"None"
            c("Celltype..major.lineage.","quantile_CSS_Group")]
treat_table <- table(treat$Celltype..major.lineage., treat$quantile_CSS_Group)
untreat_table <- table(untreat$Celltype..major.lineage., untreat$quantile_CSS_Group)
# 使用 merge 函数根据行名和列名进行合并
IMM_treat <- merge(as.data.frame(as.table(treat_table)), 
                  as.data.frame(as.table(untreat_table)), 
                  by = c("Var1", "Var2"), 
                  all = TRUE, 
                  suffixes = c("_treat", "_untreat"))
df=IMM_treat####整合
combined_df <- rbind(df,IMM_treat)#除了第一个df，后面都是result_df
result_df <- combined_df %>%
  group_by(Var1, Var2) %>%
  summarize(
    Freq_treat = sum(Freq_treat),
    Freq_untreat = sum(Freq_untreat),
    .groups = 'drop'  # 避免返回分组后的数据框
  )
saveRDS(result_df,file = "F:/Desktop/毕业设计/inputdata/单细胞/result_df.rds")
result_df=result_df[1:14,]
result_df_long <- result_df %>%
  pivot_longer(cols = starts_with("Freq"), 
               names_to = "Treatment", 
               values_to = "Frequency") %>%
  mutate(Treatment = ifelse(Treatment == "Freq_treat", "Treated", "Untreated"))
# 计算每个 Var1 和 Var2 的比例
result_df_long <- result_df_long %>%
  group_by(Var1, Treatment) %>%
  mutate(Percent = Frequency / sum(Frequency)) %>%
  ungroup()
# 绘制饼图
ggplot(result_df_long, aes(x = "", y = Percent, fill = Var2)) +
  geom_bar(stat = "identity", width = 1) +  # 绘制柱状图
  facet_wrap(~Var1 + Treatment, ncol = 2) +  # 每个 Var1 一行，分成两列（Treat和Untreat）
  coord_polar(theta = "y") +  # 转换为极坐标，绘制饼图
  theme_void() +  # 去掉背景和坐标轴
  scale_fill_manual(values = c("High" = "salmon", "Low" = "skyblue"))  # 设置颜色
###############TNM分期
cancer=CANCER[5]
meta <- read.table(paste0("F:/Desktop/毕业设计/inputdata/单细胞/", cancer, "_meta.tsv"), header = TRUE, sep = "\t")
meta <- meta[which(meta$quantile_CSS_Group != "Medium"&
                   meta$Celltype..malignancy.=="Immune cells"&
                   meta$Source == "Tumor"), ]  # 过滤数据
table(meta$TNMstage)#T3N0MX为2期， T3N2b 为3期
two<-meta[which(meta$TNMstage %in%c("T3N0MX")),
            c("Celltype..major.lineage.","quantile_CSS_Group")]
three<-meta[which(meta$TNMstage=="T3N2b"),
              c("Celltype..major.lineage.","quantile_CSS_Group")]
two_table <- table(two$Celltype..major.lineage., two$quantile_CSS_Group)
three_table <- table(three$Celltype..major.lineage., three$quantile_CSS_Group)
# 使用 merge 函数根据行名和列名进行合并
TNM <- merge(as.data.frame(as.table(two_table)), 
                   as.data.frame(as.table(three_table)), 
                   by = c("Var1", "Var2"), 
                   all = TRUE, 
                   suffixes = c("_two", "_three"))
df=TNM
df$Freq_one<-0
########################################HNSC
cancer=CANCER[9]
meta <- read.table(paste0("F:/Desktop/毕业设计/inputdata/单细胞/", cancer, "_meta.tsv"), header = TRUE, sep = "\t")
meta <- meta[which(meta$quantile_CSS_Group != "Medium"&
                     meta$Celltype..malignancy.=="Immune cells"&
                     meta$Source == "Tumor"), ]  # 过滤数据
table(meta$TNMstage)
one<-meta[which(meta$TNMstage =="I"),
          c("Celltype..major.lineage.","quantile_CSS_Group")]
two<-meta[which(meta$TNMstage %in%c("II","III")),
          c("Celltype..major.lineage.","quantile_CSS_Group")]

one_table <- table(one$Celltype..major.lineage., one$quantile_CSS_Group)
two_table <- table(two$Celltype..major.lineage., two$quantile_CSS_Group)
# 使用 merge 函数根据行名和列名进行合并
TNM <- merge(as.data.frame(as.table(one_table)), 
             as.data.frame(as.table(two_table)), 
             by = c("Var1", "Var2"), 
             all = TRUE, 
             suffixes = c("_one","_two"))
TNM$Freq_three = 0
df=rbind(df,TNM)






result_df <- df %>%
  group_by(Var1, Var2) %>%
  summarize(
    Freq_one = sum(Freq_one),
    Freq_two = sum(Freq_two),
    Freq_three = sum(Freq_three),
    .groups = 'drop'  # 避免返回分组后的数据框
  )
result_df=na.omit(result_df)



# 数据转换：根据Var1和Freq列创建比例
result_df_long <- result_df %>%
  pivot_longer(cols = starts_with("Freq"), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  group_by(Var1, Frequency_Type) %>%
  mutate(Total = sum(Frequency)) %>%
  ungroup() %>%
  mutate(Percent = Frequency / Total)  # 计算每个类别的比例

# 绘制饼图
ggplot(result_df_long, aes(x = "", y = Percent, fill = Var2)) +
  geom_bar(stat = "identity", width = 1) +  # 绘制柱状图
  coord_polar(theta = "y") +  # 转换为极坐标，绘制饼图
  facet_wrap(~ Frequency_Type~Var1 , nrow = 3) +  # 为每个Freq类型画单独的饼图
  theme_void() +  # 去掉背景和坐标轴
  labs(title = "Proportions of High and Low Groups by Freq Type") +
  scale_fill_manual(values = c("High" = "salmon", "Low" = "skyblue"))  # 设置颜色
# #######################空间转录组
# library(Seurat)
# library(SeuratObject)
# library(patchwork)
# library(ggplot2)
# library(dplyr)
# library(future)
# library(viridis)
# library(CellTrek)
# trace(CellTrek:::traint, edit = T) #修改PCA从30到5
# library(rio)
# # install.packages("Matrix", type = "source")
# # install.packages("SeuratObject", type = "source")
# # install.packages("Seurat", type = "source")
# # install.packages("irlba", type = "source")
# 
# setwd("F:/Desktop/毕业设计/testdata")
# # GBM_st = Load10X_Spatial(data.dir = "./posterior/",slice = "posterior")
# # ##标准化
# # GBM_st <- SCTransform(GBM_st,assay = "Spatial")###生成sce@SCT表达矩阵
# # ##降维
# # GBM_st <- RunPCA(GBM_st)
# # GBM_st <- RunUMAP(GBM_st, dims = 1:30, label = T)
# # ##聚类分群
# # GBM_st <- FindNeighbors(GBM_st,dims = 1:20)
# # GBM_st <- FindClusters(GBM_st,resolution = 0.1)
# # table(GBM_st$seurat_clusters)
# # #差异分析
# # dif <- FindAllMarkers(GBM_st, assay = "SCT", only.pos = T)
# # sig.dif <- dif%>%group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
# # genes <- unique(sig.dif$gene)
# 
# 
# GBM_sc<- import("F:/Desktop/毕业设计/inputdata/单细胞_原版/GBM_GSE84465_sc.txt")
# GBM_meta<- import("F:/Desktop/毕业设计/inputdata/单细胞/GBM_GSE84465_sc.tsv")
# 
# GBM_meta$cell_type = GBM_meta$quantile_CSS_Group
# GBM_meta = GBM_meta[which(GBM_meta$Celltype..major.lineage.=='Mono/Macro'),]
# 
# rownames(GBM_sc) = GBM_sc$V1
# GBM_sc=GBM_sc[,-1]  
# GBM_sc=t(GBM_sc)#行基因列细胞 GBM_sc[1:5,1:5]
# GBM_sc = GBM_sc[,which(GBM_meta$Celltype..major.lineage.=='Mono/Macro')]
# rownames(GBM_meta) = colnames(GBM_sc)
# sc <- CreateSeuratObject(counts = as.matrix(GBM_sc), 
#                          meta.data = GBM_meta)
# sc <- NormalizeData(sc)
# sc <- FindVariableFeatures(sc,nfeatures = 6000)
# sc <- ScaleData(sc, features = rownames(GBM_sc))
# sc <- RunPCA(sc,features = VariableFeatures(object=sc))
# sc <- JackStraw(sc,num.replicate = 30)
# sc <- ScoreJackStraw(sc,dims = 1:10)
# JackStrawPlot(sc)
# sc<-RunUMAP(sc,dims = 1:10)
# sc@meta.data$orig.ident<-as.factor(sc@meta.data$quantile_CSS_Group)
# sc@active.ident<-as.factor(sc@meta.data$quantile_CSS_Group)
# 
# DimPlot(sc,reduction = "umap")
# 
# brain_st_cortex <- readRDS("brain_st_cortex.rds")
# brain_st_cortex<-RenameCells(brain_st_cortex,new.names = make.names(Cells(brain_st_cortex)))
# sc<-RenameCells(sc,new.names = make.names(Cells(sc)))
# 
# 
# SpatialDimPlot(brain_st_cortex,pt.size.factor = 1.5, alpha = 0.8)
# 
# ggplot(GBM_meta, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
#   geom_point(size=0.02) +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(color = "grey80", size = 0.5),  # 主要网格线
#         panel.grid.minor = element_line(color = "grey80", size = 0.25)  # 次要网格线
#   )+
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# 
# brain_traint <-CellTrek::traint(st_data = brain_st_cortex,
#                                 sc_data = sc,
#                                 sc_assay = "RNA",
#                                 st_assay = "Spatial",
#                                 cell_names = "quantile_CSS_Group")
# 
# DimPlot(brain_traint, group.by = "type") 
# brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint,
#                                      int_assay='traint', 
#                                      sc_data=sc, 
#                                      sc_assay = 'RNA', 
#                                      reduction='pca', 
#                                      intp=T, intp_pnt=5000, 
#                                      intp_lin=F, nPCs=5, ntree=1000, 
#                                      dist_thresh=0.55, 
#                                      top_spot=5, spot_n=5, 
#                                      repel_r=20, repel_iter=20, 
#                                      keep_model=T)$celltrek
# 
# brain_celltrek$cell_type <- factor(brain_celltrek$cell_type,
#                                    levels=sort(unique(brain_celltrek$cell_type)))
# 
# CellTrek::celltrek_vis(brain_celltrek@meta.data %>%
#                          dplyr::select(coord_x, coord_y, cell_type:id_new),
#                        brain_celltrek@images$anterior1@image, 
#                        brain_celltrek@images$anterior1@scale.factors$lowres)
# #细胞亚群的空间距离
# inp_df <- brain_celltrek@meta.data %>% dplyr::select(cell_names = dplyr::one_of('cell_type'), 
#                                                           coord_x, coord_y)
# inp_df$coord_x = 270-inp_df$coord_x
# head(inp_df)
# output <- kdist(inp_df = inp_df, 
#                 ref = "High", #目标细胞类型
#                 ref_type = 'all', 
#                 que = c("High","Low"),  #其余的细胞
#                 k = 10, 
#                 new_name = "High vs Low",
#                 keep_nn = F)
# head(output$kdist_df) 
# res = output$kdist_df
# res$barcode = row.names(res)
# inp_df$barcode = row.names(inp_df)
# res = left_join(res, inp_df)
# head(res)
# library(ggpubr)
# ggboxplot(data = res, 
#           x = "cell_names",
#           y = "High vs Low", 
#           fill = "cell_names", 
#           title = " ")+ 
#   stat_compare_means(method = "t.test") +
#   theme(plot.title = element_text(color="black",hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5), #,vjust = 0.5
#         legend.position = "none") + labs(y = "Distance")

#######11.20整合免疫拟时序
# library(Seurat)
# library(data.table)
# library(dplyr)
# library(tidyverse)
# library(stringr)
# library(harmony)
# library(cowplot)
# scRNA<-get(load("E:/InPut/scRNA.RData"))
# # exp<-as(as.matrix(scRNA@assays$RNA$counts),'sparseMatrix')
# p_data <- scRNA@meta.data
# p_data$celltype<-scRNA$orig.ident
# f_data<-data.frame(gene_short_name = rownames(scRNA), row.names = rownames(scRNA))
# pd <- new("AnnotatedDataFrame", p_data)
# fd <- new("AnnotatedDataFrame", f_data)
# cds <- newCellDataSet(exp,
#                       phenoData = pd, featureData = fd,
#                       lowerDetectionLimit = 0.5,
#                       expressionFamily = uninormal())#指定用于分析的统计分布，通常是 negbinomial 或 negbinomial.size
# cds <- estimateSizeFactors(cds)#估计因子大小
# # cds <- estimateDispersions(cds, num_threads = 4)#离散度#Error: estimateDispersions only works, and is only needed, when you're using a CellDataSet with a negbinomial or negbinomial.size expression family
# cds <- detectGenes(cds, min_expr = 0.1)#过滤低质量细胞
# # cds <- setOrderingFilter(cds, f_data$gene_short_name) 
# # plot_ordering_genes(cds)
# # diff<-differentialGeneTest(cds,fullModelFormulaStr = "~cell_type",cores = 4)
# cds <- reduceDimension(cds,  max_components = 2, norm_method = 'none')
# save(cds,file="整合scRNAcds.RData")
# load("整合scRNAcds.RData")
# cds <- orderCells(cds)
# 
# plot_cell_trajectory(cds, color_by = "quantile_CSS_Group",cell_size = 1) + 
#   theme_bw(base_rect_size = 1.5)+
#   scale_color_manual(values = c("#B04F48","#104F80"))+
#   theme(
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none",
#     axis.title = element_blank(),
#     plot.title = element_text(hjust = 0.5))
# plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1) + 
#   theme_bw(base_rect_size = 1.5)+
#   theme(
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     # legend.position = "none",
#     axis.title = element_blank(),
#     plot.title = element_text(hjust = 0.5))
