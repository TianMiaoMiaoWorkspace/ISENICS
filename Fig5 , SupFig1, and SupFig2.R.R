library(ggplot2)
library(ggridges)
library(tidyverse)
library(ggsignif)
library(patchwork)
library(rio)
library(GSVA)


###############UMAP#################
load("geneset.RData")
color_list <- c("#b2cbe6", "#9FA8DA", "#FCE4EC", "#FFB3C8", "#e6a0af", "#EF9A9A", 
                "#E57373", "#EF5350", "#B71C1C", "#EC407A", "#C2185B", "#880E4F", 
                "#CE93D8", "#cb86b5", "#9C27B0", "#7B1FA2", "#D4F69F", "#A1FFA1", 
                "#91FF91", "#76d4a4", "#81C784", "#2E7D32", "#FFFDE7", "#FFECB3", 
                "#DCE775", "#C0CA33", "#FFEE58", "#FDD835", "#F9A825", "#B2EBF2", 
                "#4DD0E1", "#00BCD4", "#0097A7", "#006064", "#D7CCC8", "#A1887F")

names(color_list) = c(
  "AC_Like_Malignant", "Malignant",
  "B", "Mast", "MonoMacro", "Neutrophils", "NK", "pDC", "Plasma",
  "CD4Tconv", "CD8T", "CD8Tex", "DC", "Promonocyte", "Tprolif", "Treg", 
  "Pericytes", "Astrocyte", "Neuron", "Oligodendrocyte", "SMC", "Vascular", 
  "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts", "Myocyte", "Acinar", "Ductal",
  "Hepatic progenitor", "HSC", "Progenitor", "OPC", "Endometrial stromal cells",
  "EryPro", "GMP"
)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")

cancer_list <- CANCER[11:23]
plot_list <- list() 
for (cancer in cancer_list) {
  meta <- read.table(paste0( cancer, "_meta.tsv"), header = TRUE, sep = "\t")
  meta <- meta[which(meta$quantile_CSS_Group != "Medium"), ]  
  # p1 - Celltype 
  p1 <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = Celltype..major.lineage.)) +
    geom_point(size = 0.02) + theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "white")) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_list) 
  # p2 - quantile_CSS_Group 
  p2 <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = quantile_CSS_Group)) +
    geom_point(size = 0.02) +theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "white")) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = c("#B04F48", "#104F80")) 
  
  plot_list[[cancer]] <- p1 - p2 + plot_annotation(title = cancer)  
}

combined_plot <- wrap_plots(plot_list, ncol = 2)

pdf("F:\\Desktop\\1.pdf",width = 8.3 ,height = 11.7) 
combined_plot
dev.off()
pdf("F:\\Desktop\\2.pdf",width = 8.3 ,height = 11.7) 
combined_plot
dev.off()

meta <- read.table(paste0(cancer, "_meta.tsv"), header = TRUE, sep = "\t")

#############################################R5############################################
library(dplyr)
folder_path <- "single cell"
files <- list.files(path = folder_path, pattern = "_meta.tsv$", full.names = TRUE)
merged_data <- data.frame()
for (file in files) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  selected_columns <- df %>%
    select(UMAP_1, UMAP_2, quantile_CSS_Group, Celltype..major.lineage.)%>%
    filter(quantile_CSS_Group!="Medium")
  merged_data <- bind_rows(merged_data, selected_columns)
}
head(merged_data)
write.table(merged_data, "META_combined_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

library(Seurat)
library(data.table)
library(harmony)

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
scRNAlist <-list()
for (cancer in CANCER) {
  sc<-fread(paste0( cancer, ".txt"))
  rownames(sc)=sc$V1
  sc=sc[,-1]
  sc<-data.frame(t(sc))
  sc<-as(as.matrix(sc), "dgCMatrix")
  meta <- read.table(paste0( cancer, "_meta.tsv"), header = TRUE, sep = "\t")%>%
    select(UMAP_1, UMAP_2, quantile_CSS_Group, Celltype..major.lineage.)
  rownames(meta)<-colnames(sc)
  sc = sc[,which(meta$quantile_CSS_Group!="Medium")]
  meta = meta[which(meta$quantile_CSS_Group!="Medium"),]
  scRNAlist[[cancer]] <- CreateSeuratObject(counts = sc, meta.data = meta)
}
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA<-NormalizeData(scRNA)%>%FindVariableFeatures(nfeatures = 1000)%>% ScaleData()
scRNA<-RunPCA(scRNA,verbose=F)
# ElbowPlot(scRNA, ndims = 50)
pc.num=1:20

scRNA<-FindNeighbors(scRNA,dim=pc.num)%>%FindClusters()

color_list <- c("#b2cbe6", "#9FA8DA", "#FCE4EC", "#FFB3C8", "#e6a0af", "#EF9A9A", 
                "#E57373", "#EF5350", "#B71C1C", "#EC407A", "#C2185B", "#880E4F", 
                "#CE93D8", "#cb86b5", "#9C27B0", "#7B1FA2", "#D4F69F", "#A1FFA1", 
                "#91FF91", "#76d4a4", "#81C784", "#2E7D32", "#FFFDE7", "#FFECB3", 
                "#DCE775", "#C0CA33", "#FFEE58", "#FDD835", "#F9A825", "#B2EBF2", 
                "#4DD0E1", "#00BCD4", "#0097A7", "#006064", "#D7CCC8", "#A1887F")

names(color_list) = c(
  "AC_Like_Malignant", "Malignant", 
  "B", "Mast", "MonoMacro", "Neutrophils", "NK", "pDC", "Plasma",
  "CD4Tconv", "CD8T", "CD8Tex", "DC", "Promonocyte", "Tprolif", "Treg", 
  "Pericytes", "Astrocyte", "Neuron", "Oligodendrocyte", "SMC", "Vascular", 
  "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts", "Myocyte", "Acinar", "Ductal", 
  "Hepatic progenitor", "HSC", "Progenitor", "OPC", "Endometrial stromal cells",
  "EryPro", "GMP" 
)
scRNA<-SCTransform(scRNA)
scRNA<-RunHarmony(scRNA,group.by.vars="Celltype..major.lineage.",assay.use="SCT",max.iter.harmony=20)
scRNA<-RunTSNE(scRNA,reduction = "harmony",dims=pc.num)%>%RunUMAP(reduction = "harmony",dims=pc.num)
p1=DimPlot(scRNA,group.by = "Celltype..major.lineage.",cols =color_list)
p2=DimPlot(scRNA,group.by = "quantile_CSS_Group",cols =c("#B04F48", "#104F80"))
ggsave("F:/Desktop/UMAP1.pdf",p1,width = 8,height = 6)
ggsave("F:/Desktop/UMAP2.pdf",p2,width = 8,height = 6)


celltype_counts <- table(scRNA$Celltype..major.lineage.)
celltype_df <- as.data.frame(celltype_counts)
colnames(celltype_df) <- c("Celltype", "Count")
filtered_celltypes <- celltype_df[celltype_df$Count > 100, ]
valid_celltypes <- filtered_celltypes$Celltype
scRNA <- subset(scRNA, subset = Celltype..major.lineage. %in% valid_celltypes)
scRNA <- JoinLayers(scRNA)
LayerData(scRNA, assay = "RNA", layer = "counts")



library(monocle)
library (gplots) 
library(ggthemes)
library(dplyr)

scRNA$Celltype..major.lineage. <- gsub("AC-like Malignant","AC_like_Malignant",scRNA$Celltype..major.lineage. )
tb=table(scRNA$quantile_CSS_Group, scRNA$Celltype..major.lineage.)
balloonplot(tb)
bar_data <- as.data.frame(tb)
bar_data$Var2<-gsub("Mono/Macro","MonoMacro",bar_data$Var2)
# bar_data$Var2<-gsub("AC-like Malignant","AC like Malignant",bar_data$Var2)

bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)

ggplot(bar_per, aes(y = percent, x = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity")  +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=color_list)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




######cellchat#######
#weight
interaction=import("cellchat inferred interactions strength.txt")
interaction$diff  = interaction$CSS_high -interaction$CSS_low
interaction$height = 0.5
p1=ggplot(interaction, aes(x = cancer, y = height, color = diff, size = height)) +
  geom_point(alpha = 1) +  
  scale_color_gradientn(
    colors = c("#104F80", "white", "#B04F48"), 
    values = scales::rescale(c(0,
                               -min(interaction$diff) / (max(interaction$diff) - min(interaction$diff)),
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
interaction=import("cellchat inferred interactions.txt")
interaction$diff  = interaction$CSS_high -interaction$CSS_low
interaction$height = 0.5
# summary(interaction$diff)
p2 = ggplot(interaction, aes(x = cancer, y = height, color = diff, size = height)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colors = c("#388E3C", "white", "#B04F48"),
    values = scales::rescale(c(0,
                               -min(interaction$diff)/(max(interaction$diff)-min(interaction$diff)),
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


 p1 / p2 

###########heatmap

list = list.files("cellchatRDATA")
df = data.frame()
for(cancer in CANCER){
  cellchat = get(load(list[grep(cancer, list)]))
  a=data.frame(cellchat@net[["CSS_high"]][["count"]]-
                 cellchat@net[["CSS_low"]][["count"]] )
  if ("Malignant" %in% colnames(a)) {
    a = select(a, "Malignant") 
    a$cancer = cancer
    a$cell = rownames(a)
    df = rbind(df, a) 
  } else {
    next 
  }}

unique_cells <- unique(df$cell)

all_combinations <- expand.grid(cancer = CANCER, cell = unique_cells)
all_combinations$Malignant = 0
combined_df <- all_combinations %>%
  left_join(df, by = c("cancer", "cell")) %>%

  mutate(Malignant = coalesce(Malignant.y, Malignant.x, 0)) %>%
  select(cancer, cell, Malignant)

filtered_df <- combined_df %>%
  group_by(cell) %>%  
  filter(sum(Malignant == 0) < 21) %>%
  ungroup() 

p2 = ggplot(filtered_df, aes(x = cancer, y = cell, fill = Malignant)) +
  geom_tile(na.rm = TRUE,color = "grey30") +                            
  scale_fill_gradientn(
    colors = c("#104F80","white","#B04F48"),
    values = scales::rescale(c(0,
                               -min(filtered_df$Malignant)/(max(filtered_df$Malignant)-min(filtered_df$Malignant)),
                               1)
    ) )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),  
    axis.ticks = element_blank(),     
    axis.line = element_blank()     
  )
p2


Malignant <- filtered_df %>%
  group_by(cancer) %>%            
  summarise(total_Malignant = sum(Malignant, na.rm = TRUE))%>%
  mutate(height = 0.5)
p3 = ggplot(Malignant, aes(x = cancer, y = height,fill = total_Malignant)) +
  geom_col(color = "grey50") +
  scale_fill_gradientn(
    colors = c("#104F80","white","#B04F48"),
    values = scales::rescale(c(0,
                               -min(Malignant$total_Malignant)/(max(Malignant$total_Malignant)-min(Malignant$total_Malignant)),
                               1)
                        ) )+
  theme_classic()+theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),   
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),  
    axis.ticks = element_blank(),    
    axis.line = element_blank()     
  )
p3


######################################cellchat######################################
library(GSVA)
library(Seurat)
library(dplyr)
library(CellChat)
library(ggpubr)
library(data.table)
scRNA<-get(load("scRNA.RData"))
imm <- c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono/Macro","NK","Plasma","Tprolif","Treg")
scRNA_immsubset <- subset(scRNA, subset = Celltype..major.lineage. %in% imm)
gene=rownames(scRNA_immsubset@assays$RNA$scale.data)
scRNA_immsubset <- subset(scRNA_immsubset, features = gene)
scRNA_immsubset <- JoinLayers(scRNA_immsubset)

meta=scRNA_immsubset@meta.data
high <- scRNA_immsubset[,which(meta$quantile_CSS_Group=="High")]
low  <- scRNA_immsubset[,which(meta$quantile_CSS_Group=="Low")]
##################high
expr_data <- GetAssayData(high, slot = "data")
cellchat <- createCellChat(object = high,
                           meta = high@meta.data,
                           group.by = "Celltype..major.lineage.")
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
high = cellchat
##################################low
expr_data <- GetAssayData(low, slot = "data") 
cellchat <- createCellChat(object = low,
                           meta = low@meta.data,
                           group.by = "Celltype..major.lineage.")
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)#pathway level
cellchat <- aggregateNet(cellchat)
low = cellchat

#https://www.jianshu.com/p/49a0a0b50987
cellchat <- mergeCellChat(list(CSS_low = low,CSS_high = high),  add.names = names(list(CSS_low = low,CSS_high = high)))
# save(cellchat,file = "scRNA_imm_cellchat.RData")
################################
netVisual_heatmap(cellchat, 
                        comparison = c(1,2),
                        color.heatmap = c("#104F80", "#B04F48"),
                        font.size = 15,
                        measure = "weight",
                        color.use = c("#FCE4EC","#EC407A","#C2185B", "#880E4F", "#CE93D8","#FFB3C8", "#e6a0af","#E57373","#B71C1C", "#9C27B0", "#7B1FA2"),
                        cluster.cols = T)
###############################
weight.max <- getMaxWeight(list(CSS_low = low,CSS_high = high),attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(list(CSS_low = low,CSS_high = high))) {
    netVisual_circle(list(CSS_low = low,CSS_high = high)[[i]]@net$count, 
                     weight.scale = T, label.edge= F, 
                     edge.weight.max = weight.max[2], 
                     edge.width.max = 12, 
                     title.name = paste0("Number of interactions - ", 
                                         names(list(CSS_low = low,CSS_high = high))[i]))
  }

#############pathway
rankNet(cellchat, mode = "comparison",
                 stacked = T, 
                 comparison = c(1, 2),
                 color.use = c("#104F80","#B04F48"), 
                 do.stat = TRUE)

library(ComplexHeatmap)
library(ggpubr)
object.list=list(CSS_low = low,CSS_high = high)
for(i in 1:2){
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
}
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


gg1=netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:15),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CSS_high", angle.x = 90, remove.isolate = T)
gg2=netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:15),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in CSS_high", angle.x = 90, remove.isolate = T)


ggarrange(gg1,gg2,common.legend = T, legend="right")



###############################monocol####################################
library(rio)
library(ggplot2)
library(monocle)
library(tidyverse)
library(Biobase)
library(ggpubr)
library(dplyr)
library(Seurat)
library(GSVA)
library(CellChat)

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC",
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
for (cancer in CANCER[21]) {
  list <- list.files()
  list<-list[grep(paste0(cancer), list)]
  data<-import(list[grep(".txt", list)])
  meta<-import(list[grep(".tsv", list)])
  rownames(data) = data$V1
  data=data[,-1]  
  data=t(data)
  pd <- new("AnnotatedDataFrame", as.data.frame(meta))
  sampleNames(pd) <-colnames(data)
  fd <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
  fd <- new("AnnotatedDataFrame", fd)
  cds <- newCellDataSet(as(as.matrix(data),"sparseMatrix"),
                        phenoData = pd, featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
  pheno_data <- pData(cds)
  filtered_samples <- pheno_data[pheno_data$quantile_CSS_Group!="Medium"&pheno_data$Celltype..malignancy.=="Immune cells", ]
  sample_names <- rownames(filtered_samples)
  if(nrow(filtered_samples) >0){
    cds_filtered <- cds[, sample_names]
    cds_filtered <- estimateSizeFactors(cds_filtered)
    cds_filtered <- estimateDispersions(cds_filtered)
    cds_filtered <- detectGenes(cds_filtered, min_expr = 0.01)
    markers <- differentialGeneTest(cds_filtered, 
                                    fullModelFormulaStr = "~quantile_CSS_Group",
                                    cores = 4 ,relative_expr = T)
    ordering_genes <-rownames(subset(markers , qval< 0.01 ))
    cds_filtered <- setOrderingFilter(cds_filtered, ordering_genes = ordering_genes)
    # plot_ordering_genes(cds_filtered)
    cds_filtered <- reduceDimension(cds_filtered, 
                                    method = 'PCA', 
                                    max_components = 2, 
                                    verbose = TRUE)
    cds_filtered <-orderCells(cds_filtered)
    saveRDS(cds_filtered, file = paste0(cancer,"_cds.rds"))
  }
}


cds_filtered= readRDS(paste0(cancer,"_cds.rds"))
plot_cell_trajectory(cds_filtered, color_by = "quantile_CSS_Group",cell_size = 1) + 
  theme_bw(base_rect_size = 1.5)+
  labs(title=cancer)+
  scale_color_manual(values = c("#B04F48","#104F80"))+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
plot_cell_trajectory(cds_filtered, color_by = "Pseudotime",cell_size = 1) + 
  theme_bw(base_rect_size = 1.5)+
  labs(title=cancer)+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
colnames(pData(cds_filtered))
expressed_genes=row.names(subset(fData(cds_filtered),num_cells_expressed>=10))
pseudotime_de <- differentialGeneTest(cds_filtered[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(cds_filtered[expressed_genes,],fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

pseudotime_diff<-pseudotime_de[pseudotime_de$use_for_ordering=="TRUE"&pseudotime_de$pval<=0.05,]
states_diff<-states_de[states_de$use_for_ordering=="TRUE"&states_de$pval<=0.05,]
load("F:/Desktop/毕业设计/outputdata/alldiff.RData")
GBM_diff<-data[which(data$info=="GBM_GBM_GSE84465_MonoMacro.csv_diffsig.txt"),]
a=intersect(GBM_diff$V1,pseudotime_diff$gene_short_name)##48
a=intersect(GBM_diff$V1,states_de$gene_short_name)##478
gene=intersect(a,rownames(cds_filtered))
plot_pseudotime_heatmap(cds_filtered[1:100,], 
                             num_cluster = 2,
                             show_rownames = T
                        )

#######sup2###########
device_size <- dev.size()
width <- device_size[1]  
height <- device_size[2] 
plot_list <- list() 
for (cancer in CANCER[1:23]) {
  cds_filtered= readRDS(paste0("F:/Desktop/毕业设计/inputdata/单细胞/",cancer,"_cds.rds"))
  p <- plot_cell_trajectory(cds_filtered, color_by = "quantile_CSS_Group",cell_size = 1) + 
    theme_bw(base_rect_size = 1.5)+
    labs(title=cancer)+
    scale_color_manual(values = c("#B04F48","#104F80"))+
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5))
  plot_list[[cancer]] <- p + plot_annotation(title = cancer) 
}
library(patchwork)
combined_plot <- wrap_plots(plot_list, ncol = 3)
pdf(file="Sup2.pdf", width=width, height=height)
combined_plot
dev.off()


plot_list <- list()  
for (cancer in CANCER[1:23]) {
  cds_filtered= readRDS(paste0(cancer,"_cds.rds"))
  p <- plot_cell_trajectory(cds_filtered, color_by = "Pseudotime",cell_size = 1) + 
    theme_bw(base_rect_size = 1.5)+
    labs(title=cancer)+
    # scale_color_manual(values = c("#B04F48","#104F80"))+
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
  plot_list[[cancer]] <- p + plot_annotation(title = cancer) 
}
library(patchwork)
combined_plot <- wrap_plots(plot_list, ncol = 3)
combined_plot

#######################meta#################################
library(ggplot2)
library(reshape2)
library(dplyr)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC",
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
cancer=CANCER[5]
meta <- meta[which(meta$quantile_CSS_Group != "Medium"&
                   meta$Celltype..malignancy.=="Immune cells"&
                   meta$Source == "Tumor"), ] 
table(meta$TNMstage)
two<-meta[which(meta$TNMstage %in%c("T3N0MX")),
            c("Celltype..major.lineage.","quantile_CSS_Group")]
three<-meta[which(meta$TNMstage=="T3N2b"),
              c("Celltype..major.lineage.","quantile_CSS_Group")]
two_table <- table(two$Celltype..major.lineage., two$quantile_CSS_Group)
three_table <- table(three$Celltype..major.lineage., three$quantile_CSS_Group)

TNM <- merge(as.data.frame(as.table(two_table)), 
                   as.data.frame(as.table(three_table)), 
                   by = c("Var1", "Var2"), 
                   all = TRUE, 
                   suffixes = c("_two", "_three"))
df=TNM
df$Freq_one<-0
########################################HNSC
cancer=CANCER[9]

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
    .groups = 'drop' 
  )
result_df=na.omit(result_df)




result_df_long <- result_df %>%
  pivot_longer(cols = starts_with("Freq"), 
               names_to = "Frequency_Type", 
               values_to = "Frequency") %>%
  group_by(Var1, Frequency_Type) %>%
  mutate(Total = sum(Frequency)) %>%
  ungroup() %>%
  mutate(Percent = Frequency / Total) 


ggplot(result_df_long, aes(x = "", y = Percent, fill = Var2)) +
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") +
  facet_wrap(~ Frequency_Type~Var1 , nrow = 3) + 
  theme_void() + 
  labs(title = "Proportions of High and Low Groups by Freq Type") +
  scale_fill_manual(values = c("High" = "salmon", "Low" = "skyblue")) 
