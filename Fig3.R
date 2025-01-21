#####################Bar graph of the percentage of differentially expressed genes##############################
diff<-read.table("sig_118files_DEG.txt",sep = "\t")
library(dplyr)
library(purrr)
library(tidyr)
result <- diff %>%
  group_by(cell, sig) %>%
  summarise(count = n(), .groups = 'drop')
total_counts <- result %>% group_by(cell) %>% summarise(total = sum(count))
percentages <- result %>%
  left_join(total_counts, by = "cell") %>%
  mutate(percentage = (count / total) * 100 )

p1= ggplot(percentages, aes(x = cell, y = percentage, fill = sig)) +
  geom_bar(stat = "identity", position = "stack",color="white") +
  scale_fill_manual(values = c("#5390b5", "#d56e5e"))+
  labs(title = "Percentage of Up and Down by File Category",
       x = "",
       y = "Percentage",
       fill = " ") +
  theme_classic()+
  geom_text(aes(label = count), 
            y = ifelse(percentages$sig == "Up", 
                       10,  
                       90), 
            color = "black",cex=3) 

########################################volcano#######################################
library(ggplot2)
library(ggrepel) 
library(data.table)
FC = 1.2
p = 0.05
data <- fread("ESCA_GSE173950_MonoMacro.csv", header = TRUE, sep = ",")
data=data.frame(t(data))
sample = data[1,]
gene  =  rownames(data)[-1]
data=data[-1,]
data_numeric <- as.data.frame(lapply(data, as.numeric))
colnames(data_numeric) = sample
rownames(data_numeric) = gene
data1 = data_numeric[rowSums(data_numeric) > 0, ]
gene_means <- rowMeans(data1)
gene_percentiles <- rank(gene_means) / length(gene_means)
threshold <- 0.25  
CSS_Group <- IMMfjjcss[IMMfjjcss$cancer =="ESCA" & IMMfjjcss$cell=="MonoMacro" ,c("sample","quantile_CSS_Group")]
filtered_genes <- data1[gene_percentiles > threshold, CSS_Group$sample]
Group = factor(unlist(CSS_Group$quantile_CSS_Group), levels = c("High","Low")) 
design <- model.matrix(~0+Group)
colnames(design)=levels(Group)
rownames(design)=CSS_Group$sample
dge <- DGEList(counts=filtered_genes) %>% calcNormFactors()
keep <- rowSums(cpm(dge)>0.5) >= 2
dge<-dge[keep,,keep.lib.sizes=FALSE]
dge<-calcNormFactors(dge)
fit <- voom(dge, design, normalize="quantile") %>% lmFit(design)
constrasts = paste(rev(levels(Group)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()
DEG = topTable(fit2, coef=constrasts, n=Inf) %>%  na.omit()
DEG$symbol = rownames(DEG)
DEG$sig = ifelse(DEG$logFC >=log2(FC),"Up",ifelse(DEG$logFC <=-log2(FC),"Down","Notsig"))

DEG=DEG[-1*log10(DEG$P.Value)<=30,]
# DEG$label=ifelse(-1*log10(DEG$P.Value) >= 30 | abs(DEG$logFC)>=1.5 ,
#                       as.character(DEG$symbol), '')

ggplot(DEG, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = sig)) +
  labs(x = "log2(FC)", y = "-log10(p)") +
  theme_classic() +
  scale_color_manual(values = c("#5390b5", "grey", "#d56e5e"))+
  theme_classic()+labs(title = "ESCA_GSE173950_MonoMacro")

 

######bubble #######
# symbol (3430genes)
diff<-read.table("sig_118files_DEG.txt",sep = "\t")
symbol_counts <- diff %>%
  group_by(symbol) %>%
  summarise(cell_count = n_distinct(cell)) %>%
  arrange(desc(cell_count))
mydiff <- diff %>%left_join(symbol_counts, by = "symbol") %>%
  filter(cell_count >= 3 & abs(logFC)>=0.4) %>% arrange(desc(cell_count)) %>%
  mutate(logP = -log(adj.P.Val))

p3=ggplot(mydiff, aes(x = symbol, y = cell, fill = logFC, size = logP)) +
  geom_point(alpha = 0.9, shape = 21,stroke = 0,color = "white") +  
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#E53935", midpoint = 0) + 
  labs( x = "", y = "",fill = "Log Fold Change", size = "-log10(P Value)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

p3 + p1 + patchwork::plot_layout(ncol = 1, heights = c(4,1))


#########################venn############################
diff <- read.table("F:/Desktop/毕业设计/outputdata/sig_118files_DEG.txt",sep = "\t")
alldiff <- read.table("F:/Desktop/毕业设计/outputdata/all_118files_DEG.txt",sep = "\t")
back=alldiff$symbol[!duplicated(alldiff$symbol)]
gene=diff$symbol[!duplicated(diff$symbol)]

IMM = read.table("IMMgene.txt",header=T,sep="\t")
CSSgene = unlist(get(load("geneset.RData")))
ICP =read.table("final_179_ICP.txt",header=F,sep="\t")
SASP=read.table("SASP.csv",header=T,sep=",",row.names = 1)


library(ggVennDiagram)
library(ggplot2)

ggVennDiagram(x=list(gene,IMM$Symbol))
# phyper(215-1, 1509, 28533-1509, 3430)  #0.9958466
fisher.test(matrix(c(215,1294,3215,0),nrow = 2,byrow = T))
p-value < 2.2e-16
ggVennDiagram(x=list(gene,ICP$V1))
# phyper(45-1, 178, 28533-178, 3430)  #0.9999992
fisher.test(matrix(c(45,133,3385,0),nrow = 2,byrow = T))
p-value < 2.2e-16


########################bulk without medium###############
 CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
 imm<- c("B","Mast", "MonoMacro","Neutrophils", "NK", "pDC","Plasma","CD4Tconv", "CD8T", "CD8Tex", "DC",  "Promonocyte","Tprolif", "Treg")
 cancer_dfs <- lapply(cancer_dfs, function(df) {
   existing_columns <- imm[imm %in% names(df)]
   select(df, all_of(existing_columns))
 })
 for (cancer in CANCER) {
   data = cancer_dfs[[cancer]]
   data$Low = rowSums(data == "Low")
   data$High = rowSums(data == "High")
   data$group = case_when(
     data$High >= (ncol(data) - 3) / 2 & data$Low == 0 ~ "TIME_High",
     data$Low >= (ncol(data) - 3) / 2 & data$High == 0 ~ "TIME_Low",
     TRUE ~ "TIME_Medium")
   cancer_dfs[[cancer]]=data
 }


CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM",
          "HNSC","KICH","KIRC","LIHC","LAML","LUAD","LUSC",
          "OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
FC = 2
p = 0.05

diff=data.frame()
for (cancer in CANCER) {
  data = cancer_dfs[[cancer]]
  bulk <- get(load(paste0("D:/TCGA/TCGA-", cancer, "_RNASeq_Counts.RData")))
  bulk <- as.data.frame(bulk[, substr(colnames(bulk), 14, 15) < 10])
  bulk = log2(bulk + 1)
  bulk = bulk[, data$group != "TIME_Medium"]
  bulk1 = bulk[rowSums(bulk) > 0, ]
  gene_means <- rowMeans(bulk1)
  gene_percentiles <- rank(gene_means) / length(gene_means)
  threshold <- 0.25  
  filtered_genes <- bulk1[gene_percentiles > threshold,]
  Group = factor(unlist(data$group), levels = c("TIME_High","TIME_Low")) 
      
  design <- model.matrix(~0+Group)
  colnames(design)=levels(Group)
  rownames(design)=colnames(bulk1)
  dge <- DGEList(counts=filtered_genes) %>% calcNormFactors()
  keep <- rowSums(cpm(dge)>0.5) >= 2
  dge<-dge[keep,,keep.lib.sizes=FALSE]
  dge<-calcNormFactors(dge)
  fit <- voom(dge, design, normalize="quantile") %>% lmFit(design)
  constrasts = paste(rev(levels(Group)), collapse = "-")
  cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()
  DEG = topTable(fit2, coef=constrasts, n=Inf) %>%  na.omit()
  DEG$symbol = rownames(DEG)
    
  myDEG=DEG[DEG$P.Value<=p & abs(DEG$logFC)>=log2(FC),]
  myDEG$sig = ifelse(myDEG$logFC >=log2(FC),"Up","Down")
  if(nrow(myDEG)>0){myDEG$cancer=cancer}
  diff=rbind(diff,myDEG)

}

#####################manha####################
library(dplyr)
library(ggplot2)
library(colorspace)

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
names(color.pals)<-CANCER
diff<-read.table("bulk_TIME_DEG.txt",sep = "\t")


boxhigh<-0.7
background.dat = diff %>%
  group_by(cancer) %>%
  summarise(
    y.localup = max(logFC)-boxhigh,
    y.localdown = min(logFC)+boxhigh,
    .groups = 'drop'
  ) %>%
  distinct(cancer, .keep_all = TRUE) %>%
  mutate(
    x.local = row_number(),
    compared.group = cancer
  )

dat.plot = diff %>%
  select(logFC, P.Value, cancer, symbol, sig) %>%
  left_join(background.dat %>% select(cancer, x.local), by = "cancer") %>%
  mutate(
    compared.group = cancer,
    fdr = P.Value,
    log2FC =ifelse(logFC > 0, logFC - boxhigh, logFC + boxhigh),
    proteins = symbol,
    significance = sig )

dat.infor = background.dat%>%mutate(y.infor = 0)


dat.marked.up <- diff %>%
  left_join(background.dat %>% select(cancer, x.local), by = "cancer") %>%
  group_by(cancer) %>%
  summarize(
    log2FC = max(logFC) - boxhigh,
    fdr = P.Value[which.max(logFC)],
    compared.group = unique(cancer),
    proteins = symbol[which.max(logFC)],
    significance = sig[which.max(logFC)],
    x.local = unique(x.local),
    .groups = 'drop'  # Optional: to ungroup after summarizing
  )%>%select(-cancer)
dat.marked.down <- diff %>%
  left_join(background.dat %>% select(cancer, x.local), by = "cancer") %>%
  group_by(cancer) %>%
  summarize(
    log2FC = min(logFC) + boxhigh,
    fdr = P.Value[which.min(logFC)],
    compared.group = unique(cancer),
    proteins = symbol[which.min(logFC)],
    significance = sig[which.min(logFC)],
    x.local = unique(x.local),
    .groups = 'drop'  # Optional: to ungroup after summarizing
  )%>%select(-cancer)

max_overlaps=5


manha=ggplot()+
  geom_col(background.dat,mapping=aes(x.local,y.localup),
           fill="grey80",alpha=0.1,width=0.9,just = 0.5)+
  geom_col(background.dat,mapping=aes(x.local,y.localdown),
           fill="grey80",alpha=0.1,width=0.9,just = 0.5)+
  geom_jitter(dat.plot,
              mapping = aes(x.local, log2FC,color = significance,
                            fill = significance),  
              size = 2, width = 0.4, alpha = 0.9, 
              shape = 21, stroke = 0) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) +
  scale_color_manual(values = c("#5390b5","#d56e5e")) +
  geom_tile(dat.infor,mapping=aes(x.local,y.infor,fill=compared.group,colour=compared.group),
            height=(1-boxhigh)*2,
            color = color.pals[dat.infor$cancer],
            fill = color.pals[dat.infor$cancer],
            alpha = 0.9,
            width=0.9)+
  guides(size=guide_legend(title="Count"))+ 
  labs(x=NULL,y=NULL)+
  geom_text(dat.infor,mapping=aes(x.local,y.infor,label=compared.group),
            size=4, angle = 90)+
  
  
  ggrepel::geom_label_repel(dat.marked.up,mapping=aes(x.local,
                                                      log2FC,
                                                      label=proteins,
                                                      color=significance),
                            force = 3,size=4, 
                            max.overlaps = max_overlaps,
                            min.segment.length = 0,
                            force_pull = 2,
                            box.padding = 0.05,
                            segment.linetype = 3, 
                            segment.color = 'black', 
                            segment.alpha = 0.5, 
                            direction = "y", 
                            vjust = 0.5)+
  ggrepel::geom_label_repel(dat.marked.down,mapping=aes(x.local,
                                                        log2FC,
                                                        label=proteins,
                                                        color=significance),
                            force = 3,size=4, 
                            max.overlaps = max_overlaps,
                            min.segment.length = 0,
                            force_pull = 2,
                            box.padding = 0.05,
                            segment.linetype = 3, 
                            segment.color = 'black', 
                            segment.alpha = 0.5, 
                            direction = "y", 
                            vjust = 0.5)+
  
  theme_classic()+
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    # legend.spacing.x=unit(0.1,'cm'),
    # legend.key.width=unit(0.5,'cm'),
    # legend.key.height=unit(0.5,'cm')
    legend.position = "none"
  )



#######################################GO#########################################
library(topGO)
library(Rgraphviz)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(DOSE)
library(ggrepel)
library(GOplot)
library(rio)
# diff<-import("sig_118files_DEG.txt")
diff<-import("bulk_TIME_DEG.txt")
GO = enrichGO(gene = diff$symbol, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",  pvalueCutoff = 1,qvalueCutoff = 1, readable = T) 
GO  = data.frame(GO)
# KEGG = enrichKEGG(gene = diff$V1,  keyType = "kegg",pAdjustMethod = 'fdr', organism= "human",pvalueCutoff=1,qvalueCutoff = 1)
# KEGG  = data.frame(KEGG)

df<-data.frame()
for (go_id in unique(GO$ID)) {
  go_genes <- GO[which(GO$ID==go_id),"geneID"]
  go_genes <-unlist(strsplit(go_genes, split = "/"))
  up_genes <- sum(go_genes %in% diff[diff$sig == "Up","symbol"])
  down_genes <- sum(go_genes %in% diff[diff$sig == "Down","symbol"])
  df<-rbind(df,c( go_id,  up_genes, down_genes))
}
colnames(df) <- c("goterm","up_regulated","down_regulated")
df$termnumber=as.numeric(df$up_regulated)+as.numeric(df$down_regulated)
df$category=GO$ONTOLOGY
df$pvalue  =GO$pvalue
df$description = GO$Description# 计算 rich.factor
parts <- strsplit(GO$GeneRatio, "/")
parts <- lapply(parts, function(x) as.numeric(unlist(x)))
df$GeneRatio  <- sapply(parts, function(x) x[1] / x[2])


######################################Go circus################################################
library(circlize)
library(RColorBrewer)
data <- read.table("Diff_GO.txt",header = TRUE, sep = "\t", quote = "\"", fill = TRUE, comment.char = "")
data=df
data$"-log10Pvalue" <- -log10(data$pvalue)
data$category <- factor(data$category, levels = unique(data$category))
data$gene_num.min <- 0
data <-data[order(data$GeneRatio, decreasing = TRUE), ][1:10, ]
# a=data[,c(1,7)]
# write.table(a,"F://Desktop//a.txt",row.names = F,quote = F,sep = "\t")

circle_size = unit(1, 'snpc')
circos.par(gap.degree = 2, start.degree = 90)

plot_data <- data[c('goterm', 'gene_num.min', 'termnumber')] 

col <- c("#FFF8E1","#CCEBC5","#B3CDE3","#FFEBEE")
#"#FBB4AE" "#B3CDE3" "#CCEBC5"

category_levels <- unique(data$category)
color_mapping <- setNames(col, category_levels)

ko_color <- color_mapping[data$category]

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1) 
circos.track(
  ylim = c(0, 1), track.height = 0.12,  
  bg.border = NA,
  bg.col = ko_color,  
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  
    # circos.axis(h = 'top', labels = FALSE) 
    circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = TRUE, font = 2)  
  }
)
# circos.clear()

plot_data <- data[c('goterm', 'gene_num.min', 'termnumber', '-log10Pvalue')]  
p_max <- round(max(data$'-log10Pvalue')) +1
colorsChoice <- colorRampPalette(c('white', '#455A64'))  
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.1, bg.border = NA, stack = TRUE,  
  panel.fun = function(region, value, ...) {
   
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
 
    circos.text(xlim ,ylim ,
                labels = CELL_META$cell.xlim[2], 
                sector.index = get.current.sector.index(),
                track.index = get.current.track.index(),
                adj = c(0.5, 0.5), font=2,
                cex = 0.8,niceFacing = T)
  } 
)

# circos.clear()

plot_data_up <- data[c('goterm', 'gene_num.min', 'up_regulated')]
names(plot_data_up) <- c('id', 'start', 'end')
plot_data_up$type <- 1  

plot_data_down <- data[c('goterm', 'up_regulated', 'termnumber')]
names(plot_data_down) <- c('id', 'start', 'end')
plot_data_down$type <- 2  

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- data[,c('up_regulated', 'down_regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col =c('#d56e5e', "#5390b5"))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.1, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
    # ylim = get.cell.meta.data('ycenter')
    # xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    # sector.name = label_data[get.cell.meta.data('sector.index'),1]
    # circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = T)
    # xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    # sector.name = label_data[get.cell.meta.data('sector.index'),2]
    # circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = T)
  } )


# circos.clear()

plot_data <- data[c('goterm', 'gene_num.min', 'termnumber', 'GeneRatio')] 
circos.genomicTrack(
  plot_data, 
  ylim = range(plot_data$GeneRatio, na.rm = TRUE),
  # ylim = c(0.035, 0.45),
  track.height = 0.1, bg.col = 'gray95', bg.border = NA,  
  panel.fun = function(region, value, ...) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  
    circos.genomicRect(region, value, col = "#E0F7FA",
                       border = NA, ytop.column = 1, ybottom = 0, ...)  
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3) 
  
    gene_ratios <- value$GeneRatio
    circos.text(xlim,ylim, 
                labels = format(gene_ratios, digits = 2), 
                col = 'black', adj = c(0.5, 0.5),
                font=2,cex = 0.8,niceFacing = T)
  } )
library(grid)
library(ComplexHeatmap)
# circos.clear()
category_legend <- Legend(
  labels = unique(data$category),
  type = 'points', pch = NA, background = col, 
  labels_gp = gpar(fontsize = 8), 
  grid_height = unit(0.4, 'cm'), 
  grid_width = unit(0.4, 'cm'))

updown_legend <- Legend(
  labels = c('Up', 'Down'), 
  type = 'points', pch = NA, background =c('#d56e5e', "#5390b5"), 
  labels_gp = gpar(fontsize = 8), 
  grid_height = unit(0.4, 'cm'),
  grid_width = unit(0.4, 'cm'))

Rich_legend <- Legend(
  labels = c('GeneRatio'), 
  type = 'points', pch = NA, background = c('#E0F7FA'), 
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(0.4, 'cm'),
  grid_width = unit(0.4, 'cm'))

pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('white', '#455A64'))(6)),
  legend_height = unit(3, 'cm'), 
  labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 10), 
  title_position = 'leftcenter-rot', 
  title = '-Log10(Pvalue)')

lgd <- packLegend(category_legend, updown_legend,
                  Rich_legend, pvalue_legend)
grid.draw(lgd)

circos.clear()
