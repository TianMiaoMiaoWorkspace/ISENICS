##############################fjjimm###################################
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(limma)
library(stringr)
library(tinyarray)
library(rio)
immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono_Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
fjjcss <- read.table("E:/InPut/反卷积CSS.txt",sep = "\t",header=T)

IMMfjjcss <- fjjcss %>%
  group_by(cell, cancer) %>% mutate(
    q1 = quantile(CSS, 0.25, na.rm = TRUE),
    q3 = quantile(CSS, 0.75, na.rm = TRUE),
    dif = max(CSS) - min(CSS),
    file = paste(cancer,cell),
    quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
  )%>% filter(cell %in% immune_cell)

IMMfjjcss <-IMMfjjcss %>% filter( quantile_CSS_Group!="Medium") 
list <- list.files("E:/InPut/fjj/")
library(stringr)

filtered_list <- list %>%
  purrr::keep(~ {
    cancer <- gsub("(.*)_.+", "\\1", .x)  
    cell_type<- gsub(".*_(.+)\\..*", "\\1", .x)  
    any(IMMfjjcss$cancer == cancer & IMMfjjcss$cell == cell_type)
  })

FC = 1.2
p = 0.05
library(limma)
library(data.table)
df=data.frame()
diff=data.frame()

setwd("E:\\InPut/fjj/")
for (i in filtered_list) {
  cell_type <- gsub(".*_(.+)\\..*", "\\1", i) 
  cancer <- gsub("(.*)_.+", "\\1", i)
  data <- fread(i, header = TRUE, sep = ",")
  data = data.frame(t(data))
  sample = data[1,]
  gene  = rownames(data)[-1]
  data = data[-1,]
  data <- as.data.frame(lapply(data, as.numeric))
  colnames(data) = sample
  rownames(data) = gene
 
  data1 = data[rowSums(data > 1) >= 2, ] 
  CSS_Group <- IMMfjjcss[IMMfjjcss$cancer == cancer & IMMfjjcss$cell == cell_type, 
                         c("sample", "quantile_CSS_Group")]
  filtered_genes <- data1[, CSS_Group$sample]
  Group = factor(unlist(CSS_Group$quantile_CSS_Group), levels = c("High", "Low"))  

  design <- model.matrix(~0+Group)
  colnames(design) = levels(Group)
  rownames(design) = CSS_Group$sample
  
  fit <- lmFit(filtered_genes, design)
  contrasts = paste(rev(levels(Group)), collapse = "-")
  cont.matrix <- makeContrasts(contrasts = contrasts, levels = design)

  fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()
  DEG = topTable(fit2, coef=contrasts, n=Inf) %>% na.omit()
  DEG$symbol = rownames(DEG)
  DEG$cancer=cancer
  DEG$cell=cell_type

  myDEG = DEG[DEG$P.Value <= p & abs(DEG$logFC) >= log2(FC), ]
  myDEG$sig = ifelse(myDEG$logFC >= log2(FC), "Up", "Down")

  df=rbind(df,DEG)
  diff=rbind(diff,myDEG)
  print(paste0("Down",i,"!"))
  }

write.table(diff,"sig_113files_DEG.txt",quote = F,sep = "\t")
write.table(df,"all_113files_DEG.txt",quote = F,sep = "\t")

###################################################
diff<-read.table("sig_113files_DEG.txt",sep = "\t")
alldiff<-read.table("all_113files_DEG.txt",sep = "\t")
SASP=read.table("F:\\Desktop\\major_revision\\all_SASP.csv",header=T,sep=",")
intersect(diff$symbol,c("CDKN1A","CDKN2A","GLB1","IFNG","TNF","GZMB"))
diff[which(diff$symbol %in%c("CDKN1A","CDKN2A","GLB1","IFNG","TNF","GZMB")),]

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
p1
###############################################################################
marker = alldiff[which(alldiff$symbol %in%c("CDKN1A","CDKN2A","GLB1","IFNG","TNF","GZMB")),]
marker = marker[which(marker$cell=="CD8T"),]
library(ggplot2)
library(ggrepel) 
library(data.table)
DEG <-alldiff[alldiff$cancer=="PAAD" & alldiff$cell=="Tprolif",]
DEG$sig <- ifelse(DEG$logFC >= log2(1.2) & DEG$P.Value<= 0.05, "Up",
                  ifelse(DEG$logFC <= -log2(1.2) & DEG$P.Value <= 0.05, "Down", "NotSig"))
DEG$label <- ifelse(DEG$symbol %in% c("CDKN1A", "CDKN2A",  "IFNG", "TNF", "GZMB"), DEG$symbol, "")
intersect(DEG$symbol,c("CDKN1A", "CDKN2A",  "IFNG", "TNF", "GZMB"))
ggplot(DEG, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = sig)) +
  labs(x = "log2(FC)", y = "-log10(p)") +
  theme_classic() +
  scale_color_manual(values = c("#5390b5", "grey", "#d56e5e"))+
  theme_classic()+labs(title = "PAAD Tprolif")+
  geom_text(aes(label = label), vjust = -1, hjust = 0, size = 5)


my_diff <- diff  %>%
  group_by(cell) %>%
  slice_max(order_by = logFC, n = 2) %>%  
  bind_rows(
    diff %>%
      group_by(cell) %>%
      slice_min(order_by = logFC, n = 2)  
  ) %>%
  ungroup() %>%
  bind_rows(
    diff %>%
      filter(symbol %in% c("CDKN1A", "CDKN2A"))
  ) %>%
  distinct()  
my_diff$symbol <- reorder(my_diff$symbol, my_diff$logFC, FUN = mean)


p3=ggplot(my_diff, aes(x = symbol, y = cell, fill = logFC, size = -log10(P.Value))) +
  geom_point(shape = 21, color = "grey95") + 
  scale_size(range = c(4, 10)) +  
  # scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d7263d") +

  scale_fill_gradient2(low = "#5390b5",  mid = "white", high = "#d56e5e") + 
  labs(x = "", y = "", fill = "Log Fold Change", size = "-log10(P Value)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
p3



######################################################
gene=diff$symbol[!duplicated(diff$symbol)]

IMM = read.table("IMMgene.txt",header=T,sep="\t")
CSSgene = unlist(get(load("geneset.RData")))
ICP =read.table("final_179_ICP.txt",header=F,sep="\t")
SASP=read.table("SASP.csv",header=T,sep=",",row.names = 1)


library(ggVennDiagram)
library(ggplot2)
ggVennDiagram(x=list(gene,IMM$Symbol))
# phyper(242-1, 1267, 28533-1267, 4296)  #0.9958466
fisher.test(matrix(c(215,1294,3215,0),nrow = 2,byrow = T))
p-value < 2.2e-16
ggVennDiagram(x=list(gene,ICP$V1))
# phyper(45-1, 178, 28533-178, 3430)  #0.9999992
fisher.test(matrix(c(36,142,4502,0),nrow = 2,byrow = T))
p-value < 2.2e-16

#######################################
library(dplyr)
library(tidyr)
library(limma)
library(data.table)

immune_cell<-c("B","CD4Tconv","CD8T","CD8Tex","DC","Mast","Mono_Macro","Neutrophils","NK","Plasma","Tprolif","Treg")
df <- allCSS %>%
  group_by(cell, cancer) %>% mutate(
    q1 = quantile(CSS, 0.25, na.rm = TRUE),
    q3 = quantile(CSS, 0.75, na.rm = TRUE),
    dif = max(CSS) - min(CSS),
    file = paste(cancer,cell),
    quantile_CSS_Group = ifelse(CSS >= q3, "High",ifelse(CSS <= q1, "Low","Medium"))
  )%>% filter(cell %in% immune_cell)%>%select(c(sample,cell,cancer,quantile_CSS_Group))
cancer_dfs <- split(df, df$cancer)
cancer_dfs <- lapply(cancer_dfs, function(x) {

  matrix_data <- x %>%
    pivot_wider(names_from = cell, values_from = quantile_CSS_Group) 
    # %>%tibble::column_to_rownames("sample") 
  return(matrix_data)
})

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
for (cancer in CANCER) {
  data = cancer_dfs[[cancer]]
  print(paste(cancer,ncol(data)))
  data$Low = rowSums(data == "Low")
  data$High = rowSums(data == "High")
  data$group = case_when(
    data$Low > data$High  ~ "LICSG",
    data$High >  data$Low   ~ "HICSG",
    T ~ "MICSG"
  )
  cancer_dfs[[cancer]]=data
}
save(cancer_dfs, file = "E:\\cancer_dfs.RData")


FC = 2
p = 0.01
setwd("E:\\InPut\\TCGAbulk")
diff=data.frame()
for (cancer in CANCER) {
  group=cancer_dfs[[cancer]] %>% filter(group!= "MICSG")
  data <- fread(paste0(cancer,".txt"), header = TRUE, sep = "\t")
  data = data.frame(t(data))
  sample = data[1,]
  gene  = rownames(data)[-1]
  data = data[-1,]
  data <- as.data.frame(lapply(data, as.numeric))
  colnames(data) = sample
  rownames(data) = gene

  data1 = data[rowSums(data > 1) >= 2,group$sample] 

  Group = factor(unlist(group$group), levels = c("HICSG","LICSG")) 
  design <- model.matrix(~0+Group)
  colnames(design) = levels(Group)
  rownames(design) = group$sample
  
  fit <- lmFit(data1, design)
  contrasts = paste(rev(levels(Group)), collapse = "-")
  cont.matrix <- makeContrasts(contrasts = contrasts, levels = design)

  fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()
  DEG = topTable(fit2, coef=contrasts, n=Inf) %>% na.omit()
  DEG$symbol = rownames(DEG)
  DEG$cancer=cancer

  myDEG = DEG[DEG$P.Value <= p & abs(DEG$logFC) >= log2(FC), ]
  myDEG$sig = ifelse(myDEG$logFC >= log2(FC), "Up", "Down")
  diff=rbind(diff,myDEG)
  print(paste0("OK",cancer,"!"))

}
write.table(diff,"bulk_TIME_DEG.txt",quote = F,sep = "\t")
#########################################
library(dplyr)
library(ggplot2)
library(colorspace)
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
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

lower_quartile <- quantile(diff$logFC, 0.02)
upper_quartile <- quantile(diff$logFC, 0.98)

diff<-diff[diff$logFC <= upper_quartile &diff$logFC >= -upper_quartile,]
boxhigh<-10
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
    log2FC =ifelse(logFC > 0, logFC + boxhigh, logFC - boxhigh),
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
  # geom_col(background.dat,mapping=aes(x.local,y.localup),
  #          fill="grey80",alpha=0.1,width=0.9,just = 0.5)+
  # geom_col(background.dat,mapping=aes(x.local,y.localdown),
  #          fill="grey80",alpha=0.1,width=0.9,just = 0.5)+
  geom_jitter(dat.plot,
              mapping = aes(x.local, log2FC,color = significance,
                            fill = significance), 
              size = 2, width = 0.4, alpha = 0.9, 
              shape = 21, stroke = 0) +
  scale_fill_manual(values = c("#5390b5","#d56e5e")) +
  scale_color_manual(values = c("#5390b5","#d56e5e")) +
  geom_tile(dat.infor,mapping=aes(x.local,y.infor,fill=compared.group,colour=compared.group),
            height=boxhigh*2,###############癌症方块
            color = color.pals[dat.infor$cancer],
            fill = color.pals[dat.infor$cancer],
            alpha = 0.9,
            width=0.9)+
  guides(size=guide_legend(title="Count"))+ 
  labs(x=NULL,y=NULL)+
  geom_text(dat.infor,mapping=aes(x.local,y.infor,label=compared.group),
            size=4, angle = 90)+
  
  

  
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



################################################################################
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
diff <- read.table("sig_113files_DEG.txt",sep = "\t")
GO = enrichGO(gene = diff$symbol, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",  pvalueCutoff = 1,qvalueCutoff = 1, readable = T) 
GO  = data.frame(GO) 

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

######################################################################################
library(circlize)
library(RColorBrewer)
data <- read.table("bulk_TIME_Diff_GO.txt",
          header = TRUE, sep = "\t", quote = "\"", fill = TRUE, comment.char = "")

data$"-log10Pvalue" <- -log10(data$pvalue)
data$category <- factor(data$category, levels = unique(data$category))
data$gene_num.min <- 0
data <-data[order(data$GeneRatio, decreasing = TRUE), ][1:10, ]


circle_size = unit(1, 'snpc')
circos.par(gap.degree = 2, start.degree = 90)

plot_data <- data[c('goterm', 'gene_num.min', 'termnumber')] 

col <- c("#FFF8E1","#CCEBC5","#B3CDE3","#FFEBEE")

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
plot_data$start<-as.numeric(plot_data$start)
plot_data$end<-as.numeric(plot_data$end)
label_data <- data[,c('up_regulated', 'down_regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col =c('#d56e5e', "#5390b5"))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.1, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
  } )


# circos.clear()

plot_data <- data[c('goterm', 'gene_num.min', 'termnumber', 'GeneRatio')] 
circos.genomicTrack(
  plot_data, 
  ylim = range(plot_data$GeneRatio, na.rm = TRUE),

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
