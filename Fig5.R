library(dplyr)
library(tidyr)
library(stringr)
library(GSVA)
library(data.table)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(IOBR)
library(gridExtra)
library(rio)
library(vegan) 
library(survival)
library(survminer)
library(pROC)
library(vcd)
library(pheatmap)


cancer_dfs<-import("E:\\cancer_dfs.RData")

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")

all_cancer_results <- data.frame(variable = character(),
                                 correlation_with_css = numeric(),
                                 p_value = numeric(),
                                 cancer = character(),
                                 stringsAsFactors = FALSE)
setwd("E:\\InPut\\TCGAbulk")

for (cancer in CANCER[1:22]) {
  # cancer="LAML"
  data <- fread(paste0("E:\\InPut\\TCGAbulk\\",cancer,".txt"), header = TRUE, sep = "\t")
  data = data.frame(t(data))
  sample = data[1,]
  gene  = rownames(data)[-1]
  data = data[-1,]
  data <- as.data.frame(lapply(data, as.numeric))
  colnames(data) = sample
  rownames(data) = gene

  # group=cancer_dfs[[cancer]] %>% filter(group!= "MICSG")
  # data1 = data[rowSums(data > 1) >= 2,group$sample] 

  exp = data[rowSums(data > 1) >= 2,]
  timer <- deconvo_tme(exp, method = "timer", group_list = rep(cancer, dim(exp)[2]))
  estimate <- deconvo_tme(exp, method = "estimate")
  mcpcounter<- deconvo_tme(exp, method = "mcpcounter")
  epic  <-   deconvo_tme(exp, method = "epic")
  quantiseq  <-   deconvo_tme(exp, method = "quantiseq",scale_mrna = T)
  CSS = gsva(gsvaParam(as.matrix(exp), CSSgene, kcdf = "Gaussian"))
  ##########
  CSS = data.frame(CSS=CSS[2, ] - CSS[1, ])
  data = cbind(CSS,
    timer %>% select(-ID),epic %>% select(-ID),
    mcpcounter %>% select(-ID),estimate %>% select(-ID),
    quantiseq %>% select(-ID))
  for (var in names(data)[-which(names(data) == "CSS")]) {
    test_result <- cor.test(data$CSS, data[[var]])
    all_cancer_results <- rbind(all_cancer_results,
                                data.frame(variable = var,
                                 correlation_with_css = test_result$estimate,
                                 p_value = test_result$p.value, cancer = cancer))
    }
  }


mydf <- all_cancer_results
data_wide <- all_cancer_results %>%
  dplyr::select(-p_value) %>%  
  pivot_wider(
    names_from = variable,        
    values_from = correlation_with_css
  )
data_wide=data.frame(data_wide)
rownames(data_wide)=data_wide$cancer
cormatrix=data_wide[,-1]
cormatrix=t(cormatrix)
get_stars <- function(p) {
  if (is.na(p)) {
    return("")  
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}
all_cancer_results$Significance <- sapply(all_cancer_results$p_value, get_stars)
sig <- all_cancer_results %>%
  dplyr::select(-c(correlation_with_css, p_value)) %>% 
  pivot_wider(
    names_from = variable,
    values_from = Significance)
sig=data.frame(sig)
rownames(sig)=sig$cancer
sig=sig[,-1]
sig= t(sig)
row_info=data.frame(cell_type = 
                      c("B_cell","T_cell","T_cell","Neutrophil_cell",
                        "MonoMacro_cell","DC_cell","B_cell","other_cell",
                        "T_cell","T_cell","stromal_cell","MonoMacro_cell",
                        "NK_cell","other_cell","T_cell","T_cell","other_cell",
                        "B_cell","NK_cell","MonoMacro_cell",
                        "MonoMacro_cell","Neutrophil_cell","stromal_cell",
                        "stromal_cell","other_cell","other_cell","other_cell","other_cell",
                        "B_cell","MonoMacro_cell","MonoMacro_cell",
                        "MonoMacro_cell", "Neutrophil_cell","NK_cell",
                        "T_cell","T_cell","T_cell","DC_cell","other_cell"
                        ),
                    method = c(rep("TIMER",6),rep("EPIC",8),rep("MCPcounter",10),
                               rep("ESTIMATE",4),rep("QUANTISEQ",11)))
rownames(row_info)=rownames(cormatrix)

methodcolor <- c( "#FF6D83", "#FFB55E", "#76C9A2", "#B39DFF", "#FFDA7A")
names(methodcolor) <- c("TIMER", "EPIC","MCPcounter","ESTIMATE","QUANTISEQ")
cellcolor<-c("#ffe0e9", "#8ca1c4", "grey90", "#d07b76", "#ff9b60", "#f47983", "#fff897","#e8ac9d")
names(cellcolor) <- c("B_cell", "T_cell","other_cell","MonoMacro_cell","NK_cell","DC_cell","stromal_cell","Neutrophil_cell")
ann_colors <- list(method = methodcolor,cell_type =  cellcolor)
heatmap_colors <- colorRampPalette(c("#5390b5","white","#d56e5e"))(50)  
p=pheatmap(cormatrix, 
           cluster_cols = F,
           cluster_rows = F,
           clustering_distance_rows = "euclidean",
           clustering_method = "complete",
           cutree_rows = 2,
           show_rownames = T,
           show_colnames = T,
           annotation_row = row_info,
           annotation_colors = ann_colors,
           scale = "none",
           legend = T,
           main = "",
           gaps_row=c(6,14,24,28,39),
           border_color = "white",
           color = heatmap_colors,
           display_numbers = sig,
           fontsize_row = 8,
           fontsize_col = 8,
           fontsize_number = 8,
           angle_col = 45)



#######################################
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
res <-data.frame()
for (cancer in CANCER) {
    data <- fread(paste0("E:\\InPut\\TCGAbulk\\",cancer,".txt"), header = TRUE, sep = "\t")
    data = data.frame(t(data))
    sample = data[1,]
    gene  = rownames(data)[-1]
    data = data[-1,]
    data <- as.data.frame(lapply(data, as.numeric))
    colnames(data) = sample
    rownames(data) = gene
    exp = data
    df = data.frame(sample = colnames(exp),cancer=cancer,
                    PDL1=t(exp[which(rownames(exp)=="CD274"),]),
                    p16=t(exp[which(rownames(exp)=="CDKN2A"),]),
                    p21=t(exp[which(rownames(exp)=="CDKN1A"),]),
                    IMMgsva =t(gsva(gsvaParam(as.matrix(exp),list(IMM$Symbol), kcdf = "Gaussian")))
                    )
    df = cbind(patient = substr(df$sample, 1, 12), df)
    cli <- read.csv(paste(cancer,"_clinical.csv",sep = ""),header = TRUE)
    cli <- cli %>%mutate(patient = bcr_patient_barcode,
                         status = ifelse(vital_status == "Alive", "0", "1"),
                         age=  trunc(-days_to_birth/365),
                         time = ifelse(vital_status == "Alive", days_to_last_followup, days_to_death)) %>%
      select(patient, time, status,age)
    cli <- cli[match(df$patient, cli$patient), ]
    df = cbind(df,cli[,-1])
    res <- rbind(res,df)
}

rio::export(res,"E:\\InPut\\TCGAcli.txt")
##############IMMgsva##############
res = rio::import("E:\\InPut\\TCGAcli.txt")
load("E:/cancer_dfs.RData")
group_columns <- lapply(cancer_dfs, function(df) df$group)
res$group <- unlist(group_columns)
res=res[res$group!= "MICSG",]
# res$days<- res$age * 365.25
# p1 = ggplot(res,aes(x = group, y = days, fill = group)) + 
#   theme(legend.position = "top") + 
#   geom_boxplot(width = 0.5, outlier.shape = NA) +
#   labs(x = "", y = "age") + 
#   facet_wrap(~cancer, nrow = 1, strip.position = "bottom") + 
#   theme_classic() + 
#   scale_fill_manual(values = c("#d56e5e", "#5390b5")) + 
#   theme( axis.text.x = element_blank(),
#          strip.text = element_blank(),
#          axis.ticks = element_blank()) + 
#   stat_compare_means(aes(label = ifelse(..p.signif.. == "ns", "", ..p.signif..)),
#                      method = "wilcox.test")
ggplot(res, aes(x = group, y = IMMgsva, fill = group)) + 
  theme(legend.position = "top") + 
  geom_boxplot(width = 0.5, outlier.shape = NA) +  
  labs(x = "", y = "IMMgsva") + 
  facet_wrap(~cancer, nrow = 1, strip.position = "top") + 
  theme_classic() + 
  scale_fill_manual(values = c("#d56e5e", "#5390b5")) + 
  theme(
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(color = "black"),
    axis.ticks = element_blank()
  ) + 
  stat_compare_means(aes(label = ifelse(..p.signif.. == "ns", "", ..p.signif..)),
                     method = "wilcox.test", 
                     label.y = -0.05)

# p1/p2
########TIME########
res = res[-res$time<=0,]
res$status = as.numeric(res$status)
res$group = as.factor(res$group)
res=na.omit(res)
######
data=na.omit(res)

table(data$group)
km<-survfit(Surv(time,status)~group, data)
ggsurvplot(km,data,conf.int = T,
           pval = T,  pval.method = T,
           risk.table = F,fdr = TRUE,
           palette = c("#d56e5e","#5390b5"),
           legend  = "top")
#####11,12,19,21
cancer=CANCER[21]
data=res[res$cancer==cancer,]

table(data$group)
km<-survfit(Surv(time,status)~group, data)
ggsurvplot(km,data,conf.int = T,
           pval = T,  pval.method = F,
             risk.table = F,fdr = TRUE,
             palette = c("#d56e5e","#5390b5"),
           legend  = "top")



# adjust.pval(method = "fdr")  # FDR 

############
all_data <-read.table("F:/Desktop/毕业设计/inputdata/118免疫含有diff的25个exp的gsva.txt",header = T,sep = "\t")
data = all_data %>% filter(quantile_CSS_Group!="Medium")%>%mutate(
  Group = factor(quantile_CSS_Group,levels = c("High","Low")),
  file=paste(cancer,cell))
ggplot(data, aes(x = file, y = IMMgsva, fill = Group)) +
  geom_boxplot(notch = FALSE) +
  theme_classic() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("High" = "#ffa69f", "Low" = "#989ade")) +  
  stat_compare_means(aes(group =Group ),method = "t.test",label = "p.signif",hide.ns =T)
