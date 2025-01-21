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
load("cancer_dfs.RData")

CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC",
          "LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
imm<- c("B","Mast", "MonoMacro","Neutrophils", "NK", "pDC","Plasma","CD4Tconv", "CD8T", "CD8Tex", "DC",  "Promonocyte","Tprolif", "Treg")
CSSgene = get(load("geneset.RData"))
all_cancer_results <- data.frame(variable = character(),
                                 correlation_with_css = numeric(),
                                 p_value = numeric(),
                                 cancer = character(),
                                 stringsAsFactors = FALSE)
for (cancer in CANCER) {
  data = cancer_dfs[[cancer]]
  bulk <- get(load(paste0("D:/TCGA/TCGA-", cancer, "_RNASeq_Counts.RData")))
  bulk <- as.data.frame(bulk[, substr(colnames(bulk), 14, 15) < 10])
  bulk = log2(bulk + 1)
  # exp = bulk[, data$group != "TIME_Medium"]
  exp = bulk
  timer <- deconvo_tme(exp, method = "timer", group_list = rep(cancer, dim(exp)[2]))
  estimate <- deconvo_tme(exp, method = "estimate")
  mcpcounter<- deconvo_tme(exp, method = "mcpcounter")
  epic  <-   deconvo_tme(exp, method = "epic")
  quantiseq  <-   deconvo_tme(exp, method = "quantiseq",scale_mrna = T)
  CSS = gsva(gsvaParam(as.matrix(exp), CSSgene, kcdf = "Gaussian"))
 
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
                                 p_value = test_result$p.value, cancer = cancer))}}
save(all_cancer_results, file = "all_cancer_results.RData")


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
cellcolor<-c("#FCE4EC", "#C2185B", "#BBDEFB", "#e6a0af", "#E57373", "#CE93D8", "#FFEE58", "#EF9A9A")
names(cellcolor) <- c("B_cell", "T_cell","other_cell","MonoMacro_cell","NK_cell","DC_cell","stromal_cell","Neutrophil_cell")
ann_colors <- list(method = methodcolor,cell_type =  cellcolor)
heatmap_colors <- colorRampPalette(c("#5390b5","white","#d56e5e"))(50)  
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

##############IMMgsva box#########################
CANCER<-c("BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","THCA","UCEC","UVM")
IMM = read.table("IMMgene.txt",header=T,sep="\t")
imm<- c("B","Mast", "MonoMacro","Neutrophils", "NK", "pDC","Plasma","CD4Tconv", "CD8T", "CD8Tex", "DC",  "Promonocyte","Tprolif", "Treg")
CSSgene = get(load("geneset.RData"))

res <-data.frame()
for (cancer in CANCER) {
    data = cancer_dfs[[cancer]]
    bulk <- get(load(paste0("D:/TCGA/TCGA-", cancer, "_RNASeq_Counts.RData")))
    bulk <- as.data.frame(bulk[, substr(colnames(bulk), 14, 15) < 10])
    bulk = log2(bulk + 1)
    exp = bulk
    df = data.frame(sample = colnames(exp),cancer=cancer,
                    group = data$group,
                    PDL1=t(exp[which(rownames(exp)=="CD274"),]),
                    p16=t(exp[which(rownames(exp)=="CDKN2A"),]),
                    p21=t(exp[which(rownames(exp)=="CDKN1A"),]),
                    IMMgsva =t(gsva(gsvaParam(as.matrix(exp),list(IMM$Symbol), kcdf = "Gaussian")))
                    )
    df = cbind(patient = substr(df$sample, 1, 12), df)
    cli <- read.csv(paste("cancer,"_clinical.csv",sep = ""),header = TRUE)
    cli <- cli %>%mutate(patient = bcr_patient_barcode,
                         status = ifelse(vital_status == "Alive", "0", "1"),
                         age=  trunc(-days_to_birth/365),
                         time = ifelse(vital_status == "Alive", days_to_last_followup, days_to_death)) %>%
      select(patient, time, status,age)
    cli <- cli[match(df$patient, cli$patient), ]
    df = cbind(df,cli[,-1])
    res <- rbind(res,df)
}


res=res[res$group!= "TIME_Medium",]
p1 = ggplot(res,aes(x = group, y = age, fill = group)) + 
  theme(legend.position = "top") + 
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(x = "", y = "age") + 
  facet_wrap(~cancer, nrow = 1, strip.position = "bottom") + 
  theme_classic() + 
  scale_fill_manual(values = c("#d56e5e", "#5390b5")) + 
  theme( axis.text.x = element_blank(),
         strip.text = element_blank(),
         axis.ticks = element_blank()) + 
  stat_compare_means(aes(label = ifelse(..p.signif.. == "ns", "", ..p.signif..)),
                     method = "wilcox.test", label.y = 7)
p2=ggplot(res, aes(x = group, y = IMMgsva, fill = group)) + 
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

p1/p2
##############
res=res[res$group!= "TIME_Medium",]
ggplot(res[res$cancer%in%c("COAD","GBM","LIHC","THCA"),],
       aes(x = group, y = CD274, fill = group)) + 
  theme(legend.position = "top") + 
  geom_violin(width = 0.3,trim = F) +  
  labs(x = "", y = " ") + 
  facet_wrap(~cancer, nrow = 1, strip.position = "bottom") + 
  theme_classic() + 
  scale_fill_manual(values = c("#d56e5e", "#5390b5")) + 
  theme( axis.text.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.ticks = element_blank()) + 
  stat_compare_means(aes(label = ifelse(..p.signif.. == "ns", "", ..p.signif..)),
                     method = "t.test", label.y = 15)

########

res = res[-res$time<=0,]
res$status = as.numeric(res$status)
res$group = as.factor(res$group)
# res$time = trunc(res$time/365)
cancer=CANCER[7]
if(T){
  data=res[res$cancer==cancer,c(4,9,10)]
  data=na.omit(data)
  km<-survfit(Surv(time,status)~group, data)
  ggsurvplot(km,data,pval = T,conf.int = T,title=cancer,
             # risk.table = F,
             palette = c("#B14E53","#49548A","grey"),legend  = "top")
}


############

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
