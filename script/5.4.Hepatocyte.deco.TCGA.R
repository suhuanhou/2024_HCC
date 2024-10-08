# subtype & Immune infiltration of TCGA
################################################################################
##### setting
library(dplyr)
library(ggplot2)
library(ggpubr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
# dir_dataset = file.path(dir_main, "5.Hepatocyte/deco"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/deco"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
#### Immune infiltration of TCGA
df_immu <- read.delim(file.path(dir_result, "20240717 TCGA.LIHC.cibersort.txt"))
choo_feature <- colnames(df_immu[,2:23])

df_hep <- readRDS(file.path(dir_main, "20.bisque/result/TCGA_bisque.rds"))
df_hep <- t(df_hep) %>% data.frame()
df_hep$Sample <- rownames(df_hep)

df_immu$Mixture <- gsub("-", ".", df_immu$Mixture)
df_merge <- merge(df_hep, df_immu, by.x = "Sample", by.y = "Mixture")


df_p <- data.frame(matrix("", nrow = length(choo_feature), ncol = 10,
                          dimnames = list(choo_feature, colnames(df_hep[,1:10]))))


test_data <- function(subtype){
  # subtype = 'Hep9'
  aver_Hep <- median(df_merge[[subtype]])
  group <- ifelse(df_merge[[subtype]] >= aver_Hep, "High", "Low")
  df_merge2 <- df_merge %>% mutate(group = group)
  
  group_high <- df_merge2[df_merge2$group == "High", ]
  group_low <- df_merge2[df_merge2$group == "Low", ]
  

  results <- lapply(choo_feature, function(feature) {
    test_method <- wilcox.test(group_high[[feature]], group_low[[feature]])
    p_value <- test_method$p.value
    df_p[feature, subtype] <<- round(p_value, 5)
  })
}


for (col in colnames(df_hep[,1:10])) {test_data(col)}

plot_feature <- function(subtype, feature, tissue, feature_rename = NULL){
  # subtype = 'Hep0'; feature = 'A/G'; tissue = 'Paracancer'; feature_rename = 'Albumin/Globulin ratio'
  if (is.null(feature_rename)){feature_rename <- feature}
  
  aver_Hep <- median(df_merge[[subtype]])
  group <- ifelse(df_merge[[subtype]] >= aver_Hep, "High", "Low")
  df_merge2 <- df_merge %>% mutate(group = group)
  df_merge2$feature <- df_merge2[[feature]]
  
  my_comparisons = list(c("Low", "High"))
  
  p <- ggboxplot(df_merge2, x = "group", y = "feature", fill = "group", 
                 ylab = "Proportion of immune cells", xlab = paste0('', subtype, " in ", tissue, "")) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.5, size = 2, show.legend = FALSE) +
    scale_fill_manual(values = c("Low" = "#00BFC4", "High" = "#F8766D")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          legend.position = "none", 
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 15),  
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 10),   
    )+
    ggtitle(feature_rename) ; p
  
  feature_rename <- gsub("/", " ", feature_rename)
  png(file.path(dir_result, paste0("5.", tissue, ".", subtype, ".", feature_rename, ".png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()
}


plot_feature(subtype = 'Hep0', feature = 'T.cells.CD4.memory.resting', tissue = 'TCGA', feature_rename = 'CD4.memory.resting')
plot_feature(subtype = 'Hep0', feature = 'T.cells.regulatory..Tregs.', tissue = 'TCGA', feature_rename = 'Tregs')
plot_feature(subtype = 'Hep0', feature = 'NK.cells.resting', tissue = 'TCGA', feature_rename = 'NK.resting')
plot_feature(subtype = 'Hep0', feature = 'Monocytes', tissue = 'TCGA', feature_rename = 'Monocyte')

plot_feature(subtype = 'Hep2', feature = 'T.cells.CD4.memory.resting', tissue = 'TCGA', feature_rename = 'CD4.memory.resting')
plot_feature(subtype = 'Hep2', feature = 'T.cells.gamma.delta', tissue = 'TCGA', feature_rename = 'T.gamma.delta')

plot_feature(subtype = 'Hep9', feature = 'T.cells.regulatory..Tregs.', tissue = 'TCGA', feature_rename = 'Tregs')
plot_feature(subtype = 'Hep9', feature = 'Plasma.cells', tissue = 'TCGA', feature_rename = 'Plasma')

