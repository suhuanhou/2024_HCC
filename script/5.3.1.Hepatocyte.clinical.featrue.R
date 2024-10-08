# subtype & clinical indicators
################################################################################
##### setting
library(dplyr)
library(Seurat)
library(ggpubr)
library(openxlsx)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "5.Hepatocyte/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/clinical"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
sce_para <- sce[, sce$tissue %in% c("paracancer")]
sce_cancer <- sce[, sce$tissue %in% c("cancer")]
patient <- sort(unique(sce$patient))

# clinical data
df_meta <-  read.xlsx(file.path(dir_main, "0.data_snRNA-seq/snRNA-seq clinical.xlsx"), sheet = 'clinical')
# colnames(df_meta)
choo_feature <- colnames(df_meta)[17:35]

df_p <- data.frame(matrix("", nrow = length(choo_feature), ncol = 4,
                          dimnames = list(choo_feature, c('Hep0', 'Hep2', 'Hep9', 'Hep29'))))

## correlation
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


plot_feature <- function(subtype, feature, tissue, feature_rename = NULL){
  # subtype = 'Hep0'; feature = 'A/G'; tissue = 'Paracancer'; feature_rename = 'Albumin/Globulin ratio'
  if (is.null(feature_rename)){feature_rename <- feature}
  
  aver_Hep <- median(df_merge[[subtype]])
  group <- ifelse(df_merge[[subtype]] >= aver_Hep, "High", "Low")
  df_merge2 <- df_merge %>% mutate(group = group)
  df_merge2$feature <- df_merge2[[feature]]
  
  my_comparisons = list(c("Low", "High"))
  
  p <- ggboxplot(df_merge2, x = "group", y = "feature", fill = "group", 
                 ylab = "Clinical Detection Value", xlab = paste0('', subtype, " in ", tissue, " Tissue")) +
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
  png(file.path(dir_result, paste0("5.Proportion.", tissue, ".", subtype, ".", feature_rename, ".png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()
  
}


################################################################################
##### subtype & clinical indicators: cancer tissue
df_tmp <- as.matrix(table(sce_cancer$patient, sce_cancer$subtype))
patient_totals <- rowSums(df_tmp)

df_prop <- matrix(0, nrow = length(patient), ncol = ncol(df_tmp))
rownames(df_prop) <- rownames(df_tmp)
colnames(df_prop) <- colnames(df_tmp)

# proportion of each Hep subtype
for (i in 1:nrow(df_prop)) {df_prop[i, ] <- df_tmp[i, ] / patient_totals[i]}
df_prop <- as.data.frame(df_prop)
df_prop$Hep29 <- df_prop$Hep2 + df_prop$Hep9
df_merge <- cbind(df_prop, df_meta)


test_data('Hep2')
test_data('Hep9')
test_data('Hep29')
df_p

# visualization
plot_feature(subtype = 'Hep9', feature = 'AFP', tissue = 'Cancer')
saveRDS(df_prop, 'df_prop_cancer.rds') 


#### visualization:Hep9 & AFP
if(T){
  subtype = 'Hep9'; feature = 'AFP'; tissue = 'Cancer'
  feature_rename <- feature

  aver_Hep <- median(df_merge[[subtype]])
  group <- ifelse(df_merge[[subtype]] >= aver_Hep, "Hep9-High", "Hep9-Low")
  df_tmp <- df_merge %>% mutate(group = group)
  df_tmp$AFP <- log10(df_tmp$AFP)
  df_tmp$feature <- df_tmp[[feature]]
  df_tmp$group <- factor(df_tmp$group, levels = c("Hep9-Low", "Hep9-High"))
  my_comparisons = list(c("Hep9-Low", "Hep9-High"))
  
  p <- ggboxplot(df_tmp, x = "group", y = "feature", fill = "group", 
                 # ylab = "Clinical Detection Value", xlab = paste0('', subtype, " in ", tissue, " Tissue")) +
                 ylab = "Clinical Detection Value", xlab = paste0('')) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    # stat_boxplot(geom = "errorbar", width = 0.5)+
    geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.5, size = 2, show.legend = FALSE) +
    scale_fill_manual(values = c("Hep9-Low" = "#00BFC4", "Hep9-High" = "#F8766D")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          legend.position = "none", 
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 16),  
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 10),   
    )+
    ggtitle(feature_rename) ; p
  
  png(file.path(dir_result, paste0("5.Proportion.", tissue, ".", subtype, ".", feature_rename, ".png")), width = 1500, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()
  
}

