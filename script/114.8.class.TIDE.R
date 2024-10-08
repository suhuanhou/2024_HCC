################################################################################
##### setting
library(tinyarray)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(ggpubr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
ls_color <- c("Class_1" = "#EA8379", "Class_2" = "#7DAEE0", "Class_3" = "#B395BD")
################################################################################
#####  data
if(F){
  df_ComBat <- readRDS(file.path(dir_dataset, paste0("df_trai.rds"))) 
  selected_columns <- grep("PMID31585088|TCGA", colnames(df_ComBat), value = TRUE)
  df_ComBat <- df_ComBat[, selected_columns]
  
  Expr <- t(apply(df_ComBat, 1, function(x){x-(mean(x))}))
  Expr[1:6,1:6]
  
  write.table(Expr, file = file.path(dir_dataset, paste0("114.TIDE_data.txt")), sep = "\t", row.names = TRUE, quote = FALSE)
}

# upload: http://tide.dfci.harvard.edu/

################################################################################
##### TIDE
df_Class <- readRDS(file.path(dir_result, paste0('df_Class.rds')))
rst <- read.csv(file.path(dir_dataset, "20240903 TIDE.csv"), row.names = 1, check.names = F)
# rst[1:4,1:4]

rst <- rst[match(df_Class$Sample, rownames(rst)),]
rst = merge(rst, df_Class, by.x = "row.names", by.y = 'Sample')
# table(df_TIDE$Responder, df_TIDE$Class)

unique_cluster <- unique(as.character(rst$Class))
unique_cluster <- sort(unique_cluster)
ls_comparison <- combn(unique_cluster, 2, simplify = FALSE)

################################################################################
##### visualization
choo_col <- c('Class', 'CAF', 'Exclusion', 'TIDE')
data <- rst[, colnames(rst) %in% choo_col]
data <- data[choo_col]
data_long <- data %>% pivot_longer(cols = 2:ncol(data), names_to = "Variable", values_to = "Value")


p <- ggboxplot(data_long, x = "Class", y = "Value",
               color = "Class", palette = ls_color,
               add = "jitter",
               facet.by = "Variable", short.panel.labs = T) + 
  stat_compare_means(comparisons = ls_comparison, 
                     method = "wilcox.test",
                     size = 5,
                     vjust = 1.5         
                     ) +
  theme(strip.text = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_blank(),             
        axis.title.x = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "none"
        ) +            
  facet_wrap(~ Variable, scales = "free_y"); p


png(file.path(dir_result, paste0("114.TIDE.png")), width = 3450, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()

