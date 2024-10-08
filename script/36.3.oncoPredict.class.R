# oncoPredict & CMap
################################################################################
##### setting
library(ggplot2)
library(dplyr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
# dir_in = file.path(dir_main, "114.class")
dir_dataset = file.path(dir_main, "36.oncoPredict/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "36.oncoPredict/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result);
ls_color <- c("Class_1" = "#EA8379", "Class_2" = "#7DAEE0", "Class_3" = "#B395BD")

################################################################################
##### read data 
df_onco <- read.csv(file.path(dir_dataset, "oncoPredict.class.csv"), row.names = 1, check.names = F)
# df_onco[1:4,1:4]

df_Class <- readRDS(file.path(dir_main, '114.class/result/df_Class.rds'))
# df_Class[1:4,1:2]
df_Class <- df_Class[match(rownames(df_onco), df_Class$Sample), ]


# merge
rst = merge(df_Class, df_onco, by.x = "Sample", by.y = 'row.names')
unique_cluster <- unique(as.character(rst$Class))
unique_cluster <- sort(unique_cluster)
ls_comparison <- combn(unique_cluster, 2, simplify = FALSE)


################################################################################
##### oncoPredict
rst2 <- rst[,-1]

mean_values <- rst2 %>%
  group_by(Class) %>%
  summarise(across(everything(), mean, na.rm = TRUE))


lowest_col <- mean_values %>%
  pivot_longer(cols = -Class, names_to = "Drug", values_to = "Mean") %>%
  group_by(Class) %>%
  slice_min(order_by = Mean, n = 500)

Drug_onco <- unique(lowest_col$Drug)

################################################################################
##### oncoPredict & CMap
#### class_1
CMap_1 <- read.delim(file.path(dir_dataset, "CMap/CMap_Class_1_UP.gct"), sep = "\t", header = TRUE, skip = 2)
CMap_1 <- CMap_1[-1,-1]
CMap_1$norm_cs <- as.numeric(CMap_1$norm_cs)
CMap_1 <- CMap_1[order(CMap_1$norm_cs), ]  
CMap_1[1:4,1:4]

CMap_1_Drug <- CMap_1$pert_iname[1:500]
Drug_co_1 <- intersect(CMap_1_Drug, Drug_onco)
Drug_co_1

if(T){
  choo_col <- c('Class', 'sorafenib', Drug_co_1)
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
    )+
    facet_wrap(~ Variable, scales = "free_y") +
    ylab("IC50"); p
  
  
  png(file.path(dir_result, paste0("36.Drug_co_1.png")), width = 5000, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()
  
  Drug_co_1_filter <- c('ouabain')
}



#### class_2
CMap_2 <- read.delim(file.path(dir_dataset, "CMap/CMap_Class_2_UP.gct"), sep = "\t", header = TRUE, skip = 2)
CMap_2 <- CMap_2[-1,-1]
CMap_2$norm_cs <- as.numeric(CMap_2$norm_cs)
CMap_2 <- CMap_2[order(CMap_2$norm_cs), ]  
CMap_2[1:4,1:4]

CMap_2_Drug <- CMap_2$pert_iname[1:1500]
Drug_co_2 <- intersect(CMap_2_Drug, Drug_onco)
Drug_co_2

if(T){
  choo_col <- c('Class', 'sorafenib', Drug_co_2)
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
    )+
    facet_wrap(~ Variable, scales = "free_y") +
    ylab("IC50"); p
  
  
  png(file.path(dir_result, paste0("36.Drug_co_2.png")), width = 5000, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()

  Drug_co_2_filter <- c('teniposide') 
}



#### class_3
CMap_3 <- read.delim(file.path(dir_dataset, "CMap/CMap_Class_3_UP.gct"), sep = "\t", header = TRUE, skip = 2)
CMap_3 <- CMap_3[-1,-1]
CMap_3$norm_cs <- as.numeric(CMap_3$norm_cs)
CMap_3 <- CMap_3[order(CMap_3$norm_cs), ]  
CMap_3[1:4,1:4]

CMap_3_Drug <- CMap_3$pert_iname[1:3000]
Drug_co_3 <- intersect(CMap_3_Drug, Drug_onco)
Drug_co_3

if(T){
  choo_col <- c('Class', 'sorafenib', Drug_co_3)
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
    )+
    facet_wrap(~ Variable, scales = "free_y") +
    ylab("IC50"); p
  
  
  png(file.path(dir_result, paste0("36.Drug_co_3.png")), width = 5000, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()
  
  Drug_co_3_filter <- c('TG-101348') 
}



#### all
Drug_co_all <- c(Drug_co_1_filter, Drug_co_2_filter, Drug_co_3_filter) 

choo_col <- c('Class', 'sorafenib', Drug_co_all)
data <- rst[, colnames(rst) %in% choo_col]
data <- data[choo_col]
data_long <- data %>% pivot_longer(cols = 2:ncol(data), names_to = "Variable", values_to = "Value")
data_long$Variable <- factor(data_long$Variable, levels = choo_col[-1])


p <- ggboxplot(data_long, x = "Class", y = "Value",
               color = "Class", palette = ls_color,
               add = "jitter",
               facet.by = "Variable", short.panel.labs = T) + 
  stat_compare_means(comparisons = ls_comparison, 
                     method = "wilcox.test",
                     size = 4,
                     vjust = 1.65         
  ) +
  theme(strip.text = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_blank(),             
        axis.title.x = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14), 
        legend.position = "none"
  )+
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  ylab("IC50"); p


png(file.path(dir_result, paste0("36.Drug_co_all.png")), width = 5260, height = 1140, res = 400, bg = "transparent"); print(p); dev.off()



