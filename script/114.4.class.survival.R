# conda activate class && R
################################################################################
##### setting
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(dplyr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
# dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
df_Class <- readRDS(file.path(dir_result, paste0('df_Class.rds')))

df_clin <- readRDS(file.path(dir_dataset, paste0("clin_trai.rds"))) 
df_clin <- df_clin[df_clin$Sample %in% df_Class$Sample,]
df_clin <- df_clin[match(df_Class$Sample, df_clin$Sample),]

df_merge <- merge(df_clin, df_Class, by.x = 'Sample', by.y = 'Sample')
unique_cluster <- unique(as.character(df_merge$Class))
unique_cluster <- sort(unique_cluster)
ls_comparison <- combn(unique_cluster, 2, simplify = FALSE)
df_merge$Class <- factor(df_merge$Class, levels = sort(unique_cluster))

ls_color <- c("Class_1" = "#EA8379", "Class_2" = "#7DAEE0", "Class_3" = "#B395BD")
################################################################################
##### Survival analysis
sfit = survfit(Surv(OS, Survial)~Class, data=df_merge)

p <- ggsurvplot(
  sfit, 
  pval = TRUE, 
  data = df_merge, 
  legend.title = "Class",
  legend.labs = sort(unique_cluster),
  legend = c(0.85, 0.75),
  pval.size = 10,
  pval.coord = c(0.1, 0.15),
  palette = ls_color,
)


p$plot <- p$plot + theme(
  legend.title = element_text(size = 20, face = "bold"),  
  legend.text = element_text(size = 14),
  axis.title.x = element_text(size = 18), 
  axis.title.y = element_text(size = 18),   
  axis.text.x = element_text(size = 14), 
  axis.text.y = element_text(size = 14), 
)


png(file.path(dir_result, paste0("114.tran.survial.png")), width = 1700, height = 1900, res = 400, bg = "transparent"); print(p); dev.off()


################################################################################
##### Gene mutation
meta4 <- df_merge

ls_Mutate <- c('ARID1A_mutation', 'AXIN1_mutation', 'CTNNB1_mutation', 'KEAP1_mutation', 
               'KMT2C_mutation', 'TP53_mutation', 'TSC2_mutation')


calculate_percentage <- function(df, class_col, value_cols) {
  df %>%
    group_by(across(all_of(class_col))) %>%
    summarise(across(all_of(value_cols), 
                     list(percentage = ~mean(. == 1, na.rm = TRUE) * 100), 
                     .names = "{col}_percentage"),
              .groups = 'drop')
}


df_Mut <- calculate_percentage(meta4, "Class", ls_Mutate)
df_Mut <- data.frame(df_Mut)
print(df_Mut)
 

#### Bar chart
library(gg.gap)
library(ggplot2)
library(tidyr) 
library(scales) 

df_Mut2 <- df_Mut
df_Mut2[,-1] <- df_Mut2[,-1] / 100
rownames(df_Mut2) <- df_Mut2[,1]
df_Mut2 <- df_Mut2[,-1]
colnames(df_Mut2) <- gsub('_mutation_percentage', '', colnames(df_Mut2))


df_long <- pivot_longer(df_Mut2, cols = everything(), names_to = "Gene", values_to = "Value")
df_long$Class <- rep(rownames(df_Mut2), each = ncol(df_Mut2))
df_long <- df_long[order(df_long$Class), ]


p1 <- ggplot(df_long, aes(x = Gene, y = Value, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "", x = "", y = "Mutation") +
  scale_fill_manual(values = ls_color) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 0.65), labels = label_percent(accuracy = 1)) +
  theme(
    legend.position = 'none',
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.title.y = element_text(size = 21),  
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face = "bold") 
  ); p1


png(file.path(dir_result, paste0("114.mutation.Bar.png")), width = 2500, height = 2200, res = 400, bg = "white"); print(p1); dev.off()
  
