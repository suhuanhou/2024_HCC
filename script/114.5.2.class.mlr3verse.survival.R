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
df_Class <- readRDS(file.path(dir_result, 'pred_vali.rds'))
df_clin <- readRDS(file.path(dir_dataset, paste0("clin_vali.rds")))
# df_clin[1:4,1:4]

ls_common <- intersect(df_clin$Sample, df_Class$Sample)
df_clin <- df_clin[df_clin$Sample %in% ls_common,]
df_clin <- df_clin[match(df_clin$Sample, df_Class$Sample),]
# table(df_clin$dataset)

df_merge <- merge(df_clin, df_Class, by.x = 'Sample', by.y = 'Sample')
df_merge$Survial <- as.numeric(df_merge$Survial)

unique_cluster <- unique(as.character(df_merge$Class))
unique_cluster <- sort(unique_cluster)
ls_comparison <- combn(unique_cluster, 2, simplify = FALSE)

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


png(file.path(dir_result, paste0("114.vali1.survial.png")), width = 1700, height = 1900, res = 400, bg = "transparent"); print(p); dev.off()
