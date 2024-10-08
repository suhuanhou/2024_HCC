# subtype & Survival analysis
################################################################################
##### install
'
# conda env remove --name survival
conda create -n survival -y -c conda-forge r-base=4.2.3 r-biocmanager r-devtools python=3.9
conda activate survival 

conda install -y r-Matrix
conda install -y r-ggplot2
conda install -y -c conda-forge r-survminer
conda install -y r-openxlsx

R
options(timeout = 3000)
'

################################################################################
##### setting
library(dplyr)
library(openxlsx)
library(survival)
library(survminer)
library(gridExtra)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "5.Hepatocyte/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/clinical"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
# clinical data
df_meta <-  read.xlsx(file.path(dir_main, "0.data_snRNA-seq/snRNA-seq clinical.xlsx"), sheet = 'clinical')
df_meta$Recurrence_day <- NA
for (j in 1:nrow(df_meta)){
  if(is.na(df_meta$Death_date[j])){
    df_meta$Survival_day[j] <- as.Date(df_meta$Follow_up[j]) -  as.Date(df_meta$Surgery[j])
  }else{
    df_meta$Survival_day[j] <- as.Date(df_meta$Death_date[j]) - as.Date(df_meta$Surgery[j])
  }
  
  if(!is.na(df_meta$Recurrence_date[j])){
    df_meta$Recurrence_day[j] <- as.Date(df_meta$Recurrence_date[j]) - as.Date(df_meta$Surgery[j])
  }
}


df_meta$time = round(df_meta$Survival_day/30,2) 
df_meta$event = df_meta$Death

df_prop <- readRDS('df_prop_cancer.rds')
df_merge <- cbind(df_prop, df_meta)
df_merge <- df_merge[!is.na(df_merge$Follow_up),]


################################################################################
##### survival analysis
# colnames(df_merge)
cellRatio = colnames(df_merge)[1:10]
splots <- lapply(cellRatio, function(g){
  if (median(df_merge[, g]) == 0) {
    return()
  }
  
  df_merge$cell=ifelse(df_merge[,g]>median(df_merge[,g]),'high','low')
  sfit1=survfit(Surv(time, event)~cell, data = df_merge)
  ggsurvplot(sfit1, pval = TRUE, data = df_merge, risk.table = FALSE, 
             title = paste0(g)) 
}) 

splots <- Filter(Negate(is.null), splots)

p <- arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 5, nrow = 2, risk.table.height = 0.4)

png("5.snRNA.survival.png", width = 5000, height = 5000, res = 400, bg = "transparent"); print(p); dev.off()

