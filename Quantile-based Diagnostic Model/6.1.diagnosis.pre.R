################################################################################
##### setting
library(pROC)
library(dplyr)
library(readxl)

rm(list = ls());
dir_main = dirname(rstudioapi::getActiveDocumentContext()$path); 
# dir_main = "/share/home/shh/HCC_snRNA-seq";
dir_in = file.path(dir_main, "dataset")

dir_ = file.path(dir_main, "6.diagnosis")
dir_dataset = file.path(dir_, "dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_, "result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result); set.seed(100)
################################################################################
##### mlr3verse
xlsx <- read_xlsx(file.path(dir_main, '0.data', 'mlr3verse/importance_Hep0.xlsx'))
gene_Hep0 = xlsx$Feature[1:30]
xlsx <- read_xlsx(file.path(dir_main, '0.data', 'mlr3verse/importance_Hep2.xlsx'))
gene_Hep2 = xlsx$Feature[1:30]
xlsx <- read_xlsx(file.path(dir_main, '0.data', 'mlr3verse/importance_Hep9.xlsx'))
gene_Hep9 = xlsx$Feature[1:30]
gene_list = unique(c(gene_Hep0, gene_Hep2, gene_Hep9))


auc_df <- data.frame(matrix(ncol = 0, nrow = length(gene_list)))
rownames(auc_df) <- gene_list


##### ROC
bioRoc = function(dataset_rds){
  # dataset_rds = "GSE17548.rds"
  # load data 
  dataset <- readRDS(file.path(dir_dataset, dataset_rds))
  dataset <- dataset[, colnames(dataset) %in% c('diagnosis', gene_list)]
  levels(dataset$diagnosis) <- c("non-HCC", "HCC")

  dataset_name <- gsub("\\.rds", "", dataset_rds)
  print(dataset_name)
  
  roc_color = rainbow(ncol(dataset)-1)
  aucText=c()
  
  # try(dev.off(), silent = TRUE)    

  for (gene_choose in gene_list){
    # gene_choose = 'MEP1A'
    # print(gene_choose)
    if (gene_choose %in% colnames(dataset)){
      roc_obj <- roc(dataset$diagnosis, dataset[,gene_choose], levels = c("non-HCC", "HCC"), smooth=F)
      auc_df[gene_choose, dataset_name] <<- sprintf("%.3f", roc_obj$auc)
    }else{
      auc_df[gene_choose, dataset_name] <<- NA
    }

  }
}


library(readxl)
rds_files <- read_xlsx(file.path(dir_, 'dataset.xlsx'))
dataset_select <- rds_files$dataset[which(rds_files$verification  == 'internal')]
rds_files <- paste0(dataset_select, '.rds')


for (rds_file in rds_files){
  bioRoc(rds_file)
}


auc_df$average <- apply(auc_df, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
auc_df$average <- ifelse(is.nan(auc_df$average), NA, auc_df$average)

# filter
na_count <- rowSums(is.na(auc_df))
auc_df_filter <- auc_df[na_count <= ncol(auc_df) * 0.50, ]


gene_diagnosis <- rownames(auc_df_filter[!is.na(auc_df_filter$average) & auc_df_filter$average > 0.7, ])
saveRDS(gene_diagnosis, 'gene_diagnosis.rds')

