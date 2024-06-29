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
##### filter
gene_diagnosis <- readRDS('gene_diagnosis.rds')

df_external <- data.frame(matrix(ncol = length(gene_diagnosis) + 1, nrow = 0))
colnames(df_external) <- c('diagnosis', gene_diagnosis)


##### score
dataset_score = function(dataset_rds){
  # dataset_rds = "GSE25097.rds"
  # load data 
  dataset <- readRDS(file.path(dir_dataset, dataset_rds))
  dataset <- dataset[, colnames(dataset) %in% c('diagnosis', gene_diagnosis)]
  levels(dataset$diagnosis) <- c("non-HCC", "HCC")

  dataset_name <- gsub("\\.rds", "", dataset_rds)
  print(dataset_name)
  print(table(dataset$diagnosis))

  quantiles <- apply(dataset[,2:ncol(dataset)], 2, quantile, probs = seq(0.1, 1, 0.1), na.rm = TRUE)

  dataset_new <- dataset
  for (i in 2:ncol(dataset)) {
    for (j in 1:nrow(dataset_new)) {
      if (!is.na(dataset_new[j,i])) {
        score <- sum(quantiles[,i-1] < dataset_new[j,i])
        dataset_new[j,i] <- score + 1
      }
    }
  }


  for (col in colnames(df_external)) {
    if (!(col %in% colnames(dataset_new))) {
      dataset_new[[col]] <- 5.5
    }
  }

  
  dataset_new <- dataset_new[, colnames(df_external)]
  df_external <<- rbind(df_external, dataset_new)
}



rds_files <- read_xlsx(file.path(dir_, 'dataset.xlsx'))
dataset_select <- rds_files$dataset[which(rds_files$verification  == 'external')]
rds_files <- paste0(dataset_select, '.rds')



for (rds_file in rds_files){
  dataset_score(rds_file)
}

saveRDS(df_external, 'df_external.rds')


################################################################################
##### logistic (external)
score_reduce <- readRDS('internal_reduce.rds')


df_external$predict_prob <- predict(score_reduce, newdata = df_external, type="response")
df_external$predict <- ifelse(df_external$predict_prob > 0.5, "HCC", "non-HCC")

# ROC 
df_roc = df_external
df_roc$predict <- ifelse(df_roc$predict == "HCC", 1, 0)
df_roc$diagnosis <- ifelse(df_roc$diagnosis == "HCC", 1, 0)
roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$predict_prob))
auc <- auc(roc_obj)

cutoff <- as.numeric(coords(roc_obj, "best", ret="threshold"))  
print(cutoff)
df_roc$predict <- ifelse(df_roc$predict_prob > cutoff, "HCC", "non-HCC")


plot(roc_obj, print.auc = TRUE, col = "blue", main = "ROC Curve")
legend("bottomright", legend = paste0("AUC = ", round(auc, 2)), col = "blue", lty = 1)

saveRDS(df_roc, 'df_roc_external.rds')

################################################################################
##### AUC for all datasets
if(F){
  rds_files <- read_xlsx(file.path(dir_, 'dataset.xlsx'))
  rds_files <- rds_files[rds_files$verification %in% c('internal', 'external'),]
  dataset_select <- rds_files$dataset
  rds_all <- paste0(dataset_select, '.rds')
  
  for (rds_file in rds_all){
    dataset_rds = rds_file
    # dataset_rds = "GSE84402.rds"
    # load data 
    dataset <- readRDS(file.path(dir_dataset, dataset_rds))
    dataset <- dataset[, colnames(dataset) %in% c('diagnosis', gene_diagnosis)]
    levels(dataset$diagnosis) <- c("non-HCC", "HCC")
    
    dataset_name <- gsub("\\.rds", "", dataset_rds)
    print(dataset_name)
    
    quantiles <- apply(dataset[,2:ncol(dataset)], 2, quantile, probs = seq(0.1, 1, 0.1), na.rm = TRUE)
    
    dataset_new <- dataset
    for (i in 2:ncol(dataset)) {
      for (j in 1:nrow(dataset_new)) {
        if (!is.na(dataset_new[j,i])) {
          score <- sum(quantiles[,i-1] < dataset_new[j,i])
          dataset_new[j,i] <- score + 1
        }
      }
    }
    

    for (col in colnames(df_external)) {
      if (!(col %in% colnames(dataset_new))) {
        dataset_new[[col]] <- 5.5
      }
    }
    
    
    dataset_new <- dataset_new[, colnames(df_external)]
  
    
    external_alone_df <- dataset_new
    external_alone_df$predict_prob <- predict(score_reduce, newdata = dataset_new, type="response")
    external_alone_df$predict <- ifelse(external_alone_df$predict_prob > 0.5, "HCC", "non-HCC")
    
    df_roc = external_alone_df
    df_roc$predict <- ifelse(df_roc$predict == "HCC", 1, 0)
    df_roc$diagnosis <- ifelse(df_roc$diagnosis == "HCC", 1, 0)
    roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$predict_prob))
    auc <- auc(roc_obj)
    print(auc)
    

    cutoff <- coords(roc_obj, "best", ret="threshold")  # 最佳截断值
    print(cutoff)
  }
}


################################################################################
##### visualization
library(ggplot2)
library(dplyr)

df_roc <- readRDS('df_roc_external.rds')
gene_diagnosis_internal <- readRDS('gene_diagnosis_internal.rds')
gene_diagnosis_internal <- gene_diagnosis_internal[order(gene_diagnosis_internal)]

df_roc_muti <- df_roc
df_roc_muti$predict <- NULL
colnames(df_roc_muti)[ncol(df_roc_muti)] <- 'QDM_score'

choo_col <- c('diagnosis', 'QDM_score', gene_diagnosis_internal)
df_roc_muti <- df_roc_muti[,colnames(df_roc_muti) %in% choo_col]
df_roc_muti <- df_roc_muti[, choo_col, drop = FALSE]


column_names <- colnames(df_roc_muti)[2:ncol(df_roc_muti)]
formula_str <- paste("diagnosis ~", paste(column_names, collapse = " + "))
formula <- as.formula(formula_str)

# roc 
res <- roc(formula, data = df_roc_muti)


p <- ggroc(res, legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() +
  ggtitle("external ROC Curve")+
  theme(plot.title = element_text(hjust = 0.5,size = 25),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text=element_text(size=12,colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),  
        
  ) + labs(colour = "Gene"); p

p1 <- p + annotate("text", x=0.50, y=0.05, size = 8, label=paste("QDM_score-AUC = ", round(res$QDM_score$auc,3)))
p1 


# dev.copy(png, file.path(dir_result, "ROC_Curve_external.png"), width = 3000, height = 3000, res = 400); dev.off()
