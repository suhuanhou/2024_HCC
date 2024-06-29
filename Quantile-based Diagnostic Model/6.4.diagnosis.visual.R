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
##### read data
score_reduce <- readRDS('internal_reduce.rds')

################################################################################
##### Demonstration of regression coefficients (rose diagram)  
df_coef <- data.frame(coef(score_reduce))
df_coef <- tibble::rownames_to_column(df_coef, var = "gene")
df_coef <- df_coef[-1,]
colnames(df_coef)[2] <- 'coef'
df_coef$scale <-round(abs(df_coef$coef) * 1000)
df_coef$coef2 <-round(df_coef$coef,2)
df_coef <- df_coef[order(df_coef$scale, decreasing = F), ]


label_data<-df_coef
setDT(label_data)

label_data[,new_label:=paste0(gene, ': ',coef2)]
label_data[,id:=1:nrow(label_data)]
number_of_bar <- nrow(label_data)
label_data[,angle:=90 - 360 * (label_data$id-0.5) /number_of_bar]
label_data[,":="(hjust=ifelse(angle<90,1,0),angle1=ifelse(angle<90,angle+180,angle))]


label_data$distance <- c(200, 190, 210, 250, 
                         260, 270, 330, 340, 
                         170, 170, 170, 230, 
                         270, 290, 250, 370)


ggplot(data = df_coef,aes(x=reorder(gene,scale),y=scale,fill = ifelse(coef > 0, "Positive", "Negative")))+
  geom_bar(width = 0.8,stat = "identity")+
  coord_polar(theta = "x",start=0.05)+
  ylim(-100,500)+
  scale_fill_manual(values = c("#07B5EA", "#F7756B"), labels = c("Positive", "Negative")) +
  theme_minimal()+xlab(" ")+ylab(" ")+
  theme(
    legend.position = "none",
    text = element_text(color = "gray12"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  )+
    geom_text(aes(label = paste0(gene), 
                y = distance,angle = angle),
            data = label_data[1:8,], color = "black",
            vjust = 0.3, size = 7) +
    geom_text(aes(label = paste0(gene), 
                y = distance, angle = angle+180),
            data = label_data[9:16,],color = "black",
           vjust = 0.3, fontface = "bold",size = 7)



# dev.copy(png, file.path(dir_result, "diagnosis.coefficient.rose.png"), width = 4000, height = 3500, res = 400, bg = "transparent"); dev.off()


################################################################################
##### ROC for the specified gene in the specified dataset
rds_files <- read_xlsx(file.path(dir_, 'dataset.xlsx'))
dataset_select <- rds_files$dataset[which(rds_files$verification  == 'internal')]
rds_files <- paste0(dataset_select, '.rds')

gene_diagnosis <- readRDS('gene_diagnosis.rds')

#### score
dataset_score = function(dataset_rds){
  # dataset_rds = "GSE25097.rds"
  df_score <- data.frame(matrix(ncol = length(gene_diagnosis) + 1, nrow = 0))
  colnames(df_score) <- c('diagnosis', gene_diagnosis)
  
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
  

  for (col in colnames(df_score)) {
    if (!(col %in% colnames(dataset_new))) {
      dataset_new[[col]] <- 5.5
    }
  }
  
  
  dataset_new <- dataset_new[, colnames(df_score)]
  df_score <<- rbind(df_score, dataset_new)
}


dataset_score(rds_files[1])
df_roc = df_score
df_roc$diagnosis <- ifelse(df_roc$diagnosis == "HCC", 1, 0)

# gene_diagnosis
roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$ESR1))
par(lwd = 5)
plot(roc_obj, col = "green", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_1.png"), width = 3000, height = 3000, res = 400); dev.off()


roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$GHR))
par(lwd = 5)
plot(roc_obj, col = "green", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_4.png"), width = 3000, height = 3000, res = 400); dev.off()

roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$ALDOB))
par(lwd = 5)
plot(roc_obj, col = "red", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_2.png"), width = 3000, height = 3000, res = 400); dev.off()

roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$TAT))
par(lwd = 5)
plot(roc_obj, col = "red", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_5.png"), width = 3000, height = 3000, res = 400); dev.off()


roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$PPARGC1A))
par(lwd = 5)
plot(roc_obj, col = "#A75805", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_3.png"), width = 3000, height = 3000, res = 400); dev.off()

roc_obj <- roc(df_roc$diagnosis, as.numeric(df_roc$BHMT))
par(lwd = 5)
plot(roc_obj, col = "#A75805", main = "", lwd = 10, lwd.border = 20)
text(x = 0.4, y = 0.1, labels = paste("AUC = ", round(auc(roc_obj), 2)), cex = 5)

# dev.copy(png, file.path(dir_result, "ROC_Curve_filter_6.png"), width = 3000, height = 3000, res = 400); dev.off()

