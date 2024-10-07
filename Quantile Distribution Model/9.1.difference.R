# Analyze the variability of gene expression in different datasets  
################################################################################
##### setting
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)


rm(list = ls());
dir_main = dirname(rstudioapi::getActiveDocumentContext()$path); 
# dir_main = "/share/home/shh/HCC_snRNA-seq";
dir_in = file.path(dir_main, "6.diagnosis/dataset")
dir_ = file.path(dir_main, "9.difference")
dir_dataset = file.path(dir_, "dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_, "result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result); set.seed(100)
################################################################################
##### dataset
df_condition <- read_excel(file.path(dir_main, "6.diagnosis/dataset.xlsx"), sheet = "dataset")
df_tmp <- data.frame(df_condition)
df_tmp <- df_tmp[order(df_tmp$samples, decreasing = TRUE),]
df_tmp <- df_tmp[df_tmp$verification %in% 'external',]

rds_list <- df_tmp$dataset[1:8]

func_difference = function(rds_file){
  # rds_file = "GSE144269" 
  print(rds_file)
    
  dataset <- readRDS(file.path(dir_in, paste0(rds_file, '.rds')))
  dataset <- dataset[dataset$diagnosis %in% c('HCC', 'non-HCC'),]


  if (!choo_gene %in% colnames(dataset)) {return(NULL)}
  bool_log <- max(dataset[,-1])>50
  dataset <- dataset[,c('diagnosis', choo_gene)]
  if(bool_log){dataset[,2] <- log2(dataset[,2] + 1)}   
  
  dataset <- tibble::rownames_to_column(dataset, var = "sample")
  dataset$dataset <- rds_file
  
  dataset$Expression <- dataset[,3]
  df_FC <<- rbind(df_FC, dataset)
}

################################################################################
##### Differential analysis of individual genes
choo_gene <- "PDE7B"


df_FC <- data.frame()
for (rds_file in rds_list){
  func_difference(rds_file)
}

df_FC$dataset <- factor(df_FC$dataset, levels = unique(df_FC$dataset))
df_FC$diagnosis <- factor(df_FC$diagnosis, levels = c('non-HCC', 'HCC'))


library(ggunchained)  # devtools::install_github("JanCoUnchained/ggunchained")
library(ggpubr)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  

  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



Data_summary <- summarySE(df_FC, measurevar="Expression", groupvars=c("dataset","diagnosis"))
sample_counts <- table(df_FC$dataset) 


p <- ggplot(data = df_FC, aes(x = dataset, y = Expression, fill = diagnosis)) + 
  geom_split_violin(trim = FALSE, color = "white") + 
  geom_point(data = Data_summary, aes(x = dataset, y = Expression), pch = 19, position = position_dodge(0.9), size = 1.5) + 
  geom_errorbar(data = Data_summary, aes(ymin = Expression - ci, ymax = Expression + ci), 
                width = 0.1, position = position_dodge(0.9), color = "black", alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#00BDC2", "#FF6207")) +
  theme_bw() + 
  ggtitle(choo_gene) +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.title.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none",  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) + 
  ylab('Expression') + xlab("") + 
  coord_cartesian(ylim = c(0, 12)) +
  stat_compare_means(aes(group = diagnosis),
                     label = "p.signif",
                     method = "t.test",
                     label.y = 11,  
                     hide.ns = TRUE,
                     size = 10) +  
  annotate("text", x = 1:length(sample_counts), y = rep(1, length(sample_counts)), 
           label = sample_counts, vjust = 3, hjust = 0.5, size = 6);p
#                             

dev.copy(png, file.path(dir_result, paste0("split_violin_", choo_gene, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); dev.off()


################################################################################
##### Differential analysis of individual genes
choo_gene <- "GHR"


df_FC <- data.frame()
for (rds_file in rds_list){
  func_difference(rds_file)
}

df_FC$dataset <- factor(df_FC$dataset, levels = unique(df_FC$dataset))
df_FC$diagnosis <- factor(df_FC$diagnosis, levels = c('non-HCC', 'HCC'))


library(ggunchained)  # devtools::install_github("JanCoUnchained/ggunchained")
library(ggpubr)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



Data_summary <- summarySE(df_FC, measurevar="Expression", groupvars=c("dataset","diagnosis"))
sample_counts <- table(df_FC$dataset) 


p <- ggplot(data = df_FC, aes(x = dataset, y = Expression, fill = diagnosis)) + 
  geom_split_violin(trim = FALSE, color = "white") + 
  geom_point(data = Data_summary, aes(x = dataset, y = Expression), pch = 19, position = position_dodge(0.9), size = 1.5) + 
  geom_errorbar(data = Data_summary, aes(ymin = Expression - ci, ymax = Expression + ci), 
                width = 0.1, position = position_dodge(0.9), color = "black", alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#00BDC2", "#FF6207")) +
  theme_bw() + 
  ggtitle(choo_gene) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.title.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) + 
  ylab('Expression') + xlab("") + 
  coord_cartesian(ylim = c(-0.5, 15)) + 
  stat_compare_means(aes(group = diagnosis),
                     label = "p.signif",
                     method = "t.test",
                     label.y = 14,  
                     hide.ns = TRUE,
                     size = 10) + 
  annotate("text", x = 1:length(sample_counts), y = rep(1, length(sample_counts)), 
           label = sample_counts, vjust = 3.5, hjust = 0.5, size = 6);p

dev.copy(png, file.path(dir_result, paste0("split_violin_", choo_gene, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); dev.off()


################################################################################
##### Differential analysis of individual genes
choo_gene <- "SLCO1B3"


df_FC <- data.frame()
for (rds_file in rds_list){
  func_difference(rds_file)
}

df_FC$dataset <- factor(df_FC$dataset, levels = unique(df_FC$dataset))
df_FC$diagnosis <- factor(df_FC$diagnosis, levels = c('non-HCC', 'HCC'))


library(ggunchained)  # devtools::install_github("JanCoUnchained/ggunchained")
library(ggpubr)



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



Data_summary <- summarySE(df_FC, measurevar="Expression", groupvars=c("dataset","diagnosis"))
sample_counts <- table(df_FC$dataset) 


p <- ggplot(data = df_FC, aes(x = dataset, y = Expression, fill = diagnosis)) + 
  geom_split_violin(trim = FALSE, color = "white") + 
  geom_point(data = Data_summary, aes(x = dataset, y = Expression), pch = 19, position = position_dodge(0.9), size = 1.5) + 
  geom_errorbar(data = Data_summary, aes(ymin = Expression - ci, ymax = Expression + ci), 
                width = 0.1, position = position_dodge(0.9), color = "black", alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#00BDC2", "#FF6207")) +
  theme_bw() + 
  ggtitle(choo_gene) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.title.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none",  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) + 
  ylab('Expression') + xlab("") + 
  coord_cartesian(ylim = c(-1, 15)) +
  stat_compare_means(aes(group = diagnosis),
                     label = "p.signif",
                     method = "t.test",
                     label.y = 14, 
                     hide.ns = TRUE,
                     size = 10) +  
  annotate("text", x = 1:length(sample_counts), y = rep(1, length(sample_counts)), 
           label = sample_counts, vjust = 4, hjust = 0.5, size = 6);p


dev.copy(png, file.path(dir_result, paste0("split_violin_", choo_gene, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); dev.off()


################################################################################
##### Differential analysis of individual genes
choo_gene <- "GPC3"


df_FC <- data.frame()
for (rds_file in rds_list){
  func_difference(rds_file)
}

df_FC$dataset <- factor(df_FC$dataset, levels = unique(df_FC$dataset))
df_FC$diagnosis <- factor(df_FC$diagnosis, levels = c('non-HCC', 'HCC'))


library(ggunchained)  # devtools::install_github("JanCoUnchained/ggunchained")
library(ggpubr)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  

  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



Data_summary <- summarySE(df_FC, measurevar="Expression", groupvars=c("dataset","diagnosis"))
sample_counts <- table(df_FC$dataset)


p <- ggplot(data = df_FC, aes(x = dataset, y = Expression, fill = diagnosis)) + 
  geom_split_violin(trim = FALSE, color = "white") + 
  geom_point(data = Data_summary, aes(x = dataset, y = Expression), pch = 19, position = position_dodge(0.9), size = 1.5) + 
  geom_errorbar(data = Data_summary, aes(ymin = Expression - ci, ymax = Expression + ci), 
                width = 0.1, position = position_dodge(0.9), color = "black", alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#00BDC2", "#FF6207")) +
  theme_bw() + 
  ggtitle(choo_gene) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 28),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.title.y = element_text(family = "Times", size = 20, face = "plain"),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none",  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) + 
  ylab('Expression') + xlab("") + 
  coord_cartesian(ylim = c(-2, 17)) +
  stat_compare_means(aes(group = diagnosis),
                     label = "p.signif",
                     method = "t.test",
                     label.y = 15.75,  
                     hide.ns = TRUE,
                     size = 10) + 
  annotate("text", x = 1:length(sample_counts), y = rep(1, length(sample_counts)), 
           label = sample_counts, vjust = 5, hjust = 0.5, size = 6);p


dev.copy(png, file.path(dir_result, paste0("split_violin_", choo_gene, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); dev.off()

