################################################################################
##### setting
library(openxlsx)
# library(tidyverse)
library(dplyr)
# library(tibble)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "0.data_bulk_RNA-seq/result/SEQ-DX-PXH-202402021501_result_final")
dir_dataset = file.path(dir_main, "101.bulk/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "101.bulk/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset); 

bool_save = TRUE
################################################################################
##### read data
df_exp <- read.table(file.path(dir_in, "3.Expression/all_samples_expression_raw_count.xls"), header = TRUE, sep = "\t")
df_bulk <- df_exp[,1:(ncol(df_exp)-9)] %>% data.frame()
df_bulk$gene_id <- NULL


df_tmp <- cbind(df_bulk[,56:166], df_bulk[,1:55])
df_bulk <- df_tmp


df_bulk <- df_bulk %>% group_by(gene_name) %>% summarize_all(sum)
df_bulk <- df_bulk[rowSums(df_bulk[, -1]) != 0, ] # 保留行总和不为0的行
gene_names <- df_bulk$gene_name
df_bulk$gene_name <- NULL
rownames(df_bulk) <- gene_names

class(df_bulk)
df_bulk <- data.frame(df_bulk)
saveRDS(df_bulk, file=file.path(dir_dataset, 'df_count.rds'))

