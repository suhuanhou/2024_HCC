# Obtaining hepatocyte subtype-specific gene features with mlr3verse
################################################################################
#####  setting
library(Seurat); library(dplyr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; setwd(dir_main); set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "16.mlr3verse/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "16.mlr3verse/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
################################################################################
##### read data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))

subset_size <- 10000
num_subsets <- ceiling(ncol(sce@assays$RNA@data) / subset_size)  
data_list <- vector("list", length = num_subsets) 


for (i in 1:num_subsets) {
  start_col <- (i - 1) * subset_size + 1
  end_col <- min(i * subset_size, ncol(sce@assays$RNA@data))
  data_list[[i]] <- as.matrix(sce@assays$RNA@data[, start_col:end_col])
}

data_mtx <- do.call(cbind, data_list)
# View(data[1:20, 1:20])


#### cell state 
state <- sce$subtype %>% data.frame()  
colnames(state) <- 'state'
# table(state$state)

state <- state[, "state", drop = FALSE]  
dataset <- merge(t(data_mtx), state, by = "row.names")
saveRDS(dataset, file.path(dir_dataset, "dataset.rds"))
# dataset <- readRDS(file.path(dir_dataset, "dataset.rds"))
