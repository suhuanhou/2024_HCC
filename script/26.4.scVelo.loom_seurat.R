# seurat data
################################################################################
#####  setting
library(Seurat)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "26.scVelo/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "26.scVelo/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
#####  read data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
# Eutopic <- sce[,!sce$subtype == 'Hep0']
Eutopic <- sce
write.csv(colnames(Eutopic), file = "cellID_obs.csv", row.names = FALSE) 

cell_embeddings<-Embeddings(Eutopic, reduction = "umap")
write.csv(cell_embeddings, file = "cell_embeddings.csv")

clusters_obs<-Eutopic$subtype
write.csv(clusters_obs, file = "clusters_obs.csv")