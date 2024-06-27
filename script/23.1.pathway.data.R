# pathway
################################################################################
##### setting
library(dplyr)
library(Seurat)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "23.pathway/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "23.pathway/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session_1.txt"))
setwd(dir_result)
################################################################################
##### read data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
Idents(sce) <- sce$subtype


marker_Hep02 = FindMarkers(object = sce, ident.1 = 'Hep0', ident.2 = 'Hep2')
marker_Hep09 = FindMarkers(object = sce, ident.1 = 'Hep0', ident.2 = 'Hep9')
marker_Hep2 = FindMarkers(object = sce, ident.1 = 'Hep2', ident.2 = 'Hep0')
marker_Hep9 = FindMarkers(object = sce, ident.1 = 'Hep9', ident.2 = 'Hep0')

saveRDS(marker_Hep02, file.path(dir_dataset, "marker_Hep02.rds"))
saveRDS(marker_Hep09, file.path(dir_dataset, "marker_Hep09.rds"))
saveRDS(marker_Hep2, file.path(dir_dataset, "marker_Hep2.rds"))
saveRDS(marker_Hep9, file.path(dir_dataset, "marker_Hep9.rds"))
