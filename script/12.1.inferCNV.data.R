# inferCNV
################################################################################
##### setting
library(Seurat); library(ggplot2); library(tidyverse)
library(AnnoProbe)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
################################################################################
##### read data
sce_hep <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset", "sce_Hep.rds"))

dir_dataset = file.path(dir_main, "12.inferCNV/dataset(Endo)"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "12.inferCNV/result(Endo)"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset);

sce_all <- readRDS(file.path(dir_main, "4.merge/dataset", "sce_all.rds"))
sce_all <- subset(sce_all, subset = celltype %in% c('Endothelial cell'))

sce <- merge(sce_hep, sce_all)
# table(sce$subtype)

################################################################################
#### sample
sample_size <- 1000 
sampled_cells <- integer()  

# Group sampling based on subtype columns
for (group in unique(sce$subtype)) {
  group_cells <- which(sce$subtype == group)
  
  # If the number of cells in the group is low, add all cells in the group to the sampling result
  if (length(group_cells) < sample_size) {
    sampled_cells <- c(sampled_cells, group_cells)
  } else {  # Otherwise, a random selection of cells from the group
    sampled_indices <- sample(group_cells, size = sample_size)
    sampled_cells <- c(sampled_cells, sampled_indices)
  }
}

sampled_cells <- sort(sampled_cells)
sce <- sce[,sampled_cells]
saveRDS(sce, file.path(dir_dataset, "sce_inferCNV.rds"))

sce_count = as.data.frame(GetAssayData(sce, slot = 'counts'))
# sce_count[1:4, 1:4]
################################################################################
##### cell type
groupinfo= data.frame(cellId = colnames(sce_count), cellType= sce$subtype)


geneInfor=annoGene(rownames(sce_count),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]

### Reference genome processing
## Altered chromosome ordering
geneInfor$chr_numeric <- as.numeric(gsub("chr", "", geneInfor$chr)) 
geneInfor <- geneInfor[!is.na(geneInfor$chr_numeric), ]  # Removal of sex chromosomes
geneInfor <- geneInfor[order(geneInfor$chr_numeric), ]
geneInfor$chr_numeric <- NULL

sce_count = sce_count[rownames(sce_count) %in% geneInfor[,1],]
sce_count = sce_count[match(geneInfor[,1], rownames(sce_count)),] 


write.table(sce_count ,file = 'expFile.txt', sep = '\t',quote = F)
write.table(groupinfo,file = 'groupFiles.txt', sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file = 'geneFile.txt', sep = '\t',quote = F,col.names = F,row.names = F)
