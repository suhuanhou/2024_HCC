################################################################################
##### setting
library(Seurat); library(ggplot2); library(infercnv)
library(AnnoProbe)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_ = file.path(dir_main, "12.inferCNV"); if(!dir.exists(dir_)) dir.create(dir_, recursive = TRUE)
dir_dataset = file.path(dir_, "9.dataset(Uplot)"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.inferCNV.Uplot.txt"))
setwd(dir_dataset);
options(scipen = 100)  
bool_data = F  
################################################################################
##### exp
patients <- sprintf("IDsn%02d", 1:22)


if(bool_data){
  sce_all <- readRDS(file.path(dir_main, "4.merge/dataset", "sce_all.rds"))
  sce_all <- subset(sce_all, subset = celltype %in% c('Endothelial cell'))
  sce_hep <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset", "sce_Hep.rds"))

  
  
  for (patient in patients){
    # patient = 'IDsn05'
    sce_target <- sce_hep[,sce_hep$patient %in% patient]
    
    if (ncol(sce_target)>1000){
      sce_target2 <- subset(sce_target, cells = sample(colnames(sce_target), 1000))
    }
    
    subtype_counts <- table(sce_target2$subtype) 
    single_cell_subtypes <- names(subtype_counts[subtype_counts == 1])  
    sce_target2 <- sce_target2[, !(sce_target2$subtype %in% single_cell_subtypes)] 
    
    sce <- merge(sce_target2, sce_all)
    # table(sce$subtype)
    
    ### exp
    sce_count = as.data.frame(GetAssayData(sce, slot = 'counts'))
    # sce_count[1:4, 1:4]
    
    
    groupinfo= data.frame(cellId = colnames(sce_count), cellType= sce$subtype)
    
    geneInfor=annoGene(rownames(sce_count),"SYMBOL",'human')
    geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
    geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
    
    geneInfor$chr_numeric <- as.numeric(gsub("chr", "", geneInfor$chr))  # 性染色体无法提取编号
    geneInfor <- geneInfor[!is.na(geneInfor$chr_numeric), ]  # 去除性染色体
    geneInfor <- geneInfor[order(geneInfor$chr_numeric), ]
    geneInfor$chr_numeric <- NULL
    
    
    sce_count = sce_count[rownames(sce_count) %in% geneInfor[,1],]
    sce_count = sce_count[match(geneInfor[,1], rownames(sce_count)),] 
    

    write.table(sce_count ,file = paste0(patient, '.expFile.txt'), sep = '\t',quote = F)
    write.table(groupinfo,file = paste0(patient, '.groupFiles.txt'), sep = '\t',quote = F,col.names = F,row.names = F)
    write.table(geneInfor,file = paste0(patient, '.geneFile.txt'), sep = '\t',quote = F,col.names = F,row.names = F)
  }
}
################################################################################
##### run inferCNV
for (patient in patients){
  # patient = 'IDsn05'
  
  dir_result = file.path(dir_, "9.inferCNV(Uplot)", patient); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
  if (file.exists(file.path(dir_result, "run.final.infercnv_obj"))) {next}
  
  expFile=file.path(dir_dataset, paste0(patient, '.expFile.txt'))
  groupFiles=file.path(dir_dataset, paste0(patient, '.groupFiles.txt'))
  geneFile=file.path(dir_dataset, paste0(patient, '.geneFile.txt'))
  
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c("Endothelial cell")) 
  
  
  # cutoff=1 for Smart-seq2, cutoff=0.1 for 10x Genomics
  infercnv_obj = infercnv::run(infercnv_obj, cutoff = 0.1,  out_dir = dir_result,
                               cluster_by_groups = T,
                               analysis_mode='subclusters',
                               denoise=TRUE,
                               HMM=TRUE,
                               tumor_subcluster_partition_method = "random_trees",
                               num_threads = 20
  )
  
}
