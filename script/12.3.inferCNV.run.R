args <- commandArgs(trailingOnly = TRUE)
print(args)
# args <- c('dataset(Endo)', 'result(Endo)', 'Endothelial cell')
################################################################################
##### setting
library(Seurat); library(ggplot2); library(infercnv);

dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "4.merge/dataset")
dir_dataset = file.path(dir_main, "12.inferCNV", args[1]); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "12.inferCNV", args[2]); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset);

options(scipen = 100) 
################################################################################
expFile=file.path(dir_dataset, 'expFile.txt')
groupFiles=file.path(dir_dataset, 'groupFiles.txt')
geneFile=file.path(dir_dataset, 'geneFile.txt')


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c(args[3]))  


# cutoff=1 for Smart-seq2, cutoff=0.1 for 10x Genomics
infercnv_obj = infercnv::run(infercnv_obj, cutoff = 0.1,  out_dir = dir_result, 
                             cluster_by_groups = TRUE,  
                             write_expr_matrix=T,  # infercnv.observations.txt
                             denoise = TRUE, HMM=T, output_format = 'pdf'
)

