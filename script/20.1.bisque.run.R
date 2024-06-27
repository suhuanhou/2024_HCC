# Deconvolution of cellular proportions in bulk data from snRNA-seq datasets
################################################################################
##### install
: '
# conda env remove --name bisque
conda create -n bisque -y -c conda-forge r-base=4.2.3 r-BiocManager r-devtools
conda activate bisque 

conda install -y -c conda-forge r-Rcpp r-RcppArmadillo r-RcppEigen
conda install -y -c conda-forge r-bisquerna
conda install -y -c conda-forge r-Seurat

R
options(timeout = 3000)
BiocManager::install("tidyverse")
BiocManager::install("org.Hs.eg.db") 

devtools::install_github("YuLab-SMU/GOSemSim")
devtools::install_github("YuLab-SMU/DOSE")
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("YuLab-SMU/clusterProfiler")

'
################################################################################
##### Terminal
cd /share/home/shh/HCC_snRNA-seq/20.bisque/result
conda activate bisque 
R

################################################################################
##### setting
library(RColorBrewer); library(dplyr); library(tidyverse)
library(Biobase)
library(BisqueRNA)
library(Seurat)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_dataset = file.path(dir_main, "20.bisque/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "20.bisque/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_dataset)
################################################################################
##### bulk RNA-seq（muti~1）
if (F){
  bulk_name = 'TCGA'
  df_bulk <- read.table(file.path(dir_main, '17.TCGA/LIHC/result/LIHC.count.tumor.txt'), header = T, sep = '\t')
  rownames(df_bulk) <- df_bulk[,1]
  df_bulk <- df_bulk[,-1]
  mtx_bulk <- as.matrix(df_bulk)
  bulk_eset <- Biobase::ExpressionSet(assayData = mtx_bulk)
  # patientNames(bulk_eset)
}


if (F){
  bulk_name = 'PMID31585088'
  df_exp <- read.table(file.path(dir_dataset, "PMID31585088/dataset_raw/HCC_UQ_FPKM.tsv"), sep = '\t', header = TRUE)
  df_exp$protein <- substr(df_exp$protein,1,15)
  df_exp <- df_exp[,-ncol(df_exp)]
  

  library(readxl) # BiocManager::install('readxl')
  df_sample <- read_excel(file.path(dir_dataset, "PMID31585088/dataset_raw/About the RNA and protein Identifier match.xlsx"), sheet = 1)
  df_sample <- data.frame(df_sample)
  df_sample <- df_sample[,3:4]
  colnames(df_sample) <- c("ID_Tumor","ID_Normal")
  df_sample <- df_sample[-1,]
  
  
  # Tumor or Normal
  ls_col <- colnames(df_exp)
  ls_col <- sub("^X", "", ls_col)
  
  for (i in seq_along(ls_col)) {
    if (ls_col[i] %in% df_sample$ID_Tumor) {
      ls_col[i] <- paste0("T", ls_col[i])
    } else if (ls_col[i] %in% df_sample$ID_Normal) {
      ls_col[i] <- paste0("N", ls_col[i])
    }
  }
  
  colnames(df_exp) <- ls_col
  
  df_exp <- df_exp[, grep("^p|^T|^N", colnames(df_exp), value = TRUE)]
  rownames(df_exp) <- df_exp[,1]
  df_exp$protein <- NULL
  
  
  # FPKM2TPM
  FPKM2TPM <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  
  max(df_exp) 
  df_exp <- apply(df_exp,2,FPKM2TPM)
  df_exp <- data.frame(df_exp)
  
  #### ID conversion
  library(tidyverse)  # BiocManager::install('tidyverse')
  library(clusterProfiler)  # BiocManager::install('clusterProfiler')
  library(org.Hs.eg.db)  # BiocManager::install('org.Hs.eg.db')
  
  ls_rownames <- rownames(df_exp) 
  df_exp <- df_exp %>% rownames_to_column(var = "gene")
  name <- bitr(ls_rownames,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
  

  df_exp$symbol <- df_exp$gene
  for (i in 1:nrow(df_exp)) {
    matching_row <- which(name$ENSEMBL == df_exp$gene[i])
    if (length(matching_row) > 0) {
      df_exp$symbol[i] <- name$SYMBOL[matching_row]
    }
  }
  
  
  ids <- data.frame(
    probe_id = paste0(df_exp[,1], seq_len(nrow(df_exp))),
    symbol = df_exp$symbol
  )
  
  rownames(df_exp) <- paste0(df_exp[,1], seq_len(nrow(df_exp)))
  df_exp$gene <- NULL
  df_exp$symbol <- NULL
  

  ids = ids[match(rownames(df_exp),ids$probe_id),] 
  jimmy <- function(df_exp,ids){
    tmp = by(df_exp,
             ids$symbol,
             function(x) rownames(x)[which.max(rowMeans(x))])
    probes = as.character(tmp)
    print(dim(df_exp)) 
    df_exp = df_exp[rownames(df_exp) %in% probes,]
    
    print(dim(df_exp))
    rownames(df_exp) = ids[match(rownames(df_exp),ids$probe_id),2]
    return(df_exp)
  }
  
  new_df_exp <- jimmy(df_exp, ids)
  df_exp <- new_df_exp
  

  # filter
  df_bulk <- df_exp[, grep("^T", colnames(df_exp), value = TRUE)]
  
  mtx_bulk <- as.matrix(df_bulk)
  saveRDS(mtx_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds")))

  bulk_eset <- Biobase::ExpressionSet(assayData = mtx_bulk)
  # patientNames(bulk_eset)
}


################################################################################
##### snRNA-seq data
seurat.object <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset", "sce_Hep.rds"))
sc.pheno <- data.frame(check.names = F, check.rows = F,
                             stringsAsFactors = F,
                             row.names = colnames(seurat.object),
                             SubjectName = seurat.object$orig.ident,
                             cellType = seurat.object$subtype)
sc.meta <- data.frame(labelDescription=base::c("SubjectName", "cellType"),
                            row.names=base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.data <- as.matrix(GetAssayData(seurat.object, slot = "counts"))
# sc.data[1:1000, 1:4]
sc.eset <- Biobase::ExpressionSet(assayData = sc.data, phenoData = sc.pdata)
# table(sc.eset$cellType)


##### bisque
rst <- BisqueRNA::ReferenceBasedDecomposition(bulk_eset, sc.eset, markers = NULL, use.overlap = F)
bisque <- rst$bulk.props

saveRDS(rst, file.path(dir_dataset, paste0(bulk_name, "_bisque.rst.rds")))
saveRDS(bisque, file.path(dir_result, paste0(bulk_name, "_bisque.rds")))


