# bulk RNA-seq
################################################################################
##### read data
conda activate scenic 
R 

library(Seurat)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "101.bulk/dataset")
dir_dataset = file.path(dir_main, "24.SCENIC/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "24.SCENIC/result_bulk"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result);
################################################################################
##### data 
df_bulk <- readRDS(file.path(dir_in, "df_count.rds"))
sce <- CreateSeuratObject(counts = df_bulk)

group_list <- ifelse(grepl("bulk_P", colnames(df_bulk)), "Paracancer",
                     ifelse(grepl("bulk_C", colnames(df_bulk)), "Cancer",
                            ifelse(grepl("bulk_N", colnames(df_bulk)), "Normal", NA)))
sce$subtype <- group_list
Idents(sce) <- sce$subtype

obj <- sce
write.csv(t(as.matrix(obj@assays$RNA$counts)),file = "exp.csv", quote=F)
write.table(obj@meta.data,'metadata.xls',sep='\t',quote=F)
saveRDS(sce, 'sce_bulk.rds')

q("no")


################################################################################
##### loom
conda activate scenic 
python


import loompy as lp;
import numpy as np;
import scanpy as sc;
import os, sys
from pathlib import Path

dir_main = Path("/share/home/shh/HCC_snRNA-seq")
# dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
# dir_dataset = dir_main/"24.SCENIC/dataset"
dir_result = dir_main/"24.SCENIC/result_bulk"
os.chdir(dir_result)


x=sc.read_csv(dir_result/"exp.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create('exp.loom',x.X.transpose(),row_attrs,col_attrs)


exit()


################################################################################
##### pyscenic
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result_bulk
nohup sh /share/home/shh/HCC_snRNA-seq/script/24.7.SCENIC.pyscenic.sh -i exp.loom -n 40 > scenic.bulk.1.log 2>&1 &
  

################################################################################
##### Finding cell type-specific TFs by calculating RSS values
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result_bulk
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/24.8.SCENIC.calcRSS.R \
  -l aucell.loom \
  -m metadata.xls \
  -c subtype \
  > scenic.bulk.2.log 2>&1 &
  
  
  
################################################################################
##### plot
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result_bulk  
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/24.9.SCENIC.plot.R \
> scenic.bulk.3.log 2>&1 &

  