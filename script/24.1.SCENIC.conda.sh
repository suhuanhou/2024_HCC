# snRNA-seq
################################################################################
##### install
: '
# conda env remove --name scenic
conda create -n scenic -y python=3.9 r-base=4.2.3 r-BiocManager r-devtools r-Seurat   
conda activate scenic 

pip install scanpy loompy pyscenic
conda install -y r-SeuratDisk
conda install -c conda-forge --strict-channel-priority r-arrow
conda install -y numpy=1.22.4
pip install pandas==1.5.3 
pip install numba==0.56.4 
pip install distributed==2023.12.1
conda install -y r-ggpubr

R
BiocManager::install(c("AUCell", "RcisTarget","GENIE3"))
BiocManager::install("BiocParallel")
BiocManager::install("ComplexHeatmap")
BiocManager::install("pheatmap")
BiocManager::install("tidygraph")
BiocManager::install("ggraph")
devtools::install_github("aertslab/SCENIC")
devtools::install_github("aertslab/SCopeLoomR")

'
################################################################################
##### Resource 
# feather: https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/

## annotation
# mouse: https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
# human: https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl



################################################################################
##### setting
conda activate scenic 
R 

library(Seurat)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "24.SCENIC/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "24.SCENIC/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result);
################################################################################
##### read data 
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))


obj <- sce
Idents(obj) <-  obj$subtype


write.csv(t(as.matrix(obj@assays$RNA@counts)),file = "exp.csv", quote=F)
write.table(obj@meta.data,'metadata.xls',sep='\t',quote=F)

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
dir_dataset = dir_main/"24.SCENIC/dataset"
dir_result = dir_main/"24.SCENIC/result"
os.chdir(dir_result)


x=sc.read_csv(dir_result/"exp.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create('exp.loom',x.X.transpose(),row_attrs,col_attrs)


exit()


################################################################################
##### pyscenic
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result
nohup sh /share/home/shh/HCC_snRNA-seq/script/24.2.SCENIC.pyscenic.sh -i exp.loom -n 40 > scenic.log 2>&1 &
  

################################################################################
##### Finding cell type-specific TFs by calculating RSS values
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result  
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/24.3.SCENIC.calcRSS.R \
  -l aucell.loom \
  -m metadata.xls \
  -c subtype \
  > scenic2.log 2>&1 &
  

################################################################################
##### plot
conda activate scenic 
cd /share/home/shh/HCC_snRNA-seq/24.SCENIC/result  
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/24.4.SCENIC.plot.R \
> scenic3.log 2>&1 &

