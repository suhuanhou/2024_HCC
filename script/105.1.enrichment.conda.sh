################################################################################
##### install
:'
# conda env remove --name enrichment
conda create -n enrichment -y -c conda-forge r-base=4.2.3 r-BiocManager r-devtools
conda activate enrichment 

conda install -y -c bioconda -c conda-forge r-Seurat
# conda install -y -c bioconda -c conda-forge r-HDO.db

# dir_main=/share/home/shh/HCC_snRNA-seq
# dir_dataset=${dir_main}/110.immunedeconv/dataset
# dir_result=${dir_main}/110.immunedeconv/result
# mkdir -p ${dir_dataset}
# mkdir -p ${dir_result}

 
R
options(timeout = 300)
BiocManager::install(c("readxl", "ggplot2", "RColorBrewer"))
devtools::install_github("YuLab-SMU/GOSemSim")
devtools::install_github("YuLab-SMU/DOSE")
BiocManager::install("org.Hs.eg.db") 
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("YuLab-SMU/clusterProfiler")
BiocManager::install("MAST")

install.packages ("GOplot")
install.packages ("ggridges")

# BiocManager::install(c("WGCNA", "igraph", "GeneOverlap", "ggrepel", "UCell"))
# devtools::install_github("NightingaleHealth/ggforestplot")
# devtools::install_github("smorabit/enrichment", ref="dev")

'
################################################################################
##### run
# conda activate enrichment 
# cd /share/home/shh/HCC_snRNA-seq/34.enrichment/result
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/33.2.enrichment.run.R > enrichment.log 2>&1 &
 
