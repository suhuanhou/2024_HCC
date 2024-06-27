# METAFlux: Metabolic flux analysis
################################################################################
##### install
:'
# conda env remove --name METAFlux
conda create -n METAFlux -y -c conda-forge r-base=4.2.3 r-BiocManager r-devtools r-Seurat
conda activate METAFlux 

conda install -y -c bioconda -c conda-forge r-pheatmap
conda install -y -c bioconda -c conda-forge r-ggpubr

R
options(timeout = 3000)
# BiocManager::install(c("WGCNA", "igraph", "GeneOverlap", "ggrepel", "UCell"))
# BiocManager::install(c("ggpubr"))
install.packages("pheatmap")
devtools::install_github('KChen-lab/METAFlux')

'
################################################################################
##### METAFlux
# conda activate METAFlux 
# cd /share/home/shh/HCC_snRNA-seq/33.METAFlux/result
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/33.2.METAFlux.run.R > METAFlux.log 2>&1 &
