################################################################################
##### install
: '
# conda env remove --name class
conda create -n class -y   
conda activate class 
conda install -y -c conda-forge r-Rcpp r-kernlab
conda install -y -c conda-forge r-BiocManager 
conda install -y -c conda-forge r-dplyr r-openxlsx r-broom r-factoextra
conda install -y -c conda-forge r-ranger
conda install -y -c conda-forge r-ggplot2


R
options(timeout = 3000)
BiocManager::install("mlr3verse")
BiocManager::install("NMF")
BiocManager::install(c("survival", "survminer", "gridExtra"))
BiocManager::install("GEOquery")
BiocManager::install("ConsensusClusterPlus")
install.packages("fmsb")
install.packages("gg.gap")
'

################################################################################
##### run
conda activate class
cd /share/home/shh/HCC_snRNA-seq/114.class/result

Rscript '/share/home/shh/HCC_snRNA-seq/script/114.3.class.class.R'
Rscript '/share/home/shh/HCC_snRNA-seq/script/114.4.class.survival.R'
Rscript '/share/home/shh/HCC_snRNA-seq/script/114.5.1.class.mlr3verse.train.R'
Rscript '/share/home/shh/HCC_snRNA-seq/script/114.5.2.class.mlr3verse.survival.R' 
