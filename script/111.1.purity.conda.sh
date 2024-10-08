################################################################################
##### install
:'
# conda env remove --name purity
conda create -n purity -y -c conda-forge r-base=4.2.3 r-biocmanager r-devtools
conda activate purity 

dir_main=/share/home/shh/HCC_snRNA-seq
dir_=${dir_main}/111.purity
dir_dataset=${dir_}/dataset
dir_1=${dir_}/result_ESTIMATE
dir_2=${dir_}/result_PUREE
mkdir -p ${dir_dataset}
mkdir -p ${dir_1}
mkdir -p ${dir_2}


#### ESTIMATE
conda activate purity 
conda install -y -c conda-forge r-ggpubr r-ggsci r-ggstatsplot

R
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
'



