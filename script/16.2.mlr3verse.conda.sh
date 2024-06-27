# install mlr3verse
################################################################################
##### install(linux)
: '
# conda env remove --name mlr3verse
conda create -n mlr3verse -y r-base=4.2.3   
conda activate mlr3verse 
conda install -y -c conda-forge r-Rcpp r-kernlab
conda install -y -c conda-forge r-BiocManager 
conda install -y -c conda-forge r-dplyr r-openxlsx r-broom r-factoextra
conda install -y -c conda-forge r-ranger
R
BiocManager::install("mlr3verse")

'
################################################################################
##### run
conda activate mlr3verse
cd /share/home/shh/HCC_snRNA-seq/16.mlr3verse

Hep=Hep1
nohup Rscript '/share/home/shh/HCC_snRNA-seq/script/16.3.mlr3verse.run.R' ${Hep} > mlr3verse.${Hep}.log 2>&1 &

  
Hep=Hep2
nohup Rscript '/share/home/shh/HCC_snRNA-seq/script/16.3.mlr3verse.run.R' ${Hep} > mlr3verse.${Hep}.log 2>&1 &
  
  
Hep=Hep9
nohup Rscript '/share/home/shh/HCC_snRNA-seq/script/16.3.mlr3verse.run.R' ${Hep} > mlr3verse.${Hep}.log 2>&1 &
  
