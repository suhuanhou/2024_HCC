# install infercnv
################################################################################
##### install(linux)
: '
conda create -n inferCNV
conda activate inferCNV
conda install pkg-configure jags
conda install -c conda-forge r-rjags
R
install.packages("BiocManager")
BiocManager::install("infercnv")
BiocManager::install("Seurat")
BiocManager::install("ggplot2")
BiocManager::install("future")
BiocManager::install("openxlsx")
'
################################################################################
##### run
conda activate inferCNV

# "arg1: dataset"；"arg2: result"；"arg3: references"
dataset="dataset(Endo)"
result="result(Endo)"
control="Endothelial cell"

cd /share/home/shh/HCC_snRNA-seq/12.inferCNV/${result}
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/12.3.inferCNV.run.R  "${dataset}" "${result}"  "${control}" > inferCNV.log 2>&1 &
