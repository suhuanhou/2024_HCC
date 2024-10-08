# Drug sensitivity analysis
################################################################################
##### conda
: '
# conda env remove --name oncoPredict
conda create -n oncoPredict -y r-base=4.2.3 r-BiocManager r-devtools r-Seurat   
conda activate oncoPredict 

conda install -y CMAKE
wget -P /share/home/shh/HCC_snRNA-seq/tool/pkg https://cran.r-project.org/src/contrib/Archive/ridge/ridge_2.7.tar.gz

conda install -y gsl
wget -P /share/home/shh/HCC_snRNA-seq/tool/pkg https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
cd /share/home/shh/HCC_snRNA-seq/tool/pkg
tar -zxvf gsl-2.6.tar.gz

cd gsl-2.6
./configure --prefix=/share/home/shh/HCC_snRNA-seq/tool/pkg/gsl
make
make install

R
options(timeout = 3000)
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
BiocManager::install("sva")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
devtools::install_local("/share/home/shh/HCC_snRNA-seq/tool/pkg/ridge_2.7.tar.gz")
BiocManager::install("oncoPredict")
'

################################################################################
##### run
conda activate oncoPredict 
cd /share/home/shh/HCC_snRNA-seq/36.oncoPredict/result
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/36.2.oncoPredict.run.R > oncoPredict.log 2>&1 &
 

