################################################################################
##### install
# conda env remove --name scVelo
:'
conda create -n scVelo -y python=3.9
conda activate scVelo
pip install requests
pip install scVelo
pip install tqdm ipywidgets
pip install loompy Numba==0.59.0
pip install numpy==1.24.2
'

################################################################################
##### velocyto
conda activate scVelo
module load Analysistools/velocyto
mkdir -p /share/home/shh/HCC_snRNA-seq/26.scVelo/dataset
mkdir -p /share/home/shh/HCC_snRNA-seq/26.scVelo/result
cd /share/home/shh/HCC_snRNA-seq/26.scVelo



nohup bash -c '

rmsk_gtf=/share/home/shh/HCC_snRNA-seq/tool/UCSC/GRCh38_rmsk.gtf

# genome annotation file
ref_gtf=/share/home/shh/HCC_snRNA-seq/tool/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
  
cd /share/home/shh/HCC_snRNA-seq/1.vireo/3.snRNA_bam
for file in *cellsnp.bam; do
    sample_name=$(echo "$file" | cut -c 1-5)
    sample_folder=/share/home/shh/HCC_snRNA-seq/0.data_snRNA-seq/${sample_name}_cellranger5  

    velocyto run -b ${sample_folder}/outs/filtered_feature_bc_matrix/barcodes.tsv \
    -o /share/home/shh/HCC_snRNA-seq/26.scVelo/dataset \
    -m ${rmsk_gtf} \
    ${sample_name}_cellsnp.bam \
    ${ref_gtf}

done
' > /share/home/shh/HCC_snRNA-seq/26.scVelo/result/velocyto.log 2>&1 &
  
