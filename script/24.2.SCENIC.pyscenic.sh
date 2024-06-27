#default value
input_loom=exp.loom
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done
 

dir_main=/share/home/shh/HCC_snRNA-seq
dir_SCENIC=${dir_main}/24.SCENIC
dir_dataset=${dir_SCENIC}/dataset
dir_result=${dir_SCENIC}/result

tfs=${dir_dataset}/hs_hgnc_tfs.txt
feather=${dir_dataset}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=${dir_dataset}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
pyscenic=/share/home/shh/miniconda3/envs/scenic/bin/pyscenic

# grn
$pyscenic grn \
--num_workers $n_workers \
--output grn.scRNA.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
$pyscenic ctx \
grn.scRNA.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
ctx.csv \
--output aucell.loom \
--num_workers $n_workers
