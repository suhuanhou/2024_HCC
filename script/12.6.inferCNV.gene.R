# inferCNV: CNV regions and genes
conda activate inferCNV
R
################################################################################
##### setting
library(infercnv)
library(data.table) 
library(dplyr)
library(ggplot2)
library(scales) 


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_dataset = file.path(dir_main, "12.inferCNV/dataset(Endo)"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "12.inferCNV/result(Endo)"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
dir_analysis = file.path(dir_main, "12.inferCNV/analysis(Endo)"); if(!dir.exists(dir_analysis)) dir.create(dir_analysis, recursive = TRUE)
setwd(dir_analysis)

bool_test = TRUE
################################################################################
##### CNV regions
file_path <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"
rg <- fread(file.path(dir_result, file_path), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rg <- data.frame(rg)
rg$start_raw <- rg$start
rg$end_raw <- rg$end
# cell_group_name     cnv_name     state     chr     start     end

# scale
rg$start <- round(rg$start_raw / 1000)
rg$end <- round(rg$end_raw / 1000)


plot_Hep_chr <- function(Hep, chr){
  # Hep = 'Hep8'; chr = 'chr1'
  rg_gain <- rg[rg$state > 3,]
  rg_loss <- rg[rg$state < 3,]
  
  
  df_tmp_gain <- rg_gain[rg_gain$chr == chr,]
  df_tmp_gain <- df_tmp_gain[grepl(paste0(Hep, "_"), df_tmp_gain$cell_group_name), ]
  df_tmp_gain <- df_tmp_gain[,5:6]
  
  df_tmp_loss <- rg_loss[rg_loss$chr == chr,]
  df_tmp_loss <- df_tmp_loss[grepl(paste0(Hep, "_"), df_tmp_loss$cell_group_name), ]
  df_tmp_loss <- df_tmp_loss[,5:6]
  
  
  value_gain <- unlist(Map(seq, from = df_tmp_gain$start, to = df_tmp_gain$end))
  value_loss <- unlist(Map(seq, from = df_tmp_loss$start, to = df_tmp_loss$end))
  
  
  value_all <-  data.frame(values = c(value_gain, value_loss),
                           group = c(rep("gain", length(value_gain)),
                                     rep("loss", length(value_loss))))
  
  # head(value_all)
  
  
  p <- ggplot(value_all, aes(x = values, fill = group)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
    labs(title = paste0(Hep, "-", chr),
         # x = "Position (x1M)",
         # y = "Region count"
         ) +
    theme_minimal() +
    coord_flip() +
    scale_x_reverse(expand = c(0, 0), labels = comma_format(scale = 0.001)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 15)) + 
    scale_fill_manual(values = c("gain" = "#8D0909", "loss" = "#31319E")) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = 'bold', size = 20),  
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.x = element_text(size = 15),  
      axis.text.y = element_text(size = 15),  
      panel.border = element_rect(color = "black", fill = NA, size = 1), 
      legend.position = "none", 
      # panel.grid.major = element_blank(), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),  
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  png(file.path(dir_analysis, paste0("12.region.", Hep, "_", chr, ".png")), width = 1500, height = 6000, res = 600, bg = "transparent"); print(p); dev.off()
  
}

plot_Hep_chr('Hep1','chr12')
plot_Hep_chr('Hep2','chr3')
plot_Hep_chr('Hep9','chr6')
plot_Hep_chr('Hep9','chr11')



################################################################################
##### 指定Hep和chr中，对应Hep的差异基因表达情况和染色体位置（Top10）
file_path <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"
ge <- fread(file.path(dir_result, file_path), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ge <- data.frame(ge)
ge <- ge[ge$state > 3|ge$state < 3,]
# cell_group_name     gene_region_name     state     gene     chr     start     end

# muti~1
Hep = 'Hep1'; chr = 'chr12'
Hep = 'Hep2'; chr = 'chr3'
Hep = 'Hep9'; chr = 'chr6'
Hep = 'Hep9'; chr = 'chr11'


if(Hep == 'Hep1'){marker <- readRDS(file.path(dir_main, "23.pathway/dataset/marker_Hep02.rds"))}
if(Hep == 'Hep2'){marker <- readRDS(file.path(dir_main, "23.pathway/dataset/marker_Hep2.rds"))}
if(Hep == 'Hep9'){marker <- readRDS(file.path(dir_main, "23.pathway/dataset/marker_Hep9.rds"))}


marker <- marker[marker$p_val <0.05,]
marker$abs_FC <- abs(marker$avg_log2FC)
marker <- marker[order(-marker$abs_FC), ]

ge_choo <- ge[ge$chr == chr,]
ge_choo <- ge_choo[grepl(paste0(Hep, "_"), ge_choo$cell_group_name), ]
# head(ge_choo)

marker_choo <- marker[rownames(marker) %in% ge_choo$gene,][1:10,]
print(paste0(Hep, '-', chr))
marker_choo
