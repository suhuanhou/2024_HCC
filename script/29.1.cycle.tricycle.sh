################################################################################
##### install
# conda env remove --name tricycle 
:'

conda create --name tricycle
conda create -n enrichment -y -c conda-forge r-base=4.2.3 r-BiocManager r-devtools
conda activate tricycle

conda install -y r-igraph
# devtools::install_github("hansenlab/tricycle")
conda install cairo* libxt*
conda install -c conda-forge r-cairo

R
BiocManager::install("tricycle", dependencies = TRUE)
BiocManager::install("paletteer") 
'

################################################################################
##### setting
conda activate tricycle
R


library(Seurat)  
library(tricycle)  # BiocManager::install("tricycle")
library(cowplot)
library(ggplot2)
library(gridExtra)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "29.cycle/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "29.cycle/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
sce_obj = as.SingleCellExperiment(sce)

sce_obj <- project_cycle_space(sce_obj,
                                  species = c("human"),
                                  gname.type = "SYMBOL")


sce_obj <- estimate_cycle_position(sce_obj)


#### visualization
p <- plot_emb_circle_scale(sce_obj, dimred = 3,  # seurat对象reductions的顺序
                           point.size = 1, point.alpha = 0.9) +
                           theme_bw(base_size = 14) +
                           theme(panel.grid = element_blank(),
                           axis.line = element_blank())
legend <- circle_scale_legend(text.size = 4, alpha = 0.9)



ggsave("29.tricycle_umap.png", plot = plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4)), width = 10, height = 6.35, units = "in")
