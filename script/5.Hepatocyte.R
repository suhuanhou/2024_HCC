# Single-cell analysis pipeline for hepatocytes
################################################################################
##### setting
library(Seurat); library(ggplot2); library(tidyverse); library(harmony); library(clustree)
options(future.globals.maxSize = 10000 * 1024^2)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "4.merge/dataset")
dir_dataset = file.path(dir_main, "5.Hepatocyte/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
bool_save = TRUE 
################################################################################
##### read data
sce_all <- readRDS(file.path(dir_in, "sce_all.rds"))
sce_hep <- sce_all[, sce_all$celltype %in% c("Hepatocyte")]

##### harmony
future::plan('multisession')
sce <- NormalizeData(sce_hep) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
future::plan('multisession')
sce <- RunHarmony(sce, group.by.vars = "orig.ident")
d_dim = 30 
future::plan('multisession')
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:d_dim)


#### cluster(muti~1)
# pre-cluster at different resolutions
if(F){
  resolutions = seq(0.1, 0.5, 0.1)
  sce <- FindClusters(sce, resolution = resolutions) 
  p_cutree = clustree(sce@meta.data, prefix = "RNA_snn_res.")
  p_cutree
  dev.copy(png, file.path(dir_result, "4.cluster_tree.png"), width = 9000, height = 3500, res = 400, bg = "transparent"); dev.off()
}


resolutions = 0.2
sce <- FindClusters(sce, resolution = resolutions)


if(bool_save){saveRDS(sce, file.path(dir_dataset, "sce_1.rds"))}
# sce <- readRDS(file.path(dir_dataset, "sce_1.rds"))

##### Setting the default resolution
snn_res = "RNA_snn_res.0.2"
sce$seurat_clusters <- sce[[snn_res]]
sce$seurat_clusters <- as.integer(as.character(sce$seurat_clusters))  # 将 seurat_clusters 列的值转换为数值型
Idents(sce) <- snn_res
if(bool_save){saveRDS(sce, file.path(dir_dataset, "sce_2.rds"))}
# sce <- readRDS(file.path(dir_dataset, "sce_2.rds"))

future::plan('multisession')
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, min.dist = 0.2, n.neighbors = 500L, spread = 0.5)

if(F){
  p1 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "tissue") + labs(title = "tissue") + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p2 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "orig.ident") + labs(title = "Seq_name") + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p3 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "seurat_clusters", raster=FALSE) + ggtitle("cluster") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p1 + p2 + p3
  dev.copy(png, file.path(dir_result, "4.cluster.png"), width = 9000, height = 3500, res = 400, bg = "transparent"); dev.off()
}


# Remove clusters with low cell counts
cluster_count <- table(sce$seurat_clusters)
cluster_filter <- names(cluster_count[cluster_count > 1000])
sce <- sce[, sce$seurat_clusters %in% cluster_filter]


##### marker
future::plan('multisession')
marker <- FindAllMarkers(object = sce, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.1)
marker_top20  <- marker %>% group_by(cluster) %>% top_n(20, wt = avg_log2FC) %>% ungroup()


##### annotation
if (T){
  cluster2celltype <- c("0" = "Hep1", 
                        "1" = "Hep0", 
                        "2" = "Hep2", 
                        "3" = "Hep3", 
                        "4" = "Hep4",
                        "5" = "Hep5", 
                        "6" = "Hep6", 
                        "7" = "Hep7", 
                        "8" = "Hep8", 
                        "9" = "Hep9")  
  sce$seurat_clusters <- as.character(sce$seurat_clusters)
  sce[['subtype']] = unname(cluster2celltype[sce$seurat_clusters])
}


marker$subtype <- marker$cluster
marker$subtype = unname(cluster2celltype[marker$subtype])
saveRDS(marker, file.path(dir_dataset, "marker_Hep.rds"))
# marker <- readRDS(file.path(dir_dataset, "marker_Hep.rds"))
# write.xlsx(marker, file = file.path(dir_result, "marker.xlsx"))


Idents(sce) <- "subtype"
ls_subtype <- c("Hep0", "Hep1", "Hep2", "Hep3", "Hep4", "Hep5", "Hep6", "Hep7", "Hep8", "Hep9")
sce$subtype <- factor(sce$subtype, levels = ls_subtype)
if(bool_save){saveRDS(sce, file.path(dir_dataset, "sce_Hep.rds"))}
# sce <- readRDS(file.path(dir_dataset, "sce_Hep.rds"))

################################################################################
##### visualization
sce <- readRDS(file.path(dir_dataset, "sce_Hep.rds"))
#### UMAP
ls_color = c('#3EA94B','#3E86BD','#E51F17','#DE75DD','#EC8392',
             '#31AADF','#FDCE4F','#EC7B1A','#1505A8','#C20576')

p <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "subtype", raster = FALSE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 28), 
        legend.text = element_text(size = 24), 
        panel.border = element_rect(fill = NA, color = "black", linetype = "solid")) +
  scale_color_manual(values = ls_color); p

dev.copy(png, file.path(dir_result, "5.umap_subtype.png"), width = 4500, height = 3500, res = 400, bg = "transparent"); dev.off()




##### tissue source: Donut chart
# Calculate the ratio of cancer to paracancer
ratio_tissue <- table(sce$subtype[sce$tissue == 'cancer']) / table(sce$subtype)
df_ratio <- data.frame(subtype = names(ratio_tissue), cancer = as.vector(ratio_tissue))
df_ratio$paracancer <- 1 - df_ratio$cancer


plot_donut <- function(row){
  df_row <- df_ratio[row, c('subtype', 'cancer', 'paracancer')]
  df_donut <- data.frame(category = names(df_row)[-1], fraction = unlist(df_row)[-1])
  df_donut$ymax = cumsum(df_donut$fraction)
  df_donut$ymin = c(0, head(df_donut$ymax, n = -1))
  
  p <- ggplot(data = df_donut, aes(fill = category, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(colour = "grey30", show_guide = FALSE) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", title = "") + 
    xlim(c(0, 4)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          legend.position = "none") +
    geom_text(aes(x = 0, y = 0, label = df_ratio$subtype[row]), color = "black", size = 28)
  
  png(file.path(dir_result, paste0("5.donut_", df_ratio$subtype[row], ".png")), width = 2500, height = 2500, res = 400, bg = "transparent"); print(p); dev.off()
  
}


plot_donuts <- lapply(1:nrow(df_ratio), plot_donut)