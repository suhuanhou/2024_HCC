# merge snRNA-seq data 
################################################################################
##### setting
library(Seurat); library(ggplot2); library(tidyverse); library(harmony); library(clustree)
library(openxlsx)
options(future.globals.maxSize = 10000 * 1024^2)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "3.filter/")
dir_dataset = file.path(dir_main, "4.merge/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "4.merge/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
bool_save = TRUE
################################################################################
##### load data
sce <- readRDS('sce.rds')

##### harmony
future::plan('multisession')
sce <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F)
future::plan('multisession')
sce <- RunHarmony(sce, group.by.vars = "orig.ident")
future::plan('multisession')
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)


#### cluster(muti~1)
# pre-cluster at different resolutions
if(F){
  resolutions = seq(0.1, 0.5, 0.1)
  sce <- FindClusters(sce, resolution = resolutions)
  p_cutree = clustree(sce@meta.data, prefix = "RNA_snn_res.")
  p_cutree
  dev.copy(png, file.path(dir_result, "4.cluster_tree.png"), width = 9000, height = 3500, res = 400, bg = "transparent"); dev.off()
}


resolutions = 0.1
future::plan('multisession')
sce <- FindClusters(sce, resolution = resolutions) 


snn_res = "RNA_snn_res.0.1"
sce$seurat_clusters <- sce[[snn_res]]
sce$seurat_clusters <- as.integer(as.character(sce$seurat_clusters))  # 将 seurat_clusters 列的值转换为数值型
Idents(sce) <- snn_res
future::plan('multisession')
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, min.dist = 0.3, n.neighbors = 100L, spread = 0.4)


if(F){
  p1 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "tissue") + labs(title = "tissue") + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p2 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "orig.ident") + labs(title = "Seq_name") + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p3 <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "seurat_clusters", raster=FALSE) + ggtitle("cluster") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))
  p1 + p2 + p3
  dev.copy(png, file.path(dir_result, "4.cluster.png"), width = 9000, height = 3500, res = 400, bg = "transparent"); dev.off()
}


if(bool_save){saveRDS(sce, file.path(dir_dataset, "sce_1.rds"))}
# sce <- readRDS(file.path(dir_dataset, "sce_1.rds"))


# Remove clusters with low cell counts
cluster_count <- table(sce$seurat_clusters)
cluster_filter <- names(cluster_count[cluster_count > 1000])
sce <- sce[, sce$seurat_clusters %in% cluster_filter]


##### marker
future::plan('multisession')
marker <- FindAllMarkers(object = sce, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.1)

##### sctype:automatic annotation
if (T){
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
  source(file.path(dir_main, "tool/sctype/gene_sets_prepare.R")); source(file.path(dir_main, "tool/sctype/sctype_score_.R"))
  gs_list = gene_sets_prepare(file.path(dir_main,  "tool/sctype/ScTypeDB_short.xlsx", "Liver")) # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
  es.max = sctype_score(scRNAseqData = sce[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  cL_resutls = do.call("rbind", lapply(unique(sce@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sce@meta.data[sce@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sce@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  View(sctype_scores[,1:3])
}


##### Manual annotation
if (T){
  cluster2celltype <- c("0" = "Hepatocyte", 
                        "1" = "Hepatocyte", 
                        "2" = "Hepatocyte", 
                        "3" = "Hepatocyte", 
                        "4" = "Hepatocyte",
                        "5" = "Kupffer cell", 
                        "6" = "Liver immune cell", 
                        "7" = "Endothelial cell", 
                        "8" = "Hepatic stellate cell",
                        "9" = "Cholangiocyte")  
  sce$seurat_clusters <- as.character(sce$seurat_clusters)
  sce[['subtype']] = unname(cluster2celltype[sce$seurat_clusters])
  sce$subtype <- sce$celltype
}


##### marker2
Idents(sce) <- 'celltype'
future::plan('multisession')
marker <- FindAllMarkers(object = sce, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.1)

marker$celltype <- marker$cluster
saveRDS(marker, file.path(dir_dataset, "marker_all.rds"))
# marker <- readRDS(file.path(dir_dataset, "marker_all.rds"))

marker_top20 <- marker %>% group_by(cluster) %>% top_n(20, wt = avg_log2FC)
# write.xlsx(marker_top20, file = file.path(dir_result, "4.marker_top20.xlsx"))


Idents(sce) <- "celltype"
ls_celltype <- c("Hepatocyte", "Kupffer cell", "Liver immune cell", "Endothelial cell", "Hepatic stellate cell", "Cholangiocyte")
sce$celltype <- factor(sce$celltype, levels = ls_celltype)
if(bool_save){saveRDS(sce, file.path(dir_dataset, "sce_all.rds"))}



################################################################################
# visualization
##### umap
sce <- readRDS(file.path(dir_dataset, "sce_all.rds"))

#### celltype
ls_color = c('#00BEC3','#629CFD','#4BCE75','#F7756C','#F566E1','#B69E0E')

p <- DimPlot(sce, reduction = "umap", label = F, label.size = 10, group.by = "celltype", raster = FALSE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),  
        axis.title = element_blank(), 
        legend.text = element_text(size = 24),
        legend.position = "none",        
        panel.grid = element_blank(), 
        panel.border = element_blank(),  
        axis.line = element_blank(),  
        axis.ticks = element_blank(),  
        ) +
  scale_color_manual(values = ls_color); p


png(file.path(dir_result, "4.umap_celltype.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()



#### tissue + patient
# tissue
ls_color = c("#F7756C", "#00BEC3")
p <- DimPlot(sce, group.by = "tissue", reduction = "umap", label = FALSE, raster = FALSE) + 
  ggtitle("tissue") + 
  theme(plot.title = element_text(hjust = 0.5, size = 40),  
        axis.text = element_blank(),  
        axis.title = element_blank(), 
        legend.text = element_text(size = 24), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),  
        axis.line = element_blank(),  
        axis.ticks = element_blank(),  
        legend.position = "none"  
  ) +
  scale_color_manual(values = ls_color); p

png(file.path(dir_result, "4.umap_tissue.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()

p1 <- p + theme(legend.position = "right" )
p_legend <- cowplot::get_legend(p1)
png(file.path(dir_result, "4.umap_tissue_legend.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); grid::grid.newpage(); grid::grid.draw(p_legend); dev.off()



# patient
palette_color <- paletteer::paletteer_d("ggsci::default_ucscgb")  # palette_color <- palettes_d_names
ls_color <- as.character(palette_color[2:23])
p <- DimPlot(sce, group.by = "patient", reduction = "umap", label = F, raster = FALSE) + 
  ggtitle("patient") + 
  theme(plot.title = element_text(hjust = 0.5, size = 40),  
        axis.text = element_blank(),  
        axis.title = element_blank(), 
        legend.text = element_text(size = 24), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(),  
        legend.position = "none"  
  ) +
  scale_color_manual(values = ls_color); p

png(file.path(dir_result, "4.umap_patient.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()


p1 <- p + theme(legend.position = "right" )
p_legend <- cowplot::get_legend(p1)
png(file.path(dir_result, "4.umap_patient_legend.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); grid::grid.newpage(); grid::grid.draw(p_legend); dev.off()
 


##### visualization for marker genes
#### heatmap
library(Scillus)

sce2 <- sce
sce2 <- subset(sce, downsample = 1000)

choo_marker <- c("CYP2E1", "CYP2B6", "CYP3A4", "GLUL", "SDS", 
                 "CD163", "MSR1", 
                 "ASPM", "TOP2A", "CENPF", "MKI67", "SMC4",                  
                 "AKAP12", "VWF", "LDB2", "PTPRB", "PLEKHG1",
                 "ADAMTSL2", "PTH1R", "BGN", "ADAMTS2", 
                 "CFTR", "ANXA4")

choo_marker <- choo_marker[choo_marker %in% rownames(sce2@assays$RNA@scale.data)]


ls_color = c('#00BEC3','#629CFD','#4BCE75','#F7756C','#F566E1','#B69E0E')

palette_color <- paletteer::paletteer_d("ggsci::default_ucscgb")  # palette_color <- palettes_d_names
ls_color2 <- as.character(palette_color[2:23])

p <- plot_heatmap(dataset = sce2, 
             markers = choo_marker,
             sort_var = c("celltype", "patient"),
             anno_var = c("celltype", "patient", "tissue"),
             anno_colors = list(ls_color,                            
                                ls_color2,
                                c("#F7756C", "#00BEC3"))
); p


png(file.path(dir_result, "4.marker_all.png"), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()


####  FeaturePlot
choo_marker <- c("CYP2E1", 
                 "CD163",
                 "ASPM",                  
                 "AKAP12",
                 "ADAMTSL2",
                 "CFTR")
p <- FeaturePlot(sce, features = c(choo_marker), reduction = 'umap', ncol = 2); p
png(file.path(dir_result, "4.FeaturePlot.png"), width = 3500, height = 6000, res = 400, bg = "transparent"); print(p); dev.off()



##### Sample composition among patients (stacked bar chart)
#### Cell type composition
metadata <- sce@meta.data %>% select(patient, celltype)

# Calculate the proportion of each cell type for each patient
celltype_proportions <- metadata %>%
  group_by(patient, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(patient) %>%
  mutate(proportion = count / sum(count))


# Hepatocyte ratio
ratio <- celltype_proportions[celltype_proportions$celltype == 'Hepatocyte',]
min(ratio$proportion)
max(ratio$proportion)
mean(ratio$proportion)

ls_color = c('#00BEC3','#629CFD','#4BCE75','#F7756C','#F566E1','#B69E0E')

p <- ggplot(celltype_proportions, aes(x = patient, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "patient", y = "proportion", fill = "cell type") +
  scale_fill_manual(values = ls_color) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 3, hjust = 1.7, size = 5),
    axis.title.x = element_text(margin = margin(t = -10), size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ); p

png(file.path(dir_result, "4.composition.celltype.png"), width = 2800, height = 1500, res = 600, bg = "transparent"); print(p); dev.off()

p1 <- p + theme(
  legend.position = "right",
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12))
p_legend <- cowplot::get_legend(p1)
png(file.path(dir_result, "4.composition.celltype_legend.png"), width = 1000, height = 1000, res = 400, bg = "transparent"); grid::grid.newpage(); grid::grid.draw(p_legend); dev.off()




#### tissue composition
metadata <- sce@meta.data %>% select(patient, tissue)


# Calculate the proportion of each tissue for each patient
celltype_proportions <- metadata %>%
  group_by(patient, tissue) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(patient) %>%
  mutate(proportion = count / sum(count))


# cancer tissue ratio
ratio <- celltype_proportions[celltype_proportions$tissue == 'cancer',]
min(ratio$proportion)
max(ratio$proportion)
mean(ratio$proportion)


ls_color = c("#F7756C", "#00BEC3")

p <-ggplot(celltype_proportions, aes(x = patient, y = proportion, fill = tissue)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "patient", y = "proportion", fill = "tissue") +
  scale_fill_manual(values = ls_color) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 3, hjust = 1.7, size = 5),
    axis.title.x = element_text(margin = margin(t = -10), size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ); p


png(file.path(dir_result, "4.composition.tissue.png"), width = 2000, height = 1500, res = 600, bg = "transparent"); print(p); dev.off()


p1 <- p + theme(
  legend.position = "right",
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12))
p_legend <- cowplot::get_legend(p1)
png(file.path(dir_result, "4.composition.tissue_legend.png"), width = 1000, height = 1000, res = 400, bg = "transparent"); grid::grid.newpage(); grid::grid.draw(p_legend); dev.off()
