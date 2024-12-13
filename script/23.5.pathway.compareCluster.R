cd  /share/home/shh/HCC_snRNA-seq/23.pathway/5.compareCluster
conda activate enrichment
R
################################################################################
##### setting
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(Seurat)
library(stringr)
library(dplyr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
# dir_dataset = file.path(dir_main, "23.pathway/5.compareCluster"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "23.pathway/5.compareCluster"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result)
bool_restar = F
################################################################################
##### marker
if(bool_restar){
  sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))

  heps <- c('Hep1', 'Hep2', 'Hep3', 'Hep4', 'Hep5', 'Hep6', 'Hep7', 'Hep8', 'Hep9')
  for (hep in heps) {
    marker <- FindMarkers(object = sce, ident.1 = hep, ident.2 = 'Hep0')
    marker <- rownames_to_column(marker, var = "gene")
    saveRDS(marker, file.path(dir_dataset, paste0("marker_", hep, ".rds")))
    # write.xlsx(marker, file = file.path(dir_result, paste0("DEG_", hep, ".xlsx")))
  }
}

################################################################################
##### load data
heps <- c('Hep1', 'Hep2', 'Hep3', 'Hep4', 'Hep5', 'Hep6', 'Hep7', 'Hep8', 'Hep9')
df_marker <- data.frame()

for (hep in heps) {
  marker <- readRDS(file.path(dir_result, paste0("marker_", hep, ".rds")))
  marker$Hep <- hep
  df_marker <- rbind(df_marker, marker)
}


ids <- bitr(df_marker$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
df_marker <- merge(df_marker, ids, by.x='gene', by.y='SYMBOL')


### UP
df_marker2 <- df_marker[df_marker$avg_log2FC>0,]

df_top <- data.frame()
df_top <- rbind(df_top, 
                df_marker2[df_marker2$Hep == "Hep1", ][1:50, ],
                df_marker2[df_marker2$Hep == "Hep2", ][1:40, ],                
                df_marker2[df_marker2$Hep == "Hep3", ][1:40, ],                
                df_marker2[df_marker2$Hep == "Hep4", ][1:50, ],                
                df_marker2[df_marker2$Hep == "Hep5", ][1:150, ],                
                df_marker2[df_marker2$Hep == "Hep6", ][1:200, ],                
                df_marker2[df_marker2$Hep == "Hep7", ][1:50, ],                
                df_marker2[df_marker2$Hep == "Hep8", ][1:50, ],                
                df_marker2[df_marker2$Hep == "Hep9", ][1:200, ]               
)

ls_gene <- split(df_top$ENTREZID, df_top$Hep)
CC <- compareCluster(ls_gene, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)

p <- clusterProfiler::dotplot(CC)+ 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(),

  )
 
png("23.5.compareCluster.UP.png", width = 3500, height = 3110, res = 400, bg = "transparent"); print(p); dev.off()
print(CC@compareClusterResult$Description)


### DOWN
df_marker2 <- df_marker[df_marker$avg_log2FC<0,]

df_top <- data.frame()
df_top <- rbind(df_top, 
                df_marker2[df_marker2$Hep == "Hep1", ][1:100, ],
                df_marker2[df_marker2$Hep == "Hep2", ][1:80, ],                
                df_marker2[df_marker2$Hep == "Hep3", ][1:80, ],                
                df_marker2[df_marker2$Hep == "Hep4", ][1:100, ],                
                df_marker2[df_marker2$Hep == "Hep5", ][1:100, ],                
                df_marker2[df_marker2$Hep == "Hep6", ][1:150, ],                
                df_marker2[df_marker2$Hep == "Hep7", ][1:100, ],                
                df_marker2[df_marker2$Hep == "Hep8", ][1:100, ],                
                df_marker2[df_marker2$Hep == "Hep9", ][1:200, ]               
)


ls_gene <- split(df_top$ENTREZID, df_top$Hep)
CC <- compareCluster(ls_gene, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)


p <- clusterProfiler::dotplot(CC)+ 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(),
  )


png("23.5.compareCluster.DOWN.png", width = 3375, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()
print(CC@compareClusterResult$Description)
