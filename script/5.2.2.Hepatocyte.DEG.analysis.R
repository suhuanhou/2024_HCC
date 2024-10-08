################################################################################
##### setting
library(openxlsx)
library(VennDiagram)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
# dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "5.Hepatocyte/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/DEG"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session_1.txt"))
setwd(dir_result)
################################################################################
##### load data
#### DEG for snRNA-seq
marker_Hep2 <- readRDS(file.path(dir_dataset, "marker_Hep2.rds"))
marker_Hep9 <- readRDS(file.path(dir_dataset, "marker_Hep9.rds"))
DEG_Hep2 <-marker_Hep2[marker_Hep2$p_val < 0.05,]
DEG_Hep9 <-marker_Hep9[marker_Hep9$p_val < 0.05,]

DEG_Hep2 <- DEG_Hep2[order(DEG_Hep2$avg_log2FC),]
DEG_Hep9 <- DEG_Hep9[order(DEG_Hep9$avg_log2FC),]
write.xlsx(DEG_Hep2, file = file.path(dir_result, "DEG_Hep2.xlsx"), row.names = TRUE)
write.xlsx(DEG_Hep9, file = file.path(dir_result, "DEG_Hep9.xlsx"), row.names = TRUE)

#### DEG for bulk RNA-seq
DEG_bulk_CP <- read.xlsx(file.path(dir_main, '101.bulk/result/bulk_DEG.xlsx'), sheet = 'CP_DEG')
DEG_bulk_CN <- read.xlsx(file.path(dir_main, '101.bulk/result/bulk_DEG.xlsx'), sheet = 'CP_DEG')


################################################################################
##### visualization: Venn
#### down
set1 <- DEG_bulk_CP[DEG_bulk_CP$log2FoldChange < 0, 'gene_name']
set2 <- DEG_bulk_CN[DEG_bulk_CN$log2FoldChange < 0, 'gene_name']
set3 <- rownames(DEG_Hep2[DEG_Hep2$avg_log2FC < 1,])
set4 <- rownames(DEG_Hep9[DEG_Hep9$avg_log2FC < 1,])
ls_down <- list(A = set1, B = set2, C = set3, D = set4)

dev.off()
venn.plot <- venn.diagram(
  x = ls_down,
  # category.names = c("CP", "CN", "Hep2 VS Hep0", "Hep9 VS Hep0"),
  category.names = c("", "", "", ""),
  filename = NULL,
  fill = c('red', 'blue', 'yellow', 'green'), alpha = 0.35,
  col = 'black', cex = 1.5, fontfamily = 'Arial', 
  cat.cex = 1.5, cat.fontfamily = 'Arial'
  
); grid.draw(venn.plot)

dev.copy(png, file.path(dir_result, "5.venn.down.png"), width = 2500, height = 2000, res = 400, bg = "transparent"); dev.off()

inte_down <- Reduce(intersect, list(set1, set2, set3, set4)); inte_down
inte_Hep2 <- Reduce(intersect, list(set1, set2, set3)); inte_Hep2
inte_Hep9 <- Reduce(intersect, list(set1, set2, set4)); inte_Hep9


#### up
set1 <- DEG_bulk_CP[DEG_bulk_CP$log2FoldChange > 0, 'gene_name']
set2 <- DEG_bulk_CN[DEG_bulk_CN$log2FoldChange > 0, 'gene_name']
set3 <- rownames(DEG_Hep2[DEG_Hep2$avg_log2FC > 1,])
set4 <- rownames(DEG_Hep9[DEG_Hep9$avg_log2FC > 1,])
ls_up <- list(A = set1, B = set2, C = set3, D = set4)

dev.off()
venn.plot <- venn.diagram(
  x = ls_up,
  # category.names = c("CP", "CN", "Hep2 VS Hep0", "Hep9 VS Hep0"),
  category.names = c("", "", "", ""),
  filename = NULL,
  fill = c('red', 'blue', 'yellow', 'green'), alpha = 0.35,
  col = 'black', cex = 1.5, fontfamily = 'Arial', 
  cat.cex = 1.5, cat.fontfamily = 'Arial'
  
); grid.draw(venn.plot)

dev.copy(png, file.path(dir_result, "5.venn.up.png"), width = 2500, height = 2000, res = 400, bg = "transparent"); dev.off()

inte_up <- Reduce(intersect, list(set1, set2, set3, set4)); inte_up
inte_Hep2 <- Reduce(intersect, list(set1, set2, set3)); inte_Hep2
inte_Hep9 <- Reduce(intersect, list(set1, set2, set4)); inte_Hep9

