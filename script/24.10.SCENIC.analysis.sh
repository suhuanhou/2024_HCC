################################################################################
##### setting
conda activate scenic
R

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "24.SCENIC/dataset")
dir_scRNA = file.path(dir_main, "24.SCENIC/result")
dir_bulk = file.path(dir_main, "24.SCENIC/result_bulk")
dir_analysis = file.path(dir_main, "24.SCENIC/analysis"); if(!dir.exists(dir_analysis)) dir.create(dir_analysis, recursive = TRUE)
setwd(dir_analysis);
################################################################################
##### Find the regulon in the scRNA that matches the bulk RNA-seq data
# regulon in bulk RNA-seq data
load(file.path(dir_bulk, 'regulon_RSS.Rdata'))
ls_bulk_regulon <- names(regulons)


# RSS
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
rss <- readRDS(file.path(dir_scRNA, 'subtype_rss.rds'))
ct.col <- 'subtype'

loom <- open_loom(file.path(dir_scRNA, 'aucell.loom'))
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')


rownames(regulonAUC) <- gsub("[\\(+\\)]", "", rownames(regulonAUC))
rownames(rss) <- gsub("[\\(+\\)]", "", rownames(rss))
ls_bulk_regulon <- gsub("[\\(+\\)]", "", ls_bulk_regulon)


rssPlot <- plotRSS(rss)
regulonsToPlot <- rssPlot$rowOrder
inregulons <- regulonsToPlot


regulonsToPlot = inregulons
sce$sub_celltype <- sce@meta.data[,ct.col]
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]


mtx_auc <- sub_regulonAUC@assays@data@listData$AUC
rownames(mtx_auc) <- gsub("[\\(+\\)]", "", rownames(mtx_auc))
sce@meta.data = cbind(sce@meta.data ,t(mtx_auc[regulonsToPlot,]))
Idents(sce) <- sce$sub_celltype


meta <- sce@meta.data
celltype <- ct.col
cellsPerGroup <- split(rownames(meta),meta[,celltype])
sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

rss <- regulonActivity_byGroup

rss2 <- rss[rownames(rss) %in% ls_bulk_regulon, 
            colnames(rss) %in% c('Hep0', 'Hep2', 'Hep9')]


p <- pheatmap(rss2, scale='row', show_rownames = T,
              cluster_cols = F, 
              main='scRNA-seq', fontsize = 30,
              fontsize_row = 30, fontsize_col = 40, angle_col='45'
              )

png("24.bulk_scRNA.regulon_activity.png", width = 4000, height = 8000, res = 400, bg = "transparent"); print(p); dev.off()



################################################################################
##### TF target PDE7B
library(ggplot2)

choo_gene = 'PDE7B'
df_bulk <- read.table(file.path(dir_bulk, "grn.bulk.tsv"), header = TRUE)
df_scRNA <- read.table(file.path(dir_scRNA, "grn.scRNA.tsv"), header = TRUE)


df_bulk <- df_bulk[df_bulk$target == choo_gene,]
df_scRNA <- df_scRNA[df_scRNA$target == choo_gene,]

# filter
df_bulk2 <- df_bulk[1:20,]
df_scRNA2 <- df_scRNA[1:20,]

df_bulk2$TF <- factor(df_bulk2$TF, levels = rev(df_bulk2$TF))
df_scRNA2$TF <- factor(df_scRNA2$TF, levels = rev(df_scRNA2$TF))

ls_inte <- intersect(df_bulk2$TF, df_scRNA2$TF)

# df_bulk
p <- ggplot(df_bulk2, aes(x = importance, y = TF)) +
  geom_col(fill = "steelblue") +
  labs(x = "Importance", y = "", title = paste0("bulk RNA-seq")) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),  
        panel.background = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "gray"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 40),  
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")

        )+ 
  scale_x_continuous(expand = c(0,0))

pdf(paste0("24.Importance.bulk.", choo_gene, ".pdf"), width = 10, height = 20); print(p); dev.off()


# df_scRNA
p <- ggplot(df_scRNA2, aes(x = importance, y = TF)) +
  geom_col(fill = "steelblue") +
  labs(x = "Importance", y = "", title = paste0("scRNA-seq")) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "gray"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 40),  
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )+ 
  scale_x_continuous(expand = c(0,0)) 

pdf(paste0("24.Importance.scRNA.", choo_gene, ".pdf"), width = 10, height = 20); print(p); dev.off()

