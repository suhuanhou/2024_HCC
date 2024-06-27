################################################################################
##### conda
# conda activate enrichment  # 105.1.enrichment.conda.sh
# cd /share/home/shh/HCC_snRNA-seq/23.pathway/result

# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/23.2.pathway.plot.R  "Hep02" > pathway.Hep02.log 2>&1 &
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/23.2.pathway.plot.R  "Hep09" > pathway.Hep09.log 2>&1 &
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/23.2.pathway.plot.R  "Hep2" > pathway.Hep2.log 2>&1 &
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/23.2.pathway.plot.R  "Hep9" > pathway.Hep9.log 2>&1 &
################################################################################
##### setting
# conda activate enrichment
# R

library(dplyr)
library(stringr)
library(GOplot)
library(enrichplot)
library(clusterProfiler)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "23.pathway/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "23.pathway/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session_2.txt"))
setwd(dir_result)

GO_database <- 'org.Hs.eg.db' 
KEGG_database <- 'hsa'
################################################################################
##### target
args <- commandArgs(trailingOnly = TRUE)
print(args)

# choo_subtype = 'Hep9'
# choo_subtype = 'Hep02'
choo_subtype = args[1]


################################################################################
##### load data
df_all <- readRDS(file.path(dir_dataset, paste0("marker_", choo_subtype, ".rds")))
df_all$gene <- rownames(df_all)


geneID <- bitr(df_all$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
df_all <- merge(df_all, geneID, by.x = "gene", by.y = "SYMBOL", all = FALSE)
df_all <- df_all %>% data.frame() %>%
        arrange(desc(avg_log2FC)) %>%
        distinct(gene, .keep_all = TRUE) 
# head(df_all)
# dim(df_all)
################################################################################
##### pathway
if(choo_subtype == 'Hep2'){thre_p = 0.05; thre_FC = 0.5}
if(choo_subtype == 'Hep9'){thre_p = 0.05; thre_FC = 0.35}
if(choo_subtype == 'Hep02'){thre_p = 0.05; thre_FC = 0.5}
if(choo_subtype == 'Hep09'){thre_p = 0.05; thre_FC = 0.35}
df_GO_KEGG <- df_all[df_all$avg_log2FC > thre_FC & df_all$p_val < thre_p,]

GO <- enrichGO(df_GO_KEGG$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)


KEGG <- enrichKEGG(df_GO_KEGG$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)


GSEA_input <- df_all$avg_log2FC
names(GSEA_input) <- df_all$ENTREZID
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)


# bar chart
p <- barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
pdf(paste0("GO.barplot.", choo_subtype, ".pdf")); print(p); dev.off()


p <- dotplot(KEGG)+ theme(axis.title.x = element_text(size = 30), 
                          axis.title.y = element_text(size = 40),
                          axis.text.x = element_text(size = 30), 
                          axis.text.y = element_text(size = 30)) 


png(paste0("23.KEGG.dotplot.", choo_subtype, ".png"), width = 5500, height = 4000, res = 400, bg = "transparent"); print(p); dev.off()


#### GSEA
library(ggridges)
p <- ridgeplot(GSEA_KEGG) 
pdf(paste0("KEGG.GSEA.ridge.", choo_subtype, ".pdf"), width = 20, height = 16); print(p); dev.off()

p <- gseaplot2(GSEA_KEGG,1) 
pdf(paste0("KEGG.GSEA.single.", choo_subtype, ".pdf"), width = 20, height = 16); print(p); dev.off()

p <- gseaplot2(GSEA_KEGG, 1:nrow(GSEA_KEGG)) 
pdf(paste0("KEGG.GSEA.muti.", choo_subtype, ".pdf"), width = 20, height = 16); print(p); dev.off()

################################################################################
#### data for enrichment analysis
df_save <- data.frame(GO)
colnames(df_save)[1] <- paste0('#', colnames(df_save)[1])  
write.table(df_save, file.path(dir_result, paste0("result_GO_", choo_subtype, ".txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

df_save <- data.frame(KEGG)
colnames(df_save)[1] <- paste0('#', colnames(df_save)[1])  
write.table(df_save, file.path(dir_result, paste0("result_KEGG_", choo_subtype, ".txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

df_save <- data.frame(GSEA_KEGG)
colnames(df_save)[1] <- paste0('#', colnames(df_save)[1]) 
write.table(df_save, file.path(dir_result, paste0("result_GSEA_", choo_subtype, ".txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


df_save <- data.frame(df_all)
colnames(df_save)[1] <- paste0('#', colnames(df_save)[1]) 
write.table(df_save, file.path(dir_result, paste0("result_marker_", choo_subtype, ".txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

