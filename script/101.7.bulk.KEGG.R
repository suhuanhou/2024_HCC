################################################################################
##### conda 
# conda activate enrichment  # 105.1.enrichment.conda.sh
# cd /share/home/shh/HCC_snRNA-seq/101.bulk/result
# nohup Rscript /share/home/shh/HCC_snRNA-seq/script/101.7.bulk.KEGG.R  > KEGG.log 2>&1 &

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
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
# dir_in = file.path(dir_main, "0.data_bulk_RNA-seq/result/SEQ-DX-PXH-202402021501_result_final")
dir_dataset = file.path(dir_main, "101.bulk/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "101.bulk/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result); 

KEGG_database <- 'hsa'
GO_database <- 'org.Hs.eg.db' 
thre_p = 0.05
thre_FC = 1.5
################################################################################
##### KEGG UP
df_CP <- readxl::read_excel(file.path(dir_result, "bulk_DEG.xlsx"), sheet = 'CP')
df_CP <- df_CP %>% data.frame() %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(gene_name, .keep_all = TRUE) 
df_CP <- df_CP[df_CP$log2FoldChange > thre_FC & df_CP$pvalue < thre_p,]


df_CN <- readxl::read_excel(file.path(dir_result, "bulk_DEG.xlsx"), sheet = 'CN')
df_CN <- df_CN %>% data.frame() %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(gene_name, .keep_all = TRUE) 
df_CN <- df_CN[df_CN$log2FoldChange > thre_FC & df_CN$pvalue < thre_p,]

df_UP <- df_CP[df_CP$gene_name %in% df_CN$gene_name,]
# dim(df_UP)


##### KEGG UP
geneID <- bitr(df_UP$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
df_UP_KEGG <- merge(df_UP, geneID, by.x = "gene_name", by.y = "SYMBOL", all = FALSE)


KEGG_UP <- enrichKEGG(df_UP_KEGG$ENTREZID,
                   organism = KEGG_database,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)


KEGG_UP2 <- KEGG_UP
KEGG_UP2 <- setReadable(KEGG_UP, OrgDb =GO_database, keyType="ENTREZID")  # symbol替换ID


genelist <- df_UP_KEGG$log2FoldChange
names(genelist) <- df_UP_KEGG$gene_name

p <- cnetplot(KEGG_UP2, showCategory=5, foldChange=genelist, circular=TRUE, colorEdge=TRUE)
png("101.KEGG.cnetplot.UP.png", width = 6000, height = 5000, res = 400, bg = "transparent"); print(p); dev.off()

write.table(data.frame(KEGG_UP2), file.path(dir_result, paste0("101.result_KEGG.UP.txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



################################################################################
##### KEGG DOWN
df_CP <- readxl::read_excel(file.path(dir_result, "bulk_DEG.xlsx"), sheet = 'CP')
df_CP <- df_CP %>% data.frame() %>%
  arrange(log2FoldChange) %>%
  distinct(gene_name, .keep_all = TRUE) 
df_CP <- df_CP[df_CP$log2FoldChange < -thre_FC & df_CP$pvalue < thre_p,]


df_CN <- readxl::read_excel(file.path(dir_result, "bulk_DEG.xlsx"), sheet = 'CN')
df_CN <- df_CN %>% data.frame() %>%
  arrange(log2FoldChange) %>%
  distinct(gene_name, .keep_all = TRUE) 
df_CN <- df_CN[df_CN$log2FoldChange < -thre_FC & df_CN$pvalue < thre_p,]

df_DOWN <- df_CP[df_CP$gene_name %in% df_CN$gene_name,]
# dim(df_UP)


##### KEGG DOWN
geneID <- bitr(df_CN$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
df_DOWN_KEGG <- merge(df_DOWN, geneID, by.x = "gene_name", by.y = "SYMBOL", all = FALSE)


KEGG_DOWN <- enrichKEGG(df_DOWN_KEGG$ENTREZID,
                        organism = KEGG_database,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)


KEGG_DOWN2 <- KEGG_DOWN
KEGG_DOWN2 <- setReadable(KEGG_DOWN, OrgDb =GO_database, keyType="ENTREZID")  


genelist <- df_DOWN_KEGG$log2FoldChange
names(genelist) <- df_DOWN_KEGG$gene_name

p <- cnetplot(KEGG_DOWN2, showCategory=5, foldChange=genelist, circular=TRUE, colorEdge=TRUE)
png("101.KEGG.cnetplot.DOWN.png", width = 6000, height = 5000, res = 400, bg = "transparent"); print(p); dev.off()

write.table(data.frame(KEGG_DOWN2), file.path(dir_result, paste0("101.result_KEGG.DOWN.txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

