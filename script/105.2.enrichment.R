# GO、KEGG、GSEA
################################################################################
##### setting
# conda activate enrichment 
# R

library(readxl)
library(dplyr)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)  


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "101.bulk")
dir_dataset = file.path(dir_main, "105.enrichment/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "105.enrichment/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result); 

bool_save = TRUE
################################################################################
##### read data
# df_normalize <- readRDS(file.path(dir_in, "dataset/df_normalize.rds"))
DESeq2 <- read_excel(file.path(dir_in, "result/bulk_DEG.xlsx"), sheet = "CP")
# colnames(df_DEG)[1]<-"Gene" 


# DESeq2$symbol <- row.names(DESeq2)
DESeq2$Group <- ifelse(DESeq2$log2FoldChange > 1.5, "Up", ifelse(DESeq2$log2FoldChange<(-1.5), "Down", "none"))


deg_up <- DESeq2$Row.names[DESeq2$Group == "Up"]
deg_down <- DESeq2$Row.names[DESeq2$Group == "Down"] 
# diff <- c(up, down)
# head(diff) 

## ID
entrez_up <- bitr(geneID = deg_up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_down <- bitr(geneID = deg_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


# KEGG
KEGG_up<- enrichKEGG(gene = entrez_up$ENTREZID,
                        organism = "hsa", 
                        keyType = 'kegg',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.25)


KEGG_down <- enrichKEGG(gene = entrez_down$ENTREZID,
                     organism = "hsa", 
                     keyType = 'kegg',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.25)


# symbol
KEGG_up <- setReadable(KEGG_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG_down <- setReadable(KEGG_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# head(KEGG_diff@result)


# Rich Factor
KEGG_up <- mutate(KEGG_up, RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
KEGG_down <- mutate(KEGG_down, RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

# Fold Enrichment
KEGG_up <- mutate(KEGG_up, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_up@result$RichFactor[1:6]
KEGG_up@result$FoldEnrichment[1:6]

KEGG_down <- mutate(KEGG_down, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_down@result$RichFactor[1:6]
KEGG_down@result$FoldEnrichment[1:6]


write.csv(x = KEGG_up,file = "105.KEGG_up.csv")
write.csv(x = KEGG_down,file = "105.KEGG_down.csv")
