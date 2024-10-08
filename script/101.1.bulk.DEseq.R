# DEseq2:bulk RNA-seq
################################################################################
##### setting
library(DESeq2)
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "0.data_bulk_RNA-seq/result/SEQ-DX-PXH-202402021501_result_final")
dir_dataset = file.path(dir_main, "101.bulk/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "101.bulk/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset); 

save_ggplot <- function(gg, prefix, width = 8, height = 8, dpi = 600) {
  ggsave(paste0(prefix, ".png"), plot = gg, width = width, height = height, dpi = dpi)
}

################################################################################
##### read data
df_exp <- read.table(file.path(dir_in, "3.Expression/all_samples_expression_raw_count.xls"), header = TRUE, sep = "\t")
df_bulk <- df_exp[,1:(ncol(df_exp)-10)]
rownames(df_bulk) <- df_bulk$gene_id
df_bulk$gene_id <- NULL

df_tmp <- cbind(df_bulk[,56:165], df_bulk[,1:55])
df_bulk <- df_tmp

df_symbol <- df_exp[,c('gene_id', 'gene_name')]
saveRDS(df_symbol, file=file.path(dir_dataset, 'df_symbol.rds'))


#### colData
# colnames(df_bulk)
group <- ifelse(grepl("bulk_P", colnames(df_bulk)), "Paracancer",
                ifelse(grepl("bulk_C", colnames(df_bulk)), "Cancer",
                       ifelse(grepl("bulk_N", colnames(df_bulk)), "Normal", NA)))
group <- factor(group, levels = c('Normal', 'Paracancer', 'Cancer'))
coldata <- data.frame(row.names=colnames(df_bulk), group)


#### DESeq2
dds <- DESeqDataSetFromMatrix(countData = df_bulk, colData = coldata, design= ~group)
dds <- dds[rowSums(counts(dds) > 0) >= 3,]  
dds <- DESeq(dds)
saveRDS(dds, file=file.path(dir_dataset, 'dds.rds'))

#### normalize
df_normalize <- counts(dds, normalized=TRUE) %>% as.data.frame %>% rownames_to_column("ENSEMBL")
write.table(df_normalize, file.path(dir_dataset, 'df_normalize.txt'), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
saveRDS(df_normalize, file=file.path(dir_dataset, 'df_normalize.rds'))

#### symbol
# head(df_normalize)[1:4]
# head(probe2symbol)
ids <- df_symbol
names(ids)[1:2] <- c("probe_id","symbol")

df_tmp <- df_normalize
rownames(df_tmp) <- df_tmp$ENSEMBL
df_tmp <- df_tmp[,-1]

ids = ids[match(rownames(df_tmp),ids$probe_id),]  

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))])
  probes = as.character(tmp)
  print(dim(exprSet))  
  exprSet = exprSet[rownames(exprSet) %in% probes,]
  
  print(dim(exprSet))  
  rownames(exprSet) = ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}


df_normalize_symbol <- jimmy(df_tmp, ids)  
saveRDS(df_normalize_symbol, file=file.path(dir_dataset, 'df_normalize_symbol.rds'))


################################################################################
##### DEG
setwd(dir_result)

dds <- readRDS(file.path(dir_dataset, 'dds.rds'))
df_symbol <- readRDS(file.path(dir_dataset, 'df_symbol.rds'))

DEG_analysis <- function(group1, group2, name_sheet){
  # group1 = 'Cancer'; group2 = 'Normal'; name_sheet = 'CN'
  res <- results(dds, contrast = c('group', group1, group2))
  res_df <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  df_symbol2 <- df_symbol[df_symbol$gene_id %in% rownames(res_df), ]
  df_symbol2 <- df_symbol2[match(rownames(res_df), df_symbol2$gene_id),]
  df_cont <- merge(res_df, df_symbol2, by.x = 0, by.y = "gene_id", all.x = TRUE)
  addWorksheet(wb, name_sheet)
  writeData(wb, sheet = name_sheet, x = df_cont, startCol = 1, startRow = 1)
  
  df_cont_DEG <- df_cont[df_cont$pvalue < thre_p & abs(df_cont$log2FoldChange) > thre_FC, ]
  df_cont_DEG <- df_cont_DEG[order(-df_cont_DEG$log2FoldChange), ]
  addWorksheet(wb, paste0(name_sheet, '_DEG'))
  writeData(wb, sheet = paste0(name_sheet, '_DEG'), x = df_cont_DEG, startCol = 1, startRow = 1)
  
  tmp_res <- df_cont %>%
  mutate(change = case_when(pvalue < thre_p & log2FoldChange > thre_FC ~ "Up",
                            pvalue < thre_p & log2FoldChange < -thre_FC ~ "Down",
                            TRUE ~ "Stable"))
  
  colors = c("darkgreen","grey","darkred")
  names(colors) = c("Down","Stable","Up")
  tmp_y = max(-log10(tmp_res$pvalue))
  if (tmp_y < 10){
    y_height = 10
  } else if (tmp_y < 200){
    y_height = ceiling(tmp_y / 10) * 10
  } else {
    y_height = ceiling(tmp_y / 10) * 10 + 50
  }
  tmp_x = max(abs(tmp_res$log2FoldChange))
  if (tmp_x < 10){
    x_length = 10
  }else {
    x_length = 15
  }
  
  tmp_volcano <- ggplot(tmp_res, aes(x=log2FoldChange, y=-log10(pvalue),
                                     color=change)) +
    geom_point(alpha=0.4,size=3)+
    scale_color_manual(values=colors) +
    theme_bw() +
    ylim(c(0,y_height))+
    xlim(c(-x_length,x_length)) +
    geom_vline(xintercept = c(-thre_FC,thre_FC), lty=4, col="black", lwd=0.8)+
    geom_hline(yintercept = -log10(thre_p), lty=4, col="black", lwd=0.8) +
    labs(x = expression(log[2]~'FC'),
         y=expression(-log10(pvalue)),
         title=paste0(group1, " vs. ", group2)) +
    theme(axis.title.x = element_text(colour="black",size = 22),
          axis.title.y = element_text(colour="black",size = 22),
          axis.text.x = element_text(colour="black",size = 20),
          axis.text.y = element_text(colour="black",size = 20),
          plot.title = element_text(size=24),
          legend.text = element_text(color="black", size=20))+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="none", 
          legend.title = element_blank())
  save_ggplot(tmp_volcano, file.path(dir_result, paste0("101.volcano.", name_sheet)), width=5, height=5)
  
  return(wb)
}


thre_p = 0.05
thre_FC = 1.5

wb <- createWorkbook()
wb <- DEG_analysis('Cancer', 'Paracancer', 'CP')
wb <- DEG_analysis('Cancer', 'Normal', 'CN')
wb <- DEG_analysis('Paracancer', 'Normal', 'PN')
saveWorkbook(wb, "bulk_DEG.xlsx", overwrite = TRUE)


################################################################################
##### visual
setwd(dir_result);
thre_p = 0.05
thre_FC = 1.5

df_count <- data.frame(group = '', Up = '', Down = '')
df <- read.xlsx("bulk_DEG.xlsx", sheet = "CP_DEG")  
df_tmp <- df %>%
  mutate(change = case_when(pvalue < thre_p & log2FoldChange > thre_FC ~ "Up",
                            pvalue < thre_p & log2FoldChange < -thre_FC ~ "Down",
                            TRUE ~ "Stable"))
df_count[1,] <- c('Cancer VS Paracancer', sum(df_tmp$change == "Up"), sum(df_tmp$change == "Down"))


df <- read.xlsx("bulk_DEG.xlsx", sheet = "CN_DEG")  
df_tmp <- df %>%
  mutate(change = case_when(pvalue < thre_p & log2FoldChange > thre_FC ~ "Up",
                            pvalue < thre_p & log2FoldChange < -thre_FC ~ "Down",
                            TRUE ~ "Stable"))
df_count[2,] <- c('Cancer VS Normal', sum(df_tmp$change == "Up"), sum(df_tmp$change == "Down"))


df <- read.xlsx("bulk_DEG.xlsx", sheet = "PN_DEG")  
df_tmp <- df %>%
  mutate(change = case_when(pvalue < thre_p & log2FoldChange > thre_FC ~ "Up",
                            pvalue < thre_p & log2FoldChange < -thre_FC ~ "Down",
                            TRUE ~ "Stable"))
df_count[3,] <- c('Paracancer VS Normal', sum(df_tmp$change == "Up"), sum(df_tmp$change == "Down"))



library(ggplot2)
library(reshape2)
df_long <- melt(df_count, id.vars = "group")
df_long$value <- as.numeric(df_long$value)
df_long$group <- factor(df_long$group, levels = c('Cancer VS Paracancer', 'Cancer VS Normal', 'Paracancer VS Normal'))

p <- ggplot(df_long, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Differentially expressed genes", x = "", y = "Number of genes", fill = "Change") +
  scale_fill_manual(values = c("Up" = "#F16733", "Down" = "#7D67F6")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),  
    axis.title = element_text(size = 24), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 20),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),  
    panel.grid.minor.y = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, 2500, by = 500), limits = c(0, 2500), expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +  
  coord_cartesian(clip = 'off'); p 

png(file.path(dir_result, "101.DEG.barChat.png"), width = 4500, height = 2500, res = 400, bg = "transparent"); print(p); dev.off()

