################################################################################
##### setting
library(tidyverse)
library(ggpubr)


# rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
# dir_in = file.path(dir_main, "0.data_bulk_RNA-seq/result/SEQ-DX-PXH-202402021501_result_final")
dir_dataset = file.path(dir_main, "101.bulk/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "101.bulk/result/pair_boxplot"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result); 

bool_save = TRUE 
################################################################################
##### read data
df_normalize <- readRDS(file.path(dir_dataset, "df_normalize.rds"))
ls_gene <- readRDS(file.path(dir_dataset, "gene_diagnosis_internal.rds"))  # 


df_symbol <- readRDS(file.path(dir_dataset, "df_symbol.rds"))
# dict_gene <- list(
#   ENSG00000111700 = 'SLCO1B3',
#   ENSG00000198300 = 'PEG3',
#   ENSG00000112964 = 'GHR',
#   ENSG00000147257 = 'GPC3')


# gene_id
dict_gene <- lapply(ls_gene, function(gene) {
  matching_rows <- df_symbol[df_symbol$gene_name == gene, ]
  if (nrow(matching_rows) == 1) { 
    gene_id <- matching_rows$gene_id[1] 
    return(gene_id)
  } else if (nrow(matching_rows) == 0) {
    print(paste0(gene, ' : not found'))
  } else {
    print(paste0(gene, ' : duplicate'))
  }
})


dict_gene <- setNames(ls_gene, dict_gene) 

################################################################################
##### pair_boxplot
for (ENSEMBL in names(dict_gene)) {
  # ENSEMBL = 'ENSG00000112964'
  print(dict_gene[[ENSEMBL]])

  df_gene <- df_normalize[df_normalize$ENSEMBL == ENSEMBL,]
  df_gene <- df_gene[,-1]
  # colnames(df_gene)
  
  df_gene2 <- data.frame(ID = paste0("ID", 1:55),
                         Normal = numeric(55), 
                         Paracancer = numeric(55), 
                         Cancer = numeric(55))
  

  for (i in 1:55) {
    df_gene2$Cancer[i] <- df_gene[, paste0("bulk_C_", i)]
    df_gene2$Normal[i] <- df_gene[, paste0("bulk_N_", i)]
    df_gene2$Paracancer[i] <- df_gene[, paste0("bulk_P_", i)]
  }
  

  drawdata<- df_gene2 %>% gather("Group", "value", -ID)
  
  my_comparisons <- list(c("Paracancer", "Normal"), c("Cancer", "Paracancer"),c("Cancer", "Normal"))
  color_group <- c("Normal" = "#E4F99D", "Paracancer" = "#E39C18", "Cancer" = "#FF6207")
  

  p <- ggpaired(drawdata, x = 'Group', y = 'value', id='ID',
                color = 'Group', palette = "jco",
                line.color = "gray", line.size = 0.4,
                short.panel.labs = FALSE) +
    ylab("normalized count") +
    labs(x = "") + 
    ggtitle(dict_gene[[ENSEMBL]]) + 
    theme(axis.title.x = element_text(face = "bold"),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)) +  
    stat_compare_means(
      aes(label = ..p.signif..),
      method = "t.test",
      paired = TRUE,
      size = 5,
      comparisons = my_comparisons) +
    scale_color_manual(values = color_group)
  

  png(file.path(dir_result, paste0("101.", dict_gene[[ENSEMBL]], ".png")), width = 1900, height = 1900, res = 400, bg = "transparent")
  print(p)
  dev.off()
  
}

