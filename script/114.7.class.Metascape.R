################################################################################
##### setting
library(ggplot2)
library(dplyr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
ls_color <- c("Class_1" = "#EA8379", "Class_2" = "#7DAEE0", "Class_3" = "#B395BD")
################################################################################
##### Class_1
choo_Class = "Class_1"
df_GO <- read.delim(file.path(dir_dataset, "Metascape", choo_Class, "Enrichment_GO/_FINAL_GO.csv"), sep = ',')

choo_row <- c('carboxylic acid metabolic process',
              'Biological oxidations',
              'small molecule catabolic process',
              'steroid metabolic process',
              'Nuclear receptors meta pathway',
              'Phase I - Functionalization of compounds',
              'Complement and coagulation cascades',
              'Bile secretion',
              'lipid catabolic process',
              'sulfur compound metabolic process'
)

choo_col <- c('GO', 'Description', 'LogP', 'X.GeneInGOAndHitList')

df_GO2 <- df_GO[df_GO$Description %in% choo_row, colnames(df_GO) %in% choo_col] 
colnames(df_GO2) <- c('ID', 'Item', 'Log10_P', 'Genes')
df_GO2 <- df_GO2 %>%  distinct(Item, .keep_all = TRUE)



p <- ggplot(df_GO2, aes(x = Log10_P, y = reorder(ID, Log10_P), fill = Genes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Item, x = 0), hjust = 0, size = 6, color = "black", nudge_x = 1) + 
  scale_x_continuous(
    trans = 'reverse',   
    expand = c(0, 0)    
  ) + 
  scale_fill_gradient2(low = "white", high = ls_color[choo_Class]) + 
  scale_y_discrete(limits = rev, labels = function(x) str_wrap(x, width = 30)) + 
  labs(x = "Log10(P)", y = "", fill = "Genes") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(), 
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks.x = element_line(size = 0.5, colour = "black"), 
    axis.ticks.y = element_line(size = 0.5, colour = "black"), 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(1, 0),  
    legend.justification = c(1, 0) 
  ); p


png(file.path(dir_result, paste0("114.Metascape.Class_1.png")), width = 2500, height = 2500, res = 400, bg = "transparent"); print(p); dev.off()


################################################################################
##### Class_2
choo_Class = "Class_2"
df_GO <- read.delim(file.path(dir_dataset, "Metascape", choo_Class, "Enrichment_GO/_FINAL_GO.csv"), sep = ',')

choo_row <- c('digestion',
              'regulation of Wnt signaling pathway',
              'Gastric cancer network 1',
              'regulation of MAPK cascade',
              'positive regulation of mitotic nuclear division',
              'NABA MATRISOME ASSOCIATED',
              'lncRNA in canonical Wnt signaling and colorectal cancer',
              'nucleobase-containing compound transport',
              'negative regulation of protein modification process',
              'MicroRNAs in cancer'
)

choo_col <- c('GO', 'Description', 'LogP', 'X.GeneInGOAndHitList')

df_GO2 <- df_GO[df_GO$Description %in% choo_row, colnames(df_GO) %in% choo_col] 
colnames(df_GO2) <- c('ID', 'Item', 'Log10_P', 'Genes')
df_GO2 <- df_GO2 %>%  distinct(Item, .keep_all = TRUE)


p <- ggplot(df_GO2, aes(x = Log10_P, y = reorder(ID, Log10_P), fill = Genes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Item, x = 0), hjust = 0, size = 6, color = "black", nudge_x = 0.05) + 
  scale_x_continuous(
    trans = 'reverse',   
    expand = c(0, 0)    
  ) + 
  scale_fill_gradient2(low = "white", high = ls_color[choo_Class]) + 
  scale_y_discrete(limits = rev, labels = function(x) str_wrap(x, width = 30)) + 
  labs(x = "Log10(P)", y = "", fill = "Genes") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(), 
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks.x = element_line(size = 0.5, colour = "black"), 
    axis.ticks.y = element_line(size = 0.5, colour = "black"), 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(1, 0.45),  
    legend.justification = c(1, 0) 
  ); p


png(file.path(dir_result, paste0("114.Metascape.Class_2.png")), width = 2800, height = 2800, res = 400, bg = "transparent"); print(p); dev.off()


################################################################################
##### Class_3
choo_Class = "Class_3"
df_GO <- read.delim(file.path(dir_dataset, "Metascape", choo_Class, "Enrichment_GO/_FINAL_GO.csv"), sep = ',')

choo_row <- c('Extracellular matrix organization',
              'positive regulation of locomotion',
              'Neutrophil degranulation',
              'regulation of inflammatory response',
              'extracellular matrix organization',
              'positive regulation of cell adhesion',
              'cell-cell adhesion',
              'Interleukin-4 and Interleukin-13 signaling',
              'supramolecular fiber organization',
              'tissue morphogenesis'
)

choo_col <- c('GO', 'Description', 'LogP', 'X.GeneInGOAndHitList')

df_GO2 <- df_GO[df_GO$Description %in% choo_row, colnames(df_GO) %in% choo_col] 
colnames(df_GO2) <- c('ID', 'Item', 'Log10_P', 'Genes')
df_GO2 <- df_GO2 %>%  distinct(Item, .keep_all = TRUE)


p <- ggplot(df_GO2, aes(x = Log10_P, y = reorder(ID, Log10_P), fill = Genes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Item, x = 0), hjust = 0, size = 7, color = "black", nudge_x = 0.1) + 
  scale_x_continuous(
    trans = 'reverse',   
    expand = c(0, 0)    
  ) + 
  scale_fill_gradient2(low = "white", high = ls_color[choo_Class]) + 
  scale_y_discrete(limits = rev, labels = function(x) str_wrap(x, width = 30)) + 
  labs(x = "Log10(P)", y = "", fill = "Genes") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(), 
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks.x = element_line(size = 0.5, colour = "black"), 
    axis.ticks.y = element_line(size = 0.5, colour = "black"), 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(1, 0.02),  
    legend.justification = c(1, 0) 
  ); p


png(file.path(dir_result, paste0("114.Metascape.Class_3.png")), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()
