# DEGs
################################################################################
##### setting
library(dplyr)
library(Seurat)
library(openxlsx)
# library(tibble)
library(UpSetR)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "5.Hepatocyte/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "5.Hepatocyte/5.2.DEG"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
sce <- readRDS(file.path(dir_in, "sce_Hep.rds"))
Idents(sce) <- sce$subtype

#### DEGs
heps <- c('Hep1', 'Hep2', 'Hep3', 'Hep4', 'Hep5', 'Hep6', 'Hep7', 'Hep8', 'Hep9')

for (hep in heps) {
  marker <- FindMarkers(object = sce, ident.1 = hep, ident.2 = 'Hep0')
  marker <- rownames_to_column(marker, var = "gene")
  saveRDS(marker, file.path(dir_dataset, paste0("marker_", hep, ".rds")))
  write.xlsx(marker, file = file.path(dir_result, paste0("DEG_", hep, ".xlsx")))
}


# marker_Hep2 <- readRDS(file.path(dir_dataset, "marker_Hep2.rds"))
# marker_Hep9 <- readRDS(file.path(dir_dataset, "marker_Hep9.rds"))

marker_Hep02 = FindMarkers(object = sce, ident.1 = 'Hep0', ident.2 = 'Hep2')
marker_Hep09 = FindMarkers(object = sce, ident.1 = 'Hep0', ident.2 = 'Hep9')
marker_Hep02 <- rownames_to_column(marker_Hep02, var = "gene")
marker_Hep09 <- rownames_to_column(marker_Hep09, var = "gene")
saveRDS(marker_Hep02, file.path(dir_dataset, "marker_Hep02.rds"))
saveRDS(marker_Hep09, file.path(dir_dataset, "marker_Hep09.rds"))


################################################################################
##### upset
heps <- c('Hep1', 'Hep2', 'Hep3', 'Hep4', 'Hep5', 'Hep6', 'Hep7', 'Hep8', 'Hep9')
gene_list <- list()

# markers
for (hep in heps) {
  marker_data <- read.xlsx(file.path(dir_result, paste0("DEG_", hep, ".xlsx")))
  
  top_genes <- marker_data %>%
    head(20) %>%
    pull(gene)  
  
  gene_list[[hep]] <- top_genes
}

all_genes <- unique(unlist(gene_list))

gene_matrix <- data.frame(matrix(0, nrow = length(all_genes), ncol = length(heps)))
colnames(gene_matrix) <- heps
rownames(gene_matrix) <- all_genes


for (hep in heps) {
  genes_in_hep <- gene_list[[hep]]
  gene_matrix[genes_in_hep, hep] <- 1
}

# UpSet
p <-  upset(gene_matrix, 
            sets = heps, 
            main.bar.color = "black", 
            order.by = "freq",           
            keep.order = TRUE,           
            point.size = 4,            
            line.size = 2,              
            mb.ratio = c(0.5, 0.5),     
            text.scale = c(1.5, 1.5, 0, 0, 1.5, 1.5), 
            queries = list(
              list(query = intersects, params = list("Hep1", "Hep3", "Hep4"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep1", "Hep3", "Hep4", "Hep6"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep2", "Hep5", "Hep6"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep1", "Hep2", "Hep3", "Hep4", "Hep6"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep1", "Hep2", "Hep3", "Hep4", "Hep7"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep1", "Hep2", "Hep6"), color = "red", active = TRUE),
              list(query = intersects, params = list("Hep1", "Hep2", "Hep4"), color = "red", active = TRUE)
            )
);p

png(file.path(dir_result, paste0("5.2.DEG.UpSet.png")), width = 4500, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()


## overlap
df_view <- gene_matrix
df_view$count = rowSums(df_view)

# Hep1/3/4
df_view[df_view$Hep1 == 1 & df_view$Hep3 == 1 & df_view$Hep4 == 1 & df_view$count == 3,]
# Hep1/3/4/6
df_view[df_view$Hep1 == 1 & df_view$Hep3 == 1 & df_view$Hep4 == 1 & df_view$Hep6 == 1 & df_view$count == 4,]
# Hep2/5/6
df_view[df_view$Hep2 == 1 & df_view$Hep5 == 1 & df_view$Hep6 == 1 & df_view$count == 3,]
# Hep1/2/6
df_view[df_view$Hep1 == 1 & df_view$Hep2 == 1 & df_view$Hep6 == 1 & df_view$count == 3,]
# Hep1/2/4
df_view[df_view$Hep1 == 1 & df_view$Hep2 == 1 & df_view$Hep4 == 1 & df_view$count == 3,]
# Hep1/2/3/4/7
df_view[df_view$Hep1 == 1 & df_view$Hep2 == 1 & df_view$Hep3 == 1 & df_view$Hep4 == 1 & df_view$Hep7 == 1 & df_view$count == 5,]
# Hep1/2/3/4/6
df_view[df_view$Hep1 == 1 & df_view$Hep2 == 1 & df_view$Hep3 == 1 & df_view$Hep4 == 1 & df_view$Hep6 == 1 & df_view$count == 5,]



##### UpSet: overlap genes
ls_v <- list()
ls_v[[1]] <- c("CDA", "Hep1",  -0.343516207)
ls_v[[2]] <- c("CDA", "Hep3",  -0.339059508)
ls_v[[3]] <- c("CDA", "Hep4",  -0.588042912)
ls_v[[4]] <- c("DDI2", "Hep1", -0.343516207)
ls_v[[5]] <- c("DDI2", "Hep3", -0.300531365)
ls_v[[6]] <- c("DDI2", "Hep4",  -0.531780449)
ls_v[[7]] <- c("IFNLR1", "Hep1", -0.490266882)
ls_v[[8]] <- c("IFNLR1", "Hep3", -0.611350498)
ls_v[[9]] <- c("IFNLR1", "Hep4", -0.636425719)
ls_v[[10]] <- c("MAN1C1", "Hep1", -0.6240216)
ls_v[[11]] <- c("MAN1C1", "Hep3", -0.634036122)
ls_v[[12]] <- c("MAN1C1", "Hep4", -0.704490134)
ls_v[[13]] <- c("PEX14", "Hep1", -0.33491554)
ls_v[[14]] <- c("PEX14", "Hep3", -0.376319935)
ls_v[[15]] <- c("PEX14", "Hep4", -0.275997676)
ls_v[[16]] <- c("PHACTR4", "Hep1", -0.276338215)
ls_v[[17]] <- c("PHACTR4", "Hep3", -0.361646342)
ls_v[[18]] <- c("PHACTR4", "Hep4", -0.345764696)


ls_v[[19]] <- c("ERRF", "Hep1", -0.848068046)
ls_v[[20]] <- c("ERRF", "Hep3", -0.8238224)
ls_v[[21]] <- c("ERRF", "Hep4", -1.537668288)
ls_v[[22]] <- c("ERRF", "Hep6", -1.202598079)
ls_v[[23]] <- c("PER3", "Hep1", -0.26525677)
ls_v[[24]] <- c("PER3", "Hep3", -0.451127511)
ls_v[[25]] <- c("PER3", "Hep4", -0.797754468)
ls_v[[26]] <- c("PER3", "Hep6", -0.640066833)

ls_v[[27]] <- c("SCP2", "Hep2", -1.016707292)
ls_v[[28]] <- c("SCP2", "Hep5", -0.867644316)
ls_v[[29]] <- c("SCP2", "Hep6", -1.034397181)


ls_v[[30]] <- c("EBNA1BP2", "Hep1", -0.598679596)
ls_v[[31]] <- c("EBNA1BP2", "Hep2", -0.82402321)
ls_v[[32]] <- c("EBNA1BP2", "Hep6", -0.947992183)

ls_v[[33]] <- c("MFSD2A", "Hep1", -0.547328703)
ls_v[[34]] <- c("MFSD2A", "Hep2", -0.837356722)
ls_v[[35]] <- c("MFSD2A", "Hep3", -0.455668688)
ls_v[[36]] <- c("MFSD2A", "Hep4", -0.849627047)

ls_v[[37]] <- c("RPL11", "Hep1", 0.453403238)
ls_v[[38]] <- c("RPL11", "Hep2", 0.401444723)
ls_v[[39]] <- c("RPL11", "Hep3", 1.578870637)
ls_v[[40]] <- c("RPL11", "Hep4", 1.935718427)
ls_v[[41]] <- c("RPL11", "Hep7", 0.62332579)

ls_v[[42]] <- c("ARHGEF10L", "Hep1", -0.442696344)
ls_v[[43]] <- c("ARHGEF10L", "Hep2", -0.606165831)
ls_v[[44]] <- c("ARHGEF10L", "Hep3", -0.363730134)
ls_v[[45]] <- c("ARHGEF10L", "Hep4", -1.10240336)
ls_v[[46]] <- c("ARHGEF10L", "Hep6", -0.678063741)


df <- data.frame(t(data.frame(ls_v)))
rownames(df) <-NULL
colnames(df) <- c('gene', 'subtype', 'Log2FC')
df$Log2FC <- as.numeric(df$Log2FC)



library(ggplot2)
ls_gene <- c("CDA", "DDI2", "IFNLR1", "MAN1C1", "PEX14", "PHACTR4", "ERRF",
             "PER3", "SCP2", "EBNA1BP2", "MFSD2A", "RPL11", "ARHGEF10L")

df$gene <- factor(df$gene, levels = ls_gene)

p <- ggplot(df, aes(x = gene, y = subtype, size = abs(Log2FC), color = factor(ifelse(Log2FC > 0, "UP", "DOWN")))) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(4, 12)) + 

  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue")) + 
  theme(
    panel.background = element_rect(fill = "gray98"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = 'bold'),  
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),  
    axis.title.y = element_blank(),   
    plot.title = element_blank(),   
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 15, face = 'bold', color = 'black')
  ) +
  labs(
    color = "Log2FC",
    size = "Absolute value"    
  )+
  guides(
    color = guide_legend(order = 1),     
    size = guide_legend(order = 2)
  ); p

png(file.path(dir_result, paste0("5.2.1.Upset.gene.png")), width = 4000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()
