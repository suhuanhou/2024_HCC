# cell sycle
################################################################################
##### setting
library(ggplot2)
library(Seurat);


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "4.merge/dataset")
dir_dataset = file.path(dir_main, "19.cell_cycle/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "19.cell_cycle/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset);

################################################################################
##### 细胞周期计算
# 提取sce中细胞周期的marker基因
sce.cycle <- readRDS(file.path(dir_in, "sce_all.rds"))
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(sce.cycle))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes, match=rownames(sce.cycle)) 
sce.cycle <- CellCycleScoring(sce.cycle, g2m.features=g2m_genes, s.features=s_genes)
table(sce.cycle$Phase)



df_cycle <- data.frame(sce.cycle$Phase)
colnames(df_cycle) <- c('cell_cycle')
head(df_cycle)


# saveRDS(df_cycle, file.path(dir_dataset, "df_cycle.rds"))
# df_cycle <- readRDS(file.path(dir_dataset, "df_cycle.rds"))


################################################################################
##### 细胞周期展示
sce_hep <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset", "sce_Hep.rds"))
# df_cycle <- readRDS(file.path(dir_dataset, "df_cycle.rds"))
sce_hep$cell_cycle = df_cycle[colnames(sce_hep),]
table(sce_hep$cell_cycle)
table(sce_hep$subtype)


cell_cycle_table <- table(sce_hep$subtype, sce_hep$cell_cycle)
cell_cycle_prop <- prop.table(cell_cycle_table, margin = 1)

ggplot(data = as.data.frame(cell_cycle_prop), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "subtype", y = "Proportion", fill = "Cell Cycle") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

