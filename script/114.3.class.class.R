# conda activate class  && R
################################################################################
##### setting
library(ConsensusClusterPlus)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)

n_cluster = 3
################################################################################
##### load data
df_exp <- readRDS(file.path(dir_dataset, paste0("df_ComBat.rds"))) 
gene_feature <- readRDS(file.path(dir_dataset, paste0("gene_feature.rds"))) 
df_exp <- df_exp[rownames(df_exp) %in% gene_feature,]

selected_columns <- grep("PMID31585088|TCGA", colnames(df_exp), value = TRUE)
df_exp2 <- df_exp[, selected_columns]
df_exp2 <- as.matrix(df_exp2)


results <- ConsensusClusterPlus(df_exp2, 
                            maxK=10, 
                            reps=100, 
                            pItem=0.8, 
                            pFeature=1,
                            title="CCP", 
                            clusterAlg="km", 
                            distance="euclidean", 
                            seed=100, 
                            plot="png", 
                            writeTable=TRUE)



df_Class <- data.frame(
  Sample = colnames(df_exp2),
  Class = paste0("Class_", results[[n_cluster]][['consensusClass']])
  )

head(df_Class)
table(df_Class$Class)

saveRDS(df_Class, file.path(dir_result, paste0('df_Class.rds')))