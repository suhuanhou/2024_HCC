################################################################################
##### setting
# conda activate METAFlux 
# R

library(METAFlux)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "101.bulk")
dir_dataset = file.path(dir_main, "35.METAFlux/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "35.METAFlux/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
dir_reaction = file.path(dir_result, "reaction"); if(!dir.exists(dir_reaction)) dir.create(dir_reaction, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result); 
################################################################################
##### bulk data
# data("bulk_test_example")
data("human_blood")  # medium file for human derived samples
bulk_test_example <- readRDS(file.path(dir_in, "dataset/df_RPM.rds"))

#Calculate mras for human sample data
scores<-calculate_reaction_score(bulk_test_example)

# Calculate flux for human sample data
flux_bulk <-compute_flux(mras=scores,medium=human_blood) 

# optional: flux scores cubic root normalization
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux_bulk = cbrt(flux_bulk)
saveRDS(flux_bulk, 'flux_bulk.rds')


##### scRNA data
sce <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset/sce_Hep.rds"))


# https://github.com/KChen-lab/METAFlux/issues/5
generate_boots <- function(celltype, n) {
  dt <- data.frame(cluster = celltype, id = 1:length(celltype))
  index <- do.call(cbind, sapply(1:n, function(x) {
    splits <- dt %>%
      group_by(cluster) %>%
      sample_n(dplyr::n(), replace = TRUE) %>%
      ungroup() %>%
      dplyr::select("id")
  }))
  return(index)
}

get_ave_exp <- function(i, myseurat, samples,myident) {
  meta.data=myseurat@meta.data[samples[,i],]
  sample <-myseurat@assays$RNA@counts[,samples[,i]]
  name <- colnames(sample)
  for (j in 1:length(name)) {
    name[j] <- paste0(name[j], "_", j)
  }
  colnames(sample) <- name
  rownames(meta.data) <- name
  SeuratObject<-suppressWarnings(
    CreateSeuratObject(count=sample, meta.data = meta.data))
  SeuratObject<-NormalizeData(SeuratObject,verbose = TRUE)
  ave<-GetAssayData(AverageExpression(SeuratObject,group.by = myident,return.seurat = T), assay = "RNA") %>% as.data.frame()
  return(ave)
}

edit_calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
  set.seed(seed)
  samples=generate_boots(myseurat@meta.data[,myident],n_bootstrap)
  exp <- lapply(1:n_bootstrap,get_ave_exp,myseurat,samples,myident)
  exp <- do.call(cbind, exp)
  return(exp)
}

# mean_exp = edit_calculate_avg_exp(myseurat = sc, myident = 'celltype', n_bootstrap = 50, seed = 1)
mean_exp = edit_calculate_avg_exp(myseurat = sce, myident = 'subtype', n_bootstrap=50, seed=1)  # 

#calculate metabolic reaction scores
scores <- calculate_reaction_score(data=mean_exp)

#calculate the fractions of celltype/clusters
subtype_ratio <- table(sce$subtype)/nrow(sce@meta.data) %>% as.vector()


#num_cell=number of cell types/clusters, here we have 4 cell types, thus input is 4. Make sure the sum of cell percentage is 1.The order of fraction should be the same as that of "Mean_exp" data
flux_scRNA = compute_sc_flux(num_cell = length(unique(sce$subtype)), fraction =subtype_ratio, fluxscore=scores, medium = human_blood)


#optional: flux scores normalization
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux_scRNA = cbrt(flux_scRNA)
saveRDS(flux_scRNA, 'flux_scRNA.rds')


################################################################################
##### visualization
flux_bulk <- readRDS('flux_bulk.rds')
flux_scRNA <- readRDS('flux_scRNA.rds')

data("human_gem")  # 通路涉及的化学反应式
ls_color <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)



#### bulk
#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux_bulk)){
    activity_score[d]<-mean(abs(flux_bulk[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

score_all<-as.data.frame(do.call(rbind,pathway_score))
colnames(score_all) <- colnames(flux_bulk)
# score_all[1:4, 1:4]


group_list <- ifelse(grepl("bulk_P", colnames(score_all)), "Paracancer",
                     ifelse(grepl("bulk_C", colnames(score_all)), "Cancer",
                            ifelse(grepl("bulk_N", colnames(score_all)), "Normal", NA))) %>% 
  factor(.,levels = c("Normal", "Paracancer", "Cancer"))


mean_bulk <- t(apply(score_all, 1, function(x) tapply(x, group_list, mean)))
mean_bulk <- data.frame(mean_bulk)

# foldchange
mean_bulk$CvsP <- mean_bulk$Cancer / mean_bulk$Normal

# order
mean_bulk <- mean_bulk[order(mean_bulk$CvsP, decreasing = TRUE), ]


filter_row <- c(rownames(mean_bulk)[1:20])
score_all2 <- score_all[rownames(score_all) %in% filter_row,]
ls_reaction <- rownames(score_all2)
saveRDS(ls_reaction, 'ls_reaction.rds')

### heatmap 
p <- pheatmap::pheatmap(score_all2, cluster_cols = F, color = rev(ls_color), scale = "row")
png("35.bulk.heatmap.png", width = 5000, height = 5000, res = 400, bg = "transparent"); print(p); dev.off()


### boxplot
df_bulk_boxplot <- score_all2 %>% t() %>% data.frame()
df_bulk_boxplot$group <- as.character(group_list)
df_long <- reshape2::melt(df_bulk_boxplot, id.vars = "group", variable.name = "reaction", value.name = "value")
my_comparisons = list(c("Normal", "Paracancer"), c("Paracancer", "Cancer"), c("Normal", "Cancer"))

p <- ggboxplot(df_long, x = "group", y = "value", facet.by = "reaction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + stat_compare_means(comparisons = my_comparisons)
png("35.bulk.boxplot.png", width = 5000, height = 5000, res = 400, bg = "transparent"); print(p); dev.off()



#### flux_scRNA
df_scRNA <- flux_scRNA
df_scRNA$rownames <- rownames(df_scRNA)
df_scRNA$subtype <- ""
df_scRNA$reaction <- ""


df_scRNA <- subset(df_scRNA, !grepl("internal_medium|external_medium", rownames(df_scRNA)))
df_scRNA <- subset(df_scRNA, grepl("celltype 1 |celltype 3 |celltype 10 ", rownames(df_scRNA)))


for (i in 1:nrow(df_scRNA)) {
  if (startsWith(df_scRNA$rownames[i], "celltype")) {
    split_name <- strsplit(df_scRNA$rownames[i], " ")[[1]]
    df_scRNA$subtype[i] <- paste(split_name[1], split_name[2], sep = " ")
    df_scRNA$reaction[i] <- split_name[3]
  } else if (startsWith(df_scRNA$rownames[i], "external_medium")) {
    split_name <- strsplit(df_scRNA$rownames[i], " ")[[1]]
    df_scRNA$subtype[i] <- split_name[1]
    df_scRNA$reaction[i] <- split_name[2]
  }
}


df_scRNA$rownames <- NULL


replacement <- c("celltype 1" = "Hep0", "celltype 2" = "Hep1", 
                 "celltype 3" = "Hep2", "celltype 4" = "Hep3", 
                 "celltype 5" = "Hep4", "celltype 6" = "Hep5", 
                 "celltype 7" = "Hep6", "celltype 8" = "Hep7", 
                 "celltype 9" = "Hep8", "celltype 10" = "Hep9", "Hep00" = "Hep9")

df_scRNA$subtype <- stringr::str_replace_all(df_scRNA$subtype, replacement)
# table(df_scRNA$subtype)

match_index <- match(df_scRNA$reaction, human_gem$ID)
df_scRNA$reaction <- human_gem$SUBSYSTEM[match_index]


df_scRNA <- df_scRNA %>%
  group_by(subtype, reaction) %>%
  summarise_all(.funs = mean, na.rm = TRUE) %>% 
  data.frame()




ls_reaction <- readRDS('ls_reaction.rds')
i = 0  # 计数
for (choo_reaction in ls_reaction){
  # choo_reaction = 'Bile acid recycling'
  i = i + 1
  choo_reaction2 = gsub("/", " ", choo_reaction) 
  df_choo <- subset(df_scRNA, reaction == choo_reaction)
  my_comparisons = list(c("Hep0", "Hep2"), c("Hep0", "Hep9"))
  

  data <- df_choo[, c("subtype", paste0("V", 1:50))]
  data_long <- reshape2::melt(data, id.vars = "subtype")
  
  p <- ggboxplot(data_long, x = "subtype", y = "value", fill = "subtype", ylab = "Flux") +
    stat_compare_means(comparisons = my_comparisons) +
    geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.5, size = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("Hep0" = "#3EA94B", "Hep2" = "#E51F17", "Hep9" = "#C20576")) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.position = "none", 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),  
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 10),   
          
          
    )+
    ggtitle(choo_reaction) 
  
  
  png(file.path(dir_reaction, paste0("35.subtype.", i, ".", choo_reaction2, ".png")), width = 1500, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()
}


### pathway in bulk RNA-seq
# unique(df_long$reaction)
# choo_reaction = 'Nucleotide.metabolism'

# 1\5\15\16
i = 1
choo_reaction = unique(df_long$reaction)[i]
choo_reaction2 = gsub("/", " ", choo_reaction)
df_long2 <- df_long[df_long$reaction == choo_reaction,]
my_comparisons = list(c("Normal", "Paracancer"), c("Paracancer", "Cancer"), c("Normal", "Cancer"))


p <- ggboxplot(df_long2, x = "group", y = "value", fill = "group", ylab = "Flux") +
  stat_compare_means(comparisons = my_comparisons) +
  geom_point(position = position_jitter(width = 0.2), color = "black", alpha = 0.5, size = 1, show.legend = FALSE) +
  scale_fill_manual(values = c("Normal" = "#E5FA9C", "Paracancer" = "#E59E0F", "Cancer" = "#FF6300")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),  
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 10),   
  )+
  ggtitle(choo_reaction) 


png(file.path(paste0("35.bulk.", i, ".", choo_reaction2, ".png")), width = 1500, height = 1500, res = 400, bg = "transparent"); print(p); dev.off()




