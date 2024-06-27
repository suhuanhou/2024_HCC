################################################################################
##### setting
# conda activate METAFlux 
# R

library(METAFlux)
library(dplyr)
library(Seurat)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "101.bulk")
dir_dataset = file.path(dir_main, "35.METAFlux/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "35.METAFlux/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result); 

bool_save = TRUE
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
flux_scRNA=compute_sc_flux(num_cell = length(unique(sce$subtype)), fraction =subtype_ratio, fluxscore=scores, medium = human_blood)


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
data("human_gem")  

#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux)){
    activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))
#heatmap 
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
p <- pheatmap::pheatmap(all_pathway_score,cluster_cols = F,color = rev(mapal),scale = "row")

png("35.pathway.all.png", width = 4000, height = 3500, res = 400, bg = "transparent"); print(p); dev.off()


