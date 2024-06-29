# Cell proportion and survival analysis
################################################################################
##### setting
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)


rm(list = ls());
dir_main = dirname(rstudioapi::getActiveDocumentContext()$path); set.seed(100)
# dir_main = "/share/home/shh/HCC_snRNA-seq";

dir_in = file.path(dirname(dir_main), "LIHC")
dir_dataset = file.path(dir_main, "1.bisque/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_dataset_raw = file.path(dir_main, "1.bisque/dataset_raw"); if(!dir.exists(dir_dataset_raw)) dir.create(dir_dataset_raw, recursive = TRUE)
dir_result = file.path(dir_main, "1.bisque/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset)
################################################################################
##### read data 
#### bulk RNA-seq
df_bulk <- readRDS(file.path(dir_dataset, "PMID31585088_bulk.rds"))
clinical <- readxl::read_excel(file.path(dir_dataset_raw, "PMID31585088/supplemental/1-s2.0-S0092867419310037-mmc1.xlsx"), sheet = "1. Clinical and molecular data")
ID_sample <- colnames(df_bulk)
colnames(clinical)[1] <- "submitter_id" 
clinical <- clinical[clinical$submitter_id %in%ID_sample,]



# colnames(clinical)
meta = data.frame(clinical)
meta=meta[,colnames(meta) %in% c('submitter_id', "Age",
                                 "Gender", "Survial...1..dead..0..alive.",
                                 "Overall.survial..month.", "Recurrence.free.survival..month.")]
# rownames(meta) <- gsub("-", ".", meta$submitter_id)
rownames(meta) <- meta$submitter_id
meta=meta[match(colnames(df_bulk), rownames(meta)),]

colnames(meta) <- c('submitter_id','gender','Age', 'RFS', 'days_to_death', 'vital_status')



meta$time <- meta$days_to_death

# meta$event=ifelse(meta$vital_status=='Alive',0,1)
meta$event <- meta$vital_status
table(meta$event)


meta$Age[is.na(meta$Age)] <- 0
meta$Age=as.numeric(meta$Age)
meta$age_group=ifelse(meta$Age>median(meta$Age),'older','younger')
table(meta$age_group)
table(meta$gender)



cell_ratio <- readRDS(file.path(dir_dataset, "PMID31585088_bisque.rds"))
re <- t(cell_ratio)
# rownames(re) <- gsub("-", ".", substr(rownames(re),1,12))

meta <- cbind(meta, re[rownames(meta), ])


################################################################################
##### visualization
cellRatio = colnames(meta)[10:19]
splots <- lapply(cellRatio, function(g){
  if (median(meta[, g]) == 0) {
    return()  
  }
  
  meta$cell=ifelse(meta[,g]>median(meta[,g]),'high','low')
  sfit1=survfit(Surv(time, event)~cell, data = meta)
  ggsurvplot(sfit1, pval = TRUE, data = meta, risk.table = FALSE, 
             title = paste0(g))
}) 

splots <- Filter(Negate(is.null), splots)
dev.off()
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 5, nrow = 2, risk.table.height = 0.4)


# dev.copy(png, file.path(dir_result, "cell_ratio.survival.png"), width = 6000, height = 5000, res = 400)


## Hep
plot_subtype <- function(subtype){
  # subtype = 'Hep0'
  meta$cell=ifelse(meta[,subtype]>median(meta[,subtype]),'high','low')
  sfit1=survfit(Surv(time, event)~cell, data = meta)
  p_plot <- ggsurvplot(sfit1,  data = meta, risk.table = FALSE,
             title = paste0(subtype),
             xlab = "Time(month)", 
             font.main = c(25), 
             font.x = c(20),
             font.y = c(20),
             font.tickslab = c(16),  
             font.legend = c(20),             
             size = 1.5,
             pval = TRUE, pval.size = 8, pval.coord = c(0,0.05),
             legend = c(0.2,0.3), legend.title = 'Cell ratio', legend.labs = c('High', 'Low'), 
             # surv.median.line = 'hv',
             ggtheme = theme(
               plot.title = element_text(hjust = 0.5),
               panel.background = element_blank(),  
               axis.line = element_line(color = "black"),  
               panel.grid.major = element_blank(),  
               panel.grid.minor = element_blank(), 

               )
             )

  return(p_plot)
}

p1 <- plot_subtype('Hep0'); png(file.path(dir_result, paste0("bisque_survival_PMID31585088_Hep0.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p1); dev.off()
p2 <- plot_subtype('Hep2'); png(file.path(dir_result, paste0("bisque_survival_PMID31585088_Hep2.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p2); dev.off()
p3 <- plot_subtype('Hep9'); png(file.path(dir_result, paste0("bisque_survival_PMID31585088_Hep9.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p3); dev.off()

