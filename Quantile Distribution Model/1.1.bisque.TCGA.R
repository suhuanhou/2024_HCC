# Cell proportion and survival analysis
################################################################################
##### setting
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)

rm(list = ls());
dir_main = dirname(rstudioapi::getActiveDocumentContext()$path); 
# dir_main = "/share/home/shh/HCC_snRNA-seq";

dir_in = "Path to TCGA/LIHC"
dir_dataset = file.path(dir_main, "1.bisque/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "1.bisque/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_dataset); set.seed(100)
################################################################################
##### read data 
#### bulk RNA-seq
df_bulk <- read.table(file.path(dir_in, "result/LIHC.tpm.tumor.txt"), sep = '\t', header = TRUE)
rownames(df_bulk) <- df_bulk[,1]
df_bulk <- df_bulk[,-1]
colnames(df_bulk) <- substr(colnames(df_bulk),1,12)
exprSet <- as.matrix(df_bulk)


clinical <- readRDS(file.path(dir_in, "dataset/LIHC_clinical.rds"))


patient_df <- read.table(file.path(dir_in, "result/LIHC.patient.txt"), sep = '\t', header = TRUE)
ID_tumor <- subset(patient_df, tissue_type == "tumor")$ID
clinical <- clinical[clinical$submitter_id %in%ID_tumor,]


# colnames(clinical)
meta = clinical
meta=meta[,colnames(meta) %in% c('submitter_id', "days_to_last_follow_up",
                                 "gender", "vital_status",
                                 "age_at_index", "days_to_death")]
rownames(meta) <- gsub("-", ".", meta$submitter_id)
meta=meta[match(colnames(df_bulk), rownames(meta)),]



meta$days_to_death[is.na(meta$days_to_death)] <- 0   
meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] <- 0
meta$days = as.numeric(meta$days_to_last_follow_up)+as.numeric(meta$days_to_death)
meta$time = round(meta$days/30,2) 


meta$event=ifelse(meta$vital_status=='Alive',0,1)
table(meta$event)


meta$age_at_index[is.na(meta$age_at_index)] <- 0
meta$age_at_index=as.numeric(meta$age_at_index)
meta$age_group=ifelse(meta$age_at_index>median(meta$age_at_index),'older','younger')
table(meta$age_group)
table(meta$gender)


cell_ratio <- readRDS(file.path(dir_dataset, "TCGA_bisque.rds"))
re <- t(cell_ratio)
rownames(re) <- gsub("-", ".", substr(rownames(re),1,12))
meta <- cbind(meta, re[rownames(meta), ])



################################################################################
##### visualization
cellRatio = colnames(meta)[11:20]
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

p1 <- plot_subtype('Hep0'); png(file.path(dir_result, paste0("bisque_survival_TCGA_Hep0.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p1); dev.off()
p2 <- plot_subtype('Hep2'); png(file.path(dir_result, paste0("bisque_survival_TCGA_Hep2.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p2); dev.off()
p3 <- plot_subtype('Hep9'); png(file.path(dir_result, paste0("bisque_survival_TCGA_Hep9.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p3); dev.off()

