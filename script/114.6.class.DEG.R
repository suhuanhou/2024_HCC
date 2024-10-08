################################################################################
##### setting
library(sva)
library(limma)
library(ggplot2)
library(dplyr)
library(openxlsx)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
ls_color <- c("Class_1" = "#EA8379", "Class_2" = "#7DAEE0", "Class_3" = "#B395BD")
################################################################################
##### load data
bulk_data <- list(
  readRDS(file.path(dir_dataset, "PMID31585088_bulk.rds")),
  readRDS(file.path(dir_dataset, "TCGA_bulk.rds"))
)

row_names_list <- lapply(bulk_data, rownames)
common_rownames <- Reduce(intersect, row_names_list)


subset_data <- lapply(bulk_data, function(df) df[common_rownames, , drop = FALSE])
df_merge <- do.call(cbind, subset_data)
df_merge <- log2(df_merge + 1)  


batch <- sapply(strsplit(colnames(df_merge), "_"), `[`, 1) 

df_ComBat <- ComBat(df_merge, batch = batch)
df_ComBat <- data.frame(df_ComBat)

saveRDS(df_ComBat, file.path(dir_dataset, paste0("df_trai.rds"))) 

################################################################################
##### DEG
df_Class <- readRDS(file.path(dir_result, paste0('df_Class.rds')))
df_Class <- df_Class[match(colnames(df_ComBat), df_Class$Sample),]

save_ggplot <- function(gg, prefix, width = 8, height = 8, dpi = 600) {
  ggsave(paste0(prefix, ".png"), plot = gg, width = width, height = height, dpi = dpi)
}

thre_p = 0.05
thre_FC = 1.5


DEG_analysis <- function(group){
  # group = 'Class_1'
  group_list <- df_Class$Class
  group_list <- ifelse(group_list == group, "p", "n")

  design <- model.matrix(~0+factor(group_list)) 
  colnames(design) <- gsub("factor\\(group_list\\)", "", colnames(design))

  contrast.matrix<-makeContrasts("p-n", levels = design)
  fit <- lmFit(df_ComBat,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  tempOutput  <-  topTable(fit2, coef=1, n=Inf)
  df_cont <- na.omit(tempOutput) 
  df_cont$gene <- rownames(df_cont)
  
  addWorksheet(wb, group)
  writeData(wb, sheet = group, x = df_cont, startCol = 1, startRow = 1)
  
  df_cont_DEG <- df_cont[df_cont$P.Value < thre_p & abs(df_cont$logFC) > thre_FC, ]
  df_cont_DEG <- df_cont_DEG[order(-df_cont_DEG$logFC), ]

  addWorksheet(wb, paste0(group, '_DEG'))
  writeData(wb, sheet = paste0(group, '_DEG'), x = df_cont_DEG, startCol = 1, startRow = 1)
  
  ### visualization
  tmp_res <- df_cont %>%
    mutate(change = case_when(P.Value < thre_p & logFC > thre_FC ~ "Up",
                              P.Value < thre_p & logFC < -thre_FC ~ "Down",
                              TRUE ~ "Stable"))
  
  colors = c("darkgreen","grey","darkred")
  names(colors) = c("Down","Stable","Up")
  tmp_y = max(-log10(tmp_res$P.Value))
  if (tmp_y < 10){
    y_height = 10
  } else if (tmp_y < 200){
    y_height = ceiling(tmp_y / 10) * 10
  } else {
    y_height = ceiling(tmp_y / 10) * 10 + 50
  }
  tmp_x = max(abs(tmp_res$logFC))
  if (tmp_x < 10){
    x_length = 10
  }else {
    x_length = 15
  }
  
  tmp_volcano <- ggplot(tmp_res, aes(x=logFC, y=-log10(P.Value),
                                     color=change)) +
    geom_point(alpha=0.4,size=3)+
    scale_color_manual(values=colors) +
    theme_bw() +
    ylim(c(0,y_height))+
    xlim(c(-x_length,x_length)) +
    geom_vline(xintercept = c(-thre_FC,thre_FC), lty=4, col="black", lwd=0.8)+
    geom_hline(yintercept = -log10(thre_p), lty=4, col="black", lwd=0.8) +
    labs(x = expression(log[2]~'FC'),
         y=expression(-log10(P.Value)),
         title=paste0(group)) +
    theme(axis.title.x = element_text(colour="black",size = 22),
          axis.title.y = element_text(colour="black",size = 22),
          axis.text.x = element_text(colour="black",size = 20),
          axis.text.y = element_text(colour="black",size = 20),
          plot.title = element_text(size=24),
          legend.text = element_text(color="black", size=20))+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="none", 
          legend.title = element_blank())
  save_ggplot(tmp_volcano, file.path(dir_result, paste0("114.volcano.", group)), width=5, height=5)
  
  return(wb)
}


wb <- createWorkbook()
wb <- DEG_analysis('Class_1')
wb <- DEG_analysis('Class_2')
wb <- DEG_analysis('Class_3')

saveWorkbook(wb, "Class_DEG.xlsx", overwrite = TRUE)

