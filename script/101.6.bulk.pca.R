# 20240604 
################################################################################
##### setting
library(openxlsx)
library(pheatmap)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "0.data_bulk_RNA-seq/result/SEQ-DX-PXH-202402021501_result_final")
dir_dataset = file.path(dir_main, "101.bulk/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "101.bulk/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result); 
################################################################################
##### read data
df_exp <- readxl::read_excel(file.path(dir_in, "3.Expression/all_samples_expression_RPM.xls"), sheet=1)
df_exp <- data.frame(df_exp)
df_tmp <- cbind(df_exp[,1], df_exp[,57:166], df_exp[,2:56])
rownames(df_tmp) <- df_tmp[,1]
df_tmp <- df_tmp[,-1]
df_exp <- df_tmp


# meta
df_meta <- readxl::read_excel(file.path(dir_main, "0.data_bulk_RNA-seq/bulk_RNA-seq item.xlsx"), sheet=1)
colnames(df_meta)[1] <- "ID"
df_meta <- data.frame(df_meta)
df_meta <- df_meta[,colnames(df_meta) %in% c('ID', 'Age', 'Sex', 'HBV')]
df_meta[nrow(df_meta) + 1,1] <- "ID01"
df_meta[nrow(df_meta) + 1,1] <- "ID11"
df_meta <- df_meta[order(df_meta$ID),]


df_meta[df_meta$ID == "ID01",] <- c('ID01', '50', '+', '男')
df_meta[df_meta$ID == "ID02",] <- c('ID02', '50', '+', '男')
df_meta[df_meta$ID == "ID05",] <- c('ID02', '50', '+', '男')
df_meta[df_meta$ID == "ID11",] <- c('ID11', '50', '+', '男')


df_meta$Sex <- ifelse(df_meta$Sex == "男", "M", ifelse(df_meta$Sex == "女", "F", df_meta$Sex))
df_meta$Age <- as.numeric(df_meta$Age)
df_meta$Age <- cut(df_meta$Age, breaks = c(0, 50, 60, 70, Inf),
                   labels = c("<50", "50-60", "60-70", ">70"),
                   right = FALSE)


Tissue <- ifelse(grepl("bulk_P", colnames(df_exp)), "Paracancer",
                ifelse(grepl("bulk_C", colnames(df_exp)), "Cancer",
                       ifelse(grepl("bulk_N", colnames(df_exp)), "Normal", NA)))



df_anno_col = data.frame("Tissue" = Tissue,
                     "Age" = rep(df_meta$Age, 3),
                     "Sex" = rep(df_meta$Sex, 3),
                     "HBV" = rep(df_meta$HBV, 3)
                     )
rownames(df_anno_col) <- colnames(df_exp) 

df_anno_row = data.frame("Tissue" = Tissue)
rownames(df_anno_row) <- colnames(df_exp) 

M = cor(log(df_exp+1))
M[M < 0.8] <- 0.8


dev.off()

p <- pheatmap(M,
              color=colorRampPalette(colors = c("white","#083471"))(30),
         annotation_col = df_anno_col,
         # annotation_colors = color_col
         
         annotation_colors = list(Tissue = c("Normal" = "#E4F99D", "Paracancer" = "#E39C18", "Cancer" = "#FF6207"),
                                  Age = c("<50" = "#2F9A66", "50-60" = "#50B885", "60-70" = "#70CDA5", "60-70" = "#8FCDC2", ">70" = "#C2E1E1"),
                                  Sex = c("F" = "#F356D0", "M" = "#0797FE"),
                                  HBV = c("+" = "#BD0F67", "-" = "#A2CB65")
                                  ),
         annotation_row = df_anno_row, annotation_names_row = F,
         cluster_rows = F,cluster_cols = F,
         show_rownames = F, show_colnames = F,

); p

png("101.bulk.correlation.png", width = 4000, height = 4000, res = 400, bg = "transparent"); print(p); dev.off()

