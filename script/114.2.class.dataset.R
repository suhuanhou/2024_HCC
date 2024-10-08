################################################################################
##### setting
library(dplyr)
library(openxlsx)
library(ggplot2)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "0.data_article")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)
################################################################################
##### load data
#### PMID31585088
if (F){
  bulk_name = 'PMID31585088'
  df_exp <- read.table(file.path(dir_in, "PMID31585088/HCC_UQ_FPKM.tsv"), sep = '\t', header = TRUE)
  df_exp$protein <- substr(df_exp$protein,1,15)
  df_exp <- df_exp[,-ncol(df_exp)]  
  
  
  library(readxl)
  df_sample <- read_excel(file.path(dir_in, "PMID31585088/About the RNA and protein Identifier match.xlsx"), sheet = 1)
  df_sample <- data.frame(df_sample)
  df_sample <- df_sample[,3:4]
  colnames(df_sample) <- c("ID_Tumor","ID_Normal")
  df_sample <- df_sample[-1,]
  

  ls_col <- colnames(df_exp)
  ls_col <- sub("^X", "", ls_col)
  
  for (i in seq_along(ls_col)) {
    if (ls_col[i] %in% df_sample$ID_Tumor) {
      ls_col[i] <- paste0("T", ls_col[i])
    } else if (ls_col[i] %in% df_sample$ID_Normal) {
      ls_col[i] <- paste0("N", ls_col[i])
    }
  }
  
  colnames(df_exp) <- ls_col
  

  df_exp <- df_exp[, grep("^p|^T|^N", colnames(df_exp), value = TRUE)]
  rownames(df_exp) <- df_exp[,1]
  df_exp$protein <- NULL
  
  
  FPKM2TPM <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  
  max(df_exp)  
  df_exp <- apply(df_exp,2,FPKM2TPM)
  df_exp <- data.frame(df_exp)
  
  #### ID
  library(tidyverse)  
  library(clusterProfiler)
  library(org.Hs.eg.db) 
  
  ls_rownames <- rownames(df_exp) 
  df_exp <- df_exp %>% rownames_to_column(var = "gene")
  name <- bitr(ls_rownames,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
  
  
  df_exp$symbol <- df_exp$gene
  for (i in 1:nrow(df_exp)) {
    matching_row <- which(name$ENSEMBL == df_exp$gene[i])
    if (length(matching_row) > 0) {
      df_exp$symbol[i] <- name$SYMBOL[matching_row]
    }
  }
  
  
  ids <- data.frame(
    probe_id = paste0(df_exp[,1], seq_len(nrow(df_exp))),
    symbol = df_exp$symbol
  )
  
  rownames(df_exp) <- paste0(df_exp[,1], seq_len(nrow(df_exp)))
  df_exp$gene <- NULL
  df_exp$symbol <- NULL
  
  
  ids = ids[match(rownames(df_exp),ids$probe_id),] # ID替换和重新排序
  jimmy <- function(df_exp,ids){
    tmp = by(df_exp,
             ids$symbol,
             function(x) rownames(x)[which.max(rowMeans(x))])
    probes = as.character(tmp)
    print(dim(df_exp)) 
    df_exp = df_exp[rownames(df_exp) %in% probes,]
    
    print(dim(df_exp))
    rownames(df_exp) = ids[match(rownames(df_exp),ids$probe_id),2]
    return(df_exp)
  }
  
  new_df_exp <- jimmy(df_exp, ids)
  df_exp <- new_df_exp
  df_bulk <- df_exp[, grep("^T", colnames(df_exp), value = TRUE)]  # HCC样本筛选

  
  ## clinical
  df_clin <- read_excel(file.path(dir_in, "PMID31585088/supplemental/1-s2.0-S0092867419310037-mmc1.xlsx"), sheet = "1. Clinical and molecular data")
  df_clin <- data.frame(df_clin)

  df_bulk <- df_bulk[,df_clin$Tumor..T..sample.ID]
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  df_clin$Tumor..T..sample.ID <- paste0(bulk_name, '_', df_clin$Tumor..T..sample.ID)
  
  # colnames(df_clin)
  choo_col <- c('Tumor..T..sample.ID', 
                'Gender', 
                'Age', 
                'Overall.survial..month.', 
                'Survial...1..dead..0..alive.',
                'ARID1A.mutation',                
                'AXIN1.mutation',                
                'CTNNB1.mutation',
                'KEAP1.mutation',                
                'KMT2C.mutation',                   
                'TP53.mutation',
                'TSC2.mutation'                
                )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'Gender','Age',  'OS', 'Survial',
                          'ARID1A_mutation',                
                          'AXIN1_mutation',                
                          'CTNNB1_mutation',
                          'KEAP1_mutation',                
                          'KMT2C_mutation',                   
                          'TP53_mutation',
                          'TSC2_mutation'  )  # 1,dead
  
  df_clin2[, 6:12] <- lapply(df_clin2[, 6:12], function(x) ifelse(!is.na(x) & x != "", 1, 0))
  df_clin2$dataset <- bulk_name
  
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
}

#### TCGA
if (F){
  bulk_name = 'TCGA'
  df_exp <- read.table(file.path(dir_in, "TCGA/LIHC.tpm.tumor.txt"), sep = '\t', header = TRUE)
  rownames(df_exp) <- df_exp[,1]
  df_exp <- df_exp[,-1]
  colnames(df_exp) <- substr(colnames(df_exp),1,12)
  colnames(df_exp) <- gsub("\\.", "_", colnames(df_exp))
  df_bulk <- df_exp
  
  df_clin <- read.delim(file.path(dir_in, "TCGA/LIHC.clinical.txt"), header = TRUE)
  df_clin$days_to_death[is.na(df_clin$days_to_death)] <- 0   #缺失值标记为0
  df_clin$days_to_last_follow_up[is.na(df_clin$days_to_last_follow_up)] <- 0
  df_clin$days = as.numeric(df_clin$days_to_last_follow_up)+as.numeric(df_clin$days_to_death)
  df_clin$time = round(df_clin$days/30,2)  # 时间以月份记，保留两位小数
  
  df_clin$event=ifelse(df_clin$vital_status=='Alive',0,1)

  choo_col <- c('submitter_id',
                'gender', 
                'age_at_index',
                'time', 
                'event')
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'Gender','Age',  'OS', 'Survial')  # 1,dead
  
  df_clin2$ARID1A_mutation <- 0                
  df_clin2$AXIN1_mutation <- 0               
  df_clin2$CTNNB1_mutation <- 0
  df_clin2$KEAP1_mutation <- 0               
  df_clin2$KMT2C_mutation <- 0                  
  df_clin2$TP53_mutation <- 0
  df_clin2$TSC2_mutation <- 0
  
  df_clin2$dataset <- bulk_name
  df_clin2$Sample <- gsub("-", "_", df_clin2$Sample)

  
  ls_Mutate <- c('ARID1A', 'AXIN1', 'CTNNB1', 'KEAP1', 'KMT2C', 'TP53', 'TSC2')
  for (choo_gene in ls_Mutate) {
    # choo_gene <- 'ARID1A'
    mutation_col <- paste0(choo_gene, '_mutation')
    df_Mutate <- read.table(file.path(dir_dataset, paste0("TCGA_Mutate/", choo_gene, ".tsv")), sep = '\t', header = TRUE)
    df_Mutate$Patient.ID <- gsub("-", "_", df_Mutate$Patient.ID)
    df_clin2 <- df_clin2 %>% mutate(!!mutation_col := ifelse(Sample %in% df_Mutate$Patient.ID, 1, get(mutation_col)))
  }
  

  df_bulk <- df_bulk[,df_clin2$Sample]
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
}


#### LIRI-JP
if (F){
  bulk_name = 'LIRI'
  df_exp <- read.table(file.path(dir_in, "LIRI/LIRI.FPKM.txt"), sep = '\t', header = TRUE)
  rownames(df_exp) <- df_exp[,1]
  df_exp <- df_exp[,-1]
  df_bulk <- df_exp
  
  FPKM2TPM <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  
  max(df_bulk)  
  df_bulk <- apply(df_bulk,2,FPKM2TPM)
  df_bulk <- data.frame(df_bulk)
  

  df_clin <- read.table(file.path(dir_in, "LIRI/LIRI.clin_df.txt"), header = TRUE)
  choo_col <- c('icgc_donor_id',
                'donor_sex', 
                'donor_age_at_diagnosis',
                'month', 
                'donor_vital_status')
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'Gender','Age',  'OS', 'Survial')  # 1,dead
  

  df_clin2$ARID1A_mutation <- ''                
  df_clin2$AXIN1_mutation <- ''               
  df_clin2$CTNNB1_mutation <- ''
  df_clin2$KEAP1_mutation <- ''               
  df_clin2$KMT2C_mutation <- ''                  
  df_clin2$TP53_mutation <- ''
  df_clin2$TSC2_mutation <- ''
  
  df_clin2$dataset <- bulk_name
  df_clin2$Survial <- ifelse(df_clin2$Survial == "alive", 0, 1)
  
  df_bulk <- df_bulk[,df_clin2$Sample]
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  df_clin2$Sample <- paste0(bulk_name, '_', df_clin2$Sample)
  
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
}


#### GSE141198
if (F){
  bulk_name = 'GSE141198'
  df_exp <- read.table(file.path(dir_dataset, "GSE141198_TLCN.subset1.readcount.txt"), sep = '\t', header = TRUE)
  df_exp$GeneID <- substr(df_exp$GeneID, 1, 15)
  df_exp[1:4,1:4]
  df_gene <- readRDS(file.path(dir_main, 'tool/gene.lenth.symbol/Gene.ID.Length.rds'))
  df_gene[1:3,1:3]
  
  df_exp <- merge(df_exp, df_gene, by = "GeneID", all.x = TRUE)
  
  rownames(df_exp) <- paste0(df_exp[,1], seq_len(nrow(df_exp)))
  df_exp$GeneID <- NULL
  df_exp$Length <- NULL


  df_tmp <- aggregate(. ~ Symbol, data = df_exp, FUN = sum)
  rownames(df_tmp) <- df_tmp$Symbol
  df_tmp <- df_tmp[rowSums(df_tmp[,-1] != 0) > 0, ]
  

  df_gene2 <- aggregate(Length ~ Symbol, data = df_gene, FUN = max)
  head(df_gene2)
  
  df_tmp2 <- merge(df_tmp, df_gene2, by = "Symbol", all.x = TRUE)
  ls_Length <- df_tmp2$Length 
  df_tmp$Symbol <- NULL
  

  calculate_tpm <- function(counts, lengths) {
    rpk <- counts / (lengths / 1000)  # 每千碱基的读数
    total_rpk <- rowSums(rpk)
    tpm <- rpk / total_rpk * 1e6
    return(tpm)
  }
  
  df_bulk <- calculate_tpm(df_tmp, ls_Length)

  ## clinical
  df_clin <- readRDS(file.path(dir_dataset, paste0(bulk_name, "_clin0.rds")))

  df_bulk <- df_bulk[,df_clin$title]
  colnames(df_bulk) <- paste0(bulk_name, '_', df_clin$geo_accession)
  df_clin$geo_accession <- paste0(bulk_name, '_', df_clin$geo_accession)
  
  
  # colnames(df_clin)
  choo_col <- c('geo_accession',
                "os days:ch1",
                "os event:ch1"
  )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'OS', 'Survial')  # 1,dead
  
  df_clin2$OS = as.numeric(df_clin2$OS)
  df_clin2$OS <- round(df_clin2$OS/30,2)
  df_clin2$dataset <- bulk_name
  

  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds")))
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds")))
}


#### GSE141200
if (F){
  bulk_name = 'GSE141200'
  df_exp <- read.table(file.path(dir_dataset, "GSE141200_TLCN.subset2.readcount.txt"), sep = '\t', header = TRUE)
  colnames(df_exp)[1] <- 'GeneID'
  df_exp$GeneID <- substr(df_exp$GeneID, 1, 15)
  df_exp[1:4,1:4]
  df_gene <- readRDS(file.path(dir_main, 'tool/gene.lenth.symbol/Gene.ID.Length.rds'))
  df_gene[1:3,1:3]
  
  df_exp <- merge(df_exp, df_gene, by = "GeneID", all.x = TRUE)
  rownames(df_exp) <- paste0(df_exp[,1], seq_len(nrow(df_exp)))
  df_exp$GeneID <- NULL
  df_exp$Length <- NULL

  df_tmp <- aggregate(. ~ Symbol, data = df_exp, FUN = sum)
  rownames(df_tmp) <- df_tmp$Symbol
  df_tmp <- df_tmp[rowSums(df_tmp[,-1] != 0) > 0, ]
  

  df_gene2 <- aggregate(Length ~ Symbol, data = df_gene, FUN = max)
  head(df_gene2)
  
  df_tmp2 <- merge(df_tmp, df_gene2, by = "Symbol", all.x = TRUE)
  ls_Length <- df_tmp2$Length 
  df_tmp$Symbol <- NULL

  
  calculate_tpm <- function(counts, lengths) {
    rpk <- counts / (lengths / 1000)  # 每千碱基的读数
    total_rpk <- rowSums(rpk)
    tpm <- rpk / total_rpk * 1e6
    return(tpm)
  }
  
  df_bulk <- calculate_tpm(df_tmp, ls_Length)

  ## clinical
  df_clin <- readRDS(file.path(dir_dataset, paste0(bulk_name, "_clin0.rds")))
  
  df_bulk <- df_bulk[,df_clin$title]
  colnames(df_bulk) <- paste0(bulk_name, '_', df_clin$geo_accession)
  df_clin$geo_accession <- paste0(bulk_name, '_', df_clin$geo_accession)
  
  
  # colnames(df_clin)
  choo_col <- c('geo_accession',
                "os days:ch1",
                "os event:ch1"
  )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'OS', 'Survial')  # 1,dead
  
  df_clin2$OS = as.numeric(df_clin2$OS)
  df_clin2$OS <- round(df_clin2$OS/30,2)
  df_clin2$dataset <- bulk_name
  
  
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds")))
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds")))
}


#### GSE144269
if (F){
  bulk_name = 'GSE144269'
  df_exp <- readRDS(file.path(dir_dataset, "GSE144269/GSE144269.rds"))
  df_exp <- df_exp[df_exp$diagnosis == 'HCC',]
  df_exp$diagnosis <- NULL
  df_bulk <- t(df_exp)

  library("openxlsx")
  clin1 <- read.xlsx(file.path(dir_dataset, "GSE144269/41467_2020_18186_MOESM4_ESM.xlsx"), sheet = 1)
  colnames(clin1) <- clin1[2,]
  clin1 <- clin1[-c(1:2),]
  clin1 <- clin1[-c(77:83),]
  clin1 <- data.frame(clin1)
  
  
  eSet <- readRDS(file.path(dir_dataset, "GSE144269/GSE144269_eSet.rds"))
  clin2 = pData(eSet[[1]]) 
  ID <- lapply(strsplit(clin2$title, " "), function(x) x[2])
  clin2$ID <- unlist(ID)
  clin2 <- clin2[!grepl("non-tumor", clin2$title), ]
  
  clin1 <- clin1[clin1$Patient.ID %in% clin2$ID,]
  clin2 <- clin2[clin2$ID %in% clin1$Patient.ID,]
  clin2 <- clin2[match(clin1$Patient.ID, clin2$ID),]
  clin1$Patient.ID <- clin2$geo_accession
  
  df_clin <- merge(clin1, clin2, by.x = "Patient.ID", by.y = "geo_accession")
  
  
  # colnames(df_clin)
  choo_col <- c('Patient.ID', 
                'Time', 
                'Status'
  )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'OS', 'Survial')  # 1,dead
  
  df_clin2$OS = as.numeric(df_clin2$OS)
  df_clin2$OS <- round(df_clin2$OS/30,2)
  df_clin2$Survial <- ifelse(df_clin2$Survial == "Alive", 0, 1)
  df_clin2$dataset <- bulk_name
  
  df_clin2$Sample <-  paste0(bulk_name, '_', df_clin2$Sample)
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  
  df_bulk <- df_bulk[,colnames(df_bulk) %in% df_clin2$Sample]
  df_bulk <- df_bulk[,df_clin2$Sample]
  
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
}



#### GSE14520
if (F){
  bulk_name = 'GSE14520'
  df_exp <- readRDS(file.path(dir_dataset, "GSE14520.rds"))
  df_exp <- df_exp[df_exp$diagnosis == 'HCC',]
  df_exp$diagnosis <- NULL
  df_bulk <- t(df_exp)
  
  
  clin <- read.delim(file.path(dir_dataset, 'GSE14520_Extra_Supplement.txt'), head = T)
  meta <- data.frame(clin)
  
  meta <- meta[meta$Affy_GSM %in% colnames(df_bulk),]
  df_clin <- meta[match(colnames(df_bulk), meta$Affy_GSM), ]
  
  # colnames(df_clin)
  choo_col <- c('Affy_GSM', 
                'Survival.months', 
                'Survival.status'
  )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'OS', 'Survial')  # 1,dead
  df_clin2$dataset <- bulk_name
  
  df_clin2$Sample <-  paste0(bulk_name, '_', df_clin2$Sample)
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  
  
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
}



#### GSE16757
if (F){
  bulk_name = 'GSE16757'
  df_exp <- readRDS(file.path(dir_dataset, "GSE16757/GSE16757_bulk0.rds"))
  df_bulk <- df_exp
  
  clin <- read.delim(file.path(dir_dataset, 'GSE16757/GSE16757_Clinical.data.csv'), head = T, sep = ',')
  clin <- data.frame(clin)
  clin <- clin[clin$GEO_ID %in% colnames(df_exp),]
  df_clin <- clin[match(colnames(df_exp), clin$GEO_ID), ]
  
  # colnames(df_clin)
  choo_col <- c('GEO_ID', 
                'OS.m', 
                'Death..1.yes.'
  )
  df_clin2 <- df_clin[, colnames(df_clin) %in% choo_col]
  df_clin2 <- df_clin2[,choo_col]
  colnames(df_clin2) <- c('Sample', 'OS', 'Survial')  # 1,dead
  df_clin2$dataset <- bulk_name
  
  df_clin2$Sample <-  paste0(bulk_name, '_', df_clin2$Sample)
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  
  
  saveRDS(df_clin2, file.path(dir_dataset, paste0(bulk_name, "_clin.rds"))) 
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
}


#### bRNA-seq 
if (F){
  bulk_name = 'bRNAseq'
  df_bulk <- readRDS(file.path(dir_main, "101.bulk/dataset", "df_RPM.rds"))
  selected_columns <- grep("bulk_C", colnames(df_bulk), value = TRUE)
  df_bulk <- df_bulk[, selected_columns]
  colnames(df_bulk) <- gsub("bulk_C_", "ID", colnames(df_bulk))
  colnames(df_bulk) <- sprintf("ID%02d", as.integer(gsub("ID", "", colnames(df_bulk))))
  colnames(df_bulk) <- paste0(bulk_name, '_', colnames(df_bulk))
  saveRDS(df_bulk, file.path(dir_dataset, paste0(bulk_name, "_bulk.rds"))) 
}



#### integration
clin1 <- readRDS(file.path(dir_dataset, "PMID31585088_clin.rds"))
clin2 <- readRDS(file.path(dir_dataset, "LIRI_clin.rds"))
clin3 <- readRDS(file.path(dir_dataset, "TCGA_clin.rds"))
clin_trai <- rbind(clin1, clin2, clin3)
saveRDS(clin_trai, file.path(dir_dataset, paste0("clin_trai.rds"))) 


clin1 <- readRDS(file.path(dir_dataset, "GSE141198_clin.rds"))
clin2 <- readRDS(file.path(dir_dataset, "GSE141200_clin.rds"))
clin3 <- readRDS(file.path(dir_dataset, "GSE144269_clin.rds"))
clin4 <- readRDS(file.path(dir_dataset, "GSE14520_clin.rds"))
clin5 <- readRDS(file.path(dir_dataset, "GSE16757_clin.rds"))

# colnames(clin4)
clin_trai <- readRDS(file.path(dir_dataset, paste0("clin_trai.rds"))) 
clin_trai <- clin_trai[,c("Sample", "OS", "Survial", "dataset")]

clin_vali <- rbind(clin1, clin2, clin3, clin4, clin5,
                   clin_trai)
saveRDS(clin_vali, file.path(dir_dataset, paste0("clin_vali.rds"))) 


################################################################################
##### ComBat
bulk_data <- list(
  readRDS(file.path(dir_dataset, "PMID31585088_bulk.rds")),
  readRDS(file.path(dir_dataset, "TCGA_bulk.rds")),
  readRDS(file.path(dir_dataset, "GSE14520_bulk.rds")),
  readRDS(file.path(dir_dataset, "GSE16757_bulk.rds")),
  readRDS(file.path(dir_dataset, "GSE141198_bulk.rds")),
  readRDS(file.path(dir_dataset, "GSE141200_bulk.rds")),
  readRDS(file.path(dir_dataset, "GSE144269_bulk.rds")),
  readRDS(file.path(dir_dataset, "LIRI_bulk.rds"))
)


row_names_list <- lapply(bulk_data, rownames)
common_rownames <- Reduce(intersect, row_names_list)


xlsx <- read.xlsx(file.path(dir_main, "16.mlr3verse/importance", 'importance_Hep0.xlsx'))
ls_gene1 = xlsx$Feature[1:30]
xlsx <- read.xlsx(file.path(dir_main, "16.mlr3verse/importance", 'importance_Hep2.xlsx'))
ls_gene2 = xlsx$Feature[1:30]
xlsx <- read.xlsx(file.path(dir_main, "16.mlr3verse/importance", 'importance_Hep9.xlsx'))
ls_gene3 = xlsx$Feature[1:30]
gene_feature <- unique(c(ls_gene1, ls_gene2, ls_gene3))
gene_feature <- intersect(common_rownames, gene_feature)

saveRDS(gene_feature, file.path(dir_dataset, paste0("gene_feature.rds"))) 


subset_data <- lapply(bulk_data, function(df) df[gene_feature, , drop = FALSE])
df_merge <- do.call(cbind, subset_data)
df_merge <- log2(df_merge + 1)  # max(df_merge) # min(df_merge)


library(sva)
batch <- sapply(strsplit(colnames(df_merge), "_"), `[`, 1) 

df_ComBat <- ComBat(df_merge, batch = batch)
df_ComBat <- data.frame(df_ComBat)

saveRDS(df_ComBat, file.path(dir_dataset, paste0("df_ComBat.rds"))) 


### pca
perform_pca <- function(input_df, input_batch) {
  pca_result <- prcomp(t(input_df), scale. = TRUE) 
  pca_scores <- pca_result$x[, 1:2] 
  pca_df <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], Batch = input_batch)
  return(pca_df)
}
