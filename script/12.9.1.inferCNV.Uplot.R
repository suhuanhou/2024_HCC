# inferCNV & UPhyloplot2
################################################################################
##### re-run inferCNV
:'
conda activate inferCNV
ulimit -s 10240000  # Error: C stack usage  7969412 is too close to the limit

cd /share/home/shh/HCC_snRNA-seq/12.inferCNV/
nohup Rscript /share/home/shh/HCC_snRNA-seq/script/12.9.2.inferCNV.Uplot.inferCNV.R  > inferCNV.Uplot.20-22.log 2>&1 &
'

################################################################################
##### UPhyloplot2
:'
conda activate inferCNV
dir_=/share/home/shh/HCC_snRNA-seq/12.inferCNV

for i in {01..22}
do
    choo_target="IDsn$i"
    cd ${dir_}/9.inferCNV\(Uplot\)/${choo_target}

    sed '/^Endothelial cell/d' < 17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings > ${choo_target}.cell_groupings
    cp ${choo_target}.cell_groupings ${dir_}/Uplot2/Inputs/
done

cd ${dir_}/Uplot2
python uphyloplot2.py
cp output.svg ${dir_}/9.result\(Uplot\)/12.9.Uplot.svg
'



################################################################################
##### setting
library(tidyverse); library(tidyr)
library(dplyr); library(openxlsx)

rm(list = ls()); set.seed(100)
dir_main = "/share/home/shh/HCC_snRNA-seq"
dir_ = file.path(dir_main, "12.inferCNV")
dir_dataset = file.path(dir_main, "12.inferCNV", "9.inferCNV(Uplot)"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "12.inferCNV", "9.result(Uplot)"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result)
# bool_reRun = F

##### read data
IDsn_list <- paste0("IDsn", sprintf("%02d", 1:22))

df1 <- data.frame(IDsn = IDsn_list,
                  "cnv_1" = "", "cnv_1.1" = "", "cnv_1.2" = "","cnv_1.1.1" = "",
                  "cnv_1.1.1.1" = "", "cnv_1.1.1.2" = "", "cnv_1.1.2" = "",
                  "cnv_1.2.1" = "",
                  "cnv_1.2.1.1" = "", "cnv_1.2.1.2" = "", "cnv_1.2.2" = "",
                  "cnv_1.2.2.1" = "", "cnv_1.2.2.2" = ""
) 

df2 <- data.frame(IDsn = character(0), node = character(0), CNV = character(0)) 

for (choo_target in IDsn_list) {
  # choo_target = 'IDsn01'
  print(choo_target)
  
  
  ### chr
  chr_pq=read.table(file.path(dir_, "0.resource", "chr_pq.new.txt"),header = T,sep = "\t",stringsAsFactors = F) 
  chr_pq$chr=paste("chr",chr_pq$chr,sep = "")
  
  ### cnv
  cnv_regions=read.table(file.path(dir_dataset, choo_target, "HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"),header = T,sep = "\t",stringsAsFactors = F)
  cnv_regions <- cnv_regions[grepl("Hep", cnv_regions$cell_group_name), ]
  cnv_regions$cell_group_name <- gsub("Hep[0-9]+\\.Hep[0-9]+\\.", "", cnv_regions$cell_group_name)
  
  cnv_regions=cnv_regions%>%filter(state!=3)
  # cnv_regions$cell_group_name=str_replace(cnv_regions$cell_group_name,"all.*observations\\.","")
  cnv_regions$cnv_type=""
  
  for (i in 1:nrow(cnv_regions)) {
    tmp_chr_pq=chr_pq%>%filter(chr==cnv_regions[i,"chr"])
    if(cnv_regions[i,"start"] <= tmp_chr_pq[1,"cutoff"] & cnv_regions[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
      cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
    }else if (cnv_regions[i,"start"] >= tmp_chr_pq[2,"cutoff"] & cnv_regions[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
      cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
    }else {
      cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"chr"],"p,",cnv_regions[i,"chr"],"q",sep = "")
    }
    if (cnv_regions[i,"state"] < 3) {
      cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_loss",sep = "")
    }
    if (cnv_regions[i,"state"] > 3) {
      cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_gain",sep = "")
    }
    if (str_detect(cnv_regions[i,"cnv_type"],",")) {
      tmp1=strsplit(cnv_regions[i,"cnv_type"],",")[[1]][1]
      tmp2=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
      tmp3=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
      cnv_regions[i,"cnv_type"]=paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
      rm(list = c("tmp1","tmp2","tmp3"))
    }
  }
  
  
  ##########################################################
  ratio=0.5
  cnv_regions_copy=cnv_regions
  
  cnv_regions_part1=cnv_regions[str_detect(cnv_regions$cnv_type,","),]
  cnv_regions_part2=cnv_regions[!str_detect(cnv_regions$cnv_type,","),]
  cnv_regions_part1_new=as.data.frame(matrix(NA,ncol = ncol(cnv_regions_part1), nrow = nrow(cnv_regions_part1)*2 ))
  colnames(cnv_regions_part1_new)=colnames(cnv_regions_part1)
  for (i in 1:dim(cnv_regions_part1)[1]) {
    cnv_regions_part1_new[2*i-1,]=cnv_regions_part1[i,]
    cnv_regions_part1_new[2*i,]  =cnv_regions_part1[i,]
    
    tmp_chr_pq=chr_pq[chr_pq$chr == cnv_regions_part1[i,"chr"],]
    cnv_regions_part1_new[2*i-1,"end"] = tmp_chr_pq[1,"cutoff"]
    cnv_regions_part1_new[2*i,"start"] = tmp_chr_pq[2,"cutoff"]
    
    cnv_regions_part1_new[2*i-1,"cnv_type"] = strsplit(cnv_regions_part1_new[2*i-1,"cnv_type"],",")[[1]][1]
    cnv_regions_part1_new[2*i,"cnv_type"]   = strsplit(cnv_regions_part1_new[2*i,  "cnv_type"],",")[[1]][2]
  }
  cnv_regions=cnv_regions_part2 %>% rbind(cnv_regions_part1_new)
  cnv_regions$chr=factor(cnv_regions$chr,levels = paste("chr",1:22,sep = ""))
  cnv_regions=cnv_regions%>%arrange(cell_group_name,chr,start)
  cnv_regions$region_len=cnv_regions$end-cnv_regions$start+1
  cnv_regions=cnv_regions[,c("cell_group_name","cnv_type","region_len")]
  cnv_regions=as.data.frame(cnv_regions%>%group_by(cell_group_name,cnv_type)%>%dplyr::summarize(all_region_len=sum(region_len)))
  
  ####
  cnv_regions$final_label=""

  
  for (i in 1:dim(cnv_regions)[1]) {
    arm=strsplit(cnv_regions[i,"cnv_type"],"_")[[1]][1]
    ref_len=chr_pq[paste(chr_pq$chr,chr_pq$arm,sep = "") == arm,"arm_len"]
    
    if(cnv_regions[i,"all_region_len"] >= ref_len*ratio){
      cnv_regions[i,"final_label"]=cnv_regions[i,"cnv_type"]
    }else{
      cnv_regions[i,"final_label"]=""
    }
  }
  
  cnv_regions$cnv_type=NULL
  colnames(cnv_regions)[3]="cnv_type"
  cnv_regions=cnv_regions%>%filter(cnv_type != "")
  

  ##########################################################
  cnv_1.1.1.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.1"))[,"cnv_type"]
  cnv_1.1.1.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.1.2"))[,"cnv_type"]
  cnv_1.1.2  =as.data.frame(cnv_regions%>%filter(cell_group_name=="1.1.2"  ))[,"cnv_type"]
  cnv_1.2.1.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.1"))[,"cnv_type"]
  cnv_1.2.1.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.1.2"))[,"cnv_type"]
  cnv_1.2.2.1=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.1"))[,"cnv_type"]
  cnv_1.2.2.2=as.data.frame(cnv_regions%>%filter(cell_group_name=="1.2.2.2"))[,"cnv_type"]
  
  ###cnv_1.1.1
  cnv_1.1.1=intersect(cnv_1.1.1.1,cnv_1.1.1.2)
  cnv_1.1.1.1_uniq=setdiff(cnv_1.1.1.1,cnv_1.1.1)
  cnv_1.1.1.2_uniq=setdiff(cnv_1.1.1.2,cnv_1.1.1)
  ###cnv_1.2.1
  cnv_1.2.1=intersect(cnv_1.2.1.1,cnv_1.2.1.2)
  cnv_1.2.1.1_uniq=setdiff(cnv_1.2.1.1,cnv_1.2.1)
  cnv_1.2.1.2_uniq=setdiff(cnv_1.2.1.2,cnv_1.2.1)
  ###cnv_1.2.2
  cnv_1.2.2=intersect(cnv_1.2.2.1,cnv_1.2.2.2)
  cnv_1.2.2.1_uniq=setdiff(cnv_1.2.2.1,cnv_1.2.2)
  cnv_1.2.2.2_uniq=setdiff(cnv_1.2.2.2,cnv_1.2.2)
  ###cnv_1.1
  cnv_1.1=intersect(cnv_1.1.1,cnv_1.1.2)
  cnv_1.1.1_uniq=setdiff(cnv_1.1.1,cnv_1.1)
  cnv_1.1.2_uniq=setdiff(cnv_1.1.2,cnv_1.1)
  ###cnv_1.2
  cnv_1.2=intersect(cnv_1.2.1,cnv_1.2.2)
  cnv_1.2.1_uniq=setdiff(cnv_1.2.1,cnv_1.2)
  cnv_1.2.2_uniq=setdiff(cnv_1.2.2,cnv_1.2)
  ###cnv_1
  cnv_1=intersect(cnv_1.1,cnv_1.2)
  cnv_1.1_uniq=setdiff(cnv_1.1,cnv_1)
  cnv_1.2_uniq=setdiff(cnv_1.2,cnv_1)
  
  
 
  df1[df1$IDsn == choo_target, "cnv_1"] <- paste(cnv_1, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.1"] <- paste(cnv_1.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.2"] <- paste(cnv_1.2_uniq, collapse = ", ")
  
  df1[df1$IDsn == choo_target, "cnv_1.1.1.1"] <- paste(cnv_1.1.1.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.1.1.2"] <- paste(cnv_1.1.1.2_uniq, collapse = ", ")

  df1[df1$IDsn == choo_target, "cnv_1.2.1.1"] <- paste(cnv_1.2.1.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.2.1.2"] <- paste(cnv_1.2.1.2_uniq, collapse = ", ")

  df1[df1$IDsn == choo_target, "cnv_1.2.2.1"] <- paste(cnv_1.2.2.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.2.2.2"] <- paste(cnv_1.2.2.2_uniq, collapse = ", ")
  
  df1[df1$IDsn == choo_target, "cnv_1.1.1"] <- paste(cnv_1.1.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.1.2"] <- paste(cnv_1.1.2_uniq, collapse = ", ")

  df1[df1$IDsn == choo_target, "cnv_1.2.1"] <- paste(cnv_1.2.1_uniq, collapse = ", ")
  df1[df1$IDsn == choo_target, "cnv_1.2.2"] <- paste(cnv_1.2.2_uniq, collapse = ", ")


  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1", CNV = paste(cnv_1, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.1", CNV = paste(cnv_1.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2", CNV = paste(cnv_1.2_uniq, collapse = ", ")))
  
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.1.1", CNV = paste(cnv_1.1.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.1.1.1", CNV = paste(cnv_1.1.1.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.1.1.2", CNV = paste(cnv_1.1.1.2_uniq, collapse = ", ")))
  
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.1.2", CNV = paste(cnv_1.1.2_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.1", CNV = paste(cnv_1.2.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.1.1", CNV = paste(cnv_1.2.1.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.1.2", CNV = paste(cnv_1.2.1.2_uniq, collapse = ", ")))
  
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.2", CNV = paste(cnv_1.2.2_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.2.1", CNV = paste(cnv_1.2.2.1_uniq, collapse = ", ")))
  df2 <- rbind(df2, data.frame(IDsn = choo_target, node = "cnv_1.2.2.2", CNV = paste(cnv_1.2.2.2_uniq, collapse = ", ")))
  
}


df3 <- df2 %>% separate_rows(CNV, sep = ", ")   # 拆分
df3 <- df3[!df3$CNV == '',]

# filter 
df4 <- df3
df4 <- df4[df4$node %in% c('cnv_1', 'cnv_1.1', 'cnv_1.2'),]
df4$node = NULL
df4$freq = NULL
df4 <- df4 %>%distinct()
df4 <- df4 %>% count(CNV)
df4 <- df4 %>% rename(count = n)
df4$freq <- df4$count / 22

ls_CNV <- df4[df4$freq > 0.3, 'CNV'] %>% pull()
df5 <- df3[df3$CNV %in% ls_CNV,]
df5$CNV <- gsub("^chr", "", df5$CNV)
df5 <- df5[df5$node %in% c('cnv_1', 'cnv_1.1', 'cnv_1.2'),]


df5$chr <- sub("^([0-9]+).*", "\\1", df5$CNV)
df5$chr <- as.numeric(df5$chr)
df5 <- df5[order(df5$IDsn, df5$chr), ]
write.xlsx(df5, "df5.xlsx")



################################################################################
##### Donut chart

library(Seurat); library(ggplot2)
sce <- readRDS(file.path(dir_main, "5.Hepatocyte/dataset", "sce_Hep.rds"))


df_counts <- as.data.frame.matrix(table(sce$patient, sce$subtype))
df_counts <- df_counts / rowSums(df_counts)
df_ratio <- tibble::rownames_to_column(df_counts, var = "patient")


ls_color = c('#3EA94B','#3E86BD','#E51F17','#DE75DD','#EC8392',
             '#31AADF','#FDCE4F','#EC7B1A','#1505A8','#C20576')

plot_donut <- function(row){
  # row = 1
  df_row <- df_ratio[row,]
  df_donut <- data.frame(category = names(df_row)[-1], fraction = unlist(df_row)[-1])
  df_donut$ymax = cumsum(df_donut$fraction)
  df_donut$ymin = c(0, head(df_donut$ymax, n = -1))
  
  p <- ggplot(data = df_donut, aes(fill = category, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(colour = "grey30", show_guide = FALSE) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", title = "") + 
    xlim(c(0, 4)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = ls_color) +
    geom_text(aes(x = 0, y = 0, label = df_ratio$patient[row]), color = "black", size = 24)
  
  png(file.path(dir_result, paste0("12.9.1.donut_", df_ratio$patient[row], ".png")), width = 2500, height = 2500, res = 400, bg = "transparent"); print(p); dev.off()
  
}

plot_donuts <- lapply(1:nrow(df_ratio), plot_donut)


