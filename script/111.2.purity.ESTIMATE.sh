################################################################################
##### setting
conda activate purity 
R

library(estimate)
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
library(ggsci)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "101.bulk/dataset")
dir_dataset = file.path(dir_main, "111.purity/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "111.purity/result_ESTIMATE"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result);

bool_debug = TRUE
# bool_debug = FALSE
################################################################################
##### data
df_bulk <- readRDS(file.path(dir_in, 'df_RPM.rds'))

df_tmp <- cbind(df_bulk[,56:165], df_bulk[,1:55])
df_bulk <- df_tmp
df_bulk <- df_bulk[,grep("bulk_C|bulk_P", colnames(df_bulk))]  
head(df_bulk)[,1:5]

write.table(df_bulk, file = "df_bulk.txt", sep = "\t", quote = FALSE)

################################################################################
##### run ESTIMATE
in.file <- 'df_bulk.txt' 
outputGCT(in.file, 'ESTIMATE_input.gct') 
filterCommonGenes(input.f= in.file, output.f= 'ESTIMATE_input.gct', id="GeneSymbol")
estimateScore("ESTIMATE_input.gct", "ESTIMATE_score.gct")


#### visual
plotPurity(scores="ESTIMATE_score.gct", samples="all_samples")

ESTIMATE_score <- read.table("ESTIMATE_score.gct", skip = 2,
                             header = TRUE,row.names = 1)
ESTIMATE_score <- ESTIMATE_score[,2:ncol(ESTIMATE_score)]
ESTIMATE_score <- data.frame(t(ESTIMATE_score))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score$group = c(rep('Paracancer', 55), rep('Cancer', 55))
ESTIMATE_score$group <- factor(ESTIMATE_score$group, levels = c('Paracancer', 'Cancer'))
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]

write.table(ESTIMATE_score, file = "ESTIMATE_score.txt", quote = F,sep = "\t", row.names = F)
saveRDS(ESTIMATE_score, 'ESTIMATE_score.rds')

################################################################################
##### result
df_score <- ESTIMATE_score
# df_score <- readRDS('ESTIMATE_score.rds')

p <- ggplot(df_score, aes(x = group, y = TumorPurity, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_manual(values = c("Paracancer" = "#F2A162", "Cancer" = "#E56E51")) + 
  labs(x = "", y = 'ESTIMATE Tumor Purity') +
  stat_compare_means(comparisons = list(c("Paracancer", "Cancer"))) + 
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        
png(file.path(dir_result, paste0("111.TumorPurity.all.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()
   

p <- ggplot(df_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = c("Paracancer" = "#F2A162", "Cancer" = "#E56E51")) +  
  labs(x = "", y = 'ESTIMATE Immune Score') +
  stat_compare_means(comparisons = list(c("Paracancer", "Cancer"))) +  
  theme_bw(base_size = 16) +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        

png(file.path(dir_result, paste0("111.Immune.all.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()



p <- ggplot(df_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = c("Paracancer" = "#F2A162", "Cancer" = "#E56E51")) +  
  labs(x = "", y = 'ESTIMATE Stromal Score') +
  stat_compare_means(comparisons = list(c("Paracancer", "Cancer"))) + 
  theme_bw(base_size = 16) +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
png(file.path(dir_result, paste0("111.Stromal.all.png")), width = 2000, height = 2000, res = 400, bg = "transparent"); print(p); dev.off()