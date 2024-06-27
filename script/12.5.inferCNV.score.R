# cnv score
################################################################################
##### setting
library(ggpubr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_dataset = file.path(dir_main, "12.inferCNV/dataset(Endo)"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "12.inferCNV/result(Endo)"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result)
################################################################################
##### data
cnv_table <- read.table("infercnv.observations.txt", header=T)
tmp <- read.table("infercnv.references.txt", header=T)
all_cnv_table <- cbind(tmp,cnv_table)

# sd
down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
oneCopy=up-down
oneCopy  #0.1085806
a1= down- 2*oneCopy
a2= down- 1*oneCopy
down;up
a3= up +  1*oneCopy
a4= up + 2*oneCopy 

# Replicate the table 
all_cnv_score_table <- as.matrix(all_cnv_table)
all_cnv_score_mat <- as.matrix(all_cnv_table)

# Scoring
all_cnv_score_table[all_cnv_score_mat > 0 & all_cnv_score_mat < a2] <- "A" #complete loss. 2pts
all_cnv_score_table[all_cnv_score_mat >= a2 & all_cnv_score_mat < down] <- "B" #loss of one copy. 1pts
all_cnv_score_table[all_cnv_score_mat >= down & all_cnv_score_mat <  up ] <- "C" #Neutral. 0pts
all_cnv_score_table[all_cnv_score_mat >= up  & all_cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
all_cnv_score_table[all_cnv_score_mat > a3  & all_cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
all_cnv_score_table[all_cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts

# Check
table(all_cnv_score_table[,2])
# Replace with score 
all_cnv_score_table_pts <- all_cnv_table
rm(all_cnv_score_mat)

all_cnv_score_table_pts[all_cnv_score_table == "A"] <- 2
all_cnv_score_table_pts[all_cnv_score_table == "B"] <- 1
all_cnv_score_table_pts[all_cnv_score_table == "C"] <- 0
all_cnv_score_table_pts[all_cnv_score_table == "D"] <- 1
all_cnv_score_table_pts[all_cnv_score_table == "E"] <- 2
all_cnv_score_table_pts[all_cnv_score_table == "F"] <- 2

# Scores are stored in “cnv_score_table_pts”. Use colSums to add up scores for each cell and store as vector 
all_cell_scores_CNV <- as.data.frame(colSums(all_cnv_score_table_pts))
colnames(all_cell_scores_CNV) <- "cnv_score"
head(all_cell_scores_CNV)
write.csv(x = all_cell_scores_CNV, file = "all_cnv_scores.csv")

## meta data
pbmc_harmony2 <- readRDS(file.path(dir_dataset, "sce_inferCNV.rds"))
phe = pbmc_harmony2@meta.data
# head(rownames(phe))
# head(rownames(all_cell_scores_CNV)) 

rownames(all_cell_scores_CNV)=gsub('[.]','-',rownames(all_cell_scores_CNV))
head(rownames(all_cell_scores_CNV))

phe=phe[rownames(phe) %in% rownames(all_cell_scores_CNV),]  
phe$cnv_score  =  all_cell_scores_CNV[rownames(phe),]
# head(rownames(phe))
# head(rownames(all_cell_scores_CNV)) 


##### visualization
phe$subtype <- factor(phe$subtype, levels = c("Hep0", "Hep1", "Hep2", "Hep3", "Hep4", "Hep5", "Hep6", "Hep7", "Hep8", "Hep9", "Endothelial cell"))
# table(phe$subtype)

#### boxplot
ls_color = c('#3EA94B','#3E86BD','#E51F17','#DE75DD','#EC8392',
             '#31AADF','#FDCE4F','#EC7B1A','#1505A8','#C20576', '#F6746B')

p <- ggplot(phe, aes(x = subtype, y = cnv_score, fill = subtype)) +
  geom_violin() +
  scale_fill_manual(values = ls_color) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 12),  
    axis.title.y = element_text(size = 14),  
    panel.border = element_rect(color = "black", fill = NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank()
  ) + 
  # ylim(0, 3000) + 
  guides(fill = FALSE) +
  stat_summary(geom = "crossbar", width = 0.5, fun.data = "mean_se") +
  stat_compare_means(label = "p.signif", 
                     method = "t.test",
                     ref.group = "Endothelial cell",
                     label.y = 3500, 
                     hide.ns = TRUE,
                     size = 4
  ); p

png(file.path(dir_result, "12.cnv_score.png"), width = 4000, height = 1500, res = 600, bg = "transparent"); print(p); dev.off()


