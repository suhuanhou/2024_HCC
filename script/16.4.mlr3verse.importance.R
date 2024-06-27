# importance (github = no)
################################################################################
##### setting
# conda activate mlr3verse
# R

# library(data.table)
# library(openxlsx)
# ; library(dplyr); library(broom); library(factoextra); library(future); 
library(mlr3verse)
library(ggplot2)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
# dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_result = file.path(dir_main, "16.mlr3verse/importance"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result);
################################################################################
##### read data
# (muti~1)
subtype = 'Hep0'
subtype = 'Hep2'
subtype = 'Hep9'


learner_rf <- readRDS(file.path(dir_main, paste0("16.mlr3verse/result_",subtype, "/importance_", subtype, ".rds")))
importance = as.data.table(learner_rf$importance(), keep.rownames = TRUE)
colnames(importance) = c("Feature", "Importance")

p <- ggplot(data = importance[1:20], aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  xlab("") +
  ylab("coefficient") +
  ggtitle(subtype) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5))

png(file.path(dir_result, paste0("importance_", subtype, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()


