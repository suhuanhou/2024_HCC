# visualization
################################################################################
##### setting
library(ggplot2)


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
# dir_in = file.path(dir_main, "4.merge/dataset")
dir_dataset = file.path(dir_main, "20.bisque/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "20.bisque/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_dataset)

if (F){  # muti~1
  bulk_name = 'LIRI'
  bulk_name = 'TCGA'
  bulk_name = 'PMID31585088'
}

################################################################################
#### read data
bisque <- readRDS(file.path(dir_result, paste0(bulk_name, "_bisque.rds")))
df <- as.data.frame(bisque)

# order
df <- df[, order(df["Hep9", ])]
df$subtype <- rownames(df)


df_long <- tidyr::gather(df, key = "patient", value = "ratio", -subtype)
df_long$patient <- factor(df_long$patient, levels = colnames(df))


# colors
ls_color = c('#3EA94B','#3E86BD','#E51F17','#DE75DD','#EC8392',
             '#31AADF','#FDCE4F','#EC7B1A','#1505A8','#C20576')
my_colors <- ls_color



# Drawing Stacked Bars
p <- ggplot(df_long, aes(x = patient, y = ratio, fill = subtype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(x = paste0("patients = ", ncol(bisque)), y = "Ratio of subtype") + 
  theme(legend.position = "none",  
        axis.title.x = element_text(size = 40), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 40), 
        axis.text.y = element_text(size = 20),  
        plot.margin = margin(1, 0, 0, 0, "cm") 
  ) +
  scale_y_continuous(expand = c(0, 0)); p


png(file.path(dir_result, paste0("20.", bulk_name, "_bisque.ratio.png")), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()
