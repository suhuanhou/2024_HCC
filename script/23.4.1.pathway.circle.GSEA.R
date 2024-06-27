# GSEA Circle Chart
################################################################################
##### setting
library(RColorBrewer)
library(colorspace)
library(circlize)
library(ComplexHeatmap)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
dir_dataset = file.path(dir_main, "23.pathway/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "23.pathway/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session_4.txt"))
setwd(dir_result)

bool_debug = FALSE
################################################################################
##### target
choo_ = c('Hep02', 'Hep09', 'Hep2', 'Hep9')
choo_subtype = choo_[2]


file <- 'GSEA'
################################################################################
##### color
# ls_color1 <- rainbow_hcl(30)
ls_color2 <- brewer.pal(9, "Set1")
################################################################################
##### GSEA
#### data
# file.path(dir_dataset, paste0("GSEA_", choo_subtype, ".txt"))
# file.path(paste0("circle_", choo_subtype, ".pdf"))
dat <- read.delim(file.path(dir_dataset, paste0(file, "_", choo_subtype, ".txt")), sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
dat$id <- factor(rownames(dat), levels = rownames(dat))


## KEGG Classification
library(openxlsx)
KEGG_class <- read.xlsx(file.path(dir_main, '/tool/KEGG-pathway-classification.xlsx'), sheet = "Sheet1")
KEGG_class$Pathway.ID <- paste0('hsa', KEGG_class$Pathway.ID)

missing_element <- dat$category[!(dat$category %in% KEGG_class$Pathway.ID)]
if(length(missing_element) > 0){
  print("please update hsa")
  print("https://www.kegg.jp/kegg/pathway.html")
  print(missing_element)
}

dat$category <- KEGG_class[match(dat$id, KEGG_class$Pathway.ID), 'Pathway.Class.1']
dat <- dat[order(dat$category),]

category_freq <- data.frame(table(dat$category))

ko_color <- c()
for (i in 1:length(category_freq$Freq)){
  ko_color <- c(ko_color, rep(ls_color2[i], category_freq$Freq[i]))
}

dat$color <- ko_color

##### Draw
dev.off()
## Circle 1:hsa
# pdf(file.path(paste0("circle_", choo_subtype, ".pdf")), width = 48, height = 24)
circle_size = unit(1, 'snpc')
circos.par(gap.degree = 0.5, start.degree = 90)
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max')]


circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')#ylim、xlim
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')#sector.name
    circos.axis(h = 'top', labels.cex = 0.4, labels.niceFacing = FALSE)
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  }
)

## Circle 2: gene count + p
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.rich', '-log10Pvalue')]
label_data <- dat['gene_num.rich']
p_max <- round(max(dat$'-log10Pvalue')) + 1  
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
                       ylim = get.cell.meta.data('ycenter')  
                       xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
                       sector.name = label_data[get.cell.meta.data('sector.index'),1]
                       circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  }
)
    
## Circle 3: gene
dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('id', 'gene_num.min', 'up')]
names(plot_data_up) <- c('id', 'start', 'end')
plot_data_up$type <- 1 

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('id', 'up', 'down')]
names(plot_data_down) <- c('id', 'start', 'end')
plot_data_down$type <- 2 

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE, 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...) 
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  })


## Circle 4: enrich factor
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max', 'rich.factor')] 
label_data <- dat['id']  
color_assign <- dat$color

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  #sector.name 
    circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...) 
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'grey', lwd = 0.3) 
  } )


p <- recordPlot()  
if(bool_debug == FALSE){dev.off()}


df_category_legend <- data.frame(category = dat$category, color = ko_color)
df_category_legend <- unique(df_category_legend)

category_legend <- Legend(
  labels = df_category_legend$category,
  type = 'points', pch = NA, background = df_category_legend$color, 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'), 
  type = 'points', pch = NA, background = c('red', 'blue'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')


lgd_list_vertical <- packLegend(category_legend, pvalue_legend, updown_legend)
grid.draw(lgd_list_vertical)
circos.clear()

if(bool_debug == FALSE){
  p_legend <- recordPlot()  
  
  png(file.path(dir_result, paste0("23.circle_plot_", file, "_", choo_subtype, ".png")), width = 4000, height = 3500, res = 400, bg = 'transparent'); print(p); dev.off()
  png(file.path(dir_result, paste0("23.circle_legend_", file, "_", choo_subtype, ".png")), width = 1500, height = 3000, res = 400, bg = 'transparent'); print(p_legend); dev.off()
}


