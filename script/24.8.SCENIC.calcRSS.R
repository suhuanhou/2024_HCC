################################################################################
##### test
# opt <- list(input_loom = "aucell.loom", input_meta = "metadata.xls", celltype = "subtype")

################################################################################
##### run
library(optparse)
op_list <- list(
  make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file",metavar="rds"),
  make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object",metavar="idents"),
  make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS",metavar="lab
el")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)

loom <- open_loom(opt$input_loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

meta <- read.table(opt$input_meta,sep='\t',header=T,stringsAsFactor=F)

cellinfo <- meta[,c(opt$celltype,"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
try({
  rssPlot <- plotRSS(rss)
  save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
})

saveRDS(rss,paste0(opt$celltype,"_rss.rds"))

