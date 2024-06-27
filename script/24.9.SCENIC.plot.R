# Regulon：调控子，定义为一个转录因子及其靶基因的集合
# RAS：Regulon活性得分（activity）
# RSS：Regulon特异性得分，反应Regulon与细胞类型的对应关系
# CSI: 不同Regulon之间的关联性


# https://zhuanlan.zhihu.com/p/144459785

################################################################################
##### setting
# conda activate scenic 
# R 

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(pheatmap)

library(cowplot)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(stringr)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "5.Hepatocyte/dataset")
# dir_dataset = file.path(dir_main, "24.SCENIC/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "24.SCENIC/result_bulk"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result);
################################################################################
##### visual
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
incolor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))




plot_pyscenic <- function(inloom='aucell.loom',incolor=incolor,inrss='seurat_annotations_rss.rds',inrds='subset.rds',infun='median', ct.col='seurat_annotations',inregulons=NULL,ingrn='grn.bulk.tsv',ntop1=5,ntop2=50){
  ## test
  if(F){
    inloom='aucell.loom'
    inrss='subtype_rss.rds'
    inrds="sce_bulk.rds"
    infun='median'
    ct.col='subtype'
    inregulons=NULL
    ingrn='grn.bulk.tsv'
    ntop1=5
    ntop2=50
  }
  
  
  ###load data
  loom <- open_loom(inloom)
  
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  regulonAucThresholds <- get_regulon_thresholds(loom)
  
  embeddings <- get_embeddings(loom)
  close_loom(loom)
  
  rss <- readRDS(inrss)
  sce <- readRDS(inrds)
  
  ##calculate  RSS fc
  df = do.call(rbind,
               lapply(1:ncol(rss), function(i){
                  dat= data.frame(
                  regulon  = rownames(rss),
                  cluster =  colnames(rss)[i],
                  sd.1 = rss[,i],
                  sd.2 = apply(rss[,-i], 1, get(infun))
                 )
               }))
  
  df$fc = df$sd.1 - df$sd.2
  
  #select top regulon
  ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
  
  ntopgene <- unique(ntopg$regulon)
  write.table(ntopgene,'sd_regulon_RSS.list',sep='\t',quote=F,row.names=F,col.names=F)
  #plot rss by cluster
  
  #using plotRSS
  rssPlot <- plotRSS(rss)
  regulonsToPlot <- rssPlot$rowOrder
  rp_df <- rssPlot$df
  
  write.table(regulonsToPlot,'rss_regulon.list',sep='\t',quote=F,row.names=F,col.names=F)
  write.table(rp_df,'rssPlot_data.xls',sep='\t',quote=F)
  nlen <- length(regulonsToPlot)
  hei <- ceiling(nlen)*0.4
  blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
  lgroup <- levels(rssPlot$df$cellType)
  
  nlen2 <- length(lgroup)
  wei <- nlen2*2
  pdf(paste0('regulons_RSS_',ct.col,'_in_dotplot.pdf'),wei,hei)
  print(rssPlot$plot)
  dev.off()
  
  # sd top gene
  anrow = data.frame( group = ntopg$cluster)
  lcolor <- incolor[1:length(unique(ntopg$cluster))]
  names(lcolor) <- unique(anrow$group)
  annotation_colors <- list(group=lcolor)
  
  pn1 = rss[ntopg$regulon,]
  pn2 = rss[unique(ntopg$regulon),]
  rownames(pn1) <-  make.unique(rownames(pn1))
  rownames(anrow) <- rownames(pn1)
  scale='row'
  hei <- ceiling(length(ntopg$regulon)*0.4)
  
  # pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
  # print(pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons'))
  # print(pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons'))
  # dev.off()  
  
  
  pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
  print(pheatmap(pn1,annotation_row = anrow,
                 scale=scale,annotation_colors=annotation_colors,
                 show_rownames = T, cluster_row = F, cluster_cols = F, 
                 main='sd top regulons'))
  print(pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons', angle_col = '45'))
  dev.off()
  
  
  #plotRSS gene
  
  pn2 = rss[unique(rp_df$Topic),]
  scale='row'
  hei <- ceiling(length(unique(rp_df$Topic))*0.4)
  pdf(paste0('regulon_RSS_in_plotRSS_',ct.col,'.pdf'),wei,hei + 2)
  print(pheatmap(pn2,scale=scale,show_rownames = T, cluster_row = T, cluster_cols = F, main='plotRSS unique regulons'))
  dev.off()
  
  #all regulons
  
  hei <- ceiling(length(rownames(rss))*0.2)
  pdf(paste0('all_regulons_RSS_in_',ct.col,'.pdf'),wei,hei)
  print(pheatmap(rss,scale=scale,show_rownames = T, cluster_row = T, cluster_cols = F, main='all regulons RSS'))
  dev.off()
  
  #plot rss by all cells
  if (is.null(inregulons)){
  inregulons <- regulonsToPlot
  }else{
  inregulons <- intersect(inregulons,rownames(rss))
  regulonsToPlot <- inregulons
  
  }
  pn3=as.matrix(regulonAUC@assays@data$AUC)
  regulon <- rownames(pn3)
  #regulon <- inregulons
  pn3 <- pn3[regulon,]
  #pn3 <- pn3[,sample(1:dim(pn3)[2],500)]
  
  sce$group1=sce@meta.data[,ct.col]
  
  meta <- sce@meta.data
  meta <- meta[order(meta$group1),]
  #meta <- meta[colnames(pn3),]
  ancol = data.frame(meta[,c('group1')])
  colnames(ancol) <- c('group1')
  rownames(ancol) <- rownames(meta)
  lcolor <- incolor[1:length(unique(ntopg$cluster))]
  names(lcolor) <- unique(ntopg$cluster)
  annotation_colors <- list(group1 =lcolor)
  
  df1 <- ancol
  df1$cell <- rownames(df1)
  df1 <- df1[order(df1$group1),]
  pn3 <- pn3[,rownames(df1)]
  torange=c(-2,2)
  pn3 <- scales::rescale(pn3,to=torange)
  pn3 <- pn3[,rownames(ancol)]
  
  scale='none'
  hei <- ceiling(length(unique(regulon))*0.2)
  pdf(paste0('all_regulon_activity_in_allcells.pdf'),10,hei)
  print(
  pheatmap(pn3,annotation_col = ancol,scale=scale,annotation_colors=annotation_colors,
           cluster_row = T, cluster_cols = F, show_rownames = T,show_colnames = F)
  )
  #pheatmap(pn3,scale=scale,show_rownames = T, show_colnames = F,cluster_cols=F)
  dev.off()
  
  #plot in seurat
  regulonsToPlot = inregulons
  sce$sub_celltype <- sce@meta.data[,ct.col]
  sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
  
  # cellClusters <- data.frame(row.names = colnames(sce),
  #                            seurat_clusters = as.character(sce$seurat_clusters))
  # cellTypes <- data.frame(row.names = colnames(sce),
  #                         celltype = sce$sub_celltype)
  
  sce@meta.data = cbind(sce@meta.data ,t(sub_regulonAUC@assays@data@listData$AUC[regulonsToPlot,]))
  Idents(sce) <- sce$sub_celltype
  
  nlen <- length(regulonsToPlot)
  hei <- ceiling(nlen)*0.4
  blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
  
  nlen2 <- length(unique(sce$sub_celltype))
  wei <- nlen2*2
  pdf('regulons_activity_in_dotplot.pdf',wei,hei+4)
  print(DotPlot(sce, features = unique(regulonsToPlot)) + coord_flip()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
        scale_color_gradientn(colours = blu)
        )
  dev.off()
  
  ## umap
  # hei=ceiling(nlen/4)*4
  # pdf('regulons_activity_in_umap.pdf',16,hei)
  # print(RidgePlot(sce, features = regulonsToPlot , ncol = 4))
  # print(VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ))
  # print(FeaturePlot(sce, features = regulonsToPlot))
  # dev.off()
  
  grn <- read.table(ingrn,sep='\t',header=T,stringsAsFactors=F)
  inregulons1=gsub('[(+)]','',inregulons)
  c1 <- which(grn$TF %in% inregulons1)
  grn <- grn[c1,]
  #edge1 <- data.frame()
  #node1 <- data.frame()
  pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
  for (tf in unique(grn$TF)) {
    tmp <- subset(grn,TF==tf)
    if (dim(tmp)[1] > ntop2) {
      tmp <- tmp[order(tmp$importance,decreasing=T),]
      tmp <- tmp[1:ntop2,]
    }
  
    node2 <- data.frame(tmp$target)
    node2$node.size=1.5
    node2$node.colour <- 'black'
    colnames(node2) <- c('node','node.size','node.colour')
    df1 <- data.frame(node=tf,node.size=2,node.colour='#FFDA00')
    node2 <- rbind(df1,node2)
    
    
    edge2 <- tmp
    colnames(edge2) <- c('from','to','edge.width')
    edge2$edge.colour <- "#1B9E77"
    torange=c(0.1,1)
    edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)
    
    graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
    p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
      scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
          geom_node_label(aes(label = node,colour = node.colour),size = 3.5, repel = TRUE)
    p1 <- p1 + scale_color_manual(values=c('#FFDA00','black'))+scale_edge_color_manual(values=c("#1B9E77"))
    print(p1)
  }
  dev.off()
  
  #plot activity heatmap
  meta <- sce@meta.data
  celltype <- ct.col
  cellsPerGroup <- split(rownames(meta),meta[,celltype])
  sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  # Calculate average expression:
  regulonActivity_byGroup <- sapply(cellsPerGroup,
                                    function(cells)
                                      rowMeans(getAUC(sub_regulonAUC)[,cells]))
  scale='row'
  rss <- regulonActivity_byGroup
  hei <- ceiling(length(regulonsToPlot)*0.4)
  pn1 <- rss[regulonsToPlot,]
  pdf(paste0('regulon_activity_in_',ct.col,'.pdf'),wei,hei)
  print(
  pheatmap(pn1,scale=scale,show_rownames = T, cluster_row = F, cluster_cols = F, main='regulons activity')
  )
  dev.off()
  
  rss <- rss[,c('Normal', 'Paracancer', 'Cancer')]
  rss = rss[apply(rss, 1, function(x) sd(x)!=0),]  # 去除标准差为0的行，即全为0
  
  hei <- ceiling(length(rownames(rss))*0.2)
  # pdf(paste0('24.all_regulons_activity_in_',ct.col,'.bulk.pdf'),wei,hei)
  # print(pheatmap(rss,scale=scale,show_rownames = T, cluster_row = T, cluster_cols = F, 
  #                main='all regulons activity', fontsize = 18,
  #                fontsize_row = 8, fontsize_col = 18, angle_col='45'
  #                ))
  # dev.off()
  
  p <- pheatmap(rss, scale='row', show_rownames = T,
                cluster_cols = F, 
                main='bulk RNA-seq', fontsize = 30,
                fontsize_row = 30, fontsize_col = 40, angle_col='45'
  )
  
  png(paste0('24.all_regulons_activity_in_',ct.col,'.bulk.png'), width = 4000, height = 8000, res = 400, bg = "transparent"); print(p); dev.off()
  
  
  
  ##calculate fc
  df = do.call(rbind,
               lapply(1:ncol(rss), function(i){
                  dat= data.frame(
                  regulon  = rownames(rss),
                  cluster =  colnames(rss)[i],
                  sd.1 = rss[,i],
                  sd.2 = apply(rss[,-i], 1, get(infun))
                 )
               }))
  
  df$fc = df$sd.1 - df$sd.2
  
  #select top regulon
  ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
  
  ntopgene <- unique(ntopg$regulon)
  write.table(ntopgene,'sd_regulon_activity.list',sep='\t',quote=F,row.names=F,col.names=F)
  
  anrow = data.frame( group = ntopg$cluster)
  lcolor <- incolor[1:length(unique(ntopg$cluster))]
  names(lcolor) <- unique(anrow$group)
  annotation_colors <- list(group=lcolor)
  
  pn1 = rss[ntopg$regulon,]
  pn2 = rss[unique(ntopg$regulon),]
  rownames(pn1) <-  make.unique(rownames(pn1))
  rownames(anrow) <- rownames(pn1)
  scale='row'
  hei <- ceiling(length(ntopg$regulon)*0.4)
  pdf(paste0('regulon_activity_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
  print(
  pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,
           cluster_row = F, cluster_cols = F, show_rownames = T,main='sd top regulons')
  )
  print(
  pheatmap(pn2,scale=scale,show_rownames = T, cluster_row = F, cluster_cols = F, main='sd top unique regulons ')
  )
  dev.off()

}



plot_pyscenic(inloom='aucell.loom',
  incolor=incolor,
  inrss='subtype_rss.rds',
  inrds='sce_bulk.rds',
  infun='median', 
  ct.col='subtype',
  inregulons=NULL,ingrn='grn.bulk.tsv',ntop1=5,ntop2=50)



