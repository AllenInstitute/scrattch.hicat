plot_cl_heatmap <- function(norm.dat, 
                            cl, 
                            markers, 
                            prefix = NULL,
                            hc = NULL,
                            gene.hc = NULL,
                            centered = FALSE,
                            labels = names(cl),
                            sorted = FALSE,
                            by.cl = TRUE,
                            ColSideColors = NULL,
                            maxValue = 5,
                            min.sep = 4,
                            main = "", 
                            height = 13, 
                            width = 9)
  {
    library(matrixStats)
    blue.red <-colorRampPalette(c("blue", "white", "red"))
    select.cells=names(cl)
    tmp.dat = as.matrix(norm.dat[markers,select.cells,drop=F])
    if(!is.null(ColSideColors)){
      ColSideColors=ColSideColors[, select.cells,drop=F]
    }
    if(centered){
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      breaks=c(min(min(tmp.dat)-0.1,-maxValue),  seq(-maxValue,maxValue, length.out=99), max(max(tmp.dat)+1))
    }
    else{
      tmp.dat = tmp.dat/pmax(rowMaxs(tmp.dat), 1)
      breaks=c(0, seq(0.05, 1, length.out=100))
    }
    colnames(tmp.dat)=labels
    cexCol = min(70/ncol(tmp.dat),1)
    cexRow = min(60/nrow(tmp.dat),1)
    if(is.null(gene.hc)){
      gene.hc = hclust(dist(tmp.dat), method="ward.D")
    }
    if(is.null(hc) & !sorted & length(select.cells)< 2000){
      hc = hclust(dist(t(tmp.dat)), method="ward.D")
    }
    col = blue.red(150)[51:150]
    if(!is.null(prefix)){
      pdf(paste(prefix,"pdf",sep="."), height=height, width=width)
    }
    if(by.cl){
      if(sorted){
        ord = 1:length(cl)
      }
      else{
        if(!is.null(hc)){
          ord = order(cl, order(hc$order))
        }
        else{
          ord = order(cl)
        }
      }
      sep = cl[ord]
      sep=which(sep[-1]!=sep[-length(sep)])
      sep = c(sep[1], sep[which(sep[-1] - sep[-length(sep)] >=min.sep)+1])
      heatmap.3(tmp.dat[,ord],Rowv=as.dendrogram(gene.hc), Colv=NULL, col=col, trace="none", dendrogram="none", cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors[,ord],breaks=breaks,colsep=sep, sepcolor="black",main=main)
      cells.order=colnames(tmp.dat)[ord]
    }
    else{
      if(!is.null(hc)){
        hc = as.dendrogram(hc)
      }
      heatmap.3(tmp.dat,Rowv=as.dendrogram(gene.hc), Colv=hc, col=col, trace="none", dendrogram="none", cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors,breaks=breaks,main=main)
      cells.order=colnames(tmp.dat)[hc$order]
    }
    if(!is.null(prefix)){
      dev.off()
    }
    return(cells.order)
  }


display_cl_one_vs_others <- function(select.cl, 
                                     cl, 
                                     norm.dat, 
                                     de.genes,
                                     plot = !is.null(prefix),
                                     col = NULL, 
                                     max.cl.size = NULL,
                                     main = "",
                                     height = 13, 
                                     width = 9, 
                                     min.sep = 4, 
                                     ...)
  {
    select.cells=names(cl)        
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if(!is.null(max.cl.size)){    
      tmp.cells = sample_cells(cl,  max.cl.size)
      cl = cl[tmp.cells]
    }
    cl = as.factor(cl)
    select.cells=names(cl)
    tmp.col = setNames(jet.colors(length(unique(cl))), levels(cl))
    tmp.col[select.cl] = "black"
    cl.col = tmp.col[cl]
    
    tmp.col =t(as.matrix(cl.col, ncol=1))

    colnames(tmp.col)= select.cells
    if(!is.null(col)){
      tmp.col = rbind(tmp.col, col[,select.cells])
    }
    tmp.dat = as.matrix(norm.dat[,names(cl)])
    pairs = intersect(c(paste(select.cl, setdiff(levels(cl), select.cl), sep="_"),
      paste(setdiff(levels(cl), select.cl), select.cl, sep="_")), names(de.genes))
    markers = unique(unlist(sapply(de.genes[pairs], function(tmp){
      c(head(tmp$up.genes, n.markers), head(tmp$down.genes, 
                                            n.markers))
    }, simplify = F)))
    cells_order=NULL
    if(plot){
      cells_order=plot_cl_heatmap(tmp.dat, cl, markers, ColSideColors=tmp.col, prefix=prefix, labels=NULL, by.cl=TRUE,min.sep=min.sep,main=main, height=height, width=width)
    }
    return(list(markers=markers,cells_order= cells_order))
  }
  
#' Display cluster plot
#' 
#' @param cl 
#' @param norm.dat 
#' @param prefix 
#' @param plot 
#' @param col 
#' @param max.cl.size 
#' @param markers 
#' @param de.genes 
#' @param main 
#' @param height 
#' @param width 
#' @param min.sep 
#' @param ... 
#' 
#' @author Zizhen Yao
#' 
display_cl<- function(cl, norm.dat,prefix=NULL, plot=!is.null(prefix), col=NULL, max.cl.size=NULL,markers=NULL,de.genes=NULL, main="",height=13, width=9, min.sep=10, ...)
  {
    select.cells=names(cl)        
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if(!is.null(max.cl.size)){    
      tmp.cells = sample_cells(cl,  max.cl.size)
      cl = cl[tmp.cells]
    }
    select.cells=names(cl)
    cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
    tmp.col =t(as.matrix(cl.col, ncol=1))
    colnames(tmp.col)= select.cells
    if(!is.null(col)){
      tmp.col = rbind(tmp.col, col[,select.cells])
    }
    if(is.null(markers)){
      tmp = select_markers(norm.dat,cl, de.genes=de.genes, ...)
      markers = tmp$markers
      de.genes=tmp$de.genes
    }
    cells_order=NULL
    if(plot & !is.null(markers) & length(markers)>0){
      tmp.dat = as.matrix(norm.dat[markers, names(cl),drop=F])
      cells_order=plot_cl_heatmap(tmp.dat, cl, markers, ColSideColors=tmp.col, prefix=prefix, labels=NULL, by.cl=TRUE,min.sep=min.sep,main=main, height=height, width=width)
    }
    return(list(markers=markers,de.genes=de.genes, cells_order= cells_order))
  }


display_cl_markers_co.ratio <- function(select.cl, cl, norm.dat, co.ratio, prefix,  all.col, max.cl.size=100, markers=NULL,...)
{
  cells = names(cl)[cl %in% select.cl]
  if(is.factor(cl)){
    cl =droplevels(cl[cells])
  }
  if(!is.null(max.cl.size)){
    cells = sample_cells(cl, max.cl.size)
  }
  cl = cl[cells]
  if(is.null(markers)){
    markers= display_cl(norm.dat[, cells],cl, prefix=NULL,...)
  }
  hc = hclust(as.dist(1-as.matrix(co.ratio[cells, cells])),method="average")
  ord = order(cl, order(hc$order))
  cl=cl[ord]
  cells=names(cl)
  
  tmp=plot_cl_heatmap(norm.dat, cl, markers, prefix, sorted=TRUE,by.cl=TRUE,ColSideColors=all.col[,names(cl)])
  sep=which(cl[-1]!=cl[-length(cl)])
  pdf(paste0(prefix, ".co.pdf"))
  heatmap.3(as.matrix(co.ratio[cells, cells]), col = blue.red(100), trace="none", ColSideColors=all.col[,cells], Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
  dev.off()
  return(markers)
}

plot_cl_meta_barplot <- function(cluster, meta, col=NULL, drop=FALSE)
{
  library(ggplot2)
  meta = as.factor(meta)
  final.tbl <- table(cluster, meta)
  final.tbl = final.tbl/rowSums(final.tbl)
  if(drop){
    tb.df = droplevels(as.data.frame(final.tbl))
  }
  else{
    tb.df = as.data.frame(final.tbl)
  }
  g=ggplot(data=tb.df, aes(x=cluster,y=Freq,fill=meta))+ geom_bar(stat="identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),panel.grid.major=element_blank(),panel.background=element_blank())
  if(!is.null(col)){
    g=g + scale_fill_manual(values=col)
  }
  return(g)
}

plot_cl_cells <- function(anno)
  {
    max.mag=ceiling(max(log10(table(anno$cluster_id))))
    panel_pad <- 0.05
    n_clusters = max(anno$cluster_id)
    n_guides <- data.frame(y = seq(-4 - panel_pad * 4,-3 - panel_pad * 4,by = 1/10),
                           x = 0.5,
                           xend = n_clusters + 1,
                           label = seq(5, 0, by = -0.5)) %>%
                             mutate(yend = y)
    

    n_rects <- anno  %>%
      group_by(cluster_id, cluster_color, cluster_label) %>%
        summarise(n = n()) %>%
          ungroup() %>%
            mutate(adj_n = log10(n)) %>%
              mutate(xmin = cluster_id - 0.5,
                     xmax = cluster_id + 0.5,
                     ymin = -3 - panel_pad * 4 - adj_n / max.mag,
                     ymax = -3 - panel_pad * 4)
    
    g = ggplot(data = n_rects) + geom_rect( aes(xmin = xmin,
                 xmax = xmax,
                 ymin = ymin,
                 ymax = ymax,
                 fill = cluster_color)) +
                   geom_segment(data = n_guides,
                                aes(x = x,
                                    xend = xend,
                                    y = y,
                                    yend = yend),
                                linetype = "dashed") +
                                  geom_text(data = n_guides,
                                            aes(x = 0,
                                                y = y,
                                                label = label),
                                            size = 2,
                                            hjust = 1) +
                                              scale_color_identity() +
                                                scale_fill_identity() +
                                        #scale_y_continuous(limits = c(-n_clusters - 2,2))+      
                                                  theme_void()
    return(g)
  }

