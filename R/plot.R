plot_cl_heatmap <- function(norm.dat, cl, markers, prefix,hc=NULL, gene.hc=NULL,centered=FALSE,labels=names(cl),sorted=FALSE,by.cl=TRUE,ColSideColors=NULL,maxValue=5,min.sep=4,main="")
  {
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
    dendro = "column"
    cexCol = min(70/ncol(tmp.dat),1)
    cexRow = min(60/nrow(tmp.dat),1)
    if(is.null(gene.hc)){
      gene.hc = hclust(dist(tmp.dat), method="ward")
    }
    else{
      dendro="both"
    }
    if(is.null(hc) & !sorted & length(select.cells)< 2000){
      hc = hclust(dist(t(tmp.dat)), method="ward")
    }
    col = blue.red(150)[51:150]
    pdf(paste(prefix,"pdf",sep="."), height=13, width=9)
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
      heatmap.3(tmp.dat,Rowv=as.dendrogram(gene.hc), Colv=as.dendrogram(hc), col=col, trace="none", dendrogram=dendro, cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors,breaks=breaks,main=main)
      cells.order=colnames(tmp.dat)[hc$order]
    }
    dev.off()
    return(cells.order)
  }
    
display_cl<- function(cl, norm.dat,prefix=NULL, col=NULL, max.cl.size=NULL,markers=NULL,de.genes=NULL, main="",...)
  {
    select.cells=names(cl)        
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if(!is.null(max.cl.size)){    
      tmp.cells = sample_cells_by_genecounts(cl, norm.dat, max.cl.size=max.cl.size)
      cl = cl[tmp.cells]
    }
    select.cells=names(cl)
    cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
    tmp.col =t(as.matrix(cl.col, ncol=1))
    colnames(tmp.col)= select.cells
    if(!is.null(col)){
      tmp.col = rbind(tmp.col, col[,select.cells])
    }
    tmp.dat = as.matrix(norm.dat[,names(cl)])
    if(is.null(markers)){
      tmp = select_markers(tmp.dat,cl, de.genes=de.genes, ...)
      markers = tmp$markers
      de.genes=tmp$de.genes
    }
    cells_order=NULL
    if(!is.null(prefix) & !is.null(markers)){
      cells_order=plot_cl_heatmap(tmp.dat, cl, markers, ColSideColors=tmp.col, prefix=prefix, by.cl=TRUE,min.sep=10,main=main)
    }
    return(list(markers=markers,de.genes=de.genes, cells_order= cells_order))
  }


display_cl_markers_co.ratio <- function(select.cl, all.cl, norm.dat, co.ratio, prefix,  all.col, max.cl.size=100, markers=NULL,...)
{
  cells = names(all.cl)[all.cl %in% select.cl]
  if(is.factor(all.cl)){
    tmp.cl =droplevels(all.cl[cells])
  }
  if(!is.null(max.cl.size)){
    cells = unlist(tapply(names(tmp.cl),tmp.cl, function(x){sample(x, pmin(max.cl.size, length(x)))},simplify=F))
  }
  cl = tmp.cl[cells]
  if(is.null(markers)){
    markers= display_cl(norm.dat[, cells],tmp.cl, prefix=NULL,...)
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

plot_cl_meta_barplot <- function(cluster, meta, col=NULL)
{
  meta = as.factor(meta)
  final.tbl <- table(cluster, meta)
  final.tbl = final.tbl/rowSums(final.tbl)
  tb.df = droplevels(as.data.frame(final.tbl))
  g=ggplot(data=tb.df, aes(x=cluster,y=Freq,fill=meta))+ geom_bar(stat="identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),panel.grid.major=element_blank(),panel.background=element_blank())
  if(!is.null(col)){
    g=g + scale_fill_manual(values=col)
  }
  return(g)
}


