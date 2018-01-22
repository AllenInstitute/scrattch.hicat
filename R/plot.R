plot_cl_heatmap <- function(norm.dat, cl, markers, prefix,hc=NULL, gene.hc=NULL,centered=FALSE,labels=names(cl),sorted=FALSE,by.cl=TRUE,ColSideColors=NULL,maxValue=5,min.sep=4,main="")
  {
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
    }
    else{
      heatmap.3(tmp.dat,Rowv=as.dendrogram(gene.hc), Colv=as.dendrogram(hc), col=col, trace="none", dendrogram=dendro, cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors,breaks=breaks,main=main)
    }
    dev.off()        
  }
    
display_cl<- function(cl, norm.dat,prefix, col=NULL, max.cl.size=NULL,markers=NULL,low.th=1,de.genes=NULL, main="",method="limma", de.param = de.param())
  {
    select.cells=names(cl)        
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if(!is.null(max.cl.size)){    
      if(is.matrix(norm.dat)){
        cell.gene.counts= colSums(norm.dat[,select.cells]>0)
      }
      else{
        cell.gene.counts= Matrix::colSums(norm.dat[,select.cells]>0)
      }
      cell.weights = cell.gene.counts - min(cell.gene.counts)+500
      
      tmp.cells = unlist(tapply(names(cl),cl, function(x){
        if(length(x)>max.cl.size){
          x= sample(x, max.cl.size,prob=cell.weights[x])
        }
        x
      },simplify=FALSE))
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
      tmp = select_markers(tmp.dat,cl, de.genes=de.genes,method=method, de.param= de.param)
      markers = tmp$markers
      de.genes=tmp$de.genes
      if(!is.null(prefix) & !is.null(markers)){
        markers = plot_cl_heatmap(tmp.dat, cl, markers, ColSideColors=tmp.col, prefix=prefix, by.cl=TRUE,min.sep=10,main=main)
      }
    }
    else{
      markers = plot_cl_heatmap(tmp.dat, cl, markers, ColSideColors=tmp.col, prefix=prefix, by.cl=TRUE,min.sep=10,main=main)
    }
    return(list(markers=markers,de.genes=de.genes))
  }

plot_cl_meta_barplot <- function(cluster, meta, ord = NULL,col=NULL)
{
  final.tbl <- table(cluster, meta)
  if(!is.null(ord)){
    final.tbl= final.tbl[ord,]
  }
  final.tbl = final.tbl/rowSums(final.tbl)
  tb.df = droplevels(as.data.frame(final.tbl))
  g=ggplot(data=tb.df, aes(x=cluster,y=Freq,fill=cov.compare))+ geom_bar(stat="identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),panel.grid.major=element_blank(),panel.background=element_blank())
  if(!is.null(col)){
    g=g + scale_fill_manual(values=col)
  }
  return(g)
}



displayClGroup <- function(select.cl, all.cl, norm.dat, co.ratio, prefix, n.markers=10, all.col, default.markers=NULL,markers=NULL,de.df=NULL, max.cl.size=200, rm.eigen=NULL, rm.th=0.7,...)
{
  cells = names(all.cl)[all.cl %in% select.cl]
  tmp.cl =droplevels(all.cl[cells])
  cells = unlist(tapply(names(tmp.cl),tmp.cl, function(x){sample(x, pmin(max.cl.size, length(x)))},simplify=F))
  cl = tmp.cl[cells]
  if(is.null(markers)){
    de.df = DE.genes.pw(as.matrix(norm.dat[, cells]),cl)
    markers= displayCl(norm.dat[, cells],tmp.cl, de.df=de.df, rm.eigen=rm.eigen, rm.th=rm.th, n.markers=n.markers,default.markers=default.markers,min.sep=4,...)
  }
  hc = hclust(as.dist(1-as.matrix(co.ratio[cells, cells])),method="average")
  ord = order(cl, order(hc$order))
  
  markers= displayCl(as.matrix(norm.dat[, cells]),cl, de.df=NULL, ColSideColors=all.col[,cells], by.cl=TRUE,min.sep=4,main=paste(levels(tmp.cl), collapse=" "),markers=markers,hc=hc,prefix=prefix,...)
  ord = order(cl, order(hc$order))
  sep = cl[ord]
  sep=which(sep[-1]!=sep[-length(sep)])
  pdf(paste0(prefix, ".co.pdf"))
  heatmap.3(as.matrix(co.ratio[cells, cells])[ord,ord], col = blue.red(100), trace="none", ColSideColors=all.col[,cells][,ord], Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
  dev.off()
  if(!is.null(de.df)){
    for(x in names(de.df)){
      write.csv(de.df[[x]], file=file.path("DEX",paste0(x,".csv")),quote=F)
    }
  }
  return(list(de.df=de.df, markers=markers))
}




###meta should be discretized vector.
plotTSNEMeta <- function(tsne.df, meta, meta.col=NULL,show.legend=TRUE)
  {
    tsne.df$meta = as.factor(meta)
    if(is.null(meta.col)){
      meta.col = setNames(jet.colors(length()))
    }
    cex=0.15
    p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=meta),size=cex)
    p = p+ scale_color_manual(values=as.vector(meta.col[levels(tsne.df$meta)]))
    p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    if(!show.legend){
      p = p + theme(legend.position="none")
    }
    return(p)
  }

plotTSNEGene <- function(tsne.df, norm.dat, genes)
  {
    plots=list()
    cex=0.15
    for(g in genes){
      tsne.df$expr = norm.dat[g,row.names(tsne.df)]
      p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=expr),size=cex)
      p = p+ scale_color_gradient(low="gray",high="red") + xlab(g)
      p = p + theme(legend.position="none")
      plots[[g]]= p
      ggsave(paste0(g, ".tsne.pdf"),p)
    }
    return(plots)
  }


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
     # Make the panel
     # ncol: Number of columns of plots
     # nrow: Number of rows needed, calculated from # of cols
     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                      ncol = cols, nrow = ceiling(numPlots/cols))
   }
  if (numPlots==1) {
      print(plots[[1]])  
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                            layout.pos.col = matchidx$col))
      }
    }
}



          


