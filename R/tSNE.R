run_tSNE <- function(norm.dat, select.genes, cl, cl.df, tsne.result = NULL,...)
  {
    library(Rtsne)
    if(is.null(tsne.result)){
      tsne.result = Rtsne(t(norm.dat[select.genes,names(cl)]),...)$Y
      row.names(tsne.result)=names(cl)
    }
    tsne.df = as.data.frame(tsne.result[names(cl),])
    tsne.df$cl = cl
    tsne.df$cl_label = factor(cl.df[as.character(tsne.df$cl),"cluster_label"], levels=as.character(cl.df$cluster_label))
    tsne.df$cl_label= droplevels(tsne.df$cl_label)
    colnames(tsne.df)[1:2]=c("Lim1","Lim2")
   
    cl.center=do.call("rbind",tapply(1:nrow(tsne.df), tsne.df$cl, function(x){
      center  = c(median(tsne.df[x,1]), median(tsne.df[x,2]))
      #i = min()
      #center= x[which.max(cell.cl.co.ratio[x, as.character(i)])]
      #c(x=median(tsne.df[center,1]), y= median(tsne.df[center,2]))
    }))
    row.names(cl.center)= cl.df[row.names(cl.center), "cluster_label"]
    cex=0.15
    cl.col = setNames(as.character(cl.df$cluster_color),cl.df$cluster_label)
    shape = setNames(1:length(levels(tsne.df$cl_label)) %% 20 + 1,levels(tsne.df$cl_label))
    p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=cl_label,shape=cl_label),size=cex)
    p = p+ scale_color_manual(values=as.vector(cl.col[levels(tsne.df$cl_label)]))+ scale_shape_manual(values=as.vector(shape[levels(tsne.df$cl_label)]))
    for(i in 1:nrow(cl.center)){
      p = p +  annotate("text", label=row.names(cl.center)[i], x=cl.center[i,1], y=cl.center[i,2],size=2,color="black")
    }
    p = p +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(tsne.df$cl_label)])),ncol=5)
    p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.position="bottom")
    
    return(list(tsne.df=tsne.df, p=p))    
  }




###meta is discretized. 
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



          



