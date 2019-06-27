library(ggplot2)
plot_tsne_cl <- function(norm.dat, select.genes, cl, cl.df, tsne.df = NULL, show.legend=FALSE, cex=0.15, fn.size=2, alpha.val=1, legend.size=4 , ...)
  {
    require(Rtsne)
    if(is.null(tsne.df)){
      tsne.result = Rtsne(t(as.matrix(norm.dat[select.genes,names(cl)])),...)$Y
      row.names(tsne.result)=names(cl)
      tsne.df = as.data.frame(tsne.result[names(cl),])
      colnames(tsne.df)=c("Lim1","Lim2")
    }
    tsne.df$cl = cl[row.names(tsne.df)] 
    tsne.df$cl_label = factor(cl.df[as.character(tsne.df$cl),"cluster_label"], levels=as.character(cl.df$cluster_label))
    tsne.df$cl_label= droplevels(tsne.df$cl_label)
   
    cl.center=do.call("rbind",tapply(1:nrow(tsne.df), tsne.df$cl, function(x){
      x = sample(x, pmin(length(x),500))
      center  = c(median(tsne.df[x,1]), median(tsne.df[x,2]))
      dist = as.matrix(dist(tsne.df[x,1:2]))
      tmp= x[which.min(rowSums(dist))]
      c(x=tsne.df[tmp, 1], y= tsne.df[tmp,2])
    }))
    
    row.names(cl.center)= cl.df[row.names(cl.center), "cluster_label"]
    cl.col = setNames(as.character(cl.df$cluster_color),cl.df$cluster_label)
    shape = setNames(1:length(levels(tsne.df$cl_label)) %% 20 + 1,levels(tsne.df$cl_label))
    g=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=cl_label,shape=cl_label),size=cex)
    g = g+ scale_color_manual(values=alpha(as.vector(cl.col[levels(tsne.df$cl_label)]),alpha.val))+ scale_shape_manual(values=as.vector(shape[levels(tsne.df$cl_label)]))
    for(i in 1:nrow(cl.center)){
      g = g +  annotate("text", label=row.names(cl.center)[i], x=cl.center[i,1], y=cl.center[i,2],size=fn.size,color="black")
    }
    g = g + geom_point(data=as.data.frame(cl.center), aes(x=x, y=y), size=cex*1.5)
    g = g + theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    if(show.legend){
      g = g +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(tsne.df$cl_label)])),ncol=5)
      g <- g + guides(color = guide_legend(override.aes = list(size = legend.size)))
      g = g + theme(legend.position="bottom")
    }
    else{
      g = g + theme(legend.position="none")
    }
    return(list(tsne.df=tsne.df, g=g))    
  }




###meta is discretized. 
plot_tsne_meta <- function(tsne.df, meta, meta.col=NULL,show.legend=TRUE, cex=0.15, legend.size=5)
  {
    tsne.df$meta = meta
    p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=meta),size=cex)
    if(is.factor(meta)){
      if(is.null(meta.col)){
        meta.col = setNames(jet.colors(length(levels(meta))), levels(meta))
      }
      p = p+ scale_color_manual(values=as.vector(meta.col[levels(tsne.df$meta)]))
      p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    }
    else{
      p = p+ scale_color_gradient(low="blue",high="red")
    }
    if(!show.legend){
      p = p + theme(legend.position="none") 
    }
    else{
      p = p + guides(colour = guide_legend(override.aes = list(size=legend.size)))
    }
    return(p)
  }


plot_tsne_gene <- function(tsne.df, norm.dat, genes, cex=0.15)
  {
    plots=list()
    for(g in genes){
      tsne.df$expr = norm.dat[g,row.names(tsne.df)]
      p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=expr),size=cex)
      p = p+ scale_color_gradient(low="gray",high="red") + xlab(g)
      p = p + theme(legend.position="none")
      plots[[g]]= p
    }
    return(plots)
  }


###copy from R cookbook: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
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



          



