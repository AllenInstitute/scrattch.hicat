get_RD_cl_center <- function(rd.dat, cl)
{
  cl.center=do.call("rbind",tapply(1:nrow(rd.dat), cl[row.names(rd.dat)], function(x){
    x = sample(x, pmin(length(x),500))
    center  = c(median(rd.dat[x,1]), median(rd.dat[x,2]))
    dist = as.matrix(dist(rd.dat[x,1:2]))
    tmp= x[which.min(rowSums(dist))]
    c(x=rd.dat[tmp, 1], y= rd.dat[tmp,2])
  }))
}

plot_RD_cl <- function(rd.dat, cl, cl.color, cl.label,cex=0.15, fn.size =2, alpha.val=1,show.legend=FALSE, legend.size=2)
  {
    rd.dat$cl = cl[row.names(rd.dat)] 
    rd.dat$cl_label = droplevels(factor(cl.label[as.character(rd.dat$cl)]), levels=cl.label)
    cl.center = get_RD_cl_center(rd.dat, cl)
    shape = setNames(1:length(levels(rd.dat$cl_label)) %% 20 + 1,levels(rd.dat$cl_label))
    g=ggplot(rd.dat, aes(Lim1, Lim2)) + geom_point(aes(color=cl_label,shape=cl_label),size=cex)
    g = g+ scale_color_manual(values=alpha(as.vector(cl.color[levels(rd.dat$cl_label)]),alpha.val))+ scale_shape_manual(values=as.vector(shape[levels(rd.dat$cl_label)]))
    for(i in 1:nrow(cl.center)){
      g = g +  annotate("text", label=row.names(cl.center)[i], x=cl.center[i,1], y=cl.center[i,2],size=fn.size,color="black")
    }
    g = g + geom_point(data=as.data.frame(cl.center), aes(x=x, y=y), size=cex*1.5)
    g = g + theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    if(show.legend){
      g = g +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(rd.dat$cl_label)])),ncol=5)
      g <- g + guides(color = guide_legend(override.aes = list(size = legend.size)))
      g = g + theme(legend.position="bottom")
    }
    else{
      g = g + theme(legend.position="none")
    }
    return(g)
  }



 
###meta is discretized. 
plot_RD_meta <- function(rd.dat, meta, meta.col=NULL,show.legend=TRUE, cex=0.15, legend.size=5)
  {
    rd.dat = as.data.frame(rd.dat)
    colnames(rd.dat)[1:2] = c("Lim1","Lim2")
    library(ggplot2)
    rd.dat$meta = meta
    p=ggplot(rd.dat, aes(Lim1, Lim2)) + geom_point(aes(color=meta),size=cex)
    if(is.factor(meta)){
      if(is.null(meta.col)){
        if(length(levels(meta)) > 2){
          meta.col = setNames(jet.colors(length(levels(meta))), levels(meta))
        }
        else{
          meta.col = setNames(c("blue", "orange"), levels(meta))
        }
      }
      p = p+ scale_color_manual(values=as.vector(meta.col[levels(rd.dat$meta)]))
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


plot_RD_gene <- function(rd.dat, norm.dat, genes, cex=0.15)
  {
    library(ggplot2)
    plots=list()
    for(g in genes){
      rd.dat$expr = norm.dat[g,row.names(rd.dat)]
      p=ggplot(rd.dat, aes(Lim1, Lim2)) + geom_point(aes(color=expr),size=cex)
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



