jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#' Title
#'
#' @param rd.dat 
#' @param cl 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param rd.dat 
#' @param cl 
#' @param cl.color 
#' @param cl.label 
#' @param cex 
#' @param fn.size 
#' @param alpha.val 
#' @param show.legend 
#' @param legend.size 
#' @param label.center 
#' @param bg 
#' @param fn.color 
#'
#' @return
#' @export
#'
#' @examples
plot_RD_cl <- function(rd.dat, cl, cl.color, cl.label,cex=0.15, fn.size =2, alpha.val=NULL,show.legend=FALSE, legend.size=2, label.center=TRUE, bg="blank",fn.color="black",no.shape=TRUE,ncol=4,shift.x=0, shift.y=0)
  {
    rd.dat=as.data.frame(rd.dat)
    colnames(rd.dat) = paste0("Dim", 1:ncol(rd.dat))
    rd.dat$cl = factor(cl[row.names(rd.dat)])
    if(label.center){
      cl.center = get_RD_cl_center(rd.dat, cl)
    }
    if(!no.shape){
      shape = setNames(1:length(levels(rd.dat$cl)) %% 20 + 1,levels(rd.dat$cl))
      g=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=cl,shape=cl),size=cex)
      g = g+ scale_shape_manual(values=as.vector(shape[levels(rd.dat$cl)]))      
    }
    else{
      g=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=cl),size=cex)
    }
    if(!is.null(alpha.val)){
      col = alpha(as.vector(cl.color[levels(rd.dat$cl)]),alpha.val)
    }
    else{
      col = as.vector(cl.color[levels(rd.dat$cl)])
    }
    g = g+ scale_color_manual(values=col,labels=cl.label[levels(rd.dat$cl)])
    if(label.center){
      g = g + geom_point(data=as.data.frame(cl.center), aes(x=x, y=y), size=cex*1.5)
      for(i in 1:nrow(cl.center)){
        g = g +  annotate("text", label=cl.label[row.names(cl.center)[i]], x=cl.center[i,1]+shift.x, y=cl.center[i,2] + shift.y,size=fn.size,color=fn.color)
      }
    }
    if(bg=="blank"){
      g = g + theme_void()
      #g = g + theme(panel.background=element_blank())
      #g = g + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    }
    else{
      g = g + theme(panel.background= element_rect(fill=bg, color=NA), panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor= element_blank())
    }
    if(show.legend){
      if(no.shape){
         g = g +  guides(colour = guide_legend(override.aes = list(size = legend.size),ncol=ncol))
      }
      else{
        g = g +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(rd.dat$cl)],size = legend.size)),ncol=ncol)
      }
      g = g + theme(legend.position="bottom")
    }
    else{
      g = g + theme(legend.position="none")
    }
    g = g + coord_fixed(ratio=1)
    return(g)
  }



 
###meta is discretized. 
#' Title
#'
#' @param rd.dat 
#' @param meta 
#' @param meta.col 
#' @param show.legend 
#' @param cex 
#' @param legend.size 
#' @param alpha.val 
#'
#' @return
#' @export
#'
#' @examples
plot_RD_meta <- function(rd.dat, meta, meta.col=NULL,show.legend=TRUE, cex=0.15, legend.size=5,alpha.val=1)
  {
    rd.dat = as.data.frame(rd.dat)
    colnames(rd.dat)[1:2] = c("Dim1","Dim2")
    library(ggplot2)
    rd.dat$meta = meta
    p=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=meta),size=cex)
    if(is.factor(meta)){
      rd.dat = droplevels(rd.dat)
      if(is.null(meta.col)){
        if(length(levels(meta)) > 2){
          meta.col = setNames(jet.colors(length(levels(meta))), levels(meta))
        }
        else{
          meta.col = setNames(c("blue", "orange"), levels(meta))
        }
      }      
      p = p+ scale_color_manual(values=alpha(as.vector(meta.col[levels(rd.dat$meta)]),alpha.val))
      p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    }
    else{
      p = p+ scale_color_gradient(low="blue",high="red")
    }
    if(!show.legend){
      p = p + theme(legend.position="none") 
    }
    else{
      if(is.factor(meta)){
        p = p + guides(colour = guide_legend(override.aes = list(size=legend.size)))
      }
    }
    p = p + coord_fixed(ratio=1)
    return(p)
  }


#' Title
#'
#' @param rd.dat 
#' @param norm.dat 
#' @param genes 
#' @param cex 
#'
#' @return
#' @export
#'
#' @examples
plot_RD_gene <- function(rd.dat, norm.dat, genes, cex=0.15)
  {
    library(ggplot2)
    plots=list()
    rd.dat = as.data.frame(rd.dat)
    colnames(rd.dat)[1:2] = c("Dim1","Dim2")
    for(g in genes){
      rd.dat$expr = norm.dat[g,row.names(rd.dat)]
      p=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=expr),size=cex)
      p = p+ scale_color_gradient(low="gray80",high="red") 
      p = p + theme_void() + theme(legend.position="none")
      p = p + coord_fixed(ratio=1)
      p = p + ggtitle(g)
      plots[[g]]= p
    }
    return(plots)
  }


###copy from R cookbook: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' Title
#'
#' @param ... 
#' @param plotlist 
#' @param file 
#' @param cols 
#' @param layout 
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param df 
#' @param col 
#' @param label_col 
#' @param cex 
#' @param label.cex 
#' @param init 
#' @param prefix 
#' @param bg.col 
#'
#' @return
#' @export
#'
#' @examples
plot_3d_label <- function(df, col, label_col=NULL,cex=1, label.cex=cex, init=TRUE, prefix=NULL,bg.col="gray60")
{
  library(rgl)
  if(init){
    rgl.open()
  }
  rgl.points(df$Dim1,df$Dim2, df$Dim3, col=df[[col]])
  if(!is.null(label_col)){
    cl.center = do.call("rbind", tapply(1:nrow(df), df[[label_col]], 
                                        function(x) {
                                          x = sample(x, pmin(length(x), 500))
                                          center = c(median(df$Dim1[x]), median(df$Dim2[x]), median(df$Dim3[x]))
                                          dist = as.matrix(dist(as.matrix(df[x, c("Dim1","Dim2","Dim3")])))
                                          tmp = x[which.min(rowSums(dist))]
                                          df[tmp,c("Dim1","Dim2","Dim3")]
                                        }))
    rgl.texts(cl.center$Dim1, cl.center$Dim2, cl.center$Dim3,text = row.names(cl.center), cex=label.cex)
  }
  #rgl.viewpoint(zoom=0.6)
  if(!is.null(prefix)){
    mybgplot3d(prefix,bg.col=bg.col)
  }
  else{
    bg3d(col=bg.col)
  }
}

#' Title
#'
#' @param tt 
#' @param bg.col 
#'
#' @return
#' @export
#'
#' @examples
mybgplot3d <- function(tt, bg.col=bg.col)
{
  viewport <- par3d("viewport")
  width <- viewport["width"]
  height <- viewport["height"]
  value=NULL
  if (width > 0 && height > 0) {
    filename <- tempfile(fileext = ".png")
    png(filename = filename, width = width, height = height)
    value <- try({
      plot.new()
      title(tt)
    })
    dev.off()
  }
  result <- bg3d(texture = filename, col = bg.col, lit = FALSE)
  lowlevel(structure(result, value = value))
}


#' Title
#'
#' @param df 
#' @param val 
#' @param cex 
#' @param max.val 
#' @param init 
#' @param prefix 
#' @param bg.col 
#'
#' @return
#' @export
#'
#' @examples
plot_3d_val <- function(df, val, cex=1, max.val=quantile(val, 0.99), init=TRUE, prefix=NULL,bg.col="gray60")
{
  library(rgl)
  if(init){
    rgl.open()
  }
  col = blue.red(100)[cut(val, c(seq(min(val), max.val,length.out=100), max.val+0.001),include.lowest = TRUE)]
  rgl.points(df$Dim1,df$Dim2, df$Dim3, col=col)
  if(!is.null(prefix)){
    mybgplot3d(prefix,bg.col=bg.col)
  }
  else{
    bg3d(col=bg.col)
  }
}



#' Title
#'
#' @param df 
#' @param cols 
#' @param label_cols 
#' @param cex 
#' @param label.cex 
#' @param fn 
#' @param win.dim 
#' @param layout 
#' @param bg.col 
#' @param dir 
#'
#' @return
#' @export
#'
#' @examples
plot_3d_label_multiple <- function(df, cols, label_cols, cex=0.7, label.cex=0.7, fn = NULL,win.dim = c(20,40,1200,800), layout = NULL, bg.col="gray60", dir="./")
{
  n.win = length(cols)
  library(rgl)
  print("start plotting")
  mfrow3d(1,n.win)
  next3d()
  
  open3d()
  
  ###specify dimensions of the plots
  par3d(windowRect=win.dim)
  #bg3d(bg.col)
  if(is.null(layout)){
    layout <- matrix(1:n.win, nrow=1)
    layout = rbind(layout, layout)
  }
  layout3d(layout, sharedMouse = TRUE)
  
  for (i in 1:n.win) {
    next3d()
    col = cols[[i]]
    if(length(col)==1){
      tmp.col = paste0(col,"_color")
      if(tmp.col %in% colnames(df)){ 
        plot_3d_label(df, col=tmp.col, label_col=label_cols[[i]],init=FALSE, cex=cex, label.cex=label.cex, prefix=names(cols)[i], bg.col=bg.col)
      }
      else{
        plot_3d_val(df, val=df[[col]], init=FALSE,cex=cex,prefix=names(cols)[i], bg.col=bg.col)
      }
    }
    else{
      plot_3d_val(df, val=col,init=FALSE, cex=cex,prefix=names(cols)[i],bg.col=bg.col)
    }
  }
  if(!is.null(fn)){
    writeWebGL(dir=dir, filename=fn)
  }
}

#' Title
#'
#' @param rd.dat 
#' @param k 
#' @param th 
#'
#' @return
#' @export
#'
#' @examples
clean_outliers <- function(rd.dat, k=10, th=6)
{
  knn.result = nn2(rd.dat, k=k)
  knn.dist = knn.result[[2]][,-1]
  knn.dist.mean = rowMeans(knn.dist)
  outlier = knn.dist.mean - median(knn.dist.mean) > th * mad(knn.dist.mean)
  return(list(mean.dist=knn.dist.mean, outlier=outlier))
}

#' Title
#'
#' @param rd.dat 
#' @param select.cells 
#' @param fg.col 
#' @param bg.col 
#' @param fg.alpha 
#' @param bg.alpha 
#' @param cex 
#'
#' @return
#' @export
#'
#' @examples
plot_2d_select <- function(rd.dat, select.cells, fg.col=  "red", bg.col="gray", fg.alpha=1, bg.alpha=0.5,cex=0.15)
{
  meta = factor(row.names(rd.dat) %in% select.cells)
  levels(meta) = c("bg","fg")
  meta.col = meta.col=c(fg=fg.col, bg=bg.col)
  alpha.val=c(fg=fg.alpha,bg=bg.alpha)
  plot_RD_meta(rd.dat, meta=meta, meta.col=meta.col, alpha.val = alpha.val, show.legend=FALSE, cex=cex)
}

#' Title
#'
#' @param rd.dat 
#' @param select.cells 
#' @param fg.col 
#' @param bg.col 
#' @param fg.alpha 
#' @param bg.alpha 
#' @param cex 
#' @param web.fn 
#' @param web.dir 
#'
#' @return
#' @export
#'
#' @examples
plot_3d_select <- function(rd.dat, select.cells, fg.col=  "red", bg.col="gray", fg.alpha=1, bg.alpha=0.5,cex=0.15, web.fn=NULL, web.dir="./")                           
  {
    meta = factor(row.names(rd.dat) %in% select.cells)
    levels(meta) = c("bg","fg")
    meta.col = alpha(meta.col=c(fg=fg.col, bg=bg.col), alpha.val=c(fg=fg.alpha,bg=bg.alpha))
    df = as.data.frame(rd.dat)
    df$select = meta
    df$col  = meta.col[meta]
    rgl.open()
    rgl.points(df$Dim1,df$Dim2, df$Dim3, col=df$col)
    if(!is.null(web.fn)){
      writeWebGL(dir=dir, filename=fn)
    }
  }
  

plot_RD_cl_subset<- function(rd.dat, cl, cl.color,cl.label,select.samples,missing.color="gray85",min.size=10,fg.alpha=1,bg.alpha=0.5,...)
  {
    cl= setNames(as.character(cl), names(cl))
    cl = cl[names(cl)%in% select.samples]
    tmp = setdiff(row.names(rd.dat),names(cl))
    cl = c(cl, setNames(rep("0",length(tmp)), tmp))
    cl.size= table(cl)
    cl.small = names(cl.size)[cl.size < min.size]
    cl.label[c(cl.small,"0")] = " "
    cl.color["0"] = missing.color
    alpha.val = setNames(rep(fg.alpha, length(cl.color)),names(cl.color))
    alpha.val["0"] = bg.alpha    
    cl.color = alpha(cl.color, alpha.val)
    rd.rd = rd.dat[order(row.names(rd.dat) %in% select.samples),]
    plot_RD_cl(rd.dat, cl, cl.color, cl.label,...)    
  }
  



plot_2d_umap_anno <- function(umap.fn, anno.df, dest.d="./",meta.fields=c("platform","joint_region"),alpha=0.5)
  {
    library(data.table)
    library(dplyr)
    library(ggplot2)
    umap.df <- as.data.frame(fread(umap.fn,header=TRUE))
    umap.df = umap.df[sample(1:nrow(umap.df)),]
    colnames(umap.df) = c("sample_name","Dim1","Dim2")
    umap.df = umap.df[sample(1:nrow(umap.df)),]
    umap.df = umap.df %>% left_join(anno.df) 
    umap.2d = umap.df[,c("Dim1","Dim2")]
    row.names(umap.2d)=umap.df$sample_name
    umap.fn = basename(umap.fn)
    cl = setNames(umap.df$cl, umap.df$sample_name)
    cl.df = umap.df %>% select(cluster_id, cluster_label, cluster_color,cl) %>% unique
    cl.color = setNames(cl.df$cluster_color, cl.df$cl)
    cl.label = setNames(cl.df$cluster_label, cl.df$cl)
    g= plot_RD_cl(umap.2d, cl, cl.color = cl.color, cl.label =cl.label,alpha=alpha)
    ggsave(g, file=file.path(dest.d, gsub(".csv",".pdf",umap.fn)))
    ggsave(g, file=file.path(dest.d, gsub(".csv",".png",umap.fn)))
    g= plot_RD_cl(umap.2d, cl, cl.color = cl.color, cl.label =cl.label,alpha=alpha,label.center=FALSE)
    ggsave(g, file=file.path(dest.d, gsub(".csv",".no.label.png",umap.fn)))
    
    for(m in meta.fields){
      tmp.df = umap.df[,paste0(m, c("_id","_label","_color"))] %>% unique
      colnames(tmp.df)=c("id","label","color")
      tmp.df = tmp.df %>% arrange(id)
      tmp.color = setNames(as.character(tmp.df$color), tmp.df$label)
      g= plot_RD_meta(umap.2d, factor(umap.df[,paste0(m, "_label")], levels=names(tmp.color)),meta.col = tmp.color,alpha=0.5)
      ggsave(g, file=file.path(dest.d, gsub("csv",paste0(m,".pdf"),umap.fn)),height=7,width=8.5)
      ggsave(g, file=file.path(dest.d, gsub("csv",paste0(m,".png"),umap.fn)),height=7,width=8.5)
    }
    return(umap.2d)
  }





#dest.d = "Manuscript/common/umap_constellation/"
#for(umap.fn in dir(dest.d, pattern="umap.2d.sampled.csv")){
#  plot_2d_umap(umap.fn, dest.d)
#}
plot_3d_umap_anno <- function(umap.fn, dest.d, anno.df, cols= c("region_color","cluster_color"), label_cols=list(NULL, "cluster_label"), win.dim=c(20,40,1500, 800), cex=0.7, html.fn = NULL)
  {
    umap.3d <- as.data.frame(fread(file.path(dest.d,umap.fn),header=TRUE))
    colnames(umap.3d) = c("sample_name", paste0("Dim", 1:(ncol(umap.3d)-1)))
    umap.3d <- umap.3d %>% left_join(anno.df)
    plot_3d_label_multiple(umap.3d, cols=cols, label_cols=label_cols, cex=cex, win.dim=win.dim, fn = html.fn, dir=dest.d)
  }

