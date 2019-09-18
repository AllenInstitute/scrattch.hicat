get_knn_graph <- function(rd.dat, cl, k=15, knn.outlier.th=2, outlier.frac.th=0.5)
{
  knn.result = nn2(rd.dat,k=k)
  row.names(knn.result[[1]]) = row.names(knn.result[[2]])=row.names(rd.dat)
  knn  = knn.result[[1]]
  knn.dist = knn.result[[2]]
  cl.knn.dist.mean = tapply(names(cl),cl, function(x) mean(knn.dist[x,-1]))
  cl.knn.dist.sd = tapply(names(cl),cl, function(x) sd(knn.dist[x,-1]))
  cl.knn.dist.th = (cl.knn.dist.mean + knn.outlier.th * cl.knn.dist.sd)
  
  knn.dist.th=cl.knn.dist.th[as.character(cl[row.names(knn)])]
  outlier = apply(knn.dist, 2, function(x) x>  knn.dist.th)
  row.names(outlier)  = row.names(knn.dist)
  knn[outlier] = NA
  select.cells = row.names(outlier)[rowMeans(outlier) < outlier.frac.th]  
  pred.result = predict_knn(knn[select.cells,], row.names(rd.dat), cl)
  pred.prob = pred.result$pred.prob
  knn.cell.cl.counts = round(pred.prob * ncol(knn))
  knn.cl.cl.counts = do.call("rbind",tapply(row.names(pred.prob), cl[row.names(pred.prob)], function(x)colSums(knn.cell.cl.counts[x,])))
  knn.cl.df = as.data.frame(as.table(knn.cl.cl.counts))
  colnames(knn.cl.df)[1:2] = c("cl.from","cl.to")
  from.size = rowSums(knn.cl.cl.counts)
  to.size = colSums(knn.cl.cl.counts)
  total = sum(knn.cl.cl.counts)
  knn.cl.df$cl.from.total= from.size[as.character(knn.cl.df$cl.from)]
  knn.cl.df$cl.to.total = to.size[as.character(knn.cl.df$cl.to)]
  knn.cl.df = knn.cl.df[knn.cl.df$Freq > 0,]
  knn.cl.df$pval.log = knn.cl.df$odds  = 0
  for(i in 1:nrow(knn.cl.df)){
    q = knn.cl.df$Freq[i] - 1
    k = knn.cl.df$cl.from.total[i]
    m = knn.cl.df$cl.to.total[i]
    n = total - m
    knn.cl.df$pval.log[i]=phyper(q, m=m, n=n, k=k, lower.tail = FALSE, log.p=TRUE)
    knn.cl.df$odds[i] = (q + 1) / (k * m /total)
  }
  knn.cl.df$frac = knn.cl.df$Freq/knn.cl.df$cl.from.total
  knn.cl.df$cl.from.label = cl.df[as.character(knn.cl.df$cl.from),"cluster_label"]
  knn.cl.df$cl.to.label = cl.df[as.character(knn.cl.df$cl.to),"cluster_label"]
  return(list(knn.result=knn.result, pred.result=pred.result, knn.cl.df=knn.cl.df))
}


plot_constellation <- function(knn.cl.df, cl.center.df, cl.df, out.dir)
{
  library(gridExtra)
  library(sna)
  if(!file.exists(out.dir)){
    dir.create(out.dir)
  }
  ###==== cluster coordinates to plot
  cl.center.df <- as.data.frame(cl.center.df)
  cl.center.df <- tibble::rownames_to_column(cl.center.df, "cl")
  
  cl.center.df$cluster.id <- cl.df$cluster_id[match(cl.center.df$cl, cl.df$cl)]   
  cl.center.df$cl.size <- cl.df$size[match(cl.center.df$cl, cl.df$cl)]
  
  # select rows that show edges within cluster
  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to,] 
  
  #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
  cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)] 
  cl.center.df$color <- cl.df$cluster_color[match(cl.center.df$cl, cl.df$cl)]
  cl.center.df$label <- cl.df$cluster_label[match(cl.center.df$cl, cl.df$cl)]
  
  ###==== plot nodes
  
  p.nodes <- ggplot() + 
    geom_point(data=cl.center.df,alpha=0.4, aes(x=x, y=y, size=cl.size, color=color)) +
    scale_size_area(trans="sqrt",max_size=10,breaks = c(100,1000,10000,100000)) +
    scale_color_identity() + geom_point(data=cl.center.df, alpha=0.6,shape=21, aes(x=x, y=y, size=cl.size, color=color, stroke=c(edge.frac.within*3))) + geom_text(data=cl.center.df,aes(x=x, y=y, label=cluster.id),size = 1) #+ theme_void()
  p.nodes
  #ggsave("nodes.test.pdf", p.nodes, width = 15, height = 15, units="cm")
  
  
  ###==== extract node size/stroke width to replot later without scaling
  g <- ggplot_build(p.nodes)
  dots <-g[["data"]][[1]] #dataframe with geom_point size, color, coords
  dots.stroke <- g[["data"]][[2]]
  
  nodes <- left_join(cl.center.df, dots, by=c("x","y"))
  nodes <- left_join(nodes, select(dots.stroke, stroke, x,y), by=c("x","y"))
  write.csv(nodes, file=file.path(out.dir,"comb.nodes.csv"))
  
  
  ###==== prepare data for plotting of edges between nodes
  
  #filter all edges that are <5% of total for that cluster
  #knn.cl <- knn.cl.df[knn.cl.df$frac >0.05,] #1337 lines
  knn.cl <- knn.cl.df
  #from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),] 
  nodes$cl=as.numeric(as.character(nodes$cl))
  knn.cl.d$cl.from <- as.numeric(as.character(knn.cl.d$cl.from))
  knn.cl.d$cl.to <- as.numeric(as.character(knn.cl.d$cl.to))
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, size), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[13] <- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""
  
  
  #bidirectional 
  knn.cl.bid <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- subset(knn.cl.d[i:nrow(knn.cl.d),])
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from ),] 
    
    if (dim(r)[1] != 0) {
      line$Freq.to <- r$Freq
      line$node.pt.to <- r$node.pt.from
      line$frac.to <- r$frac
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
    print(i)
  }
  
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),] 
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
    print(i)
  }
  
  
  #min frac value = 0.03
  knn.cl.uni$node.pt.to <- nodes$size[match(knn.cl.uni$cl.to, nodes$cl)]
  knn.cl.uni$Freq.to <- 1
  knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  
  
  ###-==== create line segments
  
  line.segments <- knn.cl.lines[,1:2] #cl.from and cl.to
  nodes$cl <- as.numeric((as.character(nodes$cl)))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.from"="cl"))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.to"="cl"))
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines$frac,
                              frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  #from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  
  line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  #max fraction to max point size 
  line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.from
  
  line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
    
  
  ###=== create edges, exaggerated width
  
  library(sna)
  library(Hmisc)
  library(reshape2)
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*5),line.segments$node.size.from) #exxagerated width
  line.segments$ex.line.to <- pmin((line.segments$line.width.to*5),line.segments$node.size.to) #exxagerated width
  
  allEdges <- lapply(1:nrow(line.segments), edgeMaker, len = 500, curved = TRUE, line.segments=line.segments)
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^
  
  
  
  groups <- unique(allEdges$Group)
  
  poly.Edges <- data.frame(x=numeric(), y=numeric(), Group=character(),stringsAsFactors=FALSE)
  imax <- as.numeric(length(groups))
  
  for(i in 1:imax) { 
    #svMisc::progress(i)
    #svMisc::progress(i, progress.bar=TRUE)
    select.group <- groups[i]
    #print(select.group)
    select.edge <- allEdges[allEdges$Group %in% select.group,]
    
    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    
    N <- length(x)
    leftx <- numeric(N)
    lefty <- numeric(N)
    rightx <- numeric(N)
    righty <- numeric(N)
    
    ## Start point
    perps <- perpStart(x[1:2], y[1:2], w[1]/2)
    leftx[1] <- perps[1, 1]
    lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]
    righty[1] <- perps[2, 2]
    
    ### mid points
    for (ii in 2:(N - 1)) {
      seq <- (ii - 1):(ii + 1)
      perps <- perpMid(as.numeric(x[seq]), as.numeric(y[seq]), w[ii]/2)
      leftx[ii] <- perps[1, 1]
      lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]
      righty[ii] <- perps[2, 2]
    }
    ## Last control point
    perps <- perpEnd(x[(N-1):N], y[(N-1):N], w[N]/2)
    leftx[N] <- perps[1, 1]
    lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]
    righty[N] <- perps[2, 2]
    
    lineleft <- data.frame(x=leftx, y=lefty)
    lineright <- data.frame(x=rightx, y=righty)
    lineright <- lineright[nrow(lineright):1, ]
    lines.lr <- rbind(lineleft, lineright)
    lines.lr$Group <- select.group
    
    poly.Edges <- rbind(poly.Edges,lines.lr)
    
    Sys.sleep(0.01)
    cat("\r", i, "of", imax)
    
  }
  
  
  ####plot edges
  p.edges <- ggplot(poly.Edges, aes(group=Group))
  p.edges <- p.edges +geom_polygon(aes(x=x, y=y), alpha=0.2) + theme_void()
  p.edges
  
  
  #### plot all layers
  plot.all <-  ggplot()+geom_polygon(data=poly.Edges, alpha=0.2, aes(x=x, y=y, group=Group))+ 
    geom_point(data=cl.center.df,alpha=0.4, aes(x=x, y=y, size=cl.size, color=color)) +
    scale_size_area(trans="sqrt",max_size=10,breaks = c(100,1000,10000,100000)) +
    scale_color_identity() + geom_point(data=cl.center.df, alpha=0.6,shape=21, aes(x=x, y=y, size=cl.size, color=color, stroke=c(edge.frac.within*3))) + geom_text(data=cl.center.df,aes(x=x, y=y, label=cluster.id),size = 1) + theme_void() +theme(legend.position = "none")
  plot.all
  
  st=format(Sys.time(), "%Y%m%d_%H%M_")
  
  ggsave(file.path(out.dir,paste0(st,"comb.constellation.pdf")), plot.all, width = 25, height = 25, units="cm")
  
  
  #Legends for all
  plot.all <-  ggplot()+geom_polygon(data=poly.Edges, aes(x=x, y=y, group=Group))+ 
    geom_point(data=cl.center.df,alpha=0.4, aes(x=x, y=y, size=cl.size, color=color)) +
    scale_size_area(trans="sqrt",max_size=10,breaks = c(100,1000,10000,100000)) +
    scale_color_identity() + geom_point(data=cl.center.df, alpha=0.6,shape=21, aes(x=x, y=y, size=cl.size, color=color, stroke=c(edge.frac.within*3))) + geom_text(data=cl.center.df,aes(x=x, y=y, label=cluster.id),size = 1) + theme_void() #+theme(legend.position = "none")
  dot.size.legend <- cowplot::get_legend(plot.all)
  
  #plot(dot.size.legend)
  
  cl.center.df$cluster.label <- paste(cl.center.df$cluster.id, cl.center.df$label)
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  label.col <- setNames(cl.center.df$color, cl.center.df$cluster.label)
  
  cl.center <- ggplot(cl.center.df, aes(x=cluster.id, y=cl.size)) + geom_point(aes(color=cl.center.df$cluster.label))+scale_color_manual(values=as.vector(label.col[levels(cl.center.df$cluster.label)]))
  cl.center = cl.center +  guides(color = guide_legend(override.aes = list(size = 6), ncol=5))
  
  
  cl.center.legend <- cowplot::get_legend(cl.center)  
  #plot(cl.center.legend)
  
  layout_matrix <- rbind(c(1,2,2,2,2,2))  
  
  ggsave(file.path(out.dir,paste0(st,"comb.LEGEND.pdf")),marrangeGrob(list(dot.size.legend,cl.center.legend),nrow = 1, ncol=6, layout_matrix=layout_matrix),height=20,width=20)  
}


edgeMaker <- function(whichRow, len=100, line.segments, curved=FALSE){
  
  fromC <- unlist(line.segments[whichRow,c(3,4)])# Origin
  toC <- unlist(line.segments[whichRow,c(5,6)])# Terminus
  # Add curve:
  
  graphCenter <- colMeans(line.segments[,c(3,4)])  # Center of the overall graph
  bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
  distance1 <- sum((graphCenter - bezierMid)^2)
  if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
    bezierMid <- c(toC[1], fromC[2])
    }  # To select the best Bezier midpoint
  bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
  if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve

  edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
  
  #line.width.from in 100 steps to linewidth.to
 edge$fraction <- seq(line.segments$ex.line.from[whichRow], line.segments$ex.line.to[whichRow], length.out = len)

  
    #edge$Sequence <- 1:len  # For size and colour weighting in plot
  edge$Group <- paste(line.segments[whichRow, 1:2], collapse = ">")
  return(edge)
  }



  #vwline_utils
perpStart <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 1)
        }

avgangle <- function(x, y) {
    a1 <- angle(x[1:2], y[1:2])
    a2 <- angle(x[2:3], y[2:3])
    atan2(sin(a1) + sin(a2), cos(a1) + cos(a2))
}

perp <- function(x, y, len, a, mid) {
    dx <- len*cos(a + pi/2)
    dy <- len*sin(a + pi/2)
    upper <- c(x[mid] + dx, y[mid] + dy)
    lower <- c(x[mid] - dx, y[mid] - dy)
    rbind(upper, lower)    
        }

perpMid <- function(x, y, len) {
    ## Now determine angle at midpoint
    perp(x, y, len, avgangle(x, y), 2)
        }

perpEnd <- function(x, y, len) {
    perp(x, y, len, angle(x, y), 2)
}


## x and y are vectors of length 2
angle <- function(x, y) {
    atan2(y[2] - y[1], x[2] - x[1])
        }

