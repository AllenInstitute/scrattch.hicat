

#' @param knn.cl.df output of KNN.graph. Dataframe providing information about the cluster call of nearest neighbours of cells within a cluster. Required columns: "cl.from" = ancestor cl; "cl.to="descendant cl; "Freq"=frequency of cell pairs between cl.from and cl.to; "from.total"=frequency of cells in cl.from; "to.total"=frequency of cells in cl.to; "frac"=Freq/from.total
#' @param cl.center.df dataframe containing metadata and coordinates for plotting cluster centroids. Required columns: "x" = x coordinate, "y" = y coordinate, "cl" = unique cluster id that should match "cl.to" and "cl.from" columns in knn.cl.df, "cluster_color","size" = nr of cells in cluster 
#' @param out.dir location to write plotting files to
#' @param node.label Label to identify plotted nodes.
#' @param exxageration exxageration of edge width. Default is 1 (no exxageration)
#' @param edge.param what column to use for edge width
#' @param curved Wheter edges should be curved or not. Default is TRUE.
#' @param plot.parts output of intermediate files. default is FALSE.
#' @param plot.hull plot convex around cell type neighbourhood. Provide neighbourhood_id's that need to be plotted
#' @param node.dodge whether or not nodes are allowed to overlap. Default is false 
#' @param plot.height 
#' @param plot.width 
#' @param label.size 
#' @param max_size 
#' @param label_repel
#' 
#' @example_data:
#' load("~/scrattch.hicat/data/trajectory_plot.rda")
#'  knn.cl.df = tmp.knn.cl.df.3
#'  cl.center.df = tmp.cl.center.df.3[!is.na(tmp.cl.center.df.3$x),]
#' 
#' @usage plotting.trajectory <- plot_trajectory(knn.cl.df = knn.cl.df, cl.center.df = cl.center.df, out.dir = "data/Trajectory_example/plot", node.label="cluster_id" , color.label = "cluster_color", edge.param = "to.frac",node.dodge=TRUE,  label_repel=TRUE) 




plot_trajectory <- function(knn.cl.df, 
                               cl.center.df, 
                               out.dir, 
                               node.label="cluster_id", 
                               color.label = "cluster_color",
                               edge.param = "to.frac",
                               exxageration=2, 
                               curved = TRUE, 
                               plot.parts=FALSE, 
                               plot.hull = NULL, 
                               plot.height=25, 
                               plot.width=25, 
                               node.dodge=FALSE, 
                               label.size=5, 
                               max_size=10, 
                               label_repel=FALSE,
                               size.breaks = c(100,1000,10000,100000)) { 
  
  library(gridExtra)
  #library(sna)
  library(Hmisc)
  library(reshape2)
  #library(ggalt)
  library(ggforce)
  library(dplyr)
  #library(ggrepel)
  
  st=format(Sys.time(), "%Y%m%d_%H%M%S_")
  
  
  if(!file.exists(out.dir)){
    dir.create(out.dir)
  }
  ###==== Cluster nodes will represent both cluster.size (width of point) and edges within cluster (stroke of point)
  
  # select rows that have edges within cluster
  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to,] 
  
  if(nrow(knn.cl.same)>0){
    #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
    cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)] 
  }
   
  ###==== plot nodes
  cl.center.df$color.label = cl.center.df[,color.label]
  
  labels <- cl.center.df[[node.label]] 
  
  p.nodes <-   ggplot() +     
    geom_point(data=cl.center.df,
               shape=19, 
               aes(x=x, 
                   y=y, 
                   size=cluster_size, 
                   color=alpha(color.label, 0.8))) +
    scale_size_area(trans="sqrt",
                    max_size=max_size,
                    breaks = size.breaks) +
    scale_color_identity() +  
    geom_text(data=cl.center.df,
              aes(x=x, 
                  y=y, 
                  label=labels),
              size = label.size)
  
  
  #+ theme_void()
  #p.nodes
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"nodes.org.pos.pdf")), p.nodes, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  
  ###==== extract node size/stroke width to replot later without scaling
  g <- ggplot_build(p.nodes)
  dots <-g[["data"]][[1]] #dataframe with geom_point size, color, coords
  
  nodes <- left_join(cl.center.df, dots, by=c("x","y"))
  
  
  ###==== if node.dodge==TRUE new xy coords are calculated for overlapping nodes.
  
  if (node.dodge==TRUE){
    
    #<><><># make update here to convert units by scale. check geom_mark_hull code for oneliner
    
    # dodge nodes starting at center of plot moving outward 
    
    nodes$r<- (nodes$size/10)/2
    
    
    x.list <- c(mean(nodes$x), nodes$x )
    y.list <- c(mean(nodes$y), nodes$y)
    dist.test <- as.matrix(dist(cbind(x.list, y.list)))
    nodes$distance <- dist.test[2:nrow(dist.test), 1]
    nodes <- nodes[order(nodes$distance),]
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        print(paste(d1,d2))
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          print(paste(d1,d2))
          
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        print(paste(d1,d2))
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          print(paste(d1,d2))
          
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
  }
  
  nodes <- nodes[order(nodes$cluster_id),]
  
  
  
  ## when printing lines to pdf the line width increases slightly. This causes the edge to extend beyond the node. Prevent this by converting from R pixels to points. 
  conv.factor <- ggplot2::.pt*72.27/96
  
  
  ## line width of edge can be scaled to node point size 
  nodes$node.width <- nodes$size 
  
  
  if (plot.parts == TRUE) { 
    if (node.dodge == TRUE) {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.dodge.csv"))) }
    else {
      write.csv(nodes, file=file.path(out.dir,paste0(st,"nodes.csv")))
    }
  }
  
  ###==== prepare data for plotting of edges between nodes
  
  ##filter out all edges that are <5% of total for that cluster
  #knn.cl <- knn.cl.df[knn.cl.df$frac >0.05,] #1337 lines
  knn.cl <- knn.cl.df
  ##from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),] 
  
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, node.width), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[colnames(knn.cl.d)=="node.width"]<- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  #knn.cl.d$Freq.to <- ""
  #knn.cl.d$frac.to <- ""
  
  knn.cl.bid <- NULL
  
  ## bidirectional 
  # for (i in 1:nrow(knn.cl.d)) {
  #   line <- subset(knn.cl.d[i, ])
  #   r <- subset(knn.cl.d[i:nrow(knn.cl.d), ])
  #   r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from),
  #   ]
  #   if (dim(r)[1] != 0) {
  #     line$Freq.to <- r$Freq
  #     line$node.pt.to <- r$node.pt.from
  #     line$frac.to <- r$frac
  #     knn.cl.bid <- rbind(knn.cl.bid, line)
  #   }
  # }
  # 
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    print(i)
    line <- subset(knn.cl.d[i,])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),] 
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
    #print(i)
  }
  
  
  #min frac value = 0.01
  knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
  #knn.cl.uni$Freq.to <- 1
  #knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  
  ###==== create line segments
  
  line.segments <- knn.cl.lines %>% select(cl.from, cl.to)
  #nodes$cl <- as.numeric((as.character(nodes$cl)))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.from"="cl"))
  line.segments <- left_join(line.segments,select(nodes, x, y, cl), by=c("cl.to"="cl"))
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              #freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines[[edge.param]],
                              #frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  
    ##from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  
  line.segments$line.width.from <- line.segments$node.pt.from*line.segments$frac.from
  #line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  ##max fraction to max point size 
  line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.pt.from
  
  #line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
  
  
  ###=== create edges, exaggerated width
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  #line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*exxageration),line.segments$node.size.from) #exxagerated width
  #line.segments$ex.line.to <- pmin((line.segments$line.width.to*exxageration),line.segments$node.size.to) #exxagerated width
  
  
  line.segments <- na.omit(line.segments)
  
###### calculate new xend/yend -1r on trajectory 
  # calculate coordinates for radius +15% 
line.segments$rto <- (line.segments$node.size.to/2) + (line.segments$node.size.to)*0.25
line.segments$AC = sqrt((line.segments$x.from-line.segments$x.to)^2 + (line.segments$y.from- line.segments$y.to)^2)  
  
line.segments$x.ah = line.segments$x.to + line.segments$rto*(line.segments$x.from-line.segments$x.to)/line.segments$AC
line.segments$y.ah = line.segments$y.to + line.segments$rto*(line.segments$y.from-line.segments$y.to)/line.segments$AC
#  Xah = Xto + R*(Xfrom-Xto)/|AC|
#  Yah = Yto + R*(Yfrom-Yto)/|AC|
#  |AC| = sqrt((Xfrom-Xto)^2 + (Yfrom-Yto)^2)    
###############
  

  p.edges <- ggplot(line.segments) +
    geom_curve(aes(x=x.from ,y=y.from, xend=x.ah, yend=y.ah ), 
               cex= line.segments$line.width.from,
               curvature=0.35, 
               lineend = "round", #c('round', 'butt', 'square')
               linejoin = "mitre",#c('round', 'mitre', 'bevel')
               arrow=arrow(type="open",
                 length=(unit(scales::rescale(line.segments$line.width.from, 
                                              to=c(0.15, 0.5)), "cm"))),
                colour="grey60",
               alpha=0.8)
  
  #p.edges
  
  #############################
  ##                         ##
  ##        plotting         ##
  ##                         ##
  #############################
  
  labels <- nodes[[node.label]] 
  
  
  ####plot edges
  #p.edges <- ggplot(poly.Edges, aes(group=Group))
  #p.edges <- p.edges +geom_polygon(aes(x=x, y=y), alpha=0.2) + theme_void()
  
  
 

  if (!is.null(plot.hull)) {
    #### plot all layers
    plot.all <-  ggplot()+
      geom_curve( data= line.segments,
                  aes(x=x.from ,
                      y=y.from, 
                      xend=x.ah, 
                      yend=y.ah ),
                  cex= line.segments$line.width.from,
                  curvature=0.35, 
                  arrow=arrow(type="open",
                              length=(unit(scales::rescale(line.segments$line.width.from, 
                                                           to=c(0.15, 0.5)), "cm"))),
                  colour="grey60",
                  alpha=0.8) +
      geom_point(data=nodes,
                 alpha=0.8, 
                 shape=19,
                 aes(x=x, 
                     y=y, 
                     size=cluster_size, 
                     color=cluster_color)) +
      scale_size_area(trans="sqrt",
                      max_size=max_size,
                      breaks = c(100,1000,10000,100000)) +
                      #breaks = size.breaks)
      scale_color_identity()  + 
      theme_void()+ 
      geom_mark_hull(data=nodes,
                     concavity = 8,
                     radius = unit(5,"mm"),
                     aes(filter = nodes$clade_id %in% plot.hull,x, y, 
                         color=nodes$clade_color)) +
      theme(legend.position = "none")
    
    if(label_repel ==TRUE){
      plot.all <- plot.all +
        ggrepel::geom_text_repel(data=nodes,
                                 aes(x=x, 
                                     y=y, 
                                     label=labels),
                                 size = label.size,
                                 min.segment.length = Inf) 
    }
    else{
      plot.all <- plot.all +
        geom_text(data=nodes,
                  aes(x=x, 
                      y=y, 
                      label=labels),
                  size = label.size)  }
    
    #plot.all
  } else {
    #### plot all layers
    plot.all <-  ggplot()+
      geom_curve( data= line.segments,
                  aes(x=x.from ,
                      y=y.from, 
                      xend=x.ah, 
                      yend=y.ah ),
                  cex= line.segments$line.width.from,
                  curvature=0.35, 
                  arrow=arrow(type="open",
                              length=(unit(scales::rescale(line.segments$line.width.from, 
                                                           to=c(0.15, 0.5)), "cm"))),
                  colour="grey60",
                  alpha=0.8) +
      geom_point(data=nodes,
                 alpha=0.8, 
                 shape=19,
                 aes(x=x, 
                     y=y, 
                     size=cluster_size, 
                     color=cluster_color)) +
      scale_size_area(trans="sqrt",
                      max_size=max_size,
                      breaks = c(100,1000,10000,100000)) +
      #breaks = size.breaks)
      scale_color_identity() 
    if(label_repel ==TRUE){
      plot.all <- plot.all +
        ggrepel::geom_text_repel(data=nodes,
                        aes(x=x, 
                            y=y, 
                            label=labels),
                        size = label.size,
                        min.segment.length = Inf) + 
        theme_void() +
        theme(legend.position="none")
    }
    else{
      plot.all <- plot.all +
        geom_text(data=nodes,
                  aes(x=x, 
                      y=y, 
                      label=labels),
                  size = label.size) + 
        theme_void() +
        theme(legend.position="none") }
    #plot.all
  }
  
  
  
  segment.color = NA
  
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"comb.trajectory.pdf")), plot.all, width = plot.width, height = plot.height, units="cm",useDingbats=FALSE) }
  
  
  
  #############################
  ##                         ##
  ##      plot legends       ##
  ##                         ##
  #############################
  
  
  ### plot node size legend (1)
  plot.dot.legend <- ggplot()+
    geom_curve( data= line.segments,
                aes(x=x.from ,
                    y=y.from, 
                    xend=x.ah, 
                    yend=y.ah ),
                cex= line.segments$line.width.from,
                curvature=0.35, 
                arrow=arrow(type="open",
                            length=(unit(scales::rescale(line.segments$line.width.from, 
                                                         to=c(0.15, 0.5)), "cm"))),
                colour="grey60",
                alpha=0.8)+
    geom_point(data=nodes,
               alpha=0.8, 
               shape=19,
               aes(x=x, 
                   y=y, 
                   size=cluster_size, 
                   color=cluster_color)) +
    scale_size_area(trans="sqrt",
                    max_size=max_size,
                    breaks = c(100,1000,10000,100000)) +
    scale_color_identity() + 
    geom_text(data=nodes,
              aes(x=x, 
                  y=y, 
                  label=labels),
              size = label.size)+
    theme_void()
  dot.size.legend <- cowplot::get_legend(plot.dot.legend)
  
  ### plot cluster legend (3)
  cl.center.df$cluster.label <-  cl.center.df$cluster_label
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  label.col <- setNames(cl.center.df$cluster_color, cl.center.df$cluster.label)
  cl.center.df$cluster.label <- as.factor(cl.center.df$cluster.label)
  leg.col.nr <- min((ceiling(length(cl.center.df$cluster_id)/20)),5)
  
  cl.center <- ggplot(cl.center.df, 
                      aes(x=cluster_id, y=cluster_size)) + 
    geom_point(aes(color=cluster.label))+
    scale_color_manual(values=as.vector(label.col[levels(cl.center.df$cluster.label)]))+  
    guides(color = guide_legend(override.aes = list(size = 8), ncol=leg.col.nr))
  
  cl.center.legend <- cowplot::get_legend(cl.center)  
  #plot(cl.center.legend)
  
  
  ###plot legend line width (2)
  width.1 <- max(line.segments$frac.from,line.segments$frac.to) 
  width.05 <- width.1/2
  width.025 <- width.1/4
  
  
  edge.width.data <- tibble(node.width = c(1,1,1), x=c(2,2,2), y=c(5,3.5,2), line.width=c(1,0.5,0.25), fraction=c(100, 50, 25),frac.ex=c(width.1, width.05, width.025))
  edge.width.data$fraction.ex <- round((edge.width.data$frac.ex*100), digits = 0)
  
  poly.positions <- data.frame(id=rep(c(1,2,3), each = 4), x=c(1,1,2,2,1,1,2,2,1,1,2,2), y=c(4.9,5.1,5.5,4.5,3.4,3.6,3.75,3.25,1.9,2.1,2.125,1.875)) 
  
  if (exxageration !=1) {
    edge.width.legend <- ggplot()  +  
      geom_polygon(data=poly.positions, aes(x=x,y=y, group=id), fill="grey60")+
      geom_circle(data=edge.width.data, aes(x0=x, y0=y, r=node.width/2), fill="grey80", color="grey80", alpha=0.4)+ 
      scale_x_continuous(limits=c(0,3)) + 
      theme_void() +
      coord_fixed() + 
      geom_text(data=edge.width.data, aes(x= 2.7, y=y, label=fraction.ex, hjust=0, vjust=0.5)) + 
      annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node") } 
  else { edge.width.legend <- ggplot()  +  
    geom_polygon(data=poly.positions, aes(x=x,y=y, group=id), fill="grey60")+
    geom_circle(data=edge.width.data, aes(x0=x, y0=y, r=node.width/2), fill="grey80", color="grey80", alpha=0.4)+ 
    scale_x_continuous(limits=c(0,3)) + 
    theme_void() +coord_fixed() + 
    geom_text(data=edge.width.data, aes(x= 2.7, y=y, label=fraction, hjust=0, vjust=0.5)) + 
    annotate("text", x = 2, y = 6, label = "Fraction of edges \n to node")}
  
  
  #############################
  ##                         ##
  ##      save elements      ##
  ##                         ##
  #############################
  
  
  
  layout_legend <- rbind(c(1,3,3,3,3),c(2,3,3,3,3))  
  
  if (plot.parts == TRUE) {
    ggsave(file.path(out.dir,paste0(st,"comb.LEGEND.pdf")),gridExtra::marrangeGrob(list(dot.size.legend,edge.width.legend,cl.center.legend),nrow = 3, ncol=6, layout_matrix=layout_legend),height=20,width=20,useDingbats=FALSE)  }
  
  
  g2 <- gridExtra::arrangeGrob(grobs=list(dot.size.legend,edge.width.legend,cl.center.legend), layout_matrix=layout_legend)
  
  
  ggsave(file.path(out.dir,paste0(st,"trajectory.pdf")),marrangeGrob(list(plot.all,g2),nrow = 1, ncol=1),width = plot.width, height = plot.height, units="cm",useDingbats=FALSE)
  
  
}







