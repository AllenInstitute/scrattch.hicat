build_dend <- function(cl.dat, l.rank, l.color,collapse.height=0.02, nboot=100)
  {
    require(dendextend)
    cl.cor = cor(cl.dat)
    cl.hc = hclust(as.dist(1-cl.cor),method="average")      
  
    dend = as.dendrogram(cl.hc)
    dend=unbranch_by_length(dend, collapse.height)
    dend =reorder_dend(dend,l.rank)
    dend = collapse_branch(dend, 0.005)    
    dend = dend %>% set("labels_cex", 0.7)
    tmp=cl.hc$labels[cl.hc$order]
    dend = dend %>% set("labels_col", l.color[tmp])
    dend = dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>% set("leaves_col", l.color[tmp])
    conf=NULL
    if(nboot > 0){
      require(pvclust)
      result <- pvclust::pvclust(cl.dat, method.dist = "cor" ,method.hclust = "average", nboot=nboot)
      dend = label_dend(dend)$dend
      label=get_nodes_attr(dend,"label", include_leaves = FALSE)
      label=label[!is.na(label)]
      bp=as.integer(result$edges$bp * 100)
      au=as.integer(result$edges$au * 100)
      conf = data.frame(label,bp, au)
      dend = dend %>% pvclust_show_signif_gradient(result, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","black"))) %>% pvclust_show_signif(result, signif_type="bp", signif_value=c(2,1))
    }
    return(list(dend=dend, cl.cor=cl.cor, conf=conf))
  }


reorder_dend <- function(dend, l.rank)
  {
    tmp.dend = dend
    sc=sapply(1:length(dend), function(i){
      l = dend[[i]] %>% labels
      mean(l.rank[dend[[i]] %>% labels])
    })
    print(sc)
    ord = order(sc)
    print(ord)
    if(length(dend)>1){
      for(i in 1:length(dend)){
        if(ord[i]!=i){
          dend[[i]]= tmp.dend[[ord[i]]]
        }
        if(length(dend[[i]])>1){
          dend[[i]]=reorder_dend(dend[[i]],l.rank)
        }
      }
    }
    return(dend)
  }


unbranch_by_length <- function(dend, tol)
  {
    if(length(dend)>1){
      for(i in 1:length(dend)){
        dend[[i]]=unbranch_by_length(dend[[i]],tol)
      }
      child.h = get_childrens_heights(dend)
      h = attr(dend, "height")
      print(h)
      print(h-child.h)
      select=h - child.h < tol 
      select.children = which(select )
      if(length(select.children)>0){
        unchanged = which(!select)
        new_dend = dend[unchanged]
        names(new_dend)= unchanged
        for(i in select.children){
          if(length(dend[[i]])>1){
            for(j in 1:length(dend[[i]])){  ###make sure that no more than 100 chidren
              ind = sprintf("%02d",j)
              new_dend[[paste(i,ind,sep=".")]] = dend[[i]][[j]]
            }
          }
          else{
            new_dend[[as.character(i)]]= dend[[i]]
          }
        }
        print(names(new_dend))
        new_dend = new_dend[order(as.integer(names(new_dend)))]
        class(new_dend)= 'dendrogram'
        attr(new_dend, "height")  = attr(dend, "height")
        attr(new_dend, "members") = attr(dend, "members")
        attr(new_dend, "midpoint")= attr(dend, "midpoint")
        dend= new_dend
      }
    }
    return(dend)
  }

label_dend <- function(dend,n=1)
  {  
    if(is.null(attr(dend,"label"))){
      attr(dend, "label") =paste0("n",n)
      n= n +1
    }
    if(length(dend)>1){
      for(i in 1:length(dend)){
        tmp = label_dend(dend[[i]], n)
        dend[[i]] = tmp[[1]]
        n = tmp[[2]]
      }
    }
    return(list(dend=dend, n))
  }


get_dend_parent <- function(dend)
  {
    if(length(dend)>1){
      p = attr(dend, "label")
      c= sapply(1:length(dend), function(i)attr(dend[[i]],"label"))
      edges = data.frame(parent=rep(p, length(c)),child=c)
      tmp = do.call("rbind", sapply(1:length(dend), function(i)get_dend_parent(dend[[i]]),simplify=F))
      edges=rbind(edges, tmp)
      return(edges)
    }
    return(NULL)
  }

####get dendrogram precomputed markers
get_dend_markers <- function(dend)
  {
    if(length(dend)>1){
      m = head(names(sort(-attr(dend, "markers"))), 20)
      markers=list()   
      markers[[attr(dend, "label")]]= m
      for(i in 1:length(dend)){
        markers = c(markers, get_dend_markers(dend[[i]]))
      }
      return(markers)
    }
    return(NULL)
  }


plot_dend <- function(dend, dendro_data=NULL,node_size=1,r=c(-0.1,1))
  {
    require(dendextend)
    require(ggplot2)
    if(is.null(dendro_data)){
      dendro_data = as.ggdend(dend)
      dendro_data$nodes$label =get_nodes_attr(dend, "label")
      dendro_data$nodes = dendro_data$nodes[is.na(dendro_data$nodes$leaf),]
    }
    node_data = dendro_data$nodes
    label_data <- dendro_data$labels
    segment_data <- dendro_data$segments
    if(is.null(node_data$node_color)){
      node_data$node_color="black"
    }
    ggplot() + 
    geom_text(data = node_data, aes(x = x, y = y, label = label,color=node_color),size=node_size,vjust = 1) +
    geom_segment(data = segment_data, aes(x=x,xend=xend,y=y,yend=yend), color="gray50") +
    geom_text(data = label_data, aes(x = x, y = -0.01, label = label, color = col),size=node_size,angle = 90, hjust = 1) +
    scale_color_identity() +
    theme_dendro() +
    scale_y_continuous(limits = r)
  }


dend_list <- function(dend)
  {
    l = list()
    l[[attr(dend, "label")]]=dend
    if(length(dend)>1){
      for(i in 1:length(dend)){
        l = c(l, dend_list(dend[[i]]))
      }
    }
    return(l)
  }

dend_lca <- function(dend, l1, l2, l=rep(attr(dend,"label"),length(l1)))
  {
    node.height=setNames(get_nodes_attr(dend, "height"),get_nodes_attr(dend, "label"))
    if(length(dend)> 1){
      for(i in 1:length(dend)){
        tmp.l = attr(dend[[i]],"label")
        cat(tmp.l, i, "\n")
        labels = get_subtree_label(dend[[i]])
        select = l1 %in% labels & l2 %in% labels
        if(sum(select)>0){
          select = which(select)[node.height[l[select]] > node.height[tmp.l]]
          l[select] = tmp.l
          l = dend_lca(dend[[i]],l1, l2,l)
        }
      }
    }
    return(l)
  }


