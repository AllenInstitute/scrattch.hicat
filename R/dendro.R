#require(dendextend)
#require(dplyr)      

#' Pv clust show significant gradient
#'
#' @param dend 
#' @param pvclust_obj 
#' @param signif_type 
#' @param signif_col_fun 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
pvclust_show_signif_gradient <- function (dend, pvclust_obj, signif_type = c("bp", "au"), signif_col_fun = colorRampPalette(c("black", 
    "darkred", "red")), ...) 
{
  require(dendextend)
  require(dplyr)
  signif_type <- match.arg(signif_type)
  pvalue_per_node <- pvclust_obj$edges[[signif_type]]
  ord <- rank(get_branches_heights(dend, sort = FALSE))
  pvalue_per_node <- pvalue_per_node[ord]
  signif_col <- signif_col_fun(100)
  pvalue_by_all_nodes <- rep(NA, nnodes(dend))
  ss_leaf <- which_leaf(dend)
  pvalue_by_all_nodes[!ss_leaf] <- pvalue_per_node
  pvalue_by_all_nodes <- na_locf(pvalue_by_all_nodes)
  the_cols <- signif_col[round(pvalue_by_all_nodes * 100)]
  signif_lwd = seq(0.5,2,length.out=100)
  the_lwds = signif_lwd[round(pvalue_by_all_nodes * 100)]
  dend= dend %>% assign_values_to_branches_edgePar(the_cols, "col") %>% assign_values_to_branches_edgePar(the_lwds, "lwd") %>% assign_values_to_branches_edgePar(pvalue_by_all_nodes, "conf") 
}

#' Build dend
#'
#' @param cl.dat 
#' @param cl.cor 
#' @param l.rank 
#' @param l.color 
#' @param nboot 
#' @param ncores 
#'
#' @return
#' @export
#'
#' @examples
build_dend <- function(cl.dat, cl.cor=NULL, l.rank=NULL, l.color=NULL, nboot=100, ncores=1)
  {
    require(dendextend)
    require(dplyr)
    if(is.null(cl.cor)){
      cl.cor = cor(cl.dat)
    }
    pvclust.result=NULL
    if(nboot > 0){
      require(pvclust)
      parallel= FALSE
      if(ncores > 1){
        parallel = as.integer(ncores)
      }
      pvclust.result <- pvclust::pvclust(cl.dat, method.dist = "cor" ,method.hclust = "average", nboot=nboot, parallel=parallel)
      dend = as.dendrogram(pvclust.result$hclust)
      dend = label_dend(dend)$dend
      dend = dend %>% pvclust_show_signif_gradient(pvclust.result, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","gray","darkred","black")))
      #%>% pvclust_show_signif(pvclust.result, signif_type="bp", signif_value=c(2,1))
    }
    else{
      cl.hc = hclust(as.dist(1-cl.cor),method="average")      
      dend = as.dendrogram(cl.hc)
    }
    dend =  dendextend::set(dend,"labels_cex", 0.7)
    if(!is.null(l.color)){
      dend = dendextend::set(dend, "labels_col", l.color[labels(dend)])
    }
    dend = dend %>% dendextend::set("leaves_pch", 19) %>% dendextend::set("leaves_cex", 0.5)
    if(!is.null(l.color)){
      dend = dendextend::set(dend, "leaves_col", l.color[labels(dend)])
    }
    if(!is.null(l.rank)){
      dend =reorder_dend(dend,l.rank)
    }
    return(list(dend=dend, cl.cor=cl.cor, pvclust.result=pvclust.result))
  }


#' Unbranch by conf
#'
#' @param dend 
#' @param conf.th 
#'
#' @return
#' @export
#'
#' @examples
unbranch_by_conf  <- function(dend, conf.th)
  {
    if(length(dend)>1){
      conf = c()
      for(i in 1:length(dend)){
        if(is.null(attr(dend[[i]],"edgePar"))){
          conf[i]=1
        }
        else{
          conf[i] = attr(dend[[i]],"edgePar")$conf
        }
        dend[[i]]=unbranch_by_conf(dend[[i]],conf.th)
      }      
      select = conf < conf.th
      select.children = which(select )      
      if(length(select.children)>0){
        unchanged = which(!select)
        new_dend = dend[unchanged]
        names(new_dend)= unchanged
        for(i in select.children){
          if(length(dend[[i]])>1){
            for(j in 1:length(dend[[i]])){  ###make sure that no more than 100 chidren
              ind = sprintf("%02d",j)
              attr(dend[[i]][[j]], "edgePar") = attr(dend[[i]], "edgePar")
              new_dend[[paste(i,ind,sep=".")]] = dend[[i]][[j]]
            }
          }
          else{
            new_dend[[as.character(i)]]= dend[[i]]
          }
        }
        new_dend = new_dend[order(as.integer(names(new_dend)))]
        #cat("Unbranch", attr(dend, "label"), ":", names(new_dend),"\n")
        class(new_dend)= 'dendrogram'
        attr(new_dend, "height")  = attr(dend, "height")
        attr(new_dend, "members") = attr(dend, "members")
        attr(new_dend, "midpoint")= attr(dend, "midpoint")
        attr(new_dend, "edgePar") = attr(dend, "edgePar")
        attr(new_dend, "label") = attr(dend, "label")
        dend= new_dend
      }
    }
    return(dend)
  }


#' Prune dendrogram
#'
#' @param dend 
#' @param rm.labels 
#' @param top.level 
#'
#' @return
#' @export
#'
#' @examples
prune_dend <- function(dend, rm.labels, top.level=TRUE)
  {
    if(length(dend)>1){
      new_dend = list()
      for(i in 1:length(dend)){
        new_dend[[i]]=prune_dend(dend[[i]],rm.labels, top.level=FALSE)
      }
      new_dend = new_dend[!sapply(new_dend, is.null)]
      if(length(new_dend)>1){
        member = sum(sapply(new_dend, function(x)attr(x, "member")))
        class(new_dend)= 'dendrogram'
        attr(new_dend, "height")  = attr(dend, "height")
        attr(new_dend, "members") = member
        attr(new_dend, "edgePar") = attr(dend, "edgePar")
        attr(new_dend, "label") = attr(dend, "label")
        dend= new_dend
      }
      else if(length(new_dend)==0){
        dend=NULL
      }
      else if(length(new_dend)==1){
        dend = new_dend[[1]]
      }
    }
    else{
      if(labels(dend) %in% rm.labels){
        cat("Remove nodes",labels(dend),"\n")
        dend=NULL
      }
    }
    if(top.level & !is.null(dend)){      
      dend = collapse_branch(dend)
    }
    return(dend)  
  }
  
#' Reorder dendrogram
#'
#' @param dend 
#' @param l.rank 
#' @param top.level 
#'
#' @return
#' @export
#'
#' @examples
reorder_dend <- function(dend, l.rank, top.level=TRUE)
  {
    tmp.dend = dend
    sc=sapply(1:length(dend), function(i){
      l = dend[[i]] %>% labels
      mean(l.rank[dend[[i]] %>% labels])
    })
    print(sc)
    ord = order(sc)
    if(length(dend)>1){
      for(i in 1:length(dend)){
        if(ord[i]!=i){
          dend[[i]]= tmp.dend[[ord[i]]]
        }
        if(length(dend[[i]])>1){
          dend[[i]]=reorder_dend(dend[[i]],l.rank, top.level=FALSE)
        }
      }
    }
    if(top.level){
      dend = collapse_branch(dend, 10^-10)
    }
    return(dend)
  }


#' Unbranch by length
#'
#' @param dend 
#' @param length.th 
#'
#' @return
#' @export
#'
#' @examples
unbranch_by_length <- function(dend, length.th)
  {
    if(length(dend)>1){
      for(i in 1:length(dend)){
        dend[[i]]=unbranch_by_length(dend[[i]],length.th)
      }
      child.h = get_childrens_heights(dend)
      h = attr(dend, "height")
      select=h - child.h < length.th
      select.children = which(select )
      if(length(select.children)>0){
        unchanged = which(!select)
        new_dend = dend[unchanged]
        names(new_dend)= unchanged
        for(i in select.children){
          if(length(dend[[i]])>1){
            for(j in 1:length(dend[[i]])){  ###make sure that no more than 100 chidren
              ind = sprintf("%02d",j)
              idx= paste(i,ind,sep=".")
              new_dend[[idx]] = dend[[i]][[j]]
              attr(new_dend[[idx]], "edgePar") = attr(dend[[i]], "edgePar")
            }
          }
          else{
            new_dend[[as.character(i)]]= dend[[i]]
          }
        }
        new_dend = new_dend[order(as.integer(names(new_dend)))]
        class(new_dend)= 'dendrogram'
        attr(new_dend, "height")  = attr(dend, "height")
        attr(new_dend, "members") = attr(dend, "members")
        attr(new_dend, "midpoint")= attr(dend, "midpoint")
        attr(new_dend, "edgePar") = attr(end, "edgePar")
        dend= new_dend
      }
    }
    return(dend)
  }


#' Cut tree dendrogram
#'
#' @param dend 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
cutree_dend <- function(dend, h)
  {
    if(length(dend)>1){
      for(i in 1:length(dend)){
        dend[[i]]=cutree_dend(dend[[i]],h)
      }
      child.h = get_childrens_heights(dend)
      select= child.h < h
      
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
        new_dend = new_dend[order(as.integer(names(new_dend)))]
        class(new_dend)= 'dendrogram'
        attr(new_dend, "height")  = attr(dend, "height")
        attr(new_dend, "members") = attr(dend, "members")
        attr(new_dend, "midpoint")= attr(dend, "midpoint")
        attr(new_dend, "edgePar") = attr(end, "edgePar")
        dend= new_dend
      }
    }
    return(dend)
  }




#' Label dendrogram
#'
#' @param dend 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
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


reset_dend_label <- function(dend)
  {  
    if(length(dend)>1){
      attr(dend, "label") = NULL
      for(i in 1:length(dend)){
        dend[[i]] = reset_dend_label(dend[[i]])
      }
    }
    return(dend)
  }


#' Get dendrogram parent
#'
#' @param dend 
#'
#' @return
#' @export
#'
#' @examples
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

#' get dendrogram precomputed markers
#'
#' @param dend 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
get_dend_markers <- function(dend, n = 20)
  {
    if(length(dend)>1){
      m = head(names(sort(-attr(dend, "markers"))), n)
      markers=list()   
      markers[[attr(dend, "label")]]= m
      for(i in 1:length(dend)){
        markers = c(markers, get_dend_markers(dend[[i]]))
      }
      return(markers)
    }
    return(NULL)
  }

####get dendrogram precomputed markers
#' get dendrogram precomputed markers (direction)
#'
#' @param dend 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
get_dend_markers_direction <- function(dend, n = 20)
  {
    if(length(dend)>1){
      m = unique(unlist(lapply(attr(dend, "markers.byCl"), head, n)))
      markers=list()   
      markers[[attr(dend, "label")]]= m
      for(i in 1:length(dend)){
        markers = c(markers, get_dend_markers(dend[[i]]))
      }
      return(markers)
    }
    return(NULL)
  }



#' Plot dendrogram
#'
#' @param dend 
#' @param dendro_data 
#' @param node_size 
#' @param r 
#'
#' @return
#' @export
#'
#' @examples
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


#' Dendrogram list
#'
#' @param dend 
#'
#' @return
#' @export
#'
#' @examples
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

#' Dend lca
#'
#' @param dend 
#' @param l1 
#' @param l2 
#' @param l 
#'
#' @return
#' @export
#'
#' @examples
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





#' Dendrogram match
#'
#' @param dend.list 
#' @param cl.group 
#'
#' @return
#' @export
#'
#' @examples
dend_match <- function(dend.list, cl.group){
  dend_group=sapply(dend.list, function(d){
    table(cl.group[labels(d)])
  })
  dend_group.df = as.data.frame(as.table(dend_group))
  colnames(dend_group.df) = c("group","node","intersect")
  dend_group.df = dend_group.df %>% filter(intersect > 0)
  group_size=table(cl.group)
  dend_size = sapply(dend.list, function(x)length(labels(x)))
  dend_group.df$group_size = group_size[as.character(dend_group.df$group)]
  dend_group.df$dend_size = dend_size[as.character(dend_group.df$node)]
  dend_group.df$cl.union = sapply(1:nrow(dend_group.df), function(i){
    cl1= names(cl.group)[cl.group == dend_group.df[i,"group"]]
    cl2 = labels(dend.list[[as.character(dend_group.df[i,"node"])]])
    length(union(cl1, cl2))
  })
  dend_group.df$jaccard = dend_group.df$intersect/dend_group.df$cl.union
  group.match = with(droplevels(dend_group.df), tapply(1:nrow(dend_group.df), group, function(x)x[which.max(jaccard[x])]))
  group.dend.match = dend_group.df[group.match,]
  return(group.dend.match)
}

