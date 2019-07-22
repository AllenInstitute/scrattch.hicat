select_node_specific_markers <- function(dend, norm.dat, cl, n.markers=10,de.genes=NULL,up.gene.score=NULL, down.gene.score=NULL,top.n=50,max.num=2000)
{
  m=list()
  g1= labels(dend)
  g2=  setdiff(levels(cl), g1)
  all.genes = row.names(up.gene.score)
  pairs = do.call("rbind",strsplit(colnames(up.gene.score), "_"))
  row.names(pairs)= colnames(up.gene.score)
  up.pairs = row.names(pairs)[pairs[,1] %in% g1 & pairs[,2] %in% g2]
  down.pairs = row.names(pairs)[pairs[,1] %in% g2 & pairs[,2] %in% g1]
  if(length(up.pairs)>0 & length(down.pairs)>0){
    tmp.up.gene.score = cbind(up.gene.score[,up.pairs,drop=F], down.gene.score[,down.pairs,drop=F])
    tmp.down.gene.score = cbind(down.gene.score[,up.pairs,drop=F], up.gene.score[,down.pairs,drop=F])
    up.genes = row.names(tmp.up.gene.score)[head(order(rowSums(tmp.up.gene.score)), n.markers)]
    up.genes = up.genes[rowSums(tmp.up.gene.score[up.genes,,drop=F] < max.num) > 0]
    if(length(up.genes)>0){
      up.num = colSums(tmp.up.gene.score[up.genes,,drop=F] < max.num)
      markers=up.genes
      m[[attr(dend,"label")]]=markers
    }
  }
  if(length(dend)>1){
    for(i in 1:length(dend)){
      m=c(m,select_node_specific_markers(dend[[i]],norm.dat, cl, de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,top.n=top.n, max.num=max.num, n.markers=n.markers))
    }
  }
  return(m)
}

select_dend_markers <- function(dend, cl, de.genes,norm.dat=NULL,up.gene.score=NULL, down.gene.score=NULL,...)
{
  require(dendextend)
  require(randomForest)
  print(dend)
  if(length(dend)>1){
    cl.g = sapply(dend, labels,simplify=F)
    markers = c()
    for(i in 1:(length(cl.g)-1)){
      for(j in (i+1):length(cl.g)){
        g = select_markers_pair_group(cl=cl, cl.g[[i]],cl.g[[j]],de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)
        markers=union(markers, g)          
      }
    }
    tmp.cl = setNames(rep(1:length(cl.g),sapply(cl.g,length)),unlist(cl.g))
    select.cells = names(cl)[cl %in% names(tmp.cl)]
    if(length(markers)==0){
        next
    }
    if(!is.null(norm.dat)){
        select.cells = unlist(tapply(select.cells, droplevels(cl[select.cells]),function(x)sample(x,min(length(x), 50))))
        tmp.cl = setNames(tmp.cl[as.character(cl[select.cells])],select.cells)
        rf = randomForest(t(as.matrix(norm.dat[markers,select.cells])),factor(tmp.cl))
        w = importance(rf)
        attr(dend, "markers")=w[,1]
    }
    else{
        attr(dend, "markers")=setNames(rep(1,length(markers)),markers)
    }
    for(i in 1:length(dend)){
      dend[[i]] = select_dend_markers(dend[[i]], cl=cl, norm.dat=norm.dat, de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)
    }
  }
  return(dend)
}

select_pos_dend_markers <- function(dend,cl, norm.dat)
{
  library(matrixStats)
  if(length(dend)>1){
    cl.g = sapply(dend, labels,simplify=F)
    markers=names(sort(attr(dend, "markers"),decreasing=T))
    cl.g.mean = sapply(cl.g, function(x){
      tmp.cells = intersect(names(cl)[cl %in% x], colnames(norm.dat))
      Matrix::rowMeans(norm.dat[markers,tmp.cells])
    })
    colnames(cl.g.mean)=1:ncol(cl.g.mean)
    cl.g.mean.diff1 = cl.g.mean - rowMaxs(cl.g.mean)
    cl.g.mean.diff2 = cl.g.mean - rowMins(cl.g.mean)
    cl.g.list = lapply(colnames(cl.g.mean), function(x){row.names(cl.g.mean)[cl.g.mean.diff1[,x] > -1 & cl.g.mean.diff2[,x] > 1]})
    names(cl.g.list) = sapply(1:length(dend),function(i){attr(dend[[i]],"label")})
    attr(dend, "markers.byCl")= cl.g.list
    for(i in 1:length(dend)){
      dend[[i]] = select_pos_dend_markers(dend=dend[[i]],norm.dat=norm.dat, cl=cl)
    }
  }
  return(dend)
}


map_dend_markers <- function(dend.list, map.dat,select.cells,th=0.5)
{
  map.gene.num <- matrix(0, nrow=ncol(map.dat), ncol= length(dend.list))
  row.names(map.gene.num)= colnames(map.dat)
  colnames(map.gene.num)=names(dend.list)
  for(x in names(dend.list)){
    node = dend.list[[x]]
    if(length(node)>1){
      for(i in 1:length(node)){
        l = attr(node[[i]], "label")
        print(l)
        tmp.genes= intersect(attr(node, "markers.byCl")[[i]],row.names(map.dat))
        map.gene.num[select.cells, l] = colSums(map.dat[tmp.genes,select.cells]>th)
      }
    }
  }
  return(map.gene.num)
}



node_vs_sibling_markers <- function(dend.list, norm.dat, cl,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      dend = dend.list[[x]]
      print(labels(dend)) 
      all.cl =  droplevels(cl[cl %in% labels(dend)])
      if(length(dend)>1){
        do.call("rbind",sapply(1:length(dend),function(i){          
          cl.g =  labels(dend[[i]])
          df=group_specific_markers(cl.g, norm.dat, all.cl,...)
          if(!is.null(df)){
            df$node = attr(dend[[i]],"label")
            df$parent = attr(dend,"label")
          }
          df
        },simplify=F))
      }
      else{
        NULL
      }
    },simplify=F))
  }


node_specific_markers <- function(dend.list, norm.dat, cl,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      print(x)
      df=group_specific_markers(labels(dend.list[[x]]), norm.dat, cl,...)
      if(!is.null(df)){
        df$cl = x
      }
      df
    },simplify=F))
  }
