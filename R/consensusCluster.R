library(Matrix)
#' Title
#'
#' @param result.files 
#' @param all.cells 
#'
#' @return
#' @export
#'
#' @examples
collect_co_matrix <- function(result.files,all.cells)
  {
    subsample.cl=list()
    co.matrix= matrix(0, nrow=length(all.cells),ncol=length(all.cells))
    pr.matrix = matrix(0, nrow=length(all.cells),ncol=length(all.cells))
    row.names(co.matrix)=row.names(pr.matrix)=colnames(co.matrix)=colnames(pr.matrix)=all.cells
    for(f in result.files){
      tmp=load(f)
      cl=result$cl
      cl = cl[intersect(names(cl),all.cells)]
      pr.matrix[names(cl),names(cl)]=   pr.matrix[names(cl),names(cl)] + 1
      for(x in unique(cl)){
        y=names(cl)[cl==x]
        co.matrix[y,y]= co.matrix[y,y]+1
      }
      subsample.cl[[f]]= cl
    }
    co.ratio = co.matrix/pr.matrix
    return(list(co.ratio=co.ratio,cl.list=subsample.cl))
  }


sample_cl_list <- function(cl.list, max.cl.size=500)
  {
    select.cells=c()
    for(cl in cl.list){
      cl.size=table(cl)
      more.cl = cl[setdiff(names(cl), select.cells)]
      more.size = table(more.cl)
      add.cells= unlist(lapply(names(more.size), function(x){
        sample(names(more.cl)[more.cl==x], min(more.size[[x]],max.cl.size))
      }))
      select.cells= c(select.cells, add.cells)
    }
    return(select.cells)
  }

collect_subsample_cl_matrix <- function(norm.dat,result.files,all.cells,max.cl.size=NULL)
  {
    select.cells=c()
    cl.list=list()
    for(f in result.files){
      print(f)
      tmp=load(f)
      cl= result$cl
      test.cells = setdiff(all.cells, names(cl))
      markers=unique(result$markers)
      map.df = map_by_cor(norm.dat[markers,names(cl)],cl, norm.dat[markers,test.cells],method="means")$pred.df
      test.cl = setNames(map.df$pred.cl, row.names(map.df))
      all.cl = c(setNames(as.character(cl),names(cl)), setNames(as.character(test.cl), names(test.cl)))
      cl.list[[f]] = all.cl
    }
    if(!is.null(max.cl.size)){
      select.cells= sample_cl_list(cl.list, max.cl.size=max.cl.size)
    }
    else{
      select.cells= all.cells
    }
    cl.mat = do.call("cbind", sapply(cl.list, function(cl){
      get_cl_mat(cl[select.cells])
    },simplify=F))
    return(list(cl.list=cl.list, cl.mat = cl.mat))
  }

#' Title
#'
#' @param norm.dat 
#' @param result.files 
#' @param all.cells 
#' @param max.cl.size 
#'
#' @return
#' @export
#'
#' @examples
collect_co_matrix_sparseM <- function(norm.dat,result.files,all.cells,max.cl.size=1000)
  {
    tmp = collect_subsample_cl_matrix(norm.dat,result.files,all.cells, max.cl.size=max.cl.size)
    cl.list = tmp$cl.list
    cl.mat = tmp$cl.mat
    co.ratio = crossprod(t(cl.mat))
    co.ratio@x = co.ratio@x/length(result.files)
    return(list(co.ratio=co.ratio, cl.mat=cl.mat, cl.list=cl.list))
  }



merge_co_matrix <- function(co.ratio1, co.ratio2)
  {
    all.cells=c(row.names(co.ratio1),row.names(co.ratio2))
    tmp.co = Matrix(0,nrow= ncol(co.ratio1), ncol=ncol(co.ratio2))
    colnames(tmp.co)=colnames(co.ratio2)
    row.names(tmp.co)=row.names(co.ratio1)
    tmp.co1 = rbind(co.ratio1, t(tmp.co))
    tmp.co2 = rbind(tmp.co, co.ratio2)
    co.ratio = cbind(tmp.co1, tmp.co2)
    return(co.ratio)
  }





#' Title
#'
#' @param co.ratio cell cell co-clustering matrix
#' @param cl.list  The list of subsampled clustering results. 
#' @param norm.dat The log2 transformed normalzied expression matrix 
#' @param select.cells Cells to be clustered
#' @param all.col Color bars for plotting heatmap. Default NULL
#' @param diff.th The difference of co-clustering probablities for splitting a cluster. 
#' @param prefix Default NULL. 
#' @param method Clustering methods. Default "auto"
#' @param verbose Default FALSE
#' @param result Pre-computed clustering results used for further splitting. Default NULL. 
#' @param min.cells Minimal number of cells in a cluster. Default 4
#' @param ... Other parameters passed to merge_cl
#'
#' @return A list with cluster membership, and top pairwise marker genes. 
#' @export
#'
#' @examples
iter_consensus_clust <- function(co.ratio, cl.list, norm.dat, select.cells=colnames(co.ratio), all.col=NULL, diff.th=0.25, prefix=NULL, method=c("auto", "louvain","ward"), verbose=FALSE, de.param = de.param, result=NULL, rd.dat = NULL)
  {
    require(igraph)
    if(verbose){
      print(prefix)
    }
    if(!is.null(result)){
      markers=result$markers
      cl = setNames(as.integer(as.character(result$cl)),names(result$cl))
    }
    else{
      markers=NULL
      if(length(select.cells)  < 2 * de.param$min.cells){
        return(NULL)
      }
      if(method=="auto"){
        if (length(select.cells)> 3000){
          select.method = "louvain"
        }
        else{
          if(!is.matrix(co.ratio)){
            co.ratio = as.matrix(co.ratio[select.cells, select.cells])
          }
          select.method="ward"
        }
      }
      else{
        select.method = method
      }
      if(select.method=="ward"){
        tmp.cl = init_cut(co.ratio, select.cells, cl.list, min.cells= de.param$min.cells, th = diff.th,method=select.method)
        if(is.null(tmp.cl)){
          return(NULL)
        }
      }
      else{###louvain
        adj.mat=co.ratio[select.cells, select.cells]
        gr = graph.adjacency(adj.mat, mode="undirected",weighted=TRUE)
        comm= cluster_louvain(gr)
        rm(gr)
        gc()
        if(pass_louvain(modularity(comm), adj.mat)){
          tmp.cl = setNames(comm$membership,select.cells)
        }
        else{
          return(NULL)
        }
      }
      if(length(unique(tmp.cl))==1){
        return(NULL)
      }
      if(verbose){
        print(table(tmp.cl))
      }
      if(is.null(rd.dat)){
        rd.dat =get_cell.cl.co.ratio(tmp.cl, co.ratio)
      }
      tmp= merge_cl(norm.dat=norm.dat, cl=tmp.cl, rd.dat=rd.dat, verbose=verbose,  de.param = de.param)
      if(is.null(tmp) | !is.list(tmp)) return(NULL)
      if (length(unique(tmp$cl))==1) return(NULL)
      tmp.cl= tmp$cl
      tmp.cl=merge_cl_by_co(tmp.cl, co.ratio, diff.th)
      tmp.cl = setNames(as.integer(tmp.cl),names(tmp.cl))
      if(length(unique(tmp.cl))==1) {
        return(NULL)
      }
      cl = tmp.cl
      markers=tmp$markers
      if(verbose){
        print(table(cl))
        display_cl_markers_co.ratio(unique(cl), cl, norm.dat=norm.dat, co.ratio=co.ratio, prefix=prefix,  all.col=all.col, markers=markers)
      }
    }
    cell.cl.co.ratio = get_cl_means(co.ratio, cl)
    n.cl=max(cl)
    new.cl=cl
    for(i in sort(unique(cl))){
      tmp.prefix= paste0(prefix, ".", i)
      tmp.cells=names(cl)[cl==i]
      uncertain.cells=sum(cell.cl.co.ratio[tmp.cells, as.character(i)] < 1 - diff.th)
      if(uncertain.cells < de.param$min.cells){
        next
      }
      result= iter_consensus_clust(co.ratio=co.ratio, cl.list=cl.list, norm.dat=norm.dat, select.cells=tmp.cells,prefix=tmp.prefix, all.col=all.col, diff.th =diff.th, method=method, de.param = de.param, verbose=verbose, rd.dat=rd.dat)
      if(is.null(result)){
        next
      }
      tmp.cl= result$cl
      new.cl[names(tmp.cl)] = tmp.cl + n.cl
      n.cl = max(new.cl)
      markers=union(markers, result$markers)
    }
    cl=new.cl
    cl = setNames(as.integer(as.factor(cl)), names(cl))
    return(list(cl=cl, markers=markers))
  }

merge_cl_by_co <- function(cl, co.ratio, diff.th=0.25, verbose=0){
  cell.cl.co.ratio = get_cell.cl.co.ratio(cl, co.ratio)
  cl.co.ratio <- do.call("rbind",tapply(names(cl),cl, function(x)colMeans(cell.cl.co.ratio[x,,drop=F])))
  co.within= diag(cl.co.ratio)
  co.df <- as.data.frame(as.table(cl.co.ratio),stringsAsFactors=FALSE)
  co.df = co.df[co.df[,1]<co.df[,2]& co.df[,3]>0.1,]
  co.df$within1 = co.within[co.df[,1]] 
  co.df$within2 = co.within[co.df[,2]]
  co.df$diff = pmax(co.df$within1, co.df$within2) - co.df[,3]
  co.df = co.df[co.df$diff  < diff.th,]
  co.df = co.df[order(co.df[,1],decreasing=T),]
  if(verbose > 0){
    print(co.df)
  }
  for(i in 1:nrow(co.df)){
    cl[cl==co.df[i,2]]=co.df[i,1]
  }
  cl = setNames(as.integer(as.character(cl)), names(cl))
  return(cl)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param cl 
##' @param co.ratio 
##' @param cl.mat 
##' @return 
##' @author Zizhen Yao
get_cell.cl.co.ratio <- function(cl, co.ratio=NULL, cl.mat=NULL)
  {
    if(!is.null(co.ratio)){
      cell.cl.co.ratio=get_cl_means(co.ratio, cl)
      return(cell.cl.co.ratio)    
    }
    if(!is.null(cl.mat)){
      tmp = cl.mat %*% get_cl_sums(t(cl.mat), cl)
      tmp = tmp / Matrix::rowSums(cl.mat)
      cl.size = table(cl)
      cell.cl.co.ratio=as.matrix(t(t(tmp)/as.vector(cl.size[colnames(tmp)])))
      return(cell.cl.co.ratio)
    }
    stop("Either co.ratio or cl.mat should not be NULL")
  }

    
get_cl_co_stats <- function(cl, co.ratio=NULL, cl.mat=NULL)
  {
    cell.cl.co.ratio= get_cell.cl.co.ratio(cl, co.ratio=co.ratio, cl.mat=cl.mat)
    cl.co.ratio <- get_cl_means(t(cell.cl.co.ratio), cl)

    
    #cell.co.stats <- do.call("rbind",sapply(1:ncol(cell.cl.co.ratio),function(i){
    cell.co.stats <- sapply(1:ncol(cell.cl.co.ratio),function(i){
      select.cells=names(cl)[cl==colnames(cell.cl.co.ratio)[i]]
      cohesion = setNames(cell.cl.co.ratio[select.cells, i, drop=F], select.cells)
      best.between = rowMaxs(cell.cl.co.ratio[select.cells, -i, drop=F])
      confusion = best.between / cohesion
      separability = cohesion  - best.between
      df=data.frame(cohesion, separability, confusion)
      colnames(df) = c("cohesion", "separability", "confusion")
      df
    },simplify=F)

    cell.co.stats = do.call("rbind", cell.co.stats)
    cl.co.stats = as.data.frame(do.call("rbind",tapply(1:nrow(cell.co.stats), cl[row.names(cell.co.stats)], function(x){
      sapply(cell.co.stats[x,], median)
    })))
    
    return(list(cell.cl.co.ratio=cell.cl.co.ratio,
                cl.co.ratio=cl.co.ratio,
                cell.co.stats = cell.co.stats, 
                cl.co.stats = cl.co.stats))
  }


init_cut <- function(co.ratio, select.cells, cl.list, min.cells=4, th = 0.3,method="ward",verbose=FALSE)
{
  avg.cl.num = mean(sapply(cl.list, function(cl){
    sum(table(cl[select.cells]) >= min.cells)
  }))
  tmp.dat = co.ratio[select.cells, select.cells]
  hc=  hclust(as.dist(1-as.matrix(crossprod(tmp.dat))), method="ward.D")
  tmp.cl = cutree(hc, ceiling(avg.cl.num)+2)
  tmp.cl=refine_cl(tmp.cl, co.ratio=co.ratio, min.cells=min.cells, niter=1, confusion.th=1)$cl
  if(length(unique(tmp.cl))==1){
    return(NULL)    
  }
  tmp.cl=merge_cl_by_co(tmp.cl, tmp.dat, diff.th=th)
  if(length(unique(tmp.cl))==1){
    return(NULL)    
  }
  return(cl=tmp.cl)
}


init_cut.old <- function(co.ratio, select.cells, cl.list=NULL, min.cells=4, th = 0.3,method="ward",verbose=FALSE)
  {
    tmp.dat = as.matrix(co.ratio[select.cells, select.cells])
    hc = hclust(as.dist(1-tmp.dat),method=method)
    tmp.cl=cut_co_matrix(tmp.dat, hc$order,w=min.cells-1, th = th)
    ord = colnames(tmp.dat)[hc$order]
    sep = which(tmp.cl[-1]!=tmp.cl[-length(tmp.cl)])
    if(verbose){
      pdf("co.pdf")
      heatmap.2(tmp.dat[ord,ord], Colv=NULL,Rowv=NULL,trace="none",col=blue.red(100),colsep=sep,sepcolor="black")
      dev.off()
    }    
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    tmp.cl=merge_cl_by_co(tmp.cl, tmp.dat, diff.th=th)
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    ###reassign each cells to the clusters, if it prefer to co-cluster with another cluster much better.
    tmp.cl=refine_cl(tmp.cl, co.ratio=co.ratio, min.cells=min.cells, niter=1, tol.th=1)$cl
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    tmp <- do.call("cbind",tapply(names(tmp.cl), tmp.cl, function(x){
      rowMeans(tmp.dat[,x,drop=F])
    }))
    cl.hc = hclust(dist(t(tmp)),method="average")
    cl = setNames(factor(as.character(tmp.cl), levels=cl.hc$labels[cl.hc$order]),names(tmp.cl))
    cl = setNames(as.integer(cl),names(cl))
    return(cl)
  }
 
cut_co_matrix <- function(co.ratio, ord, w=3,th=0.25)
  {
    tmp=co.ratio[ord,ord]
    l = ncol(tmp)
    ##Compute the at every position the moving window
    co.w=sapply(1:(l-w),function(x)rowMeans(tmp[,x:(x+w)]))
    diff = colSums(abs(co.w[,1:(l-w*2-1),drop=F] - co.w[,(w+2):(l-w),drop=F]))
    diff = diff/ colSums(co.w[,1:(l-w*2-1),drop=F] + co.w[,(w+2):(l-w),drop=F])
    diff = c(rep(0,w),diff,rep(0,w))

    peak=slice(Rle(diff), th)
    sep=which.max(peak)
    diff.bin=cut(diff,breaks= c(seq(-0.05,1,length.out=100),2))
    diff.col= jet.colors(100)[diff.bin]
    cl = setNames(as.integer(cut(1:ncol(co.ratio), c(0,sep,ncol(co.ratio)+1))), colnames(co.ratio)[ord])
    return(cl)
  }
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param cl 
##' @param co.ratio 
##' @param cl.mat 
##' @param co.stats 
##' @param confusion.th 
##' @param min.cells 
##' @param niter 
##' @param tol.th 
##' @param verbose 
##' @return 
##' @author Zizhen Yao
refine_cl <- function(cl, co.ratio=NULL, cl.mat=NULL, co.stats=NULL,confusion.th=0.4,min.cells=4, niter=10, tol.th=0.02, verbose=0)
  {
    if(is.null(co.stats)){
      co.stats = get_cl_co_stats(cl, co.ratio=co.ratio, cl.mat=cl.mat)
    }
    correct = 0
    iter.num = 0
    while(iter.num < niter){
      cell.cl.co.ratio = co.stats$cell.cl.co.ratio
      cl.confusion = setNames(co.stats$cl.co.stats$confusion, row.names(co.stats$cl.co.stats))
      tmp.dat = cell.cl.co.ratio[names(cl),as.character(sort(unique(cl)))]
      pred.cl <- setNames(colnames(tmp.dat)[apply(tmp.dat, 1, which.max)], row.names(tmp.dat))      
      if(sum(cl==pred.cl) <= correct){
        break
      }
      correct = sum(cl==pred.cl)
      correct.frac= correct/length(cl)
      if(verbose){
        print(correct.frac)
      }
      if(1 - correct.frac < tol.th){
        break
      }
      tmp.cells = names(pred.cl)[pred.cl!=cl[names(pred.cl)]]
      cl[tmp.cells]=pred.cl[tmp.cells]
      cl.size = table(cl)
      ###Remove small clusters with high average confusion score, assign cells to other cluster      
      cl.small = names(cl.size)[cl.size <= 2 * min.cells]
      rm.cl = cl.small[cl.confusion[cl.small] > confusion.th | cl.size[cl.small] < min.cells]
      tmp.cells = names(cl)[cl %in% rm.cl]
      tmp.dat = cell.cl.co.ratio[tmp.cells,setdiff(unique(cl),rm.cl),drop=F]
      pred.cl = setNames(colnames(tmp.dat)[apply(tmp.dat, 1, which.max)], row.names(tmp.dat))
      cl[tmp.cells] = pred.cl[tmp.cells]
      ###Recompute co stats
      co.stats = get_cl_co_stats(cl, co.ratio=co.ratio, cl.mat=cl.mat)
      iter.num = iter.num + 1
      if(length(unique(cl))==1){
        break
      }
    }
    return(list(cl=cl, co.stats=co.stats))
  }


plot_co_matrix <- function(co.ratio, cl, max.cl.size=100, col=NULL)
{
  select.cells = names(cl)
  select.cells = sample_cells(cl, max.cl.size)
  tom  = crossprod(co.ratio[select.cells, select.cells])
  row.names(tom)=colnames(tom)=select.cells
  ###
  all.hc = hclust(as.dist(1-tom),method="average")
  ord1 = all.hc$labels[all.hc$order]
  ord1 = ord1[ord1%in% select.cells]
  ord = ord1[order(cl[ord1])]
  sep = cl[ord]
  sep=which(sep[-1]!=sep[-length(sep)])
  if(is.null(col)){
    heatmap.3(as.matrix(co.ratio[ord,ord]), col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
  }
  else{
    heatmap.3(as.matrix(co.ratio[ord,ord]), col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black", ColSideColors=col[,ord])
  }
}



