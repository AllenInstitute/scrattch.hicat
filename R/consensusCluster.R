run_iter_clust_merge <- function(norm.dat, iter, output_dir="subsample_result", all.cells=colnames(norm.dat), prefix="all", de.param=de_param(), merge.type="undirectional", sample.frac = 0.8, init.result=NULL, mc.cores=1, override=FALSE)
  {
    run <- function(i,...){
      prefix = paste("iter",i,sep=".")
      print(prefix)
      library(Matrix)
      outfile= file.path(output_dir, paste0("result.",i,".rda"))
      if(file.exists(outfile)& !override){
        return(NULL)
      }
      select.cells=sample(all.cells, round(length(all.cells)*sample.frac))
      save(select.cells, file=file.path(output_dir, paste0("cells.",i,".rda")))
      
      result <- scrattch.hicat::iter_clust_merge(norm.dat=norm.dat, select.cells=select.cells,prefix=prefix, de.param = de.param, merge.type=merge.type, result=init.result, ...)
      save(result, file=outfile)
    }
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }    
    if (mc.cores==1){
      sapply(iter, function(i){run(i,...)})
    }
    else{     
      require(doMC)
      require(foreach)
      registerDoMC(cores=mc.cores)
      foreach::foreach(i=iter,.packages=c("scrattch.hicat","Matrix"), .combine='c') %dopar% { run(i) }
    }
  }
    


#' Collect co-clustering matrix from results files
#'
#' @param result.files A directory containing results files
#' @param all.cells The cells to read from results files
#'
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

get_co_ratio <- function(cl.mat, cells, n.times)
{
  co.ratio = Matrix::crossprod(cl.mat[,cells])
  co.ratio@x = co.ratio@x/n.times
  return(co.ratio)
}

#' Iterative consensus clustering
#'
#' @param co.ratio cell cel√l co-clustering matrix
#' @param cl.list  The list of subsampled clustering results. 
#' @param norm.dat The log2 transformed normalzied expression matrix 
#' @param select.cells Cells to be clustered
#' @param de.param Differentiall expressed genes criteria for merging clusters
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately. 
#' @param diff.th The difference of co-clustering probablities for splitting a cluster. 
#' @param prefix Default NULL. 
#' @param method Clustering methods. Default "auto"
#' @param verbose Default FALSE
#' @param result Pre-computed clustering results used for further splitting. Default NULL.
#'
#' @return A list with cluster membership, and top pairwise marker genes. 
#'
iter_consensus_clust <- function(cl.list, 
                                 co.ratio=NULL,  
                                 cl.mat=NULL, 
                                 norm.dat, 
                                 select.cells=names(cl.list[[1]]), 
                                 diff.th=0.25, 
                                 prefix=NULL, 
                                 method=c("auto", "louvain","ward.D"), 
                                 verbose=FALSE, 
                                 de.param = de.param, 
                                 max.cl.size = 300, 
                                 result=NULL, 
                                 split.size = de.param$min.cells*2, 
                                 merge.type=c("undirectional", "directional"))
{
  method=method[1]
  require(igraph)
  if(verbose){
    print(prefix)
  }
  if(!is.null(result)){
    markers=result$markers
    cl = setNames(as.integer(as.character(result$cl)),names(result$cl))
    cell.cl.co.ratio= get_cell.cl.co.ratio(cl, co.ratio= co.ratio, cl.mat=cl.mat[,names(cl)])
  }
  else{
    markers=NULL
    if(length(select.cells)  < split.size){
      return(NULL)
    }
    co.ratio.sampled = FALSE
    if(is.null(co.ratio)){
      cl.size = table(cl.list[[1]][select.cells])
      graph.size= sum(cl.size^2)
      if(graph.size > 10^8){
        co.ratio.sampled=TRUE
        tmp.cl.list = lapply(cl.list, function(cl)cl[select.cells])
        sampled.cells = sample_cl_list(tmp.cl.list, max.cl.size=max.cl.size)
        
        cl.size = table(cl.list[[1]][sampled.cells])
        graph.size= sum(cl.size^2)
        if(graph.size > 10^8){
          sampled.cells = sample_cells(cl.list[[1]][sampled.cells], max.cl.size)
        }
      }
      else{
        sampled.cells=select.cells
      }
      co.ratio = get_co_ratio(cl.mat, sampled.cells, n.times=length(cl.list))
    }
    if(method=="auto"){
      if (length(select.cells)> 3000){
        select.method = "louvain"
      }
      else{
        select.method="ward.D"
      }
    }
    else{
      select.method = method
    }
    if(select.method=="ward.D"){
      if(!is.matrix(co.ratio)){
        tmp.co.ratio = as.matrix(co.ratio[select.cells, select.cells])
      }
      else{
        tmp.co.ratio = co.ratio
      }
      tmp.cl = init_cut(tmp.co.ratio, select.cells, cl.list, min.cells= de.param$min.cells, th = diff.th,method=select.method)
      rm(tmp.co.ratio)
      if(is.null(tmp.cl)){
        return(NULL)
      }
    }
    else{###louvain
      if(co.ratio.sampled){
        adj.mat = co.ratio
      }
      else{
        adj.mat = co.ratio[select.cells, select.cells]
      }
      gr = graph.adjacency(adj.mat, mode="undirected",weighted=TRUE)
      comm= cluster_louvain(gr)
      rm(gr)
      
      if(pass_louvain(modularity(comm), adj.mat)){
        tmp.cl = setNames(comm$membership,colnames(adj.mat))
        if(length(unique(tmp.cl))==1){
          return(NULL)
        }
        
      }
      else{
        return(NULL)
      }
      rm(adj.mat)
      gc()
    }
    if(verbose){
      print(table(tmp.cl))
    }
    tmp.cl=merge_cl_by_co(tmp.cl, co.ratio=co.ratio, cl.mat=cl.mat[,names(tmp.cl)],diff.th)
    cell.cl.co.ratio= get_cell.cl.co.ratio(tmp.cl, co.ratio= co.ratio, cl.mat=cl.mat[,names(tmp.cl)])    
    tmp= merge_cl(norm.dat=norm.dat, cl=tmp.cl, rd.dat=cell.cl.co.ratio, verbose=verbose,  de.param = de.param, return.markers=FALSE, max.cl.size= max.cl.size, merge.type=merge.type)
    markers=tmp$markers
    if(is.null(tmp) | !is.list(tmp)) return(NULL)
    if (length(unique(tmp$cl))==1) return(NULL)
    tmp.cl= tmp$cl
    tmp.cl = setNames(as.integer(tmp.cl),names(tmp.cl))
    if(length(unique(tmp.cl))==1) {
      return(NULL)
    }
    if(co.ratio.sampled){
      co.ratio = NULL
    }
    gc()
    cell.cl.co.ratio= get_cell.cl.co.ratio(tmp.cl, co.ratio= co.ratio, cl.mat=cl.mat[,select.cells])[select.cells,]
    cl = setNames(as.integer(colnames(cell.cl.co.ratio)[apply(cell.cl.co.ratio, 1, which.max)]), row.names(cell.cl.co.ratio))    
    if(verbose){
      cat("Total:", length(cl), "\n")
      print(table(cl))
    }
  }
  n.cl=max(cl)
  new.cl=cl
  for(i in sort(unique(cl))){
    tmp.prefix= paste0(prefix, ".", i)
    tmp.cells=names(cl)[cl==i]
    uncertain.cells=sum(cell.cl.co.ratio[tmp.cells, as.character(i)] < 1 - diff.th)
    if(uncertain.cells < de.param$min.cells){
      next
    }
    result= iter_consensus_clust(cl.list=cl.list, co.ratio=co.ratio, cl.mat = cl.mat,  norm.dat=norm.dat, select.cells=tmp.cells, prefix=tmp.prefix,  diff.th =diff.th, method=method, de.param = de.param, verbose=verbose, max.cl.size=max.cl.size, merge.type=merge.type)
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



collect_subsample_cl_matrix <- function(norm.dat,result.files,all.cells,max.cl.size=NULL,mc.cores=1)
{
  select.cells=c()
  run <- function(f){
    print(f)
    tmp=load(f)
    if(is.null(result)){
      return(NULL)
    }
    cl= result$cl    
    test.cells = setdiff(all.cells, names(cl))
    if(length(test.cells) > 0){
      if(is.null(result$test.cl)){
        markers=unique(result$markers)
        map.df = map_by_cor(norm.dat[markers,names(cl)],cl, norm.dat[markers,test.cells],method="mean")$pred.df
        test.cl = setNames(map.df$pred.cl, row.names(map.df))
      }
      else{
        test.cl = result$test.cl
      }
      all.cl = c(setNames(as.character(cl),names(cl)), setNames(as.character(test.cl), names(test.cl)))
    }
    else{
      all.cl = cl
    }
    return(all.cl[all.cells])
  }
  if (mc.cores==1){
    cl.list=sapply(result.files, function(f){run(f)},simplify=F)
    cl.list = cl.list[!sapply(cl.list,is.null)]
  }
  else{
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    cl.list= foreach(i=1:niter, .combine='c') %dopar% run(f)
  }
  if(!is.null(max.cl.size)){
    select.cells= sample_cl_list(cl.list, max.cl.size=max.cl.size)
  }
  else{
    select.cells= all.cells
  }
  cl.mat = compile_cl_mat(cl.list, select.cells)
  return(list(cl.list=cl.list, cl.mat = cl.mat))
}

compile_cl_mat <- function(cl.list, select.cells)
{
  cl.mat = do.call("cbind", sapply(names(cl.list), function(x){
    print(x)
    cl = cl.list[[x]]
    get_cl_mat(cl[select.cells])
  },simplify=F))
  cl.mat= Matrix::t(cl.mat)
}



#' Refine clusters
#' 
#' @param cl Cluster membership cluster
#' @param co.ratio cell-cell co-clustering matrix. 
#' @param cl.mat cell-cluster matrix collected from bootstrapping iterations. Either co.ratio or cl.mat should not be NULL. 
#' @param confusion.th Clusters with average confusion score greater than this threshold will be removed. Cells in this cluster will be re-distributed to other most likely clusters. 
#' @param min.cells Clusters with fewer than this many cells will be removed. Cells in this cluster will be re-distributed to other most likely clusters. 
#' @param niter maxmimal mumber of refinement iterations. 
#' @param tol.th If improvement is smaller than this threshold, terminate refinement step.  
#' @param verbose If true, print out step-by-step improvement. 
#' @return 
#' @author Zizhen Yao
refine_cl <- function(cl, 
                      co.ratio=NULL, 
                      cl.mat=NULL, 
                      confusion.th=0.6, 
                      min.cells=4, 
                      niter=50, 
                      tol.th=0.02, 
                      verbose=0)
{
  ###If cl is factor, turn in to integer vector first. 
  cl = setNames(as.integer(as.character(cl)), names(cl))
  while(TRUE){
    correct = 0
    iter.num = 0
    while(iter.num < niter){
      cell.cl.co.ratio <- get_cell.cl.co.ratio(cl, co.ratio=co.ratio, cl.mat=cl.mat)
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
      iter.num = iter.num + 1
    }
    cl.size = table(cl)
    co.stats = get_cl_co_stats(cl, co.ratio=co.ratio, cl.mat=cl.mat)
    cl.confusion = setNames(co.stats$cl.co.stats$confusion, row.names(co.stats$cl.co.stats))
    ###Remove small clusters with high average confusion score, assign cells to other cluster     
    cl.small = names(cl.size)[cl.size <  min.cells]
    rm.cl = union(names(cl.confusion)[cl.confusion > confusion.th], cl.small)
    if(length(rm.cl)==0){
      break
    }
    if(length(rm.cl) == length(cl.size)){
      cl[names(cl)] = min(cl) 
      return(list(cl=cl, co.stats=co.stats))
    }
    tmp.cells = names(cl)[cl %in% rm.cl]
    tmp.dat = cell.cl.co.ratio[tmp.cells,as.character(setdiff(unique(cl),rm.cl)),drop=F]
    pred.cl = setNames(colnames(tmp.dat)[apply(tmp.dat, 1, which.max)], row.names(tmp.dat))
    cl[tmp.cells] = pred.cl[tmp.cells]
    if(length(unique(cl))==1){
      break
    }
  }
  return(list(cl=cl, co.stats=co.stats))
}

merge_cl_by_co <- function(cl, co.ratio=NULL, cl.mat=NULL, diff.th=0.25, verbose=0){
  cell.cl.co.ratio = get_cell.cl.co.ratio(cl, co.ratio=co.ratio, cl.mat=cl.mat)
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


#' Get cell co-clustering ratios
#'  
#' @param cl Vector of cluster assignments
#' @param co.ratio coclustering ratio results
#' @param cl.mat Cluster membership matrix for all cells and all clusters from all bootstrapping iterations. 
get_cell.cl.co.ratio <- function(cl, co.ratio=NULL, cl.mat=NULL)
{
  require(Matrix)
  if(!is.null(co.ratio)){
    cell.cl.co.ratio=get_cl_means(co.ratio, cl)
  }
  else if(!is.null(cl.mat)){
    tmp1= get_cl_sums(cl.mat[,names(cl)], cl)
    tmp = Matrix::crossprod(tmp1, cl.mat)
    cl.size = table(cl)
    n.times= Matrix::colSums(cl.mat)
    tmp = tmp/ as.vector(cl.size[row.names(tmp)])
    cell.cl.co.ratio = as.matrix(Matrix::t(tmp)/ n.times)
  }
  else{
    stop("Either co.ratio or cl.mat should not be NULL")
  }
  return(cell.cl.co.ratio)
}

get_cl_co_stats <- function (cl, co.ratio = NULL, cl.mat = NULL) 
{
  require(matrixStats)
  cell.cl.co.ratio = get_cell.cl.co.ratio(cl, co.ratio = co.ratio, cl.mat = cl.mat)
  cl.co.ratio <- get_cl_means(t(cell.cl.co.ratio), cl)
  cell.co.stats <- sapply(1:ncol(cell.cl.co.ratio), function(i) {
    select.cells = names(cl)[cl == colnames(cell.cl.co.ratio)[i]]
    cohesion = setNames(cell.cl.co.ratio[select.cells, i, 
                                         drop = F], select.cells)
    best.between = rowMaxs(cell.cl.co.ratio[select.cells, 
                                            -i, drop = F])
    confusion = best.between/cohesion
    separability = cohesion - best.between
    df = data.frame(cohesion, separability, confusion)
    colnames(df) = c("cohesion", "separability", "confusion")
    df
  }, simplify = F)
  cell.co.stats = do.call("rbind", cell.co.stats)
  cl.co.stats = as.data.frame(do.call("rbind", tapply(1:nrow(cell.co.stats), 
                                                      cl[row.names(cell.co.stats)], function(x) {
                                                        sapply(cell.co.stats[x, ], median)
                                                      })))
  return(list(cell.cl.co.ratio = cell.cl.co.ratio, cl.co.ratio = cl.co.ratio, 
              cell.co.stats = cell.co.stats, cl.co.stats = cl.co.stats))
}

init_cut <- function(co.ratio, select.cells, cl.list, min.cells=4, th = 0.3,method="ward.D",verbose=FALSE)
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

plot_co_matrix <- function(co.ratio, cl, max.cl.size=100, col=NULL)
{
  blue.red <- colorRampPalette(c("blue", "white", "red"))
  select.cells = names(cl)
  select.cells = sample_cells(cl, max.cl.size)
  tom  = Matrix::crossprod(co.ratio[select.cells, select.cells])
  row.names(tom)=colnames(tom)=select.cells
  ###
  all.hc = hclust(as.dist(1-tom),method="average")
  ord1 = all.hc$labels[all.hc$order]
  ord1 = ord1[ord1%in% select.cells]
  ord = ord1[order(cl[ord1])]
  sep = cl[ord]
  sep=which(sep[-1]!=sep[-length(sep)])
  if(is.null(col)){
    heatmap.3(as.matrix(co.ratio[ord,ord]), col = blue.red(150)[50:150], trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black", labRow="")
  }
  else{
    heatmap.3(as.matrix(co.ratio[ord,ord]), col = blue.red(150)[50:150], trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black", ColSideColors=col[,ord],labRow="")
  }
}


plot_cell_cl_co_matrix <- function(co.ratio, cl, max.cl.size=100, col=NULL)
{
  blue.red <- colorRampPalette(c("blue", "white", "red"))
  select.cells = sample_cells(cl, max.cl.size)
  co.stats = get_cl_co_stats(cl, co.ratio)
  mat = co.stats$cell.cl.co.ratio
  
  tom  = Matrix::tcrossprod(mat[select.cells,])
  row.names(tom)=colnames(tom)=select.cells
  ###
  all.hc = hclust(as.dist(1-tom),method="average")
  ord1 = all.hc$labels[all.hc$order]
  ord = ord1[order(cl[ord1])]
  sep = cl[ord]
  sep=which(sep[-1]!=sep[-length(sep)])
  if(is.null(col)){
    heatmap.3(mat[ord,], col = blue.red(150)[50:150], trace="none", Rowv=NULL, Colv=NULL,rowsep=sep,sepcolor="black", dendrogram="none",labRow="")
  }
  else{
    heatmap.3(mat[ord,], col = blue.red(150)[50:150], trace="none", Rowv=NULL, Colv=NULL,rowsep=sep,sepcolor="black", ColSideColors=col[,ord],dendogram="none",labRow="")
  }
}


#' Wrapper function to repeatively run clustering on subsampled cells and infer consensus clusters
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1).
#' @param select.cells The cells to be clustered. Default: columns of norm.dat
#' @param niter The number of iteractions to run. Default 100.
#' @param sample.frac The fraction of of cells sampled per run. Default: 0.8. 
#' @param output_dir The output directory to store clutering results for each iteraction. 
#' @param mc.cores The number of cores to be used for parallel processing. 
#' @param de.param The differential gene expression threshold. See de_param() function for details. 
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately. 
#' @param override binary variable determine if the clustering results already stored in output_dir should be overriden.  
#' @param init.result The pre-set high level clusters. If set, the function will only find finer splits of the current clusters.  
#' @param ... Other parameters passed to iter_clust
#'
#' @export
#' 
run_consensus_clust <- function(norm.dat, 
                                select.cells=colnames(norm.dat), 
                                niter=100, 
                                sample.frac=0.8, 
                                co.result=NULL, 
                                output_dir="subsample_result",
                                mc.cores=1, 
                                de.param=de_param(), 
                                merge.type=c("undirectional","directional"), 
                                override=FALSE, 
                                init.result=NULL, 
                                cut.method="auto",
                                confusion.th=0.6,
                                ...)
{
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  all.cells=select.cells
  if(!is.null(init.result)){
    all.cells= intersect(all.cells, names(init.result$cl))
  }
  if(is.null(co.result)){
    tmp=run_iter_clust_merge(norm.dat=norm.dat, iter=1:niter, output_dir=output_dir, all.cells=all.cells, prefix=prefix, de.param=de.param, merge.type=merge.type, sample.frac=0.8, init.result=init.result, mc.cores=mc.cores)    
    result.files=file.path(output_dir, dir(output_dir, "result.*.rda"))
    co.result <- collect_subsample_cl_matrix(norm.dat,result.files,all.cells)
  }
  cl.size = table(co.result$cl.list[[1]])
  graph.size = sum(cl.size^2)
  if(graph.size < 10^9){
    consensus.result = iter_consensus_clust(cl.list=co.result$cl.list, cl.mat = co.result$cl.mat, norm.dat=norm.dat, select.cells=all.cells, de.param = de.param, merge.type=merge.type, method=cut.method, result= init.result)
    refine.result = refine_cl(consensus.result$cl, cl.mat = co.result$cl.mat, tol.th=0.01, confusion.th=confusion.th, min.cells= de.param$min.cells)
    markers = consensus.result$markers
  }
  else{
    result <- scrattch.hicat::iter_clust(norm.dat=norm.dat, select.cells=all.cells, de.param = de.param, merge.type= merge.type, result= init.result,...)
    
    cl=merge_cl_by_co(result$cl, co.ratio=NULL, cl.mat=co.result$cl.mat, diff.th=0.25)
    refine.result = refine_cl(cl, cl.mat = co.result$cl.mat, tol.th=0.01, confusion.th=0.6, min.cells=de.param$min.cells)
    markers=result$markers      
  }
  cl = refine.result$cl
  merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[markers,], de.param = de.param, merge.type=merge.type, return.markers=FALSE)
  return(list(co.result=co.result, cl.result=merge.result))
}
