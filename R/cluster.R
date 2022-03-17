#' Compute jaccard distances between columns of a matrix
#' 
#' @param m A matrix or sparse matrix
#' 
#' @return a sparse matrix of Jaccard distances.
#' 
jaccard <-  function(m) {
  library(Matrix)
  ## common values:                                                                                                                                              
  A <-  Matrix::tcrossprod(m)
  A <- as(A, "dgTMatrix")
  ## counts for each row                                                                                                                                         
  b <- Matrix::rowSums(m)
  ## Jacard formula: #common / (#i + #j - #common)                                                                                                               
  x = A@x / (b[A@i+1] + b[A@j+1] - A@x)
  A@x = x
  return(A)
}  

knn_jaccard <- function(knn)
{
  knn.df = data.frame(i = rep(1:nrow(knn), ncol(knn)), j=as.vector(knn))
  knn.mat = sparseMatrix(i = knn.df[[1]], j=knn.df[[2]], x=1, giveCsparse=FALSE)
  jaccard(knn.mat)
}


pass_louvain <- function(mod.sc, adj.mat)
{
  library(Matrix)
  p <- mean(Matrix::colSums(adj.mat > 0) - 1) / nrow(adj.mat)
  n <- ncol(adj.mat)
  rand.mod1 <- 0.97 * sqrt((1 - p) / (p * n))
  rand.mod2 <- (1 - 2 / sqrt(n)) * (2 / (p * n)) ^ (2 / 3)
  rand.mod.max <- max(rand.mod1, rand.mod2, na.rm = TRUE)
  #cat("Modularity:",mod.sc, "threshold:",rand.mod.max, "\n")
  return(mod.sc > rand.mod.max)
}

#' Perform Jaccard/Louvain clustering based on RANN
#'
#' @param dat A matrix of features (rows) x samples (columns)
#' @param k K nearest neighbors to use
#' 
#' @return A list object with the cluster factor object and (cl) and Jaccard/Louvain results (result)
#' 
jaccard_louvain.RANN <- function(dat, 
                                 k = 10)
{
  library(igraph)
  library(matrixStats)
  library(RANN)  
  
  knn.matrix = RANN::nn2(dat, k = k)[[1]]
  jaccard.adj  <- knn_jaccard(knn.matrix)  
  jaccard.gr <- igraph::graph.adjacency(jaccard.adj, 
                                        mode = "undirected", 
                                        weighted = TRUE)
  louvain.result <- igraph::cluster_louvain(jaccard.gr)
  mod.sc <- igraph::modularity(louvain.result)
  
  if(pass_louvain(mod.sc, jaccard.adj)) {
    
    cl <- setNames(louvain.result$membership, row.names(dat))
    
    return(list(cl = cl, result = louvain.result))
    
  } else{
    return(NULL)
  }
}


#' Perform Jaccard/Louvain clustering using Rphenograph
#' 
#' @param dat A matrix of features (rows) x samples (columns)
#' @param k K-nearest neighbors to use for clustering
#' 
#' @return A list object with the cluster factor object and (cl) and Rphenograph results (result)
#' 
jaccard_louvain <- function(dat, k = 10)
{
  suppressPackageStartupMessages(library(Rphenograph))
  
  rpheno <- Rphenograph(dat, k = k)
  
  cl <- setNames(rpheno[[2]]$membership, row.names(dat)[as.integer(rpheno[[2]]$names)])
  
  return(list(cl = cl, result = rpheno))
}

filter_RD <- function(rd.dat, rm.eigen, rm.th, verbose=FALSE)
{
  library(matrixStats)
  rm.cor=cor(rd.dat, rm.eigen[row.names(rd.dat),])
  rm.cor[is.na(rm.cor)]=0
  rm.score = rowMaxs(abs(rm.cor))
  if(verbose){
    print("rm score")
    print(tail(sort(rm.score)))
  }
  select = colSums(t(abs(rm.cor)) >= rm.th) ==0
  if(sum(!select)>0 & verbose){
    print("Remove dimension:")
    print(rm.score[!select])
  }
  if(sum(select)==0){
    return(NULL)
  }
  rd.dat = rd.dat[,select,drop=F]
}



#' Perform Jaccard/Leiden clustering 
#' 
#'
#' @param dat A matrix of samples (rows) x features (columns)
#' @param k K nearest neighbors to use
#' 
#' @return A list object with the cluster factor object and (cl) and Jaccard/Leiden results (result)
#' @param num_iter If n_iter < 0 the optimiser continues iterating until it encounters an iteration that did not improve the partition.
#' 
#' 
#load("data/Leiden_test/L4_5_6_IT.rd.dat.rda")
#

jaccard_leiden <- function(dat, k = 10, weight = NULL,
                           num_iter = 5,
                           resolution_parameter = 0.01,                           
                           random_seed = NULL,
                           verbose = FALSE, ...) {
  
  library(igraph)
  library(matrixStats)
  library(RANN)  
  library(Matrix)

  if(verbose){
    cat("Compute Jaccard distance\n")
  }
  knn.matrix = RANN::nn2(dat, k = k)[[1]]
  jaccard.adj  <- knn_jaccard(knn.matrix)
  jaccard.gr <- igraph::graph.adjacency(jaccard.adj, 
                                        mode = "undirected", 
                                        weighted = TRUE)
  
  partition_type <- 'CPMVertexPartition'
  
  #write line to determine optimal resolutoin↔
  if(verbose){
    cat("Leiden clustering\n")
  }
  cluster_result <- leidenbase::leiden_find_partition(jaccard.gr,
                                                      partition_type = partition_type,
                                                      num_iter=num_iter,
                                                      resolution_parameter=resolution_parameter,
                                                      verbose = FALSE )
  
  cl <- setNames(cluster_result$membership, row.names(dat))
  return(list(cl = cl, result = cluster_result))
}





#' One round of clustering in the iteractive clustering pipeline 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1).
#' @param select.cells The cells to be clustered. Default: columns of norm.dat
#' @param counts Raw gene counts. Default NULL, inferred from norm.dat.
#' @param method Clustering method. It can be "louvain", "hclust" and "kmeans". Default "louvain"
#' @param vg.padj.th High variance gene adjusted pvalue cut off. Default 0.5.
#' @param dim.method Dimension reduction techniques. Current options include "pca" and "WGCNA". Default "pca"
#' @param max.dim The number of top dimensions retained. Default 20. Since clustering is performed iteratively, not all relevant dimensions need to be captured in one iterations. 
#' @param rm.eigen The reduced dimensions that need to be masked and removed. Default NULL.  
#' @param rm.th The cutoff for correlation between reduced dimensions and rm.eigen. Reduced dimensions with correlatin with any rm.eigen vectors are not used for clustering. Default 0.7
#' @param de.param The differential gene expression threshold. See de_param() function for details. 
#' @param type Can either be "undirectional" or "directional". If "undirectional", the differential gene threshold de.param is applied to combined up-regulated and down-regulated genes, if "directional", then the differential gene threshold is applied to both up-regulated and down-regulated genes. 
#' @param maxGenes Only used when dim.method=="WGCNA". The maximum number of genes to calculate gene modules. 
#' @param sampleSize The number of sampled cells to compute reduced dimensions.
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.    
#' @param prefix Used to keep track of intermediate results in "verbose" mode. Default NULL.
#' @param verbose Default FALSE
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#'         
onestep_clust <- function(norm.dat, 
                          select.cells = colnames(norm.dat),
                          counts = NULL, 
                          method = c("louvain","ward.D", "leiden","kmeans"), 
                          vg.padj.th = 0.5, 
                          dim.method = c("pca","WGCNA"), 
                          max.dim = 20, 
                          rm.eigen = NULL, 
                          rm.th = 0.7, 
                          de.param = de_param(),
                          merge.type = c("undirectional", "directional"),
                          genes.allowed = row.names(norm.dat),
                          maxGenes = 3000,
                          sampleSize = 4000,
                          max.cl.size = 300,
                          k.nn = 15,
                          prefix = NULL, 
                          verbose = FALSE, 
                          regress.x=NULL,
                          max.cl=NULL)

{
  library(matrixStats)

  method <- match.arg(method)
  dim.method <- match.arg(dim.method)
  merge.type <- match.arg(merge.type)
  
  if(!is.null(regress.x)){
    print("regression")
    tmp =  lm_normalize(as.matrix(norm.dat[,select.cells]), regress.x[select.cells], R_2.th=0.1)
    norm.dat = tmp[[1]]
  }
  if(length(select.cells)>sampleSize){
    sampled.cells = sample(select.cells, pmin(length(select.cells),sampleSize))
  }
  else{
    sampled.cells = select.cells
  }
  
  ###Find high variance genes
  tmp = get_cl_present(norm.dat, setNames(rep(1, length(select.cells)),select.cells), de.param$low.th)
  select.genes = row.names(norm.dat)[which(tmp * length(select.cells) >= de.param$min.cells)]
  
  ###Find high variance genes.
  if(is.null(counts)){
    if(is.matrix(norm.dat)){
      counts = 2^(norm.dat[select.genes,sampled.cells])-1
    }
    else{
      counts = norm.dat[select.genes,sampled.cells]
      counts@x = 2^(counts@x) - 1
    }
  }
  plot_file=NULL
  if(verbose & !is.null(prefix)){
    plot_file=paste0(prefix,".vg.pdf")
  }
  vg = find_vg(counts,plot_file=plot_file)
  rm(counts)
  if(dim.method=="auto"){
    if(length(select.cells)> 1000){
      dim.method="pca"
    }
    else{
      dim.method="WGCNA"
    }
  }
  if(dim.method=="WGCNA"){
    ###Ignore vg.padj.th for WGCNA, choose top "maxGgenes" for analysis
    select.genes = as.character(vg[which(vg$loess.padj < 1),"gene"])
    select.genes = head(select.genes[order(vg[select.genes, "loess.padj"],-vg[select.genes, "z"])],maxGenes)
    rd.dat = rd_WGCNA(norm.dat, select.genes=select.genes, select.cells=select.cells, sampled.cells=sampled.cells, de.param=de.param, max.mod=max.dim, max.cl.size=max.cl.size)$rd.dat
  }
  else{
    ###If most genes are differentially expressed, then use absolute dispersion value
    select.genes = as.character(vg[which(vg$loess.padj < vg.padj.th | vg$dispersion >3),"gene"])
    select.genes = intersect(select.genes, genes.allowed)
    select.genes = head(select.genes[order(vg[select.genes, "loess.padj"],-vg[select.genes, "z"])],maxGenes)
    if(verbose){
      cat("Num high variance genes:",length(select.genes),"\n")
    }
    if(length(select.genes)< de.param$min.genes){
      return(NULL)
    }
    rd.result = rd_PCA(norm.dat,select.genes, select.cells, sampled.cells=sampled.cells, max.pca = max.dim)    
    rd.dat = rd.result$rd.dat
    if(verbose){
      cat("PCA dimensions:",ncol(rd.dat),"\n")
    }
  }
  if(is.null(rd.dat)||ncol(rd.dat)==0){
    return(NULL)
  }
  if(!is.null(rm.eigen)){
    rd.dat <- filter_RD(rd.dat, rm.eigen, rm.th, verbose=verbose)
  }
 
  if(is.null(rd.dat)||ncol(rd.dat)==0){
    return(NULL)
  }
  if(verbose){
    print(method)
  }
  if(is.null(max.cl)){
    max.cl = pmin(ncol(rd.dat)*2 + 1, round(length(select.cells)/de.param$min.cells))
  }
  if(method=="louvain"){
    k = pmin(k.nn, round(nrow(rd.dat)/2))
    tmp = jaccard_louvain(rd.dat, k)
    if(is.null(tmp)){
      return(NULL)
    }
    cl = tmp$cl
    if(length(unique(cl))>max.cl){
      tmp.means = get_cl_means(rd.dat, cl)
      tmp.hc = hclust(dist(t(tmp.means)), method="average")
      tmp.cl= cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
      cl = setNames(tmp.cl[as.character(cl)], names(cl))
    }
  }
  
  else if(method=="leiden"){
    k = pmin(k.nn, round(nrow(rd.dat)/2))
    tmp = jaccard_leiden(rd.dat, k)
    if(is.null(tmp)){
      return(NULL)
    }
    cl = tmp$cl
    if(length(unique(cl))>max.cl){
      tmp.means =do.call("cbind",tapply(names(cl),cl, function(x){
        colMeans(rd.dat[x,,drop=F])
      },simplify=F))
      tmp.hc = hclust(dist(t(tmp.means)), method="average")
      tmp.cl= cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
      cl = setNames(tmp.cl[as.character(cl)], names(cl))
    }
  }  
  
  else if(method=="ward.D"){
    hc = hclust(dist(rd.dat),method="ward.D")
    #print("Cluster cells")
    cl = cutree(hc, max.cl)
  }
  else if(method=="kmeans"){
    cl = kmeans(rd.dat, max.cl)$cluster
  }
  else{
    stop(paste("Unknown clustering method", method))
  }
  #print(table(cl))
  merge.result=merge_cl(norm.dat, cl=cl, rd.dat=rd.dat, merge.type=merge.type, de.param=de.param, max.cl.size=max.cl.size,verbose=verbose)
  
  gc()
  if(is.null(merge.result))return(NULL)
  sc = merge.result$sc
  #print(sc)
  cl = merge.result$cl
  if(length(unique(cl))<=1){
    return(NULL)
  }
  de.genes = merge.result$de.genes
  markers= merge.result$markers
  cl.dat = get_cl_means(norm.dat[markers,], cl[sample_cells(cl, max.cl.size)])
  cl.hc = hclust(dist(t(cl.dat)),method="average")
  cl = setNames(factor(as.character(cl), levels= colnames(cl.dat)[cl.hc$order]), names(cl))
  levels(cl) = 1:length(levels(cl))
  result=list(cl=cl, markers=markers)
  if(verbose){
    cat("Expand",prefix, "\n")
    cl.size=table(cl)
    print(cl.size)
  }
  return(result)
}


#' Iterative clustering algorithm for single cell RNAseq dataset
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param select.cells The cells to be clustered
#' @param prefix The character string to indicate current iteration.
#' @param split.size The minimal cluster size for further splitting
#' @param result The current clustering result as basis for further splitting.
#' @param method Clustering method. It can be "auto", "louvain", "hclust"
#' @param ... Other parameters passed to method `onestep_clust()`
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#'         
iter_clust <- function(norm.dat, 
                       select.cells = colnames(norm.dat),
                       prefix = NULL, 
                       split.size = 10, 
                       result = NULL,
                       method = "auto",
                       overwrite=TRUE,
                       verbose=verbose,
                       ...)
{
  if(!is.null(prefix)) {
    cat(prefix, length(select.cells),"\n")
  }

  if(method == "auto"){
    if(length(select.cells) > 3000){
      select.method="louvain"
    }
    else{
      select.method="ward.D"
    }
  }
  else{
    select.method=method
  }
  if(is.null(result)){
    outfile=paste0(prefix, ".rda")
    if(file.exists(outfile) & !overwrite){
      load(outfile)       
    }
    else{
      if(length(select.cells) <= 3000){
        if(!is.matrix(norm.dat)){
          norm.dat = as.matrix(norm.dat[,select.cells])
        }
      }
      result=onestep_clust(norm.dat, select.cells=select.cells, prefix=prefix,method=select.method,verbose=verbose,...)
      if(verbose){
        save(result, file=outfile)
      }
      gc()
    }
    if(is.null(result)){
      return(NULL)
    }
  }
  
  select.cells= intersect(select.cells, names(result$cl))
  cl = result$cl[select.cells]
  gene.mod = result$gene.mod
  markers=result$markers
  cl = setNames(as.integer(cl),names(cl))
  new.cl =cl
  cl.size = table(cl)
  to.split = names(cl.size)[cl.size >=split.size]
  if(length(to.split)>0){
    n.cl = 1
    for(x in sort(unique(cl))){
      tmp.cells = names(cl)[cl==x]
      if(!x %in% to.split){
        new.cl[tmp.cells]=n.cl
      }
      else{
        tmp.prefix = paste(prefix, x, sep=".")
        tmp.result=iter_clust(norm.dat=norm.dat, select.cells=tmp.cells, prefix=tmp.prefix,split.size=split.size,method= method,overwrite=overwrite,verbose=verbose,...)
        gc()
        if(is.null(tmp.result)){
          new.cl[tmp.cells]=n.cl
        }
        else{
          tmp.cl = tmp.result$cl
          if(length(unique(tmp.cl)>1)){
            new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
            markers=union(markers, tmp.result$markers)
          }
        }
      }
      n.cl = max(new.cl)+1
    }
    cl = new.cl
  }
  result=list(cl=cl, markers=markers)
  return(result)
}


#' Reorder cluster based on hiearchical clustering of clusters based on average cluster values for the input data matrix 
#'
#' @param cl A vector of cluster membership with cells as names, and cluster id as values. 
#' @param dat The data matrix with cells as columns. 
#'
#' @return Reorder cluster membership vector. The cluster id start from 1 to the number of clusters.
#' @export 
#'
reorder_cl <- function(cl, dat)
{
  cl.means = get_cl_means(dat,cl)
  cl.hc = hclust(as.dist(1 - cor(cl.means)),method = "average")
  cl = setNames(factor(as.character(cl), levels = cl.hc$labels[cl.hc$order]),names(cl))
  cl = setNames(as.integer(cl),names(cl))
}



combine_finer_split <- function(cl, finer.cl)
{
  if(is.factor(cl)){
    cl = setNames(as.integer(cl), names(cl))
  }
  if(is.factor(finer.cl)){
    finer.cl = setNames(as.integer(finer.cl), names(finer.cl))
  }
  max.cl = max(cl)
  finer.cl = finer.cl + max.cl
  cl[names(finer.cl)] = finer.cl
  return(cl)
}


iter_clust_merge <- function(norm.dat, select.cells, merge.type="undirectional", de.param = de_param(), max.cl.size = 300,...)
{
  result <- iter_clust(norm.dat=norm.dat, select.cells=select.cells, de.param = de.param, merge.type=merge.type, ...)
  result=merge_cl(norm.dat, cl=result$cl, rd.dat.t = norm.dat[result$markers,], merge.type=merge.type, de.param=de.param, max.cl.size=max.cl.size)
  return(result)
}
