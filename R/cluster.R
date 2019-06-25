# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# jaccard()
#
# knn_jaccard()
#
# pass_louvain()
#
# jaccard_louvain.RANN()
#   knn_jaccard() cluster.R
#   pass_louvain() cluster.R
#    
# jaccard_louvain()
#
# onestep_clust()
#   lm_normalize() normalize.R
#   find_vg() vg.R
#   rd_WGCNA() reduceDimension_WGCNA.R
#   rd_PCA() reduceDimension_PCA.R
#   jaccard_louvain() cluster.R
#   merge_cl() merge_cl.R
#   get_cl_means() util.R
#   display_cl() plot.R
#
# iter_clust()
#   onestep_clust() cluster.R
#   iter_clust() cluster.R
#
# reorder_cl()
#   get_cl_means() util.R
#   
# combine_finer_split()
#


#' Compute jaccard distances between columns of a matrix
#' 
#' @param m A matrix or sparse matrix
#' 
#' @return a sparse matrix of class dgTMatrix containing Jaccard distances between every pair of samples.
#' @export
#' 
jaccard <-  function(m) {
  
  # common values:                                                                                                                                              
  A <- Matrix::tcrossprod(m)
  A <- as(A, "dgTMatrix")
  
  # counts for each row                                                                                                                                         
  b <- Matrix::rowSums(m)
  
  # Jacard formula: #common / (#i + #j - #common)
  x <- A@x / (b[A@i+1] + b[A@j+1] - A@x)
  
  A@x <- x
  
  return(A)
}

#' Compute a jaccard distance matrix based on a matrix of K-Nearest Neigbors
#' 
#' @param knn a matrix or sparse matrix of k-nearest neighbors
#' 
#' @return a sparse matrix of Jaccard distances produced by the \code{jaccard()} function.
#' @export
#' 
knn_jaccard <- function(knn) {
  
  knn.df <- data.frame(i = rep(1:nrow(knn), 
                               ncol(knn)), 
                       j = as.vector(knn))
  
  knn.mat <- Matrix::sparseMatrix(i = knn.df[[1]], 
                                  j = knn.df[[2]], 
                                  x = 1)
  
  jaccard(knn.mat)
}


#' Test louvain modularity scores based on an adjacency matrix
#'
#' @param mod.sc a vector of modularity scores
#' @param adj.mat an adjacency matrix
#'
#' @return a logical vector
#' @export
#'
pass_louvain <- function(mod.sc, 
                         adj.mat) {
  
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
                                 k = 10) {
  
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
jaccard_louvain <- function(dat, 
                            k = 10) {
  
  rpheno <- Rphenograph::Rphenograph(dat, k = k)
  
  cl <- setNames(rpheno[[2]]$membership, 
                 row.names(dat)[as.integer(rpheno[[2]]$names)])
  
  return(list(cl = cl, 
              result = rpheno))
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
#' @param max.genes Only used when dim.method=="WGCNA". The maximum number of genes to calculate gene modules. 
#' @param sample.size The number of sampled cells to compute reduced dimensions.
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.    
#' @param prefix Used to keep track of intermediate results in "verbose" mode. Default NULL.
#' @param verbose Default FALSE
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#'         
onestep_clust <- function(norm.dat, 
                          select.cells = NULL, 
                          counts = NULL, 
                          method = "louvain", 
                          vg.padj.th = 0.5, 
                          dim.method = "pca", 
                          max.dim = 20, 
                          rm.eigen = NULL, 
                          rm.th = 0.7, 
                          de.param = de_param(),
                          merge.type = "undirectional", 
                          max.genes = 3000,
                          sample.size = 4000,
                          max.cl.size = 300,
                          k.nn = 15,
                          prefix = NULL, 
                          verbose = FALSE, 
                          regress.x = NULL) {
  
  if(is.null(select.cells)) {
    select.cells <- colnames(norm.dat)
  }
  
  if(is.null(counts)){
    counts <- 2 ^ (norm.dat[select.genes, sampled.cells]) - 1
  }
  
  method <- match.arg(method,
                      choices = c("louvain","ward.D","kmeans"))
  
  dim.method <- match.arg(dim.method,
                          choices = c("pca","WGCNA","auto"))
  
  merge.type <- match.arg(merge.type,
                          choices = c("undirectional","directional"))
  
  if(!is.null(regress.x)){
    if(verbose) cat("running regression\n")
    
    reg_results <- lm_normalize(as.matrix(norm.dat[,select.cells]), 
                                regress.x[select.cells], 
                                R_2.th = 0.1)
    norm.dat <- reg_results[[1]]
  }
  
  if(length(select.cells) > sample.size) {
    sampled.cells <- sample(select.cells, 
                            pmin(length(select.cells), sample.size))
  } else{
    sampled.cells <- select.cells
  }
  print(head(select.cells))
  #Filter genes based on de_params: low.th and min.cells
  if(is.matrix(norm.dat)){
    cells_gt_low.th <- rowSums(norm.dat[,select.cells] > de.param$low.th)
  } else {
    cells_gt_low.th <- Matrix::rowSums(norm.dat[,select.cells] > de.param$low.th)
  }
  genes_gt_min.cells <- cells_gt_low.th >= de.param$min.cells
  select.genes <- row.names(norm.dat)[which(genes_gt_min.cells)]
  
  #Find high variance genes.

  plot_file <- NULL
  
  if(verbose & !is.null(prefix)){
    plot_file <- paste0(prefix,".vg.pdf")
  }

  vg <- find_vg(as.matrix(counts[select.genes, sampled.cells]),
                plot_file = plot_file)
  
  if(dim.method == "auto") {
    if(length(select.cells) > 1000) {
      dim.method <- "pca"
    } else {
      dim.method <- "WGCNA"
    }
  }
  
  if(dim.method == "WGCNA"){
    ###Ignore vg.padj.th for WGCNA, choose top "maxGgenes" for analysis
    select.genes <- as.character(vg[which(vg$loess.padj < 1), "gene"])
    select.genes <- head(select.genes[order(vg[select.genes, "padj"], -vg[select.genes, "z"])],
                         max.genes)
    rd.dat <- rd_WGCNA(norm.dat = norm.dat, 
                       select.genes = select.genes, 
                       select.cells = select.cells, 
                       sampled.cells = sampled.cells, 
                       de.param = de.param, 
                       max.mod = max.dim, 
                       max.cl.size = max.cl.size)$rd.dat
  } else if(dim.method == "pca") {
    ###If most genes are differentially expressed, then use absolute dispersion value
    select.genes <- as.character(vg[which(vg$loess.padj < vg.padj.th | vg$dispersion > 3), "gene"])
    select.genes <- head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])], max.genes)
    
    if(verbose){
      cat("Num high variance genes:", length(select.genes), "\n")
    }
    if(length(select.genes) < de.param$min.genes){
      if(verbose) warning("Not enough differentially expressed genes. Returning NULL.")
      return(NULL)
    }
    rd.dat <- rd_PCA(norm.dat = norm.dat,
                     select.genes, 
                     select.cells, 
                     sampled.cells = sampled.cells, 
                     max.pca = max.dim)$rd.dat
  }
  
  if(is.null(rd.dat) || ncol(rd.dat) == 0) {
    if(verbose) warning("Dimensionality reduction returned NULL.")
    return(NULL)
  }
  
  if(!is.null(rm.eigen)){
    rm.cor <- cor(rd.dat, rm.eigen[row.names(rd.dat),])
    rm.cor[is.na(rm.cor)] <- 0
    rm.score <- rowMaxs(abs(rm.cor))
    select <- rm.score < rm.th
    if(sum(!select) > 0 & verbose) {
      if(verbose) print("Removing dimension based on rm.eigen correlation:")
      print(rm.score[!select])
    }
    if(sum(select) == 0) {
      if(verbose) warning("All reduced dimensions match rm.eigen. Returning NULL.")
      return(NULL)
    }
    rd.dat <- rd.dat[, select, drop = FALSE]
  }
  
  if(verbose) {
    print(method)
  }
  
  max.cl <- ncol(rd.dat) * 2 + 1
  
  if(method == "louvain") {
    k <- pmin(k.nn, round(nrow(rd.dat) / 2))
    jl_result <- jaccard_louvain(rd.dat, k)
    
    if(is.null(jl_result)){
      if(verbose) warning("jaccard_louvain() returned NULL.")
      return(NULL)
    }
    
    cl <- jl_result$cl
    
    if(length(unique(cl)) > max.cl) {
      # If overclustered, merge based on average hierarchical clustering of means
      # to max.cl using cutree
      jl_means <- do.call("cbind", 
                          tapply(names(cl), 
                                 cl, 
                                 function(x) {
                                   colMeans(rd.dat[x, , drop = FALSE])
                                 },
                                 simplify = FALSE))
      jl_hc <- hclust(dist(t(jl_means)), 
                      method = "average")
      jl_cut <- cutree(cl_hc, 
                       pmin(max.cl, length(unique(cl))))
      cl <- setNames(jl_cut[as.character(cl)], names(cl))
    }
  } else if(method == "ward.D") {
    ward_hc <- hclust(dist(rd.dat),
                 method = "ward.D")
    cl <- cutree(ward_hc, max.cl)
  } else if(method == "kmeans") {
    cl <- kmeans(rd.dat, max.cl)$cluster
  }
  
  rd.dat.t <- t(rd.dat)
  merge.result <- merge_cl(norm.dat, 
                           cl = cl, 
                           rd.dat.t = rd.dat.t, 
                           merge.type = merge.type, 
                           de.param = de.param, 
                           max.cl.size = max.cl.size)
  gc()
  if(is.null(merge.result)) {
    if(verbose) warning("merge_cl() returned NULL.")
    return(NULL)
  }
  
  sc <- merge.result$sc

  cl <- merge.result$cl

  # Assemble results
  if(length(unique(cl)) > 1) {
    if(verbose){
      cat("Expand", prefix, "\n")
      cl.size <- table(cl)
      print(cl.size)
      save(cl, 
           file = paste0(prefix, ".cl.rda"))
    }
    de.genes <- merge.result$de.genes
    markers <- merge.result$markers
    cl.dat <- get_cl_means(norm.dat[markers,], 
                           cl[sample_cells(cl, max.cl.size)])
    cl.hc <- hclust(dist(t(cl.dat)),method="average")
    cl <- setNames(factor(as.character(cl), levels = colnames(cl.dat)[cl.hc$order]), names(cl))
    if(verbose & !is.null(prefix)){
      display_cl(cl, 
                 norm.dat, 
                 prefix = prefix, 
                 markers = markers, 
                 max.cl.size = max.cl.size)
    }
    levels(cl) <- 1:length(levels(cl))
    result <- list(cl = cl, 
                   markers = markers)
    return(result)
  }
  return(NULL)
}


#' Iterative clustering algorithm for single cell RNAseq dataset
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param select.cells The cells to be clustered
#' @param prefix The character string to indicate current iteration.
#' @param split.size The minimal cluster size for further splitting
#' @param result The current clustering result as basis for further splitting.
#' @param method Clustering method. It can be "auto", "louvain", "hclust"
#' @param ... Other parameters passed to `onestep_clust()`
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that separate clusters     
#'         
iter_clust <- function(norm.dat, 
                       select.cells = NULL,
                       prefix = NULL, 
                       split.size = 10, 
                       result = NULL,
                       method = "auto",
                       ...) {
  
  method <- match.arg(method,
                      choices = c("auto","louvain","hclust"))
  
  if(is.null(select.cells)) {
    select.cells <- colnames(norm.dat)
  }
  
  if(!is.null(prefix)) { 
    print(prefix)
  }
  
  if(method == "auto"){
    if(length(select.cells) > 3000) {
      method <- "louvain"
    } else {
      method <- "ward.D"
    }
  }
  
  if(length(select.cells) <= 3000) {
    if(!is.matrix(norm.dat)){
      norm.dat <- as.matrix(norm.dat[,select.cells])
    }
  }
  
  if(is.null(result)){        
    result <- onestep_clust(norm.dat, 
                            select.cells = select.cells, 
                            prefix = prefix,
                            method = select.method,
                            ...)
    gc()
  }
  
  if(is.null(result)) {
    return(NULL)
  }
  
  select.cells <- intersect(select.cells, names(result$cl))
  cl <- result$cl[select.cells]
  gene.mod <- result$gene.mod
  markers <- result$markers
  
  cl <- setNames(as.integer(cl),names(cl))
  new.cl <- cl
  cl.size <- table(cl)
  to.split <- names(cl.size)[cl.size >= split.size]
  
  if(length(to.split) > 0) {
    n.cl <- 1
    for(x in sort(unique(cl))) {
      tmp.cells <- names(cl)[cl == x]
      
      if(!x %in% to.split) {
        new.cl[tmp.cells] <- n.cl
      } else {
        tmp.prefix <- paste(prefix, x, sep = ".")
        tmp.result <- iter_clust(norm.dat = norm.dat, 
                                 select.cells = tmp.cells, 
                                 prefix = tmp.prefix,
                                 split.size = split.size,
                                 method = method,
                                 ...)
        gc()
        
        if(is.null(tmp.result)) {
          new.cl[tmp.cells] <- n.cl
        } else {
          tmp.cl <- tmp.result$cl
          if(length(unique(tmp.cl) > 1)) {
            new.cl[names(tmp.cl)] <- n.cl + as.integer(tmp.cl)
            markers <- union(markers, tmp.result$markers)
          }
        }
      }
      n.cl <- max(new.cl) + 1
    }
    
    cl <- new.cl
  }
  result <- list(cl = cl, 
                 markers = markers)
  return(result)
}


#' Reorder cluster based on hiearchical clustering of mean values from the input data matrix 
#'
#' @param cl A vector of cluster membership with cells as names, and cluster id as values. 
#' @param dat The data matrix with cells as columns. 
#'
#' @return Reordered cluster membership vector. The cluster ids will start from 1 to the number of clusters.
#' @export 
#'
reorder_cl <- function(cl, 
                       dat) {
  
  cl.means <- get_cl_means(dat,
                           cl)
  
  cl.hc <- hclust(as.dist(1 - cor(cl.means)),
                  method = "average")
  
  cl <- setNames(factor(as.character(cl), 
                        levels = cl.hc$labels[cl.hc$order]),
                 names(cl))
  
  cl <- setNames(as.integer(cl),
                 names(cl))
  
  return(cl)
}



combine_finer_split <- function(cl,
                                finer.cl) {
  
  if(is.factor(cl)){
    cl <- setNames(as.integer(cl),
                   names(cl))
  }
  
  if(is.factor(finer.cl)){
    finer.cl <- setNames(as.integer(finer.cl),
                         names(finer.cl))
  }
  
  max.cl <- max(cl)
  finer.cl <- finer.cl + max.cl
  
  cl[names(finer.cl)] <- finer.cl
  
  return(cl)
}
