# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# get_pair_matrix_coor()
#
# get_pair_matrix()
#   get_pair_matrix_coor() util.R
#
# set_pair_matrix()
#   get_pair_matrix_coor() util.R
#
# get_pairs()
#
# convert_pair_matrix()
#   get_pairs() util.R
#   set_pair_matrix() util.R
#
# convert_pair_matrix_str()
#
# get_cl_mat()
#
# get_cl_sums()
#   get_cl_mat() util.R
#
# get_cl_means()
#   get_cl_sums() util.R
#
# get_cl_medians()
#
# get_cl_prop()
#   get_cl_mat() util.R
#
# sparse_cor()
#
# calc_tau()
#
# sample_cells()
#
# cpm()
#
# logCPM()
#   cpm() util.R
# 
# pair_cor()
#

#' Convert matrix row/column positions to vector position
#' 
#' @param m a matrix object
#' @param rows row positions, either as indices or matches to row names
#' @param cols row positions, either as indices or matches to row names
#' 
#' @return a vector containing the numeric position at [rows,cols]
#' 
get_pair_matrix_coor <- function(m, 
                                 rows, 
                                 cols) {
  
  if(!is.numeric(rows)){
    rows <- match(rows, rownames(m))
  }
  
  if(!is.numeric(cols)){
    cols <- match(cols, colnames(m))
  }
  
  coor <- (cols - 1L) * nrow(m) + rows
  
  return(coor)
}

#' Subset a matrix as a vector using row and column positions
#' 
#' @param m a matrix object
#' @param rows row positions, either as indices or matches to row names
#' @param cols row positions, either as indices or matches to row names
#' 
#' @return a vector with values extracted from m at [rows/cols]
#' 
get_pair_matrix <- function(m, 
                            rows, 
                            cols) {
  v <- as.vector(m)
  coor <- get_pair_matrix_coor(m, rows, cols)
  return(v[coor])
}


#' Update a matrix with values from a 1d vector using row and column positions
#' 
#' @param m a matrix object
#' @param rows row positions, either as indices or matches to row names
#' @param cols row positions, either as indices or matches to row names
#' @param vals values to insert at [rows/cols]
#' 
#' @return a matrix with updated values
#' 
set_pair_matrix <- function(m, 
                            rows, 
                            cols, 
                            vals) {
  
  coor <- get_pair_matrix_coor(m, rows, cols)
  m[coor] <- vals
  
  return(m)
}

#' Convert underscore_separated pair names to a data.frame
#' 
#' @param pairs.str a character vector with underscore_separated values
#' 
#' @return a data.frame with columns "P1" and "P2" containing the separated pair names
#' 
#' @export
#' 
get_pairs <- function(pairs.str) {
  
  pairs_split <- strsplit(pairs.str, "_")
  pairs_mat <- do.call("rbind", pairs_split)
  pairs_df <- as.data.frame(pairs_mat,
                            stringsAsFactors = FALSE)
  
  row.names(pairs_df) <- pairs.str
  colnames(pairs_df) <- c("P1", "P2")
  
  return(pairs_df)
}


#' Convert paired cluster comparison values to a matrix
#' 
#' @param pair.num a named numeric vector of values. Names correspond to compared elements separted by "_", e.g. "c1_c2", "c23_c59"
#' @param l labels for columns. Default is NULL, which will compute them from names(pair.num)
#' @param directed If FALSE (default), the first value in each pair will specify used as columns, with the second as rows. 
#' If TRUE, first values will be rows, and second will be columns.
#' 
#' @return a matrix containing values from pair.num, and named for each element separated by "_" in names(pair.num)
#'  
#' @examples
#' 
#' pair_values <- seq(1,27,3)
#' names(pair_values) <- paste(rep(letters[1:3], each = 3), rep(letters[1:3], 3), sep = "_")
#' pair_values
#' 
#' pair_matrix <- convert_pair_matrix(pair_values, directed = FALSE)
#' pair_matrix
#' 
#' pair_matrix <- convert_pair_matrix(pair_values, directed = TRUE)
#' pair_matrix
#' 
convert_pair_matrix <- function(pair.num, 
                                l = NULL,
                                directed = FALSE) {
  
  pairs <- get_pairs(names(pair.num))
  
  if(is.null(l)){
    l <- sort(unique(c(pairs[,1], pairs[,2])))
  }
  
  n.cl <- length(l)
  
  pair.num.mat <- matrix(0, nrow = n.cl, ncol = n.cl)
  rownames(pair.num.mat) <- l
  colnames(pair.num.mat) <- l
  
  pair.num.mat <- set_pair_matrix(pair.num.mat, pairs[,1], pairs[,2], pair.num)
  
  if(!directed) {
    pair.num.mat <- set_pair_matrix(pair.num.mat, pairs[,2], pairs[,1], pair.num)
  }
  
  return(pair.num.mat)
}

#' Convert paired cluster comparison values to a matrix
#' 
#' @param pair.num a vector of values with names of compared elements separted by "_", e.g. "c1_c2", "c23_c59"
#' @param l labels for columns. Default is NULL, which will compute them from names(pair.num)
#' @param directed If FALSE (default), the first value in each split will be used as columns, with the second as rows. 
#' If TRUE, first values will be rows, and second will be columns.
#' 
#' @return a matrix containing values from pair.num, and named for each element separated by "_" in names(pair.num)
#'  
#' @examples
#' 
#' pair_values <- seq(1,27,3)
#' names(pair_values) <- paste(rep(letters[1:3], each = 3), rep(letters[1:3], 3), sep = "_")
#' pair_values
#' 
#' pair_matrix <- convert_pair_matrix(pair_values, directed = FALSE)
#' pair_matrix
#' 
#' pair_matrix <- convert_pair_matrix(pair_values, directed = TRUE)
#' pair_matrix
#' 
convert_pair_matrix_str <- function(pair.str, 
                                    l = NULL,
                                    directed = FALSE) {
  
  pairs <- do.call("rbind", strsplit(names(pair.str),"_"))
  
  if(is.null(l)){
    l <- sort(unique(as.vector(pairs)))
  }
  
  n.cl <- length(l)
  
  pair.str.mat <- matrix("", nrow = n.cl, ncol = n.cl)
  rownames(pair.str.mat) <- l
  colnames(pair.str.mat) <- l
  
  for(i in 1:nrow(pairs)){
    pair.str.mat[pairs[i,1], pairs[i,2]] <- pair.str[i]
    if(!directed) {
      pair.str.mat[pairs[i,2], pairs[i,1]] <- pair.str[i]
    }
  }
  return(pair.str.mat)
}


get_weighted_means <- function(mat, w)
  {
    rowSums(Matrix::tcrossprod(mat,diag(w)))/sum(w)
  }

#' Generate a sparse matrix one-hot representation of clusters x samples
#' 
#' @param cl a cluster factor object
#' 
#' @return a sparse, one-hot matrix indicating which cluster(columns) each sample (rows) belongs to.
#' @export
#' 
get_cl_mat <- function(cl, all.cells=NULL) {
  
  if(!is.factor(cl)){
    cl <- as.factor(cl)
  }
  cl <- droplevels(cl)
  if(is.null(all.cells)){
    all.cells = names(cl)
    i = 1:length(cl)
  }
  else{
    i = match(names(cl),all.cells)
  }
  j = as.integer(cl)
  cl.mat <- Matrix::sparseMatrix(i = i,
                                 j = j,
                                 x = 1, dims=c(length(all.cells), max(j)))  
  rownames(cl.mat) <- all.cells
  colnames(cl.mat) <- levels(cl)
  return(cl.mat)
}

#' Compute cluster sums for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with sums for each cluster
#' @export
#' 
get_cl_sums_R<- function(mat, 
                        cl)
{  
  if(all(names(cl) %in% colnames(mat))){
    cl.mat <- get_cl_mat(cl, all.cells=colnames(mat))
    cl.sums <- Matrix::tcrossprod(mat, Matrix::t(cl.mat))
  }
  else{
    cl.mat <- get_cl_mat(cl, all.cells=row.names(mat))
    cl.sums <- Matrix::crossprod(mat, cl.mat)
  }
  cl.sums <- as.matrix(cl.sums)
  return(cl.sums)
}

###Get row means without subsetting the matrix
get_row_sums <- function(mat, select.row=1:nrow(mat), select.col=1:ncol(mat))
  {
    if(!is.integer(select.col)){
      select.col = match(select.col, colnames(mat))
    }
    cl = setNames(rep(1, length(select.col)), colnames(mat)[select.col])
    sums = get_cl_sums(mat, cl)
    sums[select.row,]
  }

get_row_means <- function(mat, select.row=1:nrow(mat), select.col=1:ncol(mat))
  {
    if(!is.integer(select.col)){
      select.col = match(select.col, colnames(mat))
    }
    cl = setNames(rep(1, length(select.col)), colnames(mat)[select.col])
    means = get_cl_means(mat, cl)
    means[select.row,]
  }


get_row_vars <- function(mat, select.row=1:nrow(mat), select.col=1:ncol(mat), means=NULL)
  {
    if(!is.integer(select.col)){
      select.col = match(select.col, colnames(mat))
    }
    cl = setNames(rep(1, length(select.col)), colnames(mat)[select.col])
    if(is.null(means)){
      means = get_cl_means(mat, cl)
    }
    sqr_means = get_cl_sqr_means(mat, cl)    
    vars = sqr_means - means^2
    vars[select.row,]
  }



#' Compute cluster means for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with means for each cluster
#' @export
#' 
get_cl_means_R<- function(mat, 
                         cl) {
  
  cl.sums <- get_cl_sums_R(mat, cl)
  
  cl.size <- table(cl)
  
  cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
  
  return(cl.means)
}

get_cl_means <- function(mat,cl)
{
  if(!is.factor(cl)){
    cl = setNames(factor(cl),names(cl))
  }
  result = get_cl_stats(mat, cl, stats="means")
  result[,levels(cl),drop=F]
}

get_cl_sums <- function(mat,cl)
{
  if(!is.factor(cl)){
    cl = setNames(factor(cl),names(cl))
  }
  result = get_cl_stats(mat, cl, stats="sums")
  result[,levels(cl),drop=F]
}

get_cl_sqr_means <- function(mat,cl)
{
  if(!is.factor(cl)){
    cl = setNames(factor(cl),names(cl))
  }
  result = get_cl_stats(mat, cl, stats="sqr_means")
  result[,levels(cl),drop=F]
}


get_cl_present<- function(mat, cl, low.th)
{
  if(!is.factor(cl)){
    cl = setNames(factor(cl),names(cl))
  }
  result = get_cl_stats(mat, cl, stats="present",lowth=low.th)
  result[,levels(cl),drop=F]
}

get_cl_vars <- function(mat, cl, cl.means=NULL, cl.sqr.means = NULL)
{  
  if(is.null(cl.means)){
    cl.means = get_cl_means(mat,cl)
  }
  if(is.null(cl.sqr.means)){
    cl.sqr.means = get_cl_sqr_means(mat,cl)
  }
  cl.vars = cl.sqr.means - cl.means^2
  cl.size = as.vector(table(cl)[colnames(cl.vars)])
  cl.vars  = t(t(cl.vars)  *  cl.size/(cl.size - 1))
  return(cl.vars)
}

#' Compute cluster medians for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with medians for each cluster
#' @export
#' 
get_cl_medians_R <- function(mat, cl)
{
  library(Matrix)
  library(matrixStats)
  
  cl.med <- do.call("cbind",
                    tapply(names(cl), 
                           cl, 
                           function(x){
                             matrixStats::rowMedians(as.matrix(mat[,x]))
                           }
                    )
  )
  
  rownames(cl.med) <- rownames(mat)
  
  return(cl.med)
}

get_cl_medians <- function(mat, cl)
  {
    if(!is.factor(cl)){
      cl = setNames(factor(cl),names(cl))
    }
    rcpp_get_cl_medians(mat, cl)
  }

#' Compute cluster proportions for each row in a matrix
#'
#' @param mat a gene (rows) x samples(columns) sparse matrix
#' @param cl A cluster factor object
#' @param threshold The minimum expression value used to binarize the results
#'
#' @return a matrix of genes (rows) x cluster(columns) with proportions for each cluster
#' @export
#' 
get_cl_prop <- function(mat, 
                        cl, 
                        threshold = 0) {
  
  cl.mat <- get_cl_mat(cl)
  
  cl.prop <- Matrix::tcrossprod(mat[,rownames(cl.mat)] > threshold, Matrix::t(cl.mat))
  
  cl.prop <- as.matrix(cl.prop) / Matrix::colSums(cl.mat)[col(cl.prop)]
  
  return(cl.prop)
}


#' Compute correlation scores for columns of a sparse matrix
#' 
#' @param m a sparse matrix, preferrably dgCMatrix
#' 
#' @return a matrix of correlation values between each column of m
#'  
sparse_cor <- function(m) {
  #library(Matrix)
  
  n_rows <- nrow(m)
  n_cols <- ncol(m)
  
  ii <- unique(m@i) + 1 # rows with a non-zero element
  
  Ex <- Matrix::colMeans(m)
  
  nozero <- as.vector(m[ii,]) - rep(Ex, each = length(ii))        # colmeans
  
  covmat <- (crossprod(matrix(nozero, ncol = n_cols)) +
               crossprod(t(Ex)) * (n_rows - length(ii))
  ) / (n_rows - 1)
  
  sdvec <- sqrt(diag(covmat))
  
  cormat <- covmat / crossprod(t(sdvec))
  
  return(cormat)
}

#' Calculate Tau scores for each gene
#' 
#' @param m Matrix of expression values
#' @param byRow if TRUE, treats genes as row values and samples as columns.
#' 
#' @return a vector of Tau scores
#' 
calc_tau <- function(m, 
                     byRow = TRUE) {
  if(!byRow){
    m <- t(m)
  }
  
  row_maxes <- apply(m, 1, max)
  
  m <- m / row_maxes
  tau <- Matrix::rowSums(1 - m) / (ncol(m) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}

#' Downsample cells from each cluster
#' 
#' @param cl A cluster factor object
#' @param sample.size A maximum number of cells to take from each cluster, or a named numeric object with the number of cells to be sampled and names matching cluster levels
#' @param weights A named numeric vector with weights for each cell, passed to the prob parameter of the sample() function. Default is NULL.
#' @param seed A seed value for random sampling. If NULL (default), will be randomized.
#' 
#' @return A cluster factor object containing sampled cells
#' 
sample_cells <- function(cl, 
                         sample.size, 
                         weights = NULL,
                         seed = NULL) {
  
  n_cl <- unique(cl)
  
  if(length(sample.size) == 1) {
    
    sample.size <- setNames(rep(sample.size, length(n_cl)), n_cl)
    
  }
  
  cl.cells <- split(names(cl), cl)
  
  sampled.cells <- sapply(names(cl.cells), 
                          function(x) {
                            cells <- cl.cells[[x]]
                            
                            if(sample.size[[x]] >= length(cells)){
                              return(cells)
                            }
                            
                            to.sample <- pmin(sample.size[[x]], length(cells))
                            
                            if(!is.null(weights)){
                              set.seed(seed)
                              
                              sampled <- sample(cells, to.sample, prob = weights[cells])
                              
                            } else{
                              set.seed(seed)
                              
                              sampled <- sample(cells, to.sample)
                            }
                            
                            sampled
                            
                          }, simplify = FALSE)
  
  sampled.cells <- unlist(sampled.cells)
  
  return(sampled.cells)
}

#' Convert a matrix of raw counts to a matrix of Counts per Million values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to samples, and rows to genes.
#' 
#' @param counts a matrix, dgCMatrix, or dgTMatrix of count values.
#' 
#' @return a matrix, dgCMatrix, or dgTMatrix of CPM values (matching input)
#' 
#' @export
#' 
cpm <- function(counts, sf=NULL, denom=1e6) {
  if(is.null(sf)){
    sf <- Matrix::colSums(counts)
  }
  sf = sf/denom
  if(is.matrix(counts)){
    return(sweep(counts, 2, sf, "/", check.margin=FALSE))
  } else if(class(counts) == "dgCMatrix") {
    sep <- counts@p
    sep <- sep[-1] - sep[-length(sep)]
    j <- S4Vectors::Rle(1:length(sep), sep)
    counts@x <- counts@x / sf[as.integer(j)]
  }
  else if(class(counts)=="dgTMatrix"){
    j = counts@j
    counts@x = counts@x/sf[j+1]
  }
  else{
    stop(paste("cpm function for", class(counts)[1], "not supported"))
  }
  return(counts)
}

#' Convert a matrix of raw counts to a matrix of log2(Counts per Million + 1) values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to samples, and rows to genes.
#' 
#' @param counts a matrix, dgCMatrix, or dgTMatrix of count values.
#' 
#' @return a matrix, dgCMatrix, or dgTMatrix of log2(CPM + 1) values (matching input)
#' 
#' @export
#' 
logCPM <- function(counts) {
  
  norm.dat <- cpm(counts)
  
  if(is.matrix(norm.dat)){
    norm.dat <- log2(norm.dat + 1)
  } else {
    norm.dat@x <- log2(norm.dat@x + 1)
  }
  
  norm.dat
}

#' Compute correlations between each matching row or column of two matrices
#' 
#' e.g. row 1 of mat1 will be correlated with row 1 of mat2; row 2 of mat 1 with row 2 of mat2, etc.
#' 
#' The matrices must have the same number of rows or columns
#' 
#' @param mat1 a numeric matrix for correlation
#' @param mat2 a second numeric matrix for correlation
#' @param margin 1 for rows, 2 for columns (as for the MARGIN parameter of apply())
#' 
#' @return a numeric vector with paired correlation values
#' 
#' @export
#' 
pair_cor <- function(mat1, 
                     mat2,
                     margin = 1) {
  
  if(!margin %in% c(1,2)) {
    stop("margin must be either 1 (rows) or 2 (columns).")
  }
  
  if(margin == 1) {
    
    if(nrow(mat1) != nrow(mat2)) {
      stop("mat1 and mat2 must have an equal number of rows when margin = 1.")
    }
    
    
  } else if(margin == 2) {
    
    if(ncol(mat1) != ncol(mat2)) {
      stop("mat1 and mat2 must have an equal number of columns when margin = 2.")
    }
    
    mat1 <- t(mat1)
    mat2 <- t(mat2)
  }
  
  mat1 <- mat1 - rowMeans(mat1)
  mat2 <- mat2 - rowMeans(mat2)
  sd1 <- rowSds(mat1)
  sd2 <- rowSds(mat2)
  cors <- rowSums(mat1 * mat2) / ((ncol(mat1) - 1) * sd1 * sd2)  
  return(cors)
}


filter_by_size <- function(cat, min.size)
  {
    cat.size = table(cat)
    select.cat = names(cat.size)[cat.size >= min.size]
    cat = cat[cat %in% select.cat]
    if(is.factor(cat)){
      cat = droplevels(cat)
    }
  }

l2norm <- function(X, by="column")
{
  if (by=="column") {
    l2norm <- sqrt(Matrix::colSums(X^2))
    if (!any(l2norm==0)) {
      X=sweep(X, 2, l2norm, "/", check.margin=FALSE)
    }
    else{
      warning("L2 norms of zero detected for distance='Cosine, no transformation")
    }
  } else {
    l2norm <- sqrt(Matrix::rowSums(X^2))
    if (!any(l2norm==0)) {      
      X= X/l2norm
    }
    else{
      warning("L2 norms of zero detected for distance='Cosine'")
      X = X/ pmax(l2norm,1)
    }
  }
  X
}




#' Perform matrix operation using RCPP 
#' 
#' Currently supported operations are: sums, means, median, present, sqr_means
#' 
#' 
#' @param mat a numeric matrix
#' @param cl a cluster factor object
#' @param stats operation string. Currently support: sums, means, median, present, sqr_means
#' @param sparse boolean] to indicate if the matrix is sparse or not
#' @param transpose boolean to indicate if the matrix is transposed or not
#' @param parallel boolean to indicate if running the rcpp function in parallel or not
#' @param numThreads integer to indicate the number of threads
#' 
#' @return the resulting numeric matrix after applying the operation
#' 
#' @export
#' 
#' 
#'

get_cl_stats <- function(mat, 
                         cl, 
                         stats = c("sums","means","medians","present","sqr_sums","sqr_means"),
                         low.th=1,
                         parallel = c(FALSE,TRUE),
                         mc.cores = 1,...)
{
  if(!is.factor(cl)){
    cl = as.factor(cl)
  }
  library(RcppParallel)
  mc.cores = mc.cores
  setThreadOptions(numThreads = mc.cores)
  transpose=FALSE
  sparse=TRUE
  if(!all(names(cl) %in% colnames(mat))){
    transpose=TRUE
  }
  if(is.matrix(mat)){
    sparse=FALSE
  }
  if (sparse) {
    if (transpose) {
      if (parallel) { #sparse transpose parallel
        if (stats == "sums") {
          result=rcpp_get_cl_sums_RcppParallel_transpose(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_RcppParallel_transpose(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_RcppParallel_transpose(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_RcppParallel_transpose(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_RcppParallel_transpose(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_RcppParallel_transpose(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
        
      } else { #sparse transpose non-parallel
        if (stats == "sums") {
          result=rcpp_get_cl_sums_transpose(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_transpose(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_transpose(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_transpose(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_transpose(mat, cl)
        } else if (stats == "sqr_sums") {
          result = rcpp_get_cl_sqr_sums_transpose(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
        
      }
    } else { 
      if (parallel) { #sparse non-transpose parallel
        if (stats == "sums") {
          result=rcpp_get_cl_sums_RcppParallel(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_RcppParallel(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_RcppParallel(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_RcppParallel(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_RcppParallel(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_RcppParallel(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
        
      } else { #sparse non-transpose non-parallel
        if (stats == "sums") {
          result=rcpp_get_cl_sums(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
        
      }
    }
  } else {
    if (transpose) {
      if (parallel) {#dense transpose parallel        
        if (stats == "sums") { 
          result=rcpp_get_cl_sums_RcppParallel_transpose_dense(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_RcppParallel_transpose_dense(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_RcppParallel_transpose_dense(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_RcppParallel_transpose_dense(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_RcppParallel_transpose_dense(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_RcppParallel_transpose_dense(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
        
      } else { #dense transpose non-parallel
        
        if (stats == "sums") {
          result=rcpp_get_cl_sums_transpose_dense(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_transpose_dense(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_transpose_dense(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_transpose_dense(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_transpose_dense(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_transpose_dense(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
      }
    } else {
      if (parallel) { #dense non-transpose parallel
        
        if (stats == "sums") {
          result=rcpp_get_cl_sums_RcppParallel_dense(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_RcppParallel_dense(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_RcppParallel_dense(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_RcppParallel_dense(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_RcppParallel_dense(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_RcppParallel_dense(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
      } else { #dense non-transpose non-parallel
        
        if (stats == "sums") {
          result=rcpp_get_cl_sums_dense(mat, cl)
        } else if (stats == "means") {
          result=rcpp_get_cl_means_dense(mat, cl)
        } else if (stats == "medians") {
          result=rcpp_get_cl_medians_dense(mat, cl)
        } else if (stats == "present") {
          result=rcpp_get_cl_present_dense(mat, cl,lowth=low.th)
        } else if (stats == "sqr_means") {
          result=rcpp_get_cl_sqr_means_dense(mat, cl)
        } else if (stats == "sqr_sums") {
          result=rcpp_get_cl_sqr_sums_dense(mat, cl)
        } else {
          stop("the stats value is not supported")
        }
      }
    } 
  }
  result = result[,levels(cl),drop=F]
}


standardize <- function(X, by="column")
{  
  if(by=="column"){
    mean = Matrix::colMeans(X)
    X=sweep(X, 2, mean, "-")
    sd = sqrt(Matrix::colSums(X^2)/nrow(X))
    X = sweep(X, 2, sd, "/")
  }
  else{
    mean = Matrix::rowMeans(X)
    X=sweep(X, 1, mean, "-")
    sd = sqrt(Matrix::rowSums(X^2)/ncol(X))
    X = sweep(X, 1, sd, "/")    
  }
}

