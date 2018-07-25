#' Convert a matrix of count values to CPM (Counts per Million)
#' 
#' @param counts a standard or sparse matrix
#' 
#' @return a matrix object of the same type as counts with normalized values
#' 
cpm <- function(counts)
  {
    library(Matrix)
    t(t(counts) * 10^6 / colSums(counts))
  }

#' Convert matrix row/column positions to vector position
#' 
#' @param m a matrix object
#' @param rows row positions, either as indices or matches to row names
#' @param cols row positions, either as indices or matches to row names
#' 
#' @return a vector containing the numeric position at [rows,cols]
#' 
get_pair_matrix_coor <- function(m, rows, cols)
  {
    if(!is.numeric(rows)){
      rows <- match(rows, row.names(m))
    }
  
    if(!is.numeric(cols)){
      cols <- match(cols, colnames(m))
    }
  
    coor <- (cols - 1) * nrow(m) + rows

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
get_pair_matrix <- function(m, rows, cols)
  {
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
set_pair_matrix <- function(m, rows, cols, vals)
  {
    coor <- get_pair_matrix_coor(m, rows, cols)
    m[coor] <- vals
    return(m)
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
convert_pair_matrix <- function(pair.num, 
                                l = NULL,
                                directed = FALSE)
  {
    pairs <- do.call("rbind", strsplit(names(pair.num),"_"))
    
    if(is.null(l)){
      l <- sort(unique(as.vector(pairs)))
    }
    
    n.cl <- length(l)
    
    pair.num.mat <- matrix(0, nrow = n.cl, ncol = n.cl)
    rownames(pair.num.mat) <- l
    colnames(pair.num.mat) <- l
    
    for(i in 1:nrow(pairs)){
      if(directed) {
        pair.num.mat[pairs[i,1], pairs[i,2]] <- pair.num[i]
      } else if(!directed) {
        pair.num.mat[pairs[i,2], pairs[i,1]] <- pair.num[i]
      }
    }
    
    return(pair.num.mat)
  }

#' Generate a sparse matrix one-hot representation of clusters x samples
#' 
#' @param cl a cluster factor object
#' 
#' @return a sparse, one-hot matrix indicating which cluster(columns) each sample (rows) belongs to.
#' 
get_cl_mat <- function(cl)
  {
    library(Matrix)
  
    if(!is.factor(cl)){
      cl <- as.factor(cl)
    }
  
    cl.mat <- sparseMatrix(i = 1:length(cl),  
                           j = as.integer(cl), 
                           x = 1)
    
    rownames(cl.mat) <- names(cl)
    colnames(cl.mat) <- levels(cl)
    
    return(cl.mat)
  }

#' Compute cluster sums for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with sums for each cluster
#' 
get_cl_sums <- function(mat, cl)
  {
    library(Matrix)
  
    cl.mat <- get_cl_mat(cl)
    
    cl.sums <- Matrix::tcrossprod(mat[,rownames(cl.mat)], Matrix::t(cl.mat))
    
    cl.sums <- as.matrix(cl.sums)
    
    return(cl.sums)
  }

#' Compute cluster means for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with means for each cluster
#'
get_cl_means <- function(mat, cl)
  {
    library(Matrix)
  
    cl.sums <- get_cl_sums(mat, cl)
    
    cl.size <- table(cl)
    
    cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
    
    return(cl.means)
}

#' Compute cluster medians for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with edians for each cluster
#'
get_cl_medians <- function(mat, cl)
{
  library(Matrix)
  library(matrixStats)
  
  cl.med <- do.call("cbind",
                    tapply(names(cl), 
                           cl, 
                           function(x){
                             rowMedians(as.matrix(mat[,x]))
                           }
                    )
  )
  
  rownames(cl.med) <- rownames(mat)
  
  return(cl.med)
}

#' Compute correlation scores for columns of a sparse matrix
#' 
#' @param m a sparse matrix, preferrably dgCMatrix
#' 
#' @return a matrix of correlation values between each column of m
#'  
sparse_cor <- function(m){
  library(Matrix)
  
  n_rows <- nrow(m)
  n_cols <- ncol(m)
  
  ii <- unique(m@i) + 1 # rows with a non-zero element
  
  Ex <- colMeans(m)
  
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
calc_tau <- function(m, byRow = TRUE)
{
  if(!byRow){
    m <- t(m)
  }
  m <- m / rowMaxs(m)
  tau <- rowSums(1 - m) / (ncol(m) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}

#' Downsample cells from each cluster
#' 
#' @param cl A cluster factor object
#' @param sample.size A maximum number of cells to take from each cluster, or a named numeric object with the number of cells to be sampled and names matching cluster levels
#' 
#' @return A cluster factor object containing sampled cells
#' 
sample_cells<- function(cl, sample.size, weights = NULL)
{
  
  n_cl <- unique(cl)
  
  if(length(sample.size) == 1) {
    
    sample.size <- setNames(rep(sample.size, length(n_cl)), n_cl)
  
  }
  
  cl.cells <- split(names(cl), cl)
  
  sampled.cells <- unlist(sapply(names(cl.cells), 
                                 function(x){
    cells <- cl.cells[[x]]
    if(sample.size[[x]] == length(cells)){
      return(cells)
    }
    
    to.sample <- pmin(sample.size[[x]], length(cells))
    
    if(!is.null(weights)){
      sampled <- sample(cells, to.sample, prob = weights[cells])
    }
    
    else{
      sampled <- sample(cells, to.sample)
    }
    
    sampled
    
  }, simplify = FALSE))
  
  return(sampled.cells)
}
