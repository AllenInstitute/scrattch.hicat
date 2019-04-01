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


get_pairs <- function(pairs.str) {
  pairs <- as.data.frame(do.call("rbind", strsplit(pairs.str, "_")), 
                         stringsAsFactors = FALSE)
  
  row.names(pairs) <- pairs.str
  colnames(pairs) <- c("P1", "P2")
  
  return(pairs)
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


#' Generate a sparse matrix one-hot representation of clusters x samples
#' 
#' @param cl a cluster factor object
#' 
#' @return a sparse, one-hot matrix indicating which cluster(columns) each sample (rows) belongs to.
#' @export
#' 
get_cl_mat <- function(cl) {
  
  if(!is.factor(cl)){
    cl <- as.factor(cl)
  }
  cl = droplevels(cl)
  cl.mat <- Matrix::sparseMatrix(i = 1:length(cl),  
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
#' @export
#' 
get_cl_sums <- function(mat, 
                        cl) {
  
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
#' @export
#' 
get_cl_means <- function(mat, 
                         cl) {
  
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
#' @return a matrix of genes (rows) x clusters (columns) with medians for each cluster
#' @export
#' 
get_cl_medians <- function(mat, cl)
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


#' Compute cluster proportions for each row in a matrix
#'
#' @param mat a gene (rows) x samples(columns) sparse matrix
#' @param cl A cluster factor object
#' @param thresshold The minimum expression value used to binarize the results
#'
#' @return a matrix of genes (rows) x cluster(columns) with proportions for each cluster
#' @export
#' 
get_cl_prop <- function(mat, 
                        cl, 
                        thresshold = 1) {
  
  cl.mat <- get_cl_mat(cl)
  
  cl.prop <- Matrix::tcrossprod(mat[,rownames(cl.mat)] > thresshold, Matrix::t(cl.mat))
  
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
calc_tau <- function(m, 
                     byRow = TRUE) {
  if(!byRow){
    m <- t(m)
  }
  m <- m / matrixStats::rowMaxs(m)
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
sample_cells<- function(cl, 
                        sample.size, 
                        weights = NULL) {
  
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

#' Convert a matrix of raw counts to a matrix of Counts per Million values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to samples, and rows to genes.
#' 
#' @param counts a matrix of count values.
cpm <- function(counts) {
  
  library(Matrix)
  
  sf <- Matrix::colSums(counts) / 1e6
  
  if(is.matrix(counts)){    
    return(t(t(counts) / sf))
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

logCPM <- function(counts)
  {
    norm.dat = cpm(counts)
    if(is.matrix(norm.dat)){
      norm.dat = log2(norm.dat+1)
    }
    else{
      norm.dat@x = log2(norm.dat@x + 1)
    }
    norm.dat
  }

#' Compute correlation each row of matrix1 with the corresponding row of matrix2 
#' matrix1 and matrix2 must have the same dimemsion. 
#' 
#' 
pair_cor <- function(mat1, mat2)
  {
    mat1 = mat1 - rowMeans(mat1)
    mat2 = mat2 - rowMeans(mat2)
    sd1 = rowSds(mat1)
    sd2 = rowSds(mat2)
    rowSums(mat1 * mat2)/((ncol(mat1)-1)*sd1*sd2)
  }
