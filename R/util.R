#' Convert a matrix of count values to CPM (Counts per Million)
#' 
#' @param counts a standard or sparse matrix
#' 
#' @return a matrix object of the same type as counts with normalized values
#' 
cpm <- function(counts)
  {
    library(Matrix)
    t(t(counts)*10^6/colSums(counts))
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

get_cl_mat <- function(cl)
  {
    if(!is.factor(cl)){
      cl <- factor(cl)
    }
    cl.mat <- sparseMatrix(i = 1:length(cl),  
                           j = as.integer(cl), 
                           x = 1)
    rownames(cl.mat) <- names(cl)
    colnames(cl.mat) <- levels(cl)
    return(cl.mat)
  }

get_cl_sums <- function(mat, cl)
  {
    require(Matrix)
    cl.mat <- get_cl_mat(cl)
    tmp <- Matrix::tcrossprod(mat[,rownames(cl.mat)], Matrix::t(cl.mat))
    cl.sums <- as.matrix(tmp)
    return(cl.sums)
  }

get_cl_means <- function(mat, cl)
  {
    cl.sums <- get_cl_sums(mat, cl)
    cl.size <- table(cl)
    cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
    return(cl.means)
}


get_cl_medians <- function(mat, cl)
{
  library(matrixStats)
  cl.med = do.call("cbind",tapply(names(cl), cl, function(x){
    rowMedians(as.matrix(mat[,x]))
  }))
  row.names(cl.med)=row.names(mat)
  return(cl.med)
}


sparse_cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans
  
  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
             crossprod(t(Ex))*(n-length(ii))
             )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}

calc_tau <- function(m, byRow=TRUE)
{
  if(!byRow){
    m = t(m)
  }
  m = m/rowMaxs(m)
  tau = rowSums(1 - m)/(ncol(m) - 1)
  tau[is.na(tau)]=0
  return(tau)
}


sample_cells <- function(cl,sample.size, weights=NULL)
{
  tmp = unique(cl)
  if(length(sample.size)==1){
    sample.size = setNames(rep(sample.size, length(tmp)), tmp)
  }
  cl.cells= split(names(cl),cl)
  sampled.cells = unlist(sapply(names(cl.cells), function(x){
    cells= cl.cells[[x]]
    if(sample.size[[x]]==length(cells)){
      return(cells)
    }
    to.sample = pmin(sample.size[[x]], length(cells))
    if(!is.null(weights)){
      sampled= sample(cells, to.sample, prob= weights[x])
    }
    else{
      sampled= sample(cells, to.sample)
    }
    sampled
  },simplify=FALSE))
}
  

sample_cells_by_genecounts <- function(cl, norm.dat, max.cl.size=200)
{
  select.cells=names(cl)
  if(is.matrix(norm.dat)){
    cell.gene.counts= colSums(norm.dat[,select.cells]>0)
  }
  else{
    cell.gene.counts= Matrix::colSums(norm.dat[,select.cells]>0)
  }
  cell.weights = cell.gene.counts - min(cell.gene.counts)+1
  sample_cells(cl, weights=cell.weights, max.cl.size=max.cl.size)  
}

translate_pair_name <- function(pairs, id.map)
{
  pairs.df = do.call("rbind", strsplit(pairs,"_"))
  new.pairs  = paste(id.map[as.character(pairs.df[,1])], id.map[as.character(pairs.df[,2])],sep="_")
  return(new.pairs)
}
