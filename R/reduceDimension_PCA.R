
                                        #' prcomp irlba
#'
#' @param x 
#' @param max.rank 
#' @param maxit 
#' @param tol 
#' @param center 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
prcomp.irlba <- function(x, max.rank=500, maxit=1000, tol=1e-05, center=TRUE,...)
  {
    library(irlba)
    s <- irlba::irlba(x, nv=max.rank, nu=max.rank, maxit = maxit, tol=tol, ...)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    dimnames(s$v) <- list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v)
    r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
  }


#' Title
#'
#' @param norm.dat 
#' @param select.genes 
#' @param select.cells 
#' @param sampled.cells 
#' @param max.pca 
#' @param th 
#'
#' @return
#' @export
#'
#' @examples
rd_PCA <- function(norm.dat, select.genes=row.names(norm.dat), select.cells=colnames(norm.dat),sampled.cells=select.cells, max.pca=10, th=2, verbose=FALSE, method="zscore")
{
  library(Matrix)
  library(stats)

  tmp = get_PCA(norm.dat[select.genes, sampled.cells], max.pca=max.pca, verbose=verbose,th=th,  method=method)
  if(is.null(tmp)){
    return(NULL)
  }
  rot = tmp$rot
  rd.dat = tmp$rd.dat
  pca  = tmp$pca
  if(length(sampled.cells)< length(select.cells)){
    if(verbose){
      print("Project")
    }
    tmp.dat = norm.dat[row.names(rot), select.cells,drop=F]
    rd.dat = as.matrix(Matrix::crossprod(tmp.dat, rot))  
  }
  return(list(rd.dat=rd.dat, pca=pca))
}

get_PCA <- function(dat, max.pca, verbose=FALSE, method="zscore",th=2)
  {
    library(Matrix)
    library(stats)
    
    pca = stats::prcomp(t(as.matrix(dat)),tol=0.01)
    if(method=="elbow"){
      dim.elbow = findElbowPoint(pca$sdev^2)
      if(verbose){
        cat("elbow dim:", dim.elbow, "\n")
      }      
      select = 1:min(max.pca,dim.elbow)      
    }
    else if(method=="zscore"){
      v = pca.importance[2,]
      select= which((v - mean(v))/sd(v)>th)
      select = head(select,max.pca)
    }
    else{
      Stop("Unknown method")
    }
    if(length(select)==0){
      return(NULL)
    }
    rot  =  pca$rotatio[,select,drop=F]
    rd.dat = pca$x[,select,drop=F]
    return(list(rot=rot, rd.dat = rd.dat,pca=pca))
  }



top_loading_genes <- function(rot,top.n=10)
  {
    top.genes= apply(rot, 2, function(x){
      tmp=head(order(abs(x),decreasing=T), top.n)
      split(row.names(rot)[tmp], x[tmp]>0)
    })
  }


###Taken from https://github.com/kevinblighe/PCAtools/blob/master/R/findElbowPoint.R
findElbowPoint <- function(variance) {
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  
                                        # Finding distance from each point on the curve to the diagonal.
  dy <- -diff(range(variance))
  dx <- length(variance) - 1
  l2 <- sqrt(dx^2 + dy^2)
  dx <- dx/l2
  dy <- dy/l2
  
  dy0 <- variance - variance[1]
  dx0 <- seq_along(variance) - 1
  
  parallel.l2 <- sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
  normal.x <- dx0 - dx * parallel.l2
  normal.y <- dy0 - dy * parallel.l2
  normal.l2 <- sqrt(normal.x^2 + normal.y^2)
  
  #Picking the maximum normal that lies below the line.
  #If the entire curve is above the line, we just pick the last point.
  below.line <- normal.x < 0 & normal.y < 0
  if (!any(below.line)) {
    length(variance)
  } else {
    which(below.line)[which.max(normal.l2[below.line])]
  }
}
