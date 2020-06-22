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
rd_PCA <- function(norm.dat, select.genes=row.names(norm.dat), select.cells=colnames(norm.dat),sampled.cells=select.cells, max.pca=10, th=2)
{
  library(Matrix)
  library(stats)
  
  pca = stats::prcomp(t(as.matrix(norm.dat[select.genes,sampled.cells])),tol=0.01)

  #pca.importance = summary(pca)$importance
  #v = pca.importance[2,]
  #m = mean(tail(v, length(v))/2)
  #sd = sd(tail(v, length(v))/2)
  #z= (v - m)/sd
  select = 1:min(max.pca,findElbowPoint(pca$sdev^2))
  if(length(select)==0){
    return(NULL)
  }
  if(length(sampled.cells)< length(select.cells)){
    rot  =  pca$rotatio[,select,drop=F]
    tmp.dat = norm.dat[row.names(rot), select.cells,drop=F]
    rd.dat = as.matrix(crossprod(tmp.dat, rot))
  }
  else{
    rd.dat=pca$x[,select,drop=F]
  }
  return(list(rd.dat=rd.dat, pca=pca))
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
