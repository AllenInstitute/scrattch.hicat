

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
rd_PCA <- function(norm.dat, select.genes=row.names(norm.dat), select.cells=colnames(norm.dat),sampled.cells=select.cells, max.pca=10, th=2, verbose=FALSE, method="zscore", mc.cores=1)
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
    require(parallel)
    rd.dat = parallel::pvec(select.cells, function(x){
      tmp.dat = norm.dat[row.names(rot), x,drop=F]
      rd.dat = as.matrix(Matrix::crossprod(tmp.dat, rot))
      return(list(rd.dat))
    }, mc.cores=mc.cores)
    rd.dat = do.call("rbind",rd.dat)
  }
  return(list(rd.dat=rd.dat, pca=pca))
}

get_PCA <- function(dat, max.pca, verbose=FALSE, method="zscore",th=2,fun="prcomp", rot=TRUE, init.pca = 200)
  {
    library(Matrix)
    library(stats)
    dat = as.matrix(dat)
    if(rot){
      dat = t(dat)
    }
    if(fun=="prcomp"){
      pca = stats::prcomp(dat,tol=0.01)
    }
    else{
      require("irlba")
      pca = prcomp_irlba(dat, n = min(init.pca, nrow(dat)));
    }
    if(method=="elbow"){
      dim.elbow = findElbowPoint(pca$sdev^2)
      if(verbose){
        cat("elbow dim:", dim.elbow, "\n")
      }      
      select = 1:min(max.pca,dim.elbow)      
    }
    else if(method=="zscore"){
      pca.importance = summary(pca)$importance
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
