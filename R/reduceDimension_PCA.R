prcomp.irlba <- function(x, max.rank=500, maxit=1000, tol=1e-05, center=TRUE,...)
  {
    library(irlba)
    s <- irlba(x, nv=max.rank, nu=max.rank, maxit = maxit, tol=tol, ...)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    dimnames(s$v) <- list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v)
    r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
  }


rd_PCA <- function(norm.dat, select.genes=row.names(norm.dat), select.cells=colnames(norm.dat),sampled.cells=select.cells, max.pca=10, th=2)
{
  require(Matrix)
  
  pca = prcomp(t(as.matrix(norm.dat[select.genes,sampled.cells])),tol=0.01)
  pca.importance = summary(pca)$importance
  v = pca.importance[2,]
  
  select= which((v - mean(v))/sd(v)>th) 
  tmp = head(select,max.pca)
  if(length(tmp)==0){
    return(NULL)
  }
  if(length(sampled.cells)< length(select.cells)){
    rot  =  pca$rotatio[,tmp,drop=F]
    tmp.dat = norm.dat[row.names(rot), select.cells,drop=F]
    rd.dat = as.matrix(Matrix::t(tmp.dat)  %*% rot)
  }
  else{
    rd.dat=pca$x[,tmp,drop=F]
  }
  return(list(rd.dat=rd.dat, pca=pca))
}
  
