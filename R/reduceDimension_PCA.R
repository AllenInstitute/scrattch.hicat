rd_PCA <- function(norm.dat, select.genes, select.cells,sampled.cells=select.cells, max.pca=10, w=NULL)
{
  if(is.null(w)){
    pca = prcomp(t(as.matrix(norm.dat[select.genes,sampled.cells])),tol=0.01)
    pca.importance = summary(pca)$importance
    v = pca.importance[2,]
  }
  else{
    dat <- scale(norm.dat[select.genes, select.cells], center = TRUE, scale = TRUE)
    w = w[select.genes, select.cells]
    
    for (x in 1:nsamp) {
      for (y in 1:nsamp) {
        wt1 <- w[, x] * w[, y]
        cov1 <- cov.wt(e[, c(x, y)], wt = wt1, center = FALSE)$cov
        e.cov[x, y] <- cov1[1, 2]
      }
    }
        
    nsamp <- ncol(e)
    wt <- crossprod(w)
    cov = cov.wt(dat, wt, center=FALSE)$cov
    eig1 <- eigen(cov, symmetric = TRUE)
    eig.val <- eig1$values
    eig.val[is.na(eig.val) | eig.val < 0] <- 0
    eig.vec <- eig1$vectors
    dimnames(eig.vec) <- list(colnames(e), paste0("PC", 1:ncol(eig.vec)))
    pca1 <- list()
    pca1$sdev <- sqrt(eig.val)
    pca1$rotation <- eig.vec
    pca1$x <- e %*% eig.vec
  }

  select= which((v - mean(v))/sd(v)>2) 
  tmp = head(select,max.pca)
  if(length(sampled.cells)< length(select.cells)){
    rot  =  pca$rotatio[,tmp]
    tmp.dat = norm.dat[row.names(rot), select.cells]
    rd.dat = as.matrix(t(tmp.dat)  %*% rot)
  }
  else{
    rd.dat=pca$x[,tmp,drop=F]
  }
  return(rd.dat)
}
  
