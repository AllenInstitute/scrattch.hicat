rd_PCA <- function(norm.dat, 
                   select.genes = NULL, 
                   select.cells = NULL,
                   sampled.cells = NULL, 
                   max.pca = 10, 
                   th = 2,
                   w = NULL) {
  
  if(is.null(select.genes)) {
    select.genes <- row.names(norm.dat)
  }
  if(is.null(select.cells)) {
    select.cells <- colnames(norm.dat)
  }
  if(is.null(sampled.cells)) {
    sampled.cells <- select.cells
  }
  
  if(is.null(w)){
    pca <- prcomp(t(as.matrix(norm.dat[select.genes, sampled.cells])),
                  tol = 0.01)
    pca.importance <- summary(pca)$importance
    v <- pca.importance[2,]
  } else {
    dat <- scale(norm.dat[select.genes, select.cells], 
                 center = TRUE, 
                 scale = TRUE)
    
    w <- w[select.genes, select.cells]
    
    for (x in 1:nsamp) {
      for (y in 1:nsamp) {
        wt1 <- w[, x] * w[, y]
        
        cov1 <- cov.wt(e[, c(x, y)], 
                       wt = wt1, 
                       center = FALSE)$cov
        
        e.cov[x, y] <- cov1[1, 2]
      }
    }
        
    nsamp <- ncol(e)
    wt <- crossprod(w)
    
    cov <- cov.wt(dat, 
                  wt, 
                  center=FALSE)$cov
    
    eig1 <- eigen(cov, symmetric = TRUE)
    eig.val <- eig1$values
    eig.val[is.na(eig.val) | eig.val < 0] <- 0
    eig.vec <- eig1$vectors
    dimnames(eig.vec) <- list(colnames(e), paste0("PC", 1:ncol(eig.vec)))
    
    pca1 <- list(sdev = sqrt(eig.val),
                 rotation = eig.vec,
                 x = e %*% eig.vec)
  }
  
  select <- which((v - mean(v))/sd(v)>th) 
  tmp <- head(select, max.pca)
  
  if(length(tmp) == 0) {
    return(NULL)
  }
  
  if(length(sampled.cells) < length(select.cells)) {
    rot <- pca$rotatio[,tmp , drop = FALSE]
    tmp.dat <- norm.dat[row.names(rot), select.cells, drop = FALSE]
    rd.dat <- as.matrix(Matrix::t(tmp.dat) %*% rot)
  } else {
    rd.dat <- pca$x[, tmp, drop = FALSE]
  }
  
  return(list(rd.dat = rd.dat, 
              pca = pca))
}
  
