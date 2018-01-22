# cf Habib 2016 Div-Seq - SI pp 16-17
WeightedPca <- function(e, w, e.center = FALSE, e.scale = FALSE, 
                        max.pcs = min(dim(e))) {
  
  e <- scale(e, center = e.center, scale = e.scale)
  
  # Use custom code
  nsamp <- ncol(e)
  e.cov <- matrix(NA, nsamp, nsamp)
  for (x in 1:nsamp) {
    for (y in 1:nsamp) {
      wt1 <- w[, x] * w[, y]
      cov1 <- cov.wt(e[, c(x, y)], wt = wt1, center = FALSE)$cov
      e.cov[x, y] <- cov1[1, 2]
    }
  }
  eig1 <- eigen(e.cov, symmetric = TRUE)
  eig.val <- eig1$values
  eig.val[is.na(eig.val) | eig.val < 0] <- 0
  eig.vec <- eig1$vectors
  dimnames(eig.vec) <- list(colnames(e), paste0("PC", 1:ncol(eig.vec)))

  # Create pseudo prcomp object
  pca1 <- list()
  pca1$sdev <- sqrt(eig.val)
  pca1$rotation <- eig.vec
  pca1$x <- e %*% eig.vec
  
  # Use scde
  # require(scde)
  # bwpca1 <- bwpca(as.matrix(e), as.matrix(w), npcs = max.pcs, center = FALSE)
  # colnames(bwpca1$scores) <- paste0("PC", 1:ncol(bwpca1$scores))
  # pca1 <- list()
  # pc.order <- order(bwpca1$sd, decreasing = TRUE)  # Weighted variance explained is not in order
  # pca1$sdev <- bwpca1$sd[pc.order]
  # pca1$rotation <- bwpca1$rotation[, pc.order]
  # pca1$x <- bwpca1$scores[, pc.order]
  
  return(pca1)
}


###w.mat should be row centered
cov.wt1 <- function(e.mat, w.mat)
{
  nsamp <- ncol(e.mat)
  e.cov <- matrix(NA, nsamp, nsamp)
  for (x in 1:nsamp) {
    for (y in 1:nsamp) {
      wt1 <- w.mat[,x]*w.mat[, y]
      cov1 <- cov.wt(e.mat[, c(x, y)], wt = wt1, center = FALSE)$cov
      e.cov[x, y] <- cov1[1, 2]
    }
  }
}

###w.mat should be row centered
cov.wt1 <- function(e.mat, w.mat)
{
  nsamp <- ncol(e.mat)
  e.cov <- matrix(NA, nsamp, nsamp)
  for (x in 1:nsamp) {
    for (y in 1:nsamp) {
      wt1 <- w.mat[,x]*w.mat[, y]
      cov1 <- cov.wt(e.mat[, c(x, y)], wt = wt1, center = FALSE)$cov
      e.cov[x, y] <- cov1[1, 2]
    }
  }
  return(e.cov)
}

###w.mat should be row centered
cov.wt2 <- function(e.mat, w.mat)
{
  w.mat <- t(t(w.mat)/ colSums(w.mat))
  nsamp <- ncol(e.mat)
  e.cov <- matrix(NA, nsamp, nsamp)
  e.mat = sqrt(w.mat) * e.mat
  colSums(w.mat^2)
  e.cov = crossprod(e.mat)/(1 - sum(colSums(w.mat^2)))
  return(e.cov)
}


e = matrix(runif(10*5), nrow=10)
e = e * runif(10)*10
e = e - rowMeans(e)
w = matrix(runif(10*5), nrow=10)
cov1=cov.wt1(e,w)
cov2=cov.wt2(e,w)




