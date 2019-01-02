#' Doublet detection in single-cell RNA sequencing data
#'
#' This function generaetes artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets.
#'
#' @param seu A fully-processed Seurat object (i.e. after normalization, variable gene definition,
#' scaling, PCA, and tSNE).
#' @param expected.doublets The number of doublets expected to be present in the original data.
#' This value can best be estimated from cell loading densities into the 10X/Drop-Seq device.
#' @param porportion.artificial The proportion (from 0-1) of the merged real-artificial dataset
#' that is artificial. In other words, this argument defines the total number of artificial doublets.
#' Default is set to 25%, based on optimization on PBMCs (see McGinnis, Murrow and Gartner 2018, BioRxiv).
#' @param proportion.NN The proportion (from 0-1) of the merged real-artificial dataset used to define
#' each cell's neighborhood in PC space. Default set to 1%, based on optimization on PBMCs (see McGinnis,
#' Murrow and Gartner 2018, BioRxiv).
#' @return An updated Seurat object with metadata for pANN values and doublet predictions.
#' @export
#' @examples
#' seu <- doubletFinder(seu, expected.doublets = 1000, proportion.artificial = 0.25, proportion.NN = 0.01)


doubletFinder <- function(data, select.genes, proportion.artificial = 0.25, k = ncol(data) * 0.005, sc.th=0.5) {
  library(RANN)
  if (expected.doublets == 0) {  stop("Need to set number of expected doublets...")  }

  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  real.cells <- colnames(dat)

  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1-proportion.artificial)-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  
  data_wdoublets <- cbind(data, doublets)  
  data_wdoublets@x = log2(data_wdoublets@x+1)
  sampled.cells = sample(1:ncol(data), 10000)

  print("Running PCA")
  rd.dat = rd_PCA(data_wdoublets, select.genes=select.markers, select.cells= colnames(data_wdoublets), th=0, max.pca=50, sampled.cells= sampled.cells)

  print("Initialize pANN structure") 
  knn.result = RANN::nn2(rd.dat, k=k)
  knn.idx = knn.result[[1]]
  knn.dist = knn.result[[2]]
  
  doublet.freq = knn.idx  > ncol(data)
  doublet.score = rowMeans(doublet.freq)
  doublet.score = doublet.score[1:ncol(data)]
  doublet.score = setNames(doublet.score, colnames(data))

  doublet.candidate = doublet.score > sc.th
  doublet.knn.dist =knn.dist[doublet.candidate,][as.vector(doublet.freq[doublet.candidate,])]
  th = mean(doublet.knn.dist) - sd(doublet.knn.dist)
  
  
  
  return(doublet.score)
}
