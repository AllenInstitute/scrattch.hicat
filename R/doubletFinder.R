#' Doublet detection in single-cell RNA sequencing data
#'
#' Adopted from https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0 
#' (https://github.com/chris-mcginnis-ucsf/DoubletFinder) 
#'
#' This function generates artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets.
#'
#' 
#' @param data gene x sample matrix with counts (non-normalized)
#' @param select.genes list of genes with highest variance between samples
#' @param proportion.artificial The proportion (from 0-1) of the merged real-artificial dataset
#' that is artificial. In other words, this argument defines the total number of artificial doublets.
#' Default is set to 20\%
#' @param k The number of nearest neighbours of the merged real-artificial dataset used to define
#' each cell's neighborhood in PC space. Value is the minimum of 1% of cells sampled or 100.
#'
#'  
#' @param plot 
#' 
#' @return An list of doublet.scores per samples and plots depicting the doublet scores for cells and artificial doublets.
#' 
doubletFinder <- function(data, select.genes, proportion.artificial = 0.20,
                          k = NULL, plot=FALSE) {
  
    if(is.null(k)) {
      k <- round(pmin(100, ncol(data) * 0.01))
      k <- round(pmax(10, k))
    }
      
  
  library(RANN)
  
  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  real.cells <- colnames(data)
  
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1-proportion.artificial)-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  
  data_wdoublets <- cbind(data, doublets)  
  norm.dat = scrattch.hicat::logCPM(data_wdoublets)
  
  print("Running PCA")
  if(ncol(data) > 10000){
    sampled.cells = sample(1:ncol(data), pmin(ncol(data),10000))
    rd.dat = rd_PCA(norm.dat, select.genes=select.genes, select.cells= colnames(data_wdoublets), th=0, max.pca=50, sampled.cells= sampled.cells)
  }
  else{
    rd.dat = rd_PCA(norm.dat, select.genes=select.genes, select.cells= colnames(data_wdoublets), th=0, max.pca=50)
  }
  rd.dat <-rd.dat$rd.dat
  print("Initialize pANN structure") 
  knn.result = RANN::nn2(rd.dat, k=k)
  knn.idx = knn.result[[1]]
  knn.dist = knn.result[[2]]
  
  num = ncol(data)
  knn.result1 = RANN::nn2(rd.dat[1:num,], rd.dat[(num+1):nrow(rd.dat),], k = 10)
  knn.dist1 = knn.result1[[2]]
  
  dist.th = mean(as.vector(knn.dist1)) + 1.64 * sd(as.vector(knn.dist1))
  
  doublet.freq = knn.idx  > ncol(data) & knn.dist < dist.th
  doublet.score = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)]))
  names(doublet.score) = row.names(rd.dat)
  artificial.doublet.score=doublet.score[colnames(doublets)]
  doublet.score=doublet.score[colnames(data)]
    
  if (plot == TRUE) {
    print("plotting")
    ds = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)])) 
    ds=as.data.frame(ds)
    ds$sample <- colnames(data_wdoublets)
    
    ds$group <- ""
    idx <- startsWith(ds$sample,"X")
    ds[idx, "group"] <- "artifical doublets"
    idx <- !startsWith(ds$sample,"X")
    ds[idx, "group"] <- "samples"
    
    plot.title <- gsub("^.*?-","",ds[1,2])
    
    p=ggplot2::ggplot(ds, aes(x = doublet.score, fill=group, color=group)) +geom_density(alpha=0.4)+scale_color_manual(values=c("#F9627D","#2F3061"))+scale_fill_manual(values=c("#F9627D","#2F3061")) +labs(title=plot.title)
    
    return(list(doublet.score,artifical.doublet.score, p))
    
  } else {  return(list(doublet.score,artificial.doublet.score)) }
  
}



