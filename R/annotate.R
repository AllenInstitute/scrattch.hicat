# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
#
# map_by_cor()
#   get_cl_means() util.R
#
# map_cl_summary()
#   map_by_cor() annotate.R
#
# predict_annotate_cor()
#   map_by_cor() annotate.R
#   compare_annotate() annotate.R
# 
# map_sampling()
#   map_by_cor() annotate.R
#
# map_cv()
#   map_by_cor() annotate.R
#
# compare_annotate()
#
# match_cl()
#   get_cl_means() util.R
# 
# find_low_quality_cl
#   get_de_matrix() de.genes.R
#   get_pair_matrix() util.R

#' Map samples to a training dataset by correlation
#' 
#' @param train.dat Training data matrix, usually log-transformed CPM
#' @param train.cl Training cluster factor object
#' @param test.dat Data for cells to map to the training set. Should have the same genes as train.dat.
#' @param method Which statistic to compare. "median" or "mean". Default is "median".
#' 
#' @return a list object containing two objects:
#' \itemize{
#' \item pred.df: a data.frame with two columns, pred.cl and pred.score with the predicted cluster and correlation scores.
#' \item cor.matrix: a matrix object with correlation scores for each cluster.
#' }
#' 
#' @export
#' 
map_by_cor <- function(train.dat, 
                       train.cl, 
                       test.dat,
                       method = "median") {
  
  method <- match.arg(arg = method, 
                      choices = c("mean","median"))
  
  # Get medians or means for each cluster
  if(method == "median"){
    cl.meds <- tapply(names(train.cl), 
                      train.cl, 
                      function(x) {
                        train.mat <- train.dat[, x, drop = F]
                        train.mat <- as.matrix(train.mat)
                        matrixStats::rowMedians(train.mat)
                      }
    )
    
    cl.dat <- do.call("cbind", cl.meds)
  } else {
    cl.dat <- get_cl_means(train.dat, train.cl)
  }
  row.names(cl.dat) <- row.names(train.dat)
  
  # Perform correlations
  if(!is.matrix(test.dat) & nrow(test.dat)*ncol(test.dat) > 1e8){
    test.cl.cor <- qlcMatrix::corSparse(test.dat, cl.dat)
    colnames(test.cl.cor) = colnames(cl.dat)
    row.names(test.cl.cor) = colnames(test.dat)
  } else{
    test.cl.cor <- cor(as.matrix(test.dat), cl.dat)
  }
  test.cl.cor[is.na(test.cl.cor)] <- 0
  
  # Find maximum correlation
  max.cl.cor <- apply(test.cl.cor, 1, which.max)
  pred.cl <- colnames(test.cl.cor)[max.cl.cor]
  pred.cl <- setNames(pred.cl, row.names(test.cl.cor))
  
  # Get maximum correlation values
  pred.score <- apply(test.cl.cor, 1, max)
  
  # Convert to factor if train.cl was a factor and match levels.
  if(is.factor(train.cl)){
    pred.cl <- setNames(factor(pred.cl, levels = levels(train.cl)), names(pred.cl))
  }
  
  # Output results
  pred.df <- data.frame(pred.cl = pred.cl,
                        pred.score = pred.score)
  
  out_list <- list(pred.df = pred.df,
                   cor.matrix = test.cl.cor)
  
  return(out_list)    
}

#' Map a dataset to a reference, and compare existing cluster calls to the reference comparison
#' 
#' @param ref.dat Training data matrix, usually log-transformed CPM
#' @param ref.cl Training cluster factor object
#' @param map.dat Data for cells to map to the training set. Should have the same genes as train.dat.
#' @param map.cl Cluster assignments for the training set to compare to results of mapping.
#' @param method Which statistic to compare. "median" or "mean". Default is "median".
#' 
#' @return a list object with two objects:  
#' \itemize{
#' \item map.df: A data.frame with the mapping results for each sample in map.dat to the reference
#' \item cl.map.df: A data.frame with cluster-level frequency of mapping for each cluster in map.cl to ref.cl
#' }
#' 
#' @export
#' 
map_cl_summary <- function(ref.dat, 
                           ref.cl, 
                           map.dat, 
                           map.cl,
                           method = "median") {
  
  method <- match.arg(arg = method,
                      choices = c("mean","median"))
  
  # Map the training set to the reference
  map.result <- map_by_cor(train.dat = ref.dat, 
                           train.cl = ref.cl, 
                           test.dat = map.dat,
                           method = method)
  cor.matrix <- map.result$cor.matrix
  
  map.df <- map.result$pred.df
  colnames(map.df)[1] <- "map.cl"
  map.df$org.cl <- map.cl[row.names(map.df)]
  
  # Compute the fraction of times each sample was mapped to each cluster
  cl.size <- table(map.cl)
  cl.map.df <- as.data.frame(with(map.df, table(org.cl, map.cl)))
  cl.map.df$Prob <- round(cl.map.df$Freq / cl.size[as.character(cl.map.df$org.cl)], digits = 2)
  
  # Compute the mean fraction of mapping for all samples in the training set clusters
  # to the training set clusters.
  cl.map.df$pred.score <- 0
  for(i in 1:nrow(cl.map.df)){
    select <- names(map.cl)[map.cl == as.character(cl.map.df$org.cl[i])]
    cl.map.df$pred.score[i] <- mean(cor.matrix[select, as.character(cl.map.df$map.cl[i])])
  }
  cl.map.df$pred.score <- round(cl.map.df$pred.score, digits = 2)
  
  # Remove comparisons with no mapping
  cl.map.df <- cl.map.df[cl.map.df$Freq > 0, ]
  
  # Return output
  out_list <- list(map.df = map.df,
                   cl.map.df = cl.map.df)
  
  return(out_list)
}


#' Predict annotations by cluster correlation
#' 
#' This function performs map_by_cor(), then compare_annotate().
#'
#' @param cl a cluster factor object for data to map to the reference
#' @param norm.dat a normalized data matrix for data to map to the reference
#' @param ref.markers a set of reference marker genes
#' @param ref.cl a reference cluster factor object
#' @param ref.cl.df a reference cl.df data.frame that describes the reference clusters
#' @param ref.norm.dat a reference normalized data matrix
#' @param method "median" or "mean". Default is "median".
#' @param reorder Whether or not to reorder the input clusters based on the reference.
#'
#' @return a list object with annotation results
#' 
#' @export
#'
predict_annotate_cor <- function(cl, 
                                 norm.dat, 
                                 ref.markers, 
                                 ref.cl, 
                                 ref.cl.df, 
                                 ref.norm.dat, 
                                 method = "median", 
                                 reorder = TRUE) {
  
  method <- match.arg(arg = method,
                      choices = c("mean", "median"))
  
  map_results <- map_by_cor(ref.norm.dat[ref.markers,], 
                            ref.cl, 
                            norm.dat[ref.markers, names(cl)],
                            method = method)
  
  pred.cl <- setNames(factor(as.character(map_results$pred.df$pred.cl), 
                             levels = row.names(ref.cl.df)), 
                      row.names(map_results$pred.df))
  
  map_results$annotate <- compare_annotate(cl, 
                                           pred.cl, 
                                           ref.cl.df, 
                                           reorder = reorder)
  
  return(map_results)
}

#' Perform bootstrapped mapping using a fraction of provided marker genes.
#' 
#' @param train.dat Training data matrix, usually log-transformed CPM
#' @param train.cl Training cluster factor object
#' @param test.dat Data for cells to map to the training set. Should have the same genes as train.dat.
#' @param markers A vector of marker gene symbols to use for comparisons
#' @param markers.perc The fraction of randomly sampled markers to use for each round. Default = 0.8.
#' @param iter Number of iterations to perform. Default = 100.
#' @param method Method for mapping, passed to map_by_cor(). Default = "median".
#' @param verbose Whether or not to display progress notifications.
#' 
#' @return a list object with two objects:  
#' \itemize{
#' \item map.df: A data.frame with the mapping results for each sample in test.dat to the reference
#' \item map.freq: A table with the frequency of mapping of each sample to each cluster across all iterations.
#' }
#' 
map_sampling <- function(train.dat, 
                         train.cl, 
                         test.dat, 
                         markers, 
                         markers.perc = 0.8, 
                         iter = 100, 
                         method = "median",
                         verbose = TRUE) {
  
  method <- match.arg(arg = method,
                      choices = c("mean", "median"))
  
  # Perform mapping iter times.
  map.result <- sapply(1:iter, 
                       function(i){
                         if(verbose) {
                           cat("\r", paste0("Running iteration ",i," of ",iter,".        "))
                           flush.console()
                         }
                         tmp.markers <- sample(markers, round(length(markers) * markers.perc))
                         map_by_cor(train.dat[tmp.markers,], 
                                    train.cl, 
                                    test.dat[tmp.markers,], 
                                    method = method)
                       }, simplify = F)
  
  # Extract predicted cluster assignments from each iteration
  map.cl <- sapply(map.result, 
                   function(x) {
                     x$pred.df$pred.cl
                   }
  )
  # Compute fraction of times each sample mapped to each cluster
  row.names(map.cl) <- colnames(test.dat)
  map <- as.data.frame(as.table(as.matrix(map.cl)))
  map.table <- table(map$Var1, map$Freq)
  map.freq <- unclass(map.table)
  
  # Find the most frequently mapped cluster for each sample
  max.freq <- apply(map.freq, 1, which.max)
  pred.cl <- colnames(map.freq)[max.freq]
  pred.cl <- setNames(pred.cl, row.names(map.freq))
  
  # Gather results
  map.df <- data.frame(pred.cl = pred.cl, 
                       prob = matrixStats::rowMaxs(map.freq) / iter)
  
  # output results
  out_list <- list(map.df = map.df,
                   map.freq = map.freq)
  
  return(out_list)
}

#' Run a single round of cross-validation of cluster mapping using a subset of marker genes
#' 
#' @param norm.dat a normalized data matrix for clustered cells
#' @param cl a cluster factor object for the cells in norm.dat
#' @param markers a character object with the marker genes to use for cross-validation
#' @param n.bin an integer indicating the number of bins to use
#' @param g.perc the fraction of genes to use for validation.
#' @param method Method for mapping. Must be either "median" (Default) or "mean".
#' @param verbose Whether or not to display progress notifications.
#' 
#' @return a named character object with the results of one round of cross-validation
#' 
#' @export
#' 
map_cv <- function(norm.dat, 
                   cl, 
                   markers, 
                   n.bin = 5,
                   g.perc = 1, 
                   method = "median",
                   verbose = TRUE) {
  
  method <- match.arg(arg = method,
                      choices = c("mean", "median"))
  
  bins <- tapply(names(cl), 
                 cl, 
                 function(x){
                   if(length(x) > n.bin){
                     tmp <- rep_len(1:n.bin, length(x))
                   }else{
                     tmp <- sample(1:n.bin, length(x))
                   }
                   setNames(tmp[sample(length(tmp))], x)
                 })
  
  names(bins) <- NULL
  bins <- unlist(bins)
  bins <- bins[names(cl)]
  pred.cl <- setNames(rep(NA, length(cl)), names(cl))
  
  for(i in 1:n.bin) {
    if(verbose) {
      cat("\r", paste0("Running bin ",i," of ",n.bin,".        "))
      flush.console()
    }
    train.cells <- names(cl)[bins != i]
    test.cells <- names(cl)[bins == i]
    select.markers <- sample(markers, round(length(markers) * g.perc))
    
    map.result <- map_by_cor(norm.dat[select.markers,], 
                             cl[train.cells], 
                             norm.dat[select.markers, 
                                      test.cells],
                             method=method)$pred.df
    
    pred.cl[test.cells] <- as.character(map.result[test.cells, "pred.cl"])
  }
  
  return(pred.cl)
}


###cluster annotation ref.cl.df must include "cluster_label" column

#' Compare two sets of cluster assignments for the same set of cells
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param ref.cl A cluster factor object for the reference clusters
#' @param ref.cl.df A data.frame with reference cluster annotations that includes a "cluster_label" column 
#' @param reorder Whether or not to reorder the overlaps by overlap size. Default = TRUE
#' @param rename a logical indicating whether or not to rename the results.
#'
#' @return a list with 5 objects:  
#' \itemize{
#' \item cl
#' \item cl.df
#' \item g A ggplot2 dot plot object for the comparison.
#' \item tb.df
#' \item cl.id.map
#' }
#' 
#' @export
#' 
compare_annotate <- function(cl, 
                             ref.cl, 
                             ref.cl.df=NULL, 
                             reorder = TRUE,
                             rename = TRUE) {

  if(!is.factor(cl)){
    cl <- setNames(factor(cl), names(cl))
  }
  
  if(!is.factor(ref.cl)){
    if(!is.null(ref.cl.df)){
      ref.cl <- setNames(factor(as.character(ref.cl),
                                levels = row.names(ref.cl.df)),
                         names(ref.cl))
    }
   else{
     ref.cl <- setNames(factor(ref.cl), names(ref.cl))
   }
  }
  
  common.cells <- intersect(names(cl),names(ref.cl))
  
  if(length(common.cells) == 0) {
    stop("No common names in cl and ref.cl for comparison.")
  }
  ###Find clusters not present in ref.cl
  tmp.cl <- droplevels(cl[common.cells])
  ref.cl <- droplevels(ref.cl[common.cells])
  absent.cl = setdiff(levels(cl), levels(tmp.cl))
  # compare predicted cluster member with the new clustering result 
  tb <- table(tmp.cl, ref.cl)
  cl.id.map <- NULL
  
  # Reorder clusters by size of overlap if reorder == TRUE
  if(reorder){
    tmp <- apply(tb, 1, which.max)
    cl_names <- names(cl)
    cl <- factor(as.character(cl), levels = c(row.names(tb)[order(tmp)],absent.cl))
    cl <- setNames(cl, cl_names)
    if(rename){
      cl.id.map <- data.frame(new = 1:length(levels(cl)),
                              old = levels(cl))
      levels(cl) <- 1:length(levels(cl))

    }
  }
  
  # Assign the best matching old cluster to each new cluster. 
  tb <- table(cl = cl[common.cells],ref.cl = ref.cl)
  max.ref.cl <- colnames(tb)[apply(tb, 1, which.max)]
  
  cl.df <- data.frame(ref.cl = max.ref.cl)
  cl.df <- cbind(cl.df, ref.cl.df[max.ref.cl,])
  
  cl_label <- rep("", nrow(cl.df))
  
  cl_split <- split(1:nrow(cl.df), cl.df$cluster_label)
  for(label in names(cl_split)){
    x <- cl_split[[label]]
    if(length(x) > 1){
      cl_label[x] <- paste(label, 1:length(x), sep = "_")
    } else {
      cl_label[x] <- label
    }
  }
  absent.cl = row.names(tb)[rowSums(tb)==0]
  cl.df$cluster_label <- cl_label
  row.names(cl.df) <- levels(cl)
  cl.size <- table(cl)
  cl.df$size <- cl.size[row.names(cl.df)]
  
  # Plot the mapping
  tb.df <- as.data.frame(tb)
  tb.df <- tb.df[tb.df$Freq > 0,]
  
  select.cells <- names(cl)

  cl.size <- table(cl)
  ref.cl.size <- table(ref.cl)
  
  tb.df$jaccard <- as.vector(tb.df$Freq / (cl.size[as.character(tb.df[,1])] + ref.cl.size[as.character(tb.df[,2])] - tb.df$Freq))

  # Compute Jaccard statistics for each pair of clusters
  #tb.df$jaccard <- 0
  #for(i in 1:nrow(tb.df)){
  #  n_ol <- length(union(common.cells[cl[common.cells] == as.character(tb.df[i,1])],
  #                       common.cells[ref.cl[common.cells] == as.character(tb.df[i,2])]))
  #  
  #  tb.df$jaccard[i] <- tb.df$Freq[i] / n_ol
  #}
  
  tb.df$ref.cl.label <- factor(ref.cl.df[as.character(tb.df$ref.cl),"cluster_label"], levels = ref.cl.df$cluster_label)
  
  g <- ggplot2::ggplot(tb.df, 
                       ggplot2::aes(x = cl, 
                                    y = ref.cl.label)) + 
    ggplot2::geom_point(ggplot2::aes(size = sqrt(Freq),
                                     color = jaccard)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.1,
                                                       hjust = 0.2, 
                                                       angle = 90,
                                                       size = 7),
                   axis.text.y = ggplot2::element_text(size = 6)) + 
    ggplot2::scale_color_gradient(low = "yellow", 
                                  high = "darkblue") + 
    ggplot2::scale_size(range = c(0, 3))
  
  out_list <- list(cl = cl,
                   cl.df = cl.df,
                   g = g,
                   tb.df = tb.df,
                   cl.id.map = cl.id.map, absent.cl=absent.cl)
  
  return(out_list)
}

#' Correct duplicated colors
#'
#' @param colorset a character vector of R or hex colors
#'
#' @return a character vector of hex colors with duplicated colors replaced
#' 
#' @export
#'
#' @examples
#' 
#' original_colors <- c("#00FF00","#00FF00","#FF0000","#00FF00")
#' 
#' new_colors <- adjust_color(original_colors)
#' 
adjust_color <- function(colorset) {
  
  duplicated_colors <- which(duplicated(colorset))
  
  while(length(duplicated_colors) > 0) {
    
    rgb <- col2rgb(colorset)
    
    for(x in duplicated_colors) {
      
      if(x < length(duplicated_colors)) {
        tmp <- round(rgb[,x - 1] * 0.8 + sample(20, 3) + rgb[,x + 1] * 0.2) 
      } else{
        tmp <- round(rgb[,x - 1] * 0.8 + sample(40, 3))
      }
      
      rgb[,x] <- tmp
    }
    
    rgb[rgb > 255] <- 255
    
    colorset <- rgb(rgb[1,],
                    rgb[2,],
                    rgb[3,], 
                    maxColorValue = 255)
    
    duplicated_colors <- which(duplicated(colorset))
  }
  
  return(colorset)
}

#' Generate an initial cl.df object based on cl
#' 
#' @param cl a cluster factor object
#' 
#' @return a data.frame with an id, color, and size for each cluster.
#' 
#' @export
#' 
get_cl_df <- function(cl) {
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", 
                                   "#7F0000"))
  
  cl.df <- data.frame(cluster_label = sort(unique(cl)))
  cl.size <- table(cl)
  cl.df$cluster_id <- 1:nrow(cl.df)
  cl.df$cluster_color <-jet.colors(nrow(cl.df))
  cl.df$size <- cl.size[row.names(cl.df)]
  row.names(cl.df) <- cl.df$cluster_label
  
  return(cl.df)

}

#' Compute correlations of clusters to a reference set, and get the best-correlated reference clusters
#' 
#' match_cl uses pearson correlation of cluster means.
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param dat a normalized data matrix for data to map to the reference
#' @param ref.cl A cluster factor object for the reference clusters
#' @param ref.cl.df A data.frame with reference cluster annotations that includes a "cluster_label" column 
#' @param ref.dat a reference normalized data matrix
#' @param rename a logical indicating whether or not to rename the results.
#' 
#' #' @return a list with 3 objects:  
#' \itemize{
#' \item cl
#' \item cl.df
#' \item cor
#' }
#' 
#' @export
#' 
match_cl <- function(cl, 
                     dat, 
                     ref.cl, 
                     ref.cl.df, 
                     ref.dat, 
                     rename = TRUE) {
  
  cl.means <- get_cl_means(dat, cl)
  ref.cl.means <- get_cl_means(ref.dat, ref.cl)
  mat <- cor(cl.means, ref.cl.means)
  tmp <- apply(mat, 1, which.max)
  max.ref.cl <- colnames(mat)[tmp]
  
  cl_names <- names(cl)
  cl <- factor(as.character(cl), 
               levels = c(row.names(mat)[order(tmp)]))
  cl <- setNames(cl, cl_names)
  
  cl.df <- data.frame(ref.cl = max.ref.cl)
  cl.df <- cbind(cl.df, ref.cl.df[max.ref.cl,])
  row.names(cl.df) <- row.names(mat)
  cl.df <- cl.df[levels(cl),]
  mat <- mat[levels(cl),]
  
  if(rename){
    levels(cl) <- 1:length(levels(cl))
    row.names(cl.df) <- row.names(mat) <- levels(cl)
  }
  
  return(list(cl = cl, 
              cl.df = cl.df, 
              cor = mat))
}

#' Use DEGenes to identify clusters with few differentially expressed genes
#'
#' To use this, you'll need to perform pairwise gene expression tests with de_score(), and 
#' a set of clusters with high confidence (cl.good).
#' 
#' Low-quality clusters will be from all clusters in cl.df that are not in cl.good.
#'
#' @param cl.df A data.frame with cluster annotations that includes a "cluster_label" column 
#' @param cl.good a cluster factor object for cells from high-quality clusters.
#' @param de.score.mat an optional matrix of pairwise degene scores generated by get_de_matrix().
#' @param de.genes differential gene expression results like those generated by de_score().
#'
#' @return a data.frame with the closest "good" cluster for each non-good/low cluster.
#' 
#' @export
#'
find_low_quality_cl <- function(cl.df, 
                                cl.good, 
                                de.score.mat = NULL, 
                                de.genes) {
  
  if(is.null(de.score.mat)){
    de.score.mat <- get_de_matrix(de.genes, 
                                  directed = TRUE, 
                                  field = "num")
  }
  
  diag(de.score.mat) <- max(de.score.mat)
  
  tmp.mat <- de.score.mat[, cl.good]
  low.pair <- data.frame(cl.low = row.names(tmp.mat), 
                         cl.good = colnames(tmp.mat)[apply(tmp.mat, 1, which.min)],
                        stringsAsFactors = FALSE)
  
  low.pair$low.gene.counts <- cl.df[low.pair[,1], "gene.counts"]
  low.pair$high.gene.counts <- cl.df[low.pair[,2], "gene.counts"]
  
  low.pair$low.size <- cl.df[low.pair[,1], "size"]
  low.pair$high.size <- cl.df[low.pair[,2], "size"]
  
  low.pair$up.genes <- get_pair_matrix(de.score.mat, 
                                       low.pair$cl.low, 
                                       low.pair$cl.good)
  
  low.pair$cl.low.label <- cl.df[as.character(low.pair$cl.low), "cluster_label"]
  low.pair$cl.good.label <- cl.df[as.character(low.pair$cl.good), "cluster_label"]
  
  row.names(low.pair) <- low.pair$cl.low
  
  return(low.pair)
}


#' Plot low-quality clusters
#'
#' @param norm.data normalized data matrix for clustered cells.
#' @param cl A cluster factor object for all cells.
#' @param low.df a data.frame with results from find_low_quality_cl().
#' @param nn.df 
#' @param de.genes differential gene expression results like those generated by de_score().
#' @param all.col 
#'
#' @return
#'
plot_low_qc <- function(norm.dat, 
                        cl, 
                        low.df, 
                        nn.df, 
                        de.genes, 
                        all.col) {
  
  for(i in 1:nrow(low.df)) {
    x <- low.df[i, "cl.low"]
    y <- low.df[i, "cl.good"]
    i <- nn.df[x, "nn.cl"]
    j <- nn.df[y, "nn.cl"]
    
    tmp.cl <- droplevels(cl[cl %in% c(x, y, i, j)])
    tmp.cl <- tmp.cl[names(tmp.cl) %in% colnames(norm.dat)]
    
    tmp <- display_cl(tmp.cl, 
                      norm.dat, 
                      prefix = paste(levels(tmp.cl), 
                                     collapse = "_"), 
                      col = all.col, 
                      max.cl.size = 100, 
                      de.genes = de.genes)
  }
  
}
