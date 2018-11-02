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
map_by_cor <- function(train.dat, 
                       train.cl, 
                       test.dat,
                       method = "median")
{
  # Get medians or means for each cluster
  if(method == "median"){
    library(matrixStats)
    cl.meds <- tapply(names(train.cl), 
                      train.cl, 
                      function(x) {
                        train.mat <- train.dat[, x, drop = F]
                        train.mat <- as.matrix(train.mat)
                        rowMedians(train.mat)
                      }
    )
    
    cl.dat <- do.call("cbind", cl.meds)
  } else {
    cl.dat <- get_cl_means(train.dat, train.cl)
  }
  row.names(cl.dat) <- row.names(train.dat)
  
  # Perform correlations
  test.cl.cor <- cor(as.matrix(test.dat), cl.dat)
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
#' 
#' @return a list object with two objects:  
#' \itemize{
#' \item map.df: A data.frame with the mapping results for each sample in map.dat to the reference
#' \item cl.map.df: A data.frame with cluster-level frequency of mapping for each cluster in map.cl to ref.cl
#' }
map_cl_summary <- function(ref.dat, 
                           ref.cl, 
                           map.dat, 
                           map.cl)
{
  # Map the training set to the reference
  map.result <- map_by_cor(ref.dat, ref.cl, map.dat)
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


predict_annotate_cor <- function(cl, 
                                 norm.dat, 
                                 ref.markers, 
                                 ref.cl, 
                                 ref.cl.df, 
                                 ref.norm.dat, 
                                 method = "median", 
                                 reorder = FALSE)
{
  map_results <- map_by_cor(ref.norm.dat[ref.markers,], 
                            ref.cl, 
                            norm.dat[ref.markers, names(cl)],
                            method = method)
  
  pred.cl <- setNames(factor(as.character(map_result$pred.df$pred.cl), levels = row.names(ref.cl.df)), row.names(map_results$pred.df))
  
  results <- compare_annotate(cl, 
                              pred.cl, 
                              ref.cl.df, 
                              reorder = reorder)
  
  return(results)
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
map_sampling <- function(train.dat, 
                         train.cl, 
                         test.dat, 
                         markers, 
                         markers.perc = 0.8, 
                         iter = 100, 
                         method = "median",
                         verbose = TRUE)
{
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
  map.freq <- table(map$Var1, map$Freq)
  
  # Find the most frequently mapped cluster for each sample
  max.freq <- apply(map.freq, 1, which.max)
  pred.cl <- colnames(map.freq)[max.freq]
  pred.cl <- setNames(pred.cl, row.names(map.freq))
  
  # Gather results
  map.df <- data.frame(pred.cl = pred.cl, 
                       prob = rowMaxs(map.freq) / iter)
  
  # output results
  out_list <- list(map.df = map.df,
                   map.freq = map.freq)
  
  return(out_list)
}

map_cv <- function(norm.dat, 
                   cl, 
                   markers, 
                   n.bin = 5,
                   g.perc = 1) {
  
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
  names(bins)=NULL
  bins <- unlist(bins)
  bins <- bins[names(cl)]
  pred.cl <- setNames(rep(NA, length(cl)), names(cl))
  
  for(i in 1:n.bin) {
    print(i)
    train.cells <- names(cl)[bins != i]
    test.cells <- names(cl)[bins == i]
    select.markers <- sample(markers, round(length(markers) * g.perc))
    
    map.result <- map_by_cor(norm.dat[select.markers,], 
                             cl[train.cells], 
                             norm.dat[select.markers, 
                                      test.cells])$pred.df
    
    pred.cl[test.cells] <- as.character(map.result[test.cells, "pred.cl"])
  }
  
  return(pred.cl)
}


###cluster annotation ref.cl.df must include "cluster_label" column

#' Compare two sets of cluster assignments
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param ref.cl A cluster factor object for the reference clusters
#' @param ref.cl.df A data.frame with reference cluster annotations that includes a "cluster_label" column 
#' @param reorder Whether or not to reorder the overlaps by overlap size. Default = TRUE
#'
#' @return a list with 5 objects:  
#' \itemize{
#' \item cl
#' \item cl.df
#' \item g A ggplot2 dot plot object for the comparison.
#' \item tb.df
#' \item cl.id.map
#' }
compare_annotate <- function(cl, 
                             ref.cl, 
                             ref.cl.df, 
                             reorder = TRUE)
{
  library(ggplot2)
  
  common.cells <- intersect(names(cl),names(ref.cl))
  
  # compare predicted cluster member with the new clustering result 
  tb <- table(cl[common.cells], ref.cl[common.cells])
  cl.id.map <- NULL
  
  # Reorder clusters by size of overlap if reorder == TRUE
  if(reorder){
    tmp <- apply(tb, 1, which.max)
    cl_names <- names(cl)
    cl <- factor(as.character(cl), levels = row.names(tb)[order(tmp)])
    cl <- setNames(cl, cl_names)
    cl.id.map <- data.frame(new = 1:length(levels(cl)),
                            old = levels(cl))
    levels(cl) <- 1:length(levels(cl))
  }
  
  # Assign the best matching old cluster to each new cluster. 
  tb <- table(cl = cl[common.cells],
              ref.cl = ref.cl[common.cells])
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
  
  cl.df$cluster_label <- cl_label
  row.names(cl.df) <- levels(cl)
  cl.size <- table(cl)
  cl.df$size <- cl.size[row.names(cl.df)]
  
  # Plot the mapping
  tb.df <- as.data.frame(tb)
  tb.df <- tb.df[tb.df$Freq > 0,]
  
  select.cells <- names(cl)
  
  # Compute Jaccard statistics for each pair of clusters
  tb.df$jaccard <- 0
  for(i in 1:nrow(tb.df)){
    n_ol <- length(union(names(cl)[cl == as.character(tb.df[i,1])],
                         names(ref.cl)[ref.cl == as.character(tb.df[i,2])]))
    
    tb.df$jaccard[i] <- tb.df$Freq[i] / n_ol
  }
  
  tb.df$ref.cl.label <- factor(ref.cl.df[as.character(tb.df$ref.cl),"cluster_label"], levels = ref.cl.df$cluster_label)
                               
  
  g <- ggplot(tb.df, 
              aes(x = cl, 
                  y = ref.cl.label)) + 
    geom_point(aes(size = sqrt(Freq),
                   color = jaccard)) + 
    theme(axis.text.x = element_text(vjust = 0.1,
                                     hjust = 0.2, 
                                     angle = 90,
                                     size = 7),
          axis.text.y = element_text(size = 6)) + 
    scale_color_gradient(low = "yellow", high = "darkblue") + 
    scale_size(range=c(0,3))
  
  out_list <- list(cl = cl,
                   cl.df = cl.df,
                   g = g,
                   tb.df = tb.df,
                   cl.id.map = cl.id.map)
  
  return(out_list)
}

get_color <- function(color.mapping)
{
  while(1){
    tmp.color = which(duplicated(color.mapping))
    if(length(tmp.color)==0){
      break
    }
    rgb = col2rgb(color.mapping)
    for(x in tmp.color){
      if(x < length(tmp.color)){
        tmp= round(rgb[,x-1] *0.8 + sample(20, 3) + rgb[,x+1] *0.2) 
      }
      else{
        tmp= round(rgb[,x-1] *0.8 +  sample(40, 3))
      }
      rgb[,x] = tmp
    }
    rgb[rgb > 255]=255
    color.mapping = rgb(rgb[1,],rgb[2,],rgb[3,], maxColorValue=255)
  }
  return(color.mapping)
}

