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


predict_annotate_cor <- function(cl, norm.dat, ref.markers, ref.cl, ref.cl.df, ref.norm.dat, method="median", reorder=FALSE)
  {
    tmp = map_by_cor(ref.norm.dat[ref.markers,], ref.cl, norm.dat[ref.markers,names(cl)],method=method)
    pred.cl= setNames(factor(as.character(tmp$pred.df$pred.cl),levels=row.names(ref.cl.df)), row.names(tmp$pred.df))
    compare_annotate(cl, pred.cl, ref.cl.df, reorder=reorder)
  }

map_sampling <- function(train.dat, train.cl, test.dat, markers, markers.perc=0.8, iter=100, method="median")
  {
    map.result = sapply(1:iter, function(i){
      tmp.markers=sample(markers, round(length(markers)*markers.perc))
      map_by_cor(train.dat[tmp.markers,], train.cl, test.dat[tmp.markers,], method=method)
    },simplify=F)
    map.cl = sapply(map.result, function(x)x$pred.df$pred.cl)

    row.names(map.cl) = colnames(test.dat)
    map=  as.data.frame(as.table(as.matrix(map.cl)))
    map.freq <- table(map$Var1, map$Freq)
    map.df = data.frame(pred.cl=setNames(colnames(map.freq)[apply(map.freq, 1, which.max)],row.names(map.freq)), prob=rowMaxs(map.freq)/iter)
    return(list(map.df=map.df, map.freq=map.freq))
  }

map_cv <- function(norm.dat, cl, markers, n.bin=5,g.perc=1){
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  pred.cl = setNames(rep(NA, length(cl)), names(cl))
  for(i in 1:n.bin){
    print(i)
    train.cells = names(cl)[bins!=i]
    test.cells =names(cl)[bins==i]
    select.markers=sample(markers, round(length(markers)*g.perc))
    map.result <- map_by_cor(norm.dat[select.markers,], cl[train.cells], norm.dat[select.markers, test.cells])$pred.df
    pred.cl[test.cells] = as.character(map.result[test.cells, "pred.cl"])
  }
  return(pred.cl)
}


###cluster annotation ref.cl.df must include "cluster_label" column
compare_annotate<- function(cl, ref.cl, ref.cl.df, reorder=TRUE)
{
  common.cells=intersect(names(cl),names(ref.cl))
  ###compare predicted cluster member with the new clustering result 
  tb = table(cl[common.cells], ref.cl[common.cells])
  cl.id.map=NULL
  ###Reorder clusters
  if(reorder){
    tmp = apply(tb, 1, which.max)
    cl = setNames(factor(as.character(cl), levels=row.names(tb)[order(tmp)]), names(cl))
    cl.id.map <- data.frame(new=1:length(levels(cl)),old=levels(cl))
    levels(cl)=1:length(levels(cl))
  }
  ###Assign the best matching old cluster to each new cluster. 
  tb=table(cl=cl[common.cells],ref.cl=ref.cl[common.cells])
  tmp = colnames(tb)[apply(tb, 1, which.max)]
  
  cl.df = data.frame(ref.cl=tmp)
  cl.df = cbind(cl.df, ref.cl.df[tmp,])
  
  cl_label= rep("", nrow(cl.df))
  tmp = split(1:nrow(cl.df), cl.df$cluster_label)
  for(label in names(tmp)){
    x = tmp[[label]]
    if(length(x)>1){
      cl_label[x] <- paste(label, 1:length(x),sep="_")
    }
    else{
      cl_label[x]<-label
    }
  }
  cl.df$cluster_label = cl_label
  row.names(cl.df) = levels(cl)
  cl.size = table(cl)
  cl.df$size = cl.size[row.names(cl.df)]
  ###plot the mapping
  tb.df = as.data.frame(tb)
  tb.df = tb.df[tb.df$Freq > 0,]
  library(ggplot2)
  select.cells=names(cl)
  tb.df$jaccard=0
  for(i in 1:nrow(tb.df)){
    tb.df[i,"jaccard"] = tb.df[i,"Freq"]/length(union(names(cl)[cl==as.character(tb.df[i,1])],names(ref.cl)[ref.cl==as.character(tb.df[i,2])]))
  }
  tb.df$ref.cl.label = factor(ref.cl.df[as.character(tb.df$ref.cl),"cluster_label"], levels=ref.cl.df$cluster_label)
  g= ggplot(tb.df, aes(x=cl, y=ref.cl.label)) + geom_point(aes(size=sqrt(Freq),color=jaccard))
  g = g+ theme(axis.text.x=element_text(vjust=0.1,hjust = 0.2, angle=90,size=7),axis.text.y=element_text(size=6)) + scale_color_gradient(low="yellow",high="darkblue")
  g= g+scale_size(range=c(0,3))
  return(list(cl=cl, cl.df=cl.df,g = g,tb.df=tb.df,cl.id.map=cl.id.map))
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

