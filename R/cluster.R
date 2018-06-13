library(Matrix)
library(matrixStats)
library(Rphenograph)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))

jaccard <- function(m) {
  require(Matrix)
  ## common values:
  A =  m %*% t(m)
  ## indexes for non-zero common values
  im = Matrix::which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = Matrix::rowSums(m)  
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
    )  
  return( J )
}


pass_louvain <- function(mod.sc, adj.mat)
{
  p = mean(Matrix::colSums(adj.mat>0)-1)/ nrow(adj.mat)
  n = ncol(adj.mat)
  rand.mod1 <- 0.97 * sqrt((1 - p)/(p*n))
  rand.mod2 <- (1 - 2 / sqrt(n)) * (2 / (p * n))^(2/3)
  rand.mod.max <- max(rand.mod1, rand.mod2, na.rm=TRUE)
  cat("Modularity:",mod.sc, "threshold:",rand.mod.max, "\n")
  return(mod.sc > rand.mod.max)
}

###rows are cells, columns are feathers
jaccard_louvain.FNN <- function(dat, k=10, knn.matrix=NULL)
  {
    suppressMessages(library(igraph))
    if(is.null(knn.matrix)){
      knn.result= FNN::get.knn(dat, k)
      knn.matrix = knn.result[[1]]
      suppressMessages(library(FNN))
    }
    p = as.vector(t(knn.matrix))
    edge = cbind(rep(1:nrow(knn.matrix),rep(k, nrow(knn.matrix))),p)
    edge.unique = cbind(rowMins(edge), rowMaxs(edge))
    edge.unique=unique(edge.unique)
    knn.gr = igraph::graph(t(edge))
    knn.matrix = igraph::get.adjacency(knn.gr)    
    jaccard.adj  = jaccard(knn.matrix)
    jaccard.gr = igraph::graph.adjacency(jaccard.adj, mode="undirected",weighted=TRUE)
    louvain.result= igraph::cluster_louvain(jaccard.gr)
    mod.sc = igraph::modularity(louvain.result)
    if(pass_louvain(mod.sc,jaccard.adj)){
      cl  = setNames(louvain.result$membership,row.names(dat))
      return(list(cl=cl, result=louvain.result))
    }
    else{
      return(NULL)
    }
  }

jaccard_louvain <- function(dat, k=10)
{
  suppressMessages(library(Rphenograph))
  rpheno <- Rphenograph(dat, k = k)
  cl  = setNames(rpheno[[2]]$membership, row.names(dat)[as.integer(rpheno[[2]]$names)])
  return(list(cl = cl, result=rpheno))
}

#' One round of clustering in the iteractive clustering pipeline 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1).
#' @param select.cells The cells to be clustered. Default: columns of norm.dat
#' @param counts Raw gene counts. Default NULL, inferred from norm.dat.
#' @param method Clustering method. It can be "louvain", "hclust" and "kmeans". Default "louvain"
#' @param vg.padj.th High variance gene adjusted pvalue cut off. Default 0.5.
#' @param dim.method Dimension reduction techniques. Current options include "pca" and "WGCNA". Default "pca"
#' @param max.dim The number of top dimensions retained. Default 20. Since clustering is performed iteratively, not all relevant dimensions need to be captured in one iterations. 
#' @param rm.eigen The reduced dimensions that need to be masked and removed. Default NULL.  
#' @param rm.th The cutoff for correlation between reduced dimensions and rm.eigen. Reduced dimensions with correlatin with any rm.eigen vectors are not used for clustering. Default 0.7
#' @param de.param The differential gene expression threshold. See de_param() function for details. 
#' @param min.genes The minimal number of high variance and differentially expressed genes genes. Default 5. 
#' @param type Can either be "undirectional" or "directional". If "undirectional", the differential gene threshold de.param is applied to combined up-regulated and down-regulated genes, if "directional", then the differential gene threshold is applied to both up-regulated and down-regulated genes. 
#' @param maxGenes Only used when dim.method=="WGCNA". The maximum number of genes to calculate gene modules. 
#' @param sampleSize The number of sampled cells to compute reduced dimensions.
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.    
#' @param prefix Used to keep track of intermediate results in "verbose" mode. Default NULL.
#' @param verbose Default FALSE
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#' @export
#'
#' @examples
onestep_clust <- function(norm.dat, select.cells=colnames(norm.dat), counts=NULL, method=c("louvain","ward", "kmeans"), vg.padj.th=0.5, dim.method=c("pca","WGCNA"), max.dim=20, rm.eigen=NULL, rm.th=0.7, de.param = de_param(),min.genes=5, type=c("undirectional", "directional"), maxGenes=3000,sampleSize=4000,max.cl.size=300, prefix=NULL, verbose=FALSE)
                          
  {
    method=method[1]
    dim.method=dim.method[1]
    type=type[1]
    
    if(length(select.cells)>sampleSize){
      sampled.cells = sample(select.cells, pmin(length(select.cells),sampleSize))
    }
    else{
      sampled.cells = select.cells
    }
    ###Find high variance genes
    if(is.matrix(norm.dat)){
      select.genes = row.names(norm.dat)[which(rowSums(norm.dat[,select.cells] > de.param$low.th) >= de.param$min.cells)]
    }
    else{
      select.genes = row.names(norm.dat)[which(Matrix::rowSums(norm.dat[,select.cells] > de.param$low.th) >= de.param$min.cells)]
    }
    ###Find high variance genes.
    if(is.null(counts)){
      counts = 2^(norm.dat[select.genes, sampled.cells])-1
    }
    plot.fig=NULL
    if(verbose & !is.null(prefix)){
      plot.fig=paste0(prefix,".vg.pdf")
    }
    vg = findVG(as.matrix(counts[select.genes,sampled.cells]),plot.fig=plot.fig)
    if(dim.method=="auto"){
      if(length(select.cells)> 1000){
        dim.method="pca"
      }
      else{
        dim.method="WGCNA"
      }
    }
    if(dim.method=="WGCNA"){
      ###Ignore vg.padj.th for WGCNA, choose top "maxGgenes" for analysis
      select.genes = row.names(vg)[which(vg$loess.padj < 1)]
      select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
      rd.dat = rd_WGCNA(norm.dat, select.genes=select.genes, select.cells=select.cells, sampled.cells=sampled.cells, de.param=de.param, max.mod=max.dim, max.cl.size=max.cl.size)
    }
    else{
       ###If most genes are differentially expressed, then use absolute dispersion value
      select.genes = row.names(vg)[which(vg$loess.padj < vg.padj.th | vg$dispersion >3)]
      select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
      if(verbose){
        cat("Num high variance genes:",length(select.genes),"\n")
      }
      if(length(select.genes)< min.genes){
        return(NULL)
      }
      rd.dat = rd_PCA(norm.dat,select.genes, select.cells, sampled.cells=sampled.cells, max.pca = max.dim)
    }
    if(is.null(rd.dat)||ncol(rd.dat)==0){
      return(NULL)
    }
    if(!is.null(rm.eigen)){
      rm.cor=cor(rd.dat, rm.eigen[row.names(rd.dat),])
      rm.cor[is.na(rm.cor)]=0
      rm.score = rowMaxs(abs(rm.cor))
      select = rm.score < rm.th
      if(sum(!select)>0){
        print("Remove dimension:")
        print(rm.score[!select])
      }
      if(sum(select)==0){
        return(NULL)
      }
      rd.dat = rd.dat[,select,drop=F]
    }
    print(method)
    max.cl = ncol(rd.dat)*2 + 1
    if(method=="louvain"){
      tmp = jaccard_louvain(rd.dat, 15)
      if(is.null(tmp)){
        return(NULL)
      }
      cl = tmp$cl
      if(length(unique(cl))>max.cl){
        tmp.means =do.call("cbind",tapply(names(cl),cl, function(x){
          colMeans(rd.dat[x,,drop=F])
        },simplify=F))
        tmp.hc = hclust(dist(t(tmp.means)), method="average")
        tmp.cl= cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
        cl = setNames(tmp.cl[as.character(cl)], names(cl))
      }
    }
    else if(method=="ward"){
      hc = hclust(dist(rd.dat),method="ward.D")
      print("Cluster cells")
      cl = cutree(hc, max.cl)
    }
    else if(method=="kmeans"){
      cl = kmeans(rd.dat, max.cl)$cluster
    }
    else{
      stop(paste("Unknown clustering method", method))
    }
    #print(table(cl))
    merge.result=merge_cl(norm.dat, cl=cl, rd.dat=rd.dat, type=type, de.param=de.param, max.cl.size=max.cl.size)
    gc()
    if(is.null(merge.result))return(NULL)
    sc = merge.result$sc
    #print(sc)
    cl = merge.result$cl
    print(table(cl))
    if(verbose){
      save(cl, file=paste0(prefix, ".cl.rda"))
    }
    if(length(unique(cl))>1){
      de.genes = merge.result$de.genes
      markers= merge.result$markers
      cl.dat = get_cl_means(norm.dat[markers,], cl[sample_cells(cl, max.cl.size)])
      cl.hc = hclust(dist(t(cl.dat)),method="average")
      cl = setNames(factor(as.character(cl), levels= colnames(cl.dat)[cl.hc$order]), names(cl))
      if(verbose & !is.null(prefix)){
        display_cl(norm.dat, cl, prefix=prefix, markers=markers, max.cl.size=max.cl.size)
      }
      levels(cl) = 1:length(levels(cl))
      result=list(cl=cl, markers=markers)
      return(result)
    }
    return(NULL)
  }


#' Iterative clustering algorithm for single cell RNAseq dataset
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param select.cells The cells to be clustered
#' @param prefix The character string to indicate current iteration.
#' @param split.size The minimal cluster size for further splitting
#' @param result The current clustering result as basis for further splitting.
#' @param method Clustering method. It can be "auto", "louvain", "hclust"
#' @param ... Other parameters passed to method "onestep_clust"
#'
#' @return Clustering result is returned as a list with two elements: 
#'         cl: cluster membership for each cell
#'         markers: top markers that seperate clusters     
#' @export
#'
#' @examples clust.result = iter_clust(norm.dat)
#'           clust.result = iter_clust(norm.dat, de.param = de_param(q1.th=0.5, de.score.th=100))
iter_clust <- function(norm.dat, select.cells=colnames(norm.dat),prefix=NULL, split.size = 10, result=NULL,method="auto",...)
  {
    print(prefix)
    if(method=="auto"){
      if(length(select.cells)>3000){
        select.method="louvain"
      }
      else{
        select.method="ward"
      }
    }
    else{
      select.method=method
    }
    if(length(select.cells) <= 3000){
      if(!is.matrix(norm.dat)){
        norm.dat = as.matrix(norm.dat[,select.cells])
      }
    }
    if(is.null(result)){        
      result=onestep_clust(norm.dat, select.cells=select.cells, prefix=prefix,method=select.method,...)
      gc()
    }
    if(!is.null(result)){
      #save(result, file=paste0(prefix,".rda"))
      cl = result$cl[select.cells]
      gene.mod = result$gene.mod
      markers=result$markers
      cl = setNames(as.integer(cl),names(cl))
      new.cl =cl
      cl.size = table(cl)
      to.split = names(cl.size)[cl.size >=split.size]
      if(length(to.split)>0){
        n.cl = 1
        for(x in sort(unique(cl))){
          tmp.cells = names(cl)[cl==x]
          if(!x %in% to.split){
            new.cl[tmp.cells]=n.cl
          }
          else{
            tmp.prefix = paste(prefix, x, sep=".")
            tmp.result=iter_clust(norm.dat=norm.dat, select.cells=tmp.cells, prefix=tmp.prefix,split.size=split.size,method= method,...)
            gc()
            if(is.null(tmp.result)){
              new.cl[tmp.cells]=n.cl
            }
            else{
              tmp.cl = tmp.result$cl
              if(length(unique(tmp.cl)>1)){
                cat("Expand",tmp.prefix, "\n")
                print(table(tmp.cl))
                new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
                markers=union(markers, tmp.result$markers)
              }
            }
          }
          n.cl = max(new.cl)+1
        }
        cl = new.cl
      }
      result=list(cl=cl, markers=markers)
      return(result)
    }
    return(NULL)
  }


#' Reorder cluster based on hiearchical clustering of clusters based on average cluster values for the input data matrix 
#'
#' @param cl A vector of cluster membership with cells as names, and cluster id as values. 
#' @param dat The data matrix with cells as columns. 
#'
#' @return Reorder cluster membership vector. The cluster id start from 1 to the number of clusters.
#' @export 
#'
#' @examples 
reorder_cl <- function(cl, dat)
{
  cl.means = get_cl_means(dat,cl)
  cl.hc = hclust(as.dist(1- cor(cl.means)),method="average")
  cl = setNames(factor(as.character(cl), levels=cl.hc$labels[cl.hc$order]),names(cl))
  cl = setNames(as.integer(cl),names(cl))
}

#' Merge clusters based on pairwise differential expressed genes. 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param cl A vector of cluster membership with cells as names, and cluster id as values. 
#' @param rd.dat Reduced dimensions for cells. Used to determine which clusters are close to each other. Clusters are merged among nearest neighbors first. 
#' @param de.param The DE gene criteria. See de_param for details. 
#' @param type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately. 
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.   
#' @param de.method Use limma by default. We are still testing "chisq" mode. 
#' @return A list with cl (cluster membership), de.genes (differentially expressed genes), sc (cluster pairwise de.score), markers (top cluster pairwise markers)
#' @export
#'
#' @examples
merge_cl<- function(norm.dat, cl, rd.dat, de.param = de_param(), type=c("undirectional","directional"), max.cl.size=300,de.method="limma",de.genes=NULL, return.markers=TRUE, verbose=0)
  {
    cl = setNames(as.integer(as.character(cl)), names(cl))
    de.df=list()
    select.cells = names(cl)
    pairs=NULL
    if(!is.null(de.genes)){
      pairs=do.call("rbind",strsplit(names(de.genes), "_"))
      row.names(pairs)=names(de.genes)
    }
    if(is.matrix(norm.dat)){
      cell.gene.counts= colSums(norm.dat[,select.cells]>0)
    }
    else{
      cell.gene.counts= Matrix::colSums(norm.dat[,select.cells]>0)
    }
    cell.weights = cell.gene.counts - min(cell.gene.counts)+200
    while(length(unique(cl)) > 1){
      cl.size = table(cl)
      ##Down sample cells for efficiency   
      tmp.cells = sample_cells(cl, weights=cell.weights, max.cl.size=max.cl.size)
      tmp.dat = as.matrix(norm.dat[,tmp.cells])
      ###Merge small clusters with the closest neighbors first.
      cl.rd = t(get_cl_means(t(rd.dat),cl))
      while(TRUE){
        cl.size = table(cl)
        ##Compute cluster similary on reduced dimension
        if(ncol(cl.rd)>2){
          cl.sim = cor(t(cl.rd))
        }
        else{
          cl.diff=as.matrix(dist(cl.rd/as.vector(cl.size[row.names(cl.rd)])))
          cl.sim = 1 - cl.diff/max(cl.diff)
        }
        cl.small = names(cl.size)[cl.size < de.param$min.cells]
        ###Merge small clusters with its closest neighbors.
        if(length(cl.small)>0){
          tmp=as.data.frame(as.table(cl.sim[cl.small,,drop=F]))
          tmp[,1]=as.integer(as.character(tmp[,1]))
          tmp[,2]=as.integer(as.character(tmp[,2]))
          tmp = tmp[tmp[,1]!=tmp[,2],,drop=F]
          closest.pair = which.max(tmp$Freq)
          x = tmp[closest.pair,1]
          y=  tmp[closest.pair,2]
          if(verbose > 0){
            cat("Merge: ", x,y, "sim:", tmp[closest.pair,3],"\n")
          }
          cl[cl==x]= y
          cl.rd[as.character(y),]= cl.rd[as.character(y),] + cl.rd[as.character(x),]
          cl.rd = cl.rd[row.names(cl.rd)!=x,,drop=F]
        }		
        else{
          break
        }	
      }
      if(length(cl.size)==1){
        return(NULL)
      }
      
      ####if the number of clusters is fewer than 10, compute all pairs. If no clusters are significant, quit
      if(length(cl.size) < 10 & length(de.genes)==0){
        de.genes = de_score(tmp.dat, cl=cl[tmp.cells], de.param= de.param, method=de.method, de.genes=de.genes)
        sc = sapply(de.genes, function(x){
          if(length(x)>0){x$score}
          else{0}
        })
        merge.pairs=do.call("rbind",strsplit(names(de.genes), "_"))
        row.names(merge.pairs)=names(de.genes)
        pairs = merge.pairs
      }
      else{
        ####Choose top 2 nearest neighbor based on correlation matrix.
        if(nrow(cl.rd)>2){
          k = 3
          knn.matrix=t(sapply(1:nrow(cl.sim), function(i){colnames(cl.sim)[-i][head(order(cl.sim[i,-i],decreasing=T), k)]}))
          row.names(knn.matrix)=row.names(cl.sim)
          merge.pairs = do.call("rbind",apply(knn.matrix, 2,function(x)data.frame(c1=row.names(knn.matrix),c2=x)))
          merge.pairs[,1] = as.integer(as.character(merge.pairs[,1]))
          merge.pairs[,2] = as.integer(as.character(merge.pairs[,2]))
          tmp1 =pmin(merge.pairs[,1],merge.pairs[,2])
          tmp2 =pmax(merge.pairs[,1],merge.pairs[,2])
          merge.pairs = cbind(tmp1,tmp2)        
        }
        else{
          merge.pairs = matrix(as.integer(row.names(cl.rd)), nrow=1)
        }
        
        row.names(merge.pairs) = paste(merge.pairs[,1],merge.pairs[,2],sep="_")
        merge.pairs= merge.pairs[!duplicated(row.names(merge.pairs)),,drop=F]
        ###Determine the de score for these pairs
        if(nrow(merge.pairs)==0){
          break
        }
        #####Check pairs already known but not yet merged yet.
        new.pairs = setdiff(row.names(merge.pairs),names(de.genes))
        pairs = rbind(pairs, merge.pairs[new.pairs,,drop=F])
        tmp.de.genes =de_score_pairs(tmp.dat, cl=cl[tmp.cells], pairs=merge.pairs[new.pairs,,drop=F], de.param= de.param, method=de.method)$de.genes
        de.genes[names(tmp.de.genes)] = tmp.de.genes
        sc = sapply(de.genes[row.names(merge.pairs)], function(x){
          if(length(x)>0){x$score}
          else{0}
        })
      }
      sc = sort(sc)
      #print(head(sc,10))      
      to.merge = sapply(names(sc), function(p){
        x = de.genes[[p]]
        if(length(x)==0){
          to.merge = TRUE
        }
        else if(type=="undirectional"){
          to.merge=x$score < de.param$de.score.th    
        }
        else{
          to.merge=x$up.score < de.param$de.score.th | x$down.score < de.param$de.score.th 
        }
        to.merge
      })
      if(sum(to.merge)==0){
        break
      }
      sc = sc[to.merge]
      to.merge= merge.pairs[names(sc),,drop=FALSE]
      merged =c()
      ###The first pair in to.merge always merge. For the remaining pairs, if both clusters have already enough cells,
      ###or independent of previus merging, then they can be directly merged as well, without re-assessing DE genes. 
      for(i in 1:nrow(to.merge)){
        p = to.merge[i,]
        if(i == 1 | sc[i] < de.param$de.score.th /2  & length(intersect(p, merged))==0){
          if(verbose > 0){
            cat("Merge ",p[1], p[2], sc[i],"\n")
          }
          cl[cl==p[2]] = p[1]
          rm.pairs = row.names(pairs)[pairs[,1]%in% p | pairs[,2]%in% p]
          de.genes = de.genes[setdiff(names(de.genes),rm.pairs)]
        }
        merged = c(merged, p)
      }
      pairs = pairs[names(de.genes),,drop=F]
    }
    if(length(unique(cl))<2){
      return(NULL)
    }
    if(verbose > 0){
      print(table(cl))
    }
    markers = NULL
    if(return.markers){
      de.genes = de_score(as.matrix(norm.dat[,tmp.cells]), cl[tmp.cells], de.genes=de.genes, de.param=de.param)
      markers = select_markers(norm.dat, cl, de.genes=de.genes, n.markers=50)$markers
    }
    sc = sapply(de.genes, function(x){
      if(length(x)>0){x$score}
      else{0}
    })
    return(list(cl=cl, de.genes=de.genes,sc=sc, markers=markers))
  }


