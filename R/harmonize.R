###comb.dat include the following elements
###dat.list a list of data matrix
###ref.de.param.list the DE gene criteria for each reference dataset (optional)
###meta.df merged meta data for all datasets. 
###cl.list clusters for each dataset (optional)
###cl.df.list cluster annotations for each dataset (optional) 
prepare_harmonize<- function(dat.list, meta.df=NULL, cl.list=NULL, cl.df.list = NULL, de.param.list=NULL, de.genes.list=NULL, rename=TRUE)
  {
    common.genes = row.names(dat.list[[1]])
    for(x in 2:length(dat.list)){
      common.genes= intersect(common.genes, row.names(dat.list[[x]]))
    }
    if(rename){
      for(x in names(dat.list)){
        colnames(dat.list[[x]]) = paste(x, colnames(dat.list[[x]]), sep=".")
      }
      if(!is.null(cl.list)){
        for(x in names(cl.list)){
          names(cl.list[[x]]) = paste(x, names(cl.list[[x]]), sep=".")
        }
      }
    }
    platform = do.call("c",lapply(names(dat.list), function(p){
      dat = dat.list[[p]]
      setNames(rep(p, ncol(dat)), colnames(dat))
    }))
    platform = factor(platform, names(dat.list))
    gene.counts <- do.call("c",lapply(names(dat.list), function(p){
      dat = dat.list[[p]]
      setNames(Matrix::colSums(dat > 0), colnames(dat))
    }))
    df = data.frame(platform, gene.counts)
    if(!is.null(meta.df)){
      meta.df = cbind(meta.df, df[row.names(meta.df),])
    }
    else{
      meta.df = df
    }
    all.cells = unlist(lapply(dat.list, colnames))
    comb.dat = list(dat.list=dat.list, meta.df = meta.df, cl.list=cl.list, cl.df.list = cl.df.list, de.genes.list = de.genes.list, de.param.list= de.param.list, common.genes=common.genes, all.cells= all.cells)
  }



#' Test knn
#'
#' @param knn 
#' @param cl 
#' @param reference 
#' @param ref.cl 
#'
#' @return
#' @export
#'
#' @examples
test_knn <- function(knn, cl, reference, ref.cl, plot=FALSE)
  {
    library(reshape)
    library(ggplot2)
    cl=  cl[row.names(knn)]
    if(is.factor(cl)){
      cl = droplevels(cl)
    }
    ref.cl =ref.cl[reference]
    if(is.factor(ref.cl)){
      ref.cl = droplevels(ref.cl)
    }
    if(length(unique(cl)) <=1 | length(unique(ref.cl)) <= 1){
      return(NULL)
    }
    pred.result = predict_knn(knn, reference, ref.cl)
    pred.prob = pred.result$pred.prob    
    if(length(unique(pred.prob$nn.cl))){
      return(NULL)
    }
    pred.prob$cl = cl[pred.prob$query]
    pred.prob$n = NULL
    cl.pred.prob=pred.prob %>% group_by(cl,nn.cl) %>% summarize(n=sum(freq))
    cl.size = as.data.frame(table(cl))
    cl.pred.prob = cl.pred.prob %>% left_join(cl.size) %>% mutate(prob=n/Freq)
    match.df = cl.pred.prob %>% group_by(cl) %>% summarize(match.cl = nn.cl[which.max(prob)], 
      match.score = max(prob), Freq=max(Freq))
    
    cl.score = mean(match.df$match.score)
    cell.score = with(match.df, sum(match.score * Freq)/sum(Freq))
        
    g <- ggplot(cl.pred.prob, 
                aes(x = cl, y = nn.cl)) + 
                  geom_point(aes(color = prob)) + 
                    theme(axis.text.x = element_text(vjust = 0.1,
                            hjust = 0.2, 
                            angle = 90,
                            size = 7),
                          axis.text.y = element_text(size = 6)) + 
                            scale_color_gradient(low = "white", high = "darkblue") + scale_size(range=c(0,3))
    result=list(cl.score=cl.score, cell.score= cell.score, cell.pred.prob = pred.prob, cl.pred.prob = cl.pred.prob, g=g)
    return(result)
  }


#' Sample sets lists
#'
#' @param cells.list 
#' @param cl.list 
#' @param cl.sample.size 
#' @param sample.size 
#'
#' @return
#' @export
#'
#' @examples
sample_sets_list <- function(cells.list, cl.list, cl.sample.size=100, sample.size=5000)
  {
    for(x in names(cells.list)){
      if(length(cells.list[[x]]) > sample.size){
        if(is.null(cl.list[[x]])){
          cells.list[[x]] = sample(cells.list[[x]], sample.size)
        }
        else{
          tmp.cl = cl.list[[x]][cells.list[[x]]]
          if(is.factor(tmp.cl)){
            tmp.cl = droplevels(tmp.cl)
          }
          cl.size = table(tmp.cl)
          cells.list[[x]] = sample_cells(tmp.cl, max(cl.sample.size,round(sample.size/sum(cl.size >= 4))))
        }
      }
    }
    return(cells.list)
  }



#' Batch process
#'
#' @param x 
#' @param batch.size 
#' @param FUN 
#' @param mc.cores 
#' @param .combine 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
batch_process <- function(x, batch.size, FUN, mc.cores=1, .combine="c",...)
  {
    require(foreach)
    require(doMC)
    if (mc.cores == 1) {
      registerDoSEQ()
    }
    else {
      registerDoMC(cores=mc.cores)
      #on.exit(parallel::stopCluster(), add = TRUE)
    }
    bins = split(x, floor((1:length(x))/batch.size))
    results= foreach(i=1:length(bins), .combine=.combine) %dopar% FUN(bins[[i]],...)
    return(results)
  }
 
pvec_no_combine <- function (v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE, 
    mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE) 
{
    if (!is.vector(v)) 
        stop("'v' must be a vector")
    cores <- as.integer(mc.cores)
    if (cores < 1L) 
        stop("'mc.cores' must be >= 1")
    if (cores == 1L) 
        return(FUN(v, ...))
    if (mc.set.seed) 
        mc.reset.stream()
    n <- length(v)
    l <- if (n <= cores) 
        as.list(v)
    else {
        il <- as.integer(n/cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) 
            sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[si[ix]:se[ix]])
    }
    jobs <- NULL
    prepareCleanup()
    on.exit(cleanup())
    FUN <- match.fun(FUN)
    jobs <- lapply(seq_len(min(n, cores)), function(i) mcparallel(FUN(l[[i]], 
        ...), name = i, mc.set.seed = mc.set.seed, silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    return(res)
  }

knn_combine <- function(results)
{
  knn.index = do.call("rbind", lapply(results, function(x)x[[1]]))
  knn.distance = do.call("rbind", lapply(results, function(x)x[[2]]))
  return(list(knn.index, knn.distance))
}

#' get knn batch
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#' @param batch.size 
#' @param mc.cores 
#'
#' @return
#' @export
#'
#' @examples
get_knn_batch <- function(dat, ref.dat, k, method="cor", dim=NULL, batch.size, mc.cores=1,return.distance=FALSE,...)
  {
    if(return.distance){
      fun = "knn_combine"
    }
    else{
      fun = "rbind"
    }
    results <- batch_process(x=1:ncol(dat), batch.size=batch.size, mc.cores=mc.cores, .combine=fun, FUN=function(bin){
      get_knn(dat=dat[,bin,drop=F], ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance, ...)
    })
    return(results)
  }
    



#' Get KNN
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#'
#' @return
#' @export
#'
#' @examples
get_knn <- function(dat, ref.dat, k, method ="cor", dim=NULL,index=NULL, build.index=FALSE, transposed=TRUE, return.distance=FALSE)
  {
    if(transposed){
      cell.id = colnames(dat)
    }
    else{
      cell.id= row.names(dat)
    }
    
    if(method=="cor"){
      if(transposed){
        knn.result = knn_cor(ref.dat, dat,k=k)
      }
      else{
        ref.dat = Matrix::t(ref.dat)
        dat = Matrix::t(dat)
        knn.result = knn_cor(ref.dat, dat,k=k)
      }
    }
    else{
      if(transposed){
        if(is.null(index)){
          ref.dat = Matrix::t(ref.dat)
        }
        dat = Matrix::t(dat)
      }                  
      if(method=="RANN"){
        knn.result = RANN::nn2(ref.dat, dat, k=k)
      }
      else if(method %in% c("Annoy.Euclidean", "Annoy.Cosine")){
        library(BiocNeighbors)      
        if(is.null(index)){
          if (method=="Annoy.Cosine"){
            ref.dat = l2norm(ref.dat,by = "row")
          }
          if(build.index){
            index= buildAnnoy(ref.dat)
          }
        }
        if (method=="Annoy.Cosine"){       
          dat = l2norm(dat,by="row")
        }
        knn.result = queryAnnoy(X= ref.dat, query=dat, k=k, precomputed = index)
      }
      else if(method == "CCA"){
        mat3 = crossprod(ref.dat, dat)
        cca.svd <- irlba(mat3, dim=dim)
        knn.result = knn_cor(cca.svd$u, cca.svd$v,  k=k)
      }
      else{
        stop(paste(method, "method unknown"))
      }
    }
    knn.index= knn.result[[1]]
    knn.distance = knn.result[[2]]
    row.names(knn.index) = row.names(knn.distance)=cell.id
    if(!return.distance){
      return(knn.index)
    }
    else{
      list(index=knn.index, distance=knn.distance)
    }
  }


#' Select joint genes
#'
#' @param comb.dat 
#' @param ref.list 
#' @param select.cells 
#' @param maxGenes 
#' @param vg.padj.th 
#' @param max.dim 
#' @param use.markers 
#' @param top.n 
#' @param rm.eigen 
#' @param rm.th 
#'
#' @return
#' @export
#'
#' @examples
select_joint_genes  <-  function(comb.dat, ref.list, select.cells = comb.dat$all.cells, maxGenes=2000, vg.padj.th=0.5, max.dim=20,use.markers=TRUE, top.n=100,rm.eigen=NULL, rm.th=rep(0.7, ncol(rm.eigen)))
  {
    library(matrixStats)
    select.genes.list = list()
    for(ref.set in names(ref.list)){
      #print(ref.set)
      ref.cells = intersect(ref.list[[ref.set]], select.cells)
      ref.dat = comb.dat$dat.list[[ref.set]][,ref.cells]
      ###if cluster membership is available, use cluster DE genes
      if(use.markers & !is.null(comb.dat$de.genes.list[[ref.set]])){
        cl = droplevels(comb.dat$cl.list[[ref.set]][ref.cells])
        cl.size = table(cl)
        cl = droplevels(cl[cl %in% names(cl.size)[cl.size > comb.dat$de.param.list[[ref.set]]$min.cells]])
        if(length(levels(cl)) <= 1){
          return(NULL)
        }
        de.genes = comb.dat$de.genes.list[[ref.set]]
        #print(length(de.genes.list[[ref.set]]))
        select.genes = display_cl(cl, norm.dat=ref.dat, max.cl.size = 200, n.markers=20, de.genes= de.genes)$markers
        select.genes = intersect(select.genes, comb.dat$common.genes)
      }
####if cluster membership is not available, use high variance genes and genes with top PCA loading
      else{
        tmp.dat = ref.dat[Matrix::rowSums(ref.dat >= 1) >=comb.dat$de.param.list[[ref.set]]$min.cells, ]
        tmp.dat@x = 2^tmp.dat@x - 1
        vg = find_vg(tmp.dat)
        rm(tmp.dat)        
        #gc()
        select.genes = intersect(row.names(vg)[which(vg$loess.padj < vg.padj.th | vg$dispersion >3)],comb.dat$common.genes)
        
        if(length(select.genes) < 5){
          return(NULL)
        }
        select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
        rd = rd_PCA(norm.dat=ref.dat,select.genes, ref.cells, max.pca = max.dim)
        if(is.null(rd)){
          return(NULL)
        }
        rd.dat = rd$rd.dat
        rot = t(rd$pca$rotation[,1:ncol(rd$rd.dat)])
        if(!is.null(rm.eigen)){
          rm.cor=abs(cor(rd.dat, rm.eigen[row.names(rd.dat),]))
          rm.cor[is.na(rm.cor)]=0
          select = t(t(rm.cor) < rm.th)
          select = apply(select, 1, all)          
          if(sum(select)==0){
            return(NULL)
          }
          #print(rm.cor)
          #if(sum(!select)>0){
          #  print(rm.cor[!select,,drop=F])
          #}
          rot = rot[,select,drop=FALSE]
        }
        if(is.null(rot)){
          return(NULL)
        }
        rot.scaled = (rot  - rowMeans(rot))/rowSds(rot)
        gene.rank = t(apply(-abs(rot), 1, rank))
        select = gene.rank <= top.n & abs(rot.scaled ) > 2
        select.genes = colnames(select)[colSums(select)>0]
      }
      rm(ref.dat)
      select.genes.list[[ref.set]] = select.genes
    }
    gene.score = table(unlist(select.genes.list))
    if(length(gene.score)==0){
      return(NULL)
    }
    select.genes= names(head(sort(gene.score, decreasing=T), maxGenes))
    return(select.genes)
  }


cleanAnnoyIndex <- function(index)
  {
    unlink(index@path)
    rm(index)
    gc()
  }


#' compute knn
#'
#' @param comb.dat 
#' @param select.genes 
#' @param ref.list 
#' @param select.sets 
#' @param select.cells 
#' @param k 
#' @param method 
#' @param self.method 
#' @param batch.size 
#' @param mc.cores 
#'
#' @return
#' @export
#'
#' @examples
compute_knn <- function(comb.dat, select.genes, ref.list, select.sets=names(comb.dat$dat.list), select.cells=comb.dat$all.cells, k=15, cross.knn.method=c("Annoy.Cosine","cor"), self.knn.method=c("Annoy.Euclidean","RANN"), batch.size=10000, mc.cores=1)
  {
    
    cat("Number of select genes", length(select.genes), "\n")
    cat("Get knn\n")
    dat.list = comb.dat$dat.list
###index is the index of knn from all the cells
    knn.list = list()
    for(ref.set in names(ref.list)){
      cat("Ref ", ref.set, "\n")
      if(length(ref.list[[ref.set]]) <= k) {
        ##Not enough reference points to compute k
        next
      }
      k.tmp = k
      if(length(ref.list[[ref.set]]) <= k*2) {
        k.tmp = round(k/2)
      }
      dat = dat.list[[ref.set]]
      ref.cells = ref.list[[ref.set]]
      ref.dat = dat[select.genes,ref.cells]
      ref.cells = colnames(ref.dat)[Matrix::colSums(ref.dat) > 0]
      ref.dat = ref.dat[,ref.cells]
      map.cells=  intersect(select.cells, colnames(dat))
      if(length(map.cells)==0){
        next
      }
      tmp.cores = mc.cores
      if(length(map.cells)< batch.size){
        tmp.cores = 1
      }                
      rd.dat = rd_PCA(dat,select.genes=select.genes, select.cells=map.cells, max.pca = 50, sampled.cells=ref.cells, th=1, mc.cores=tmp.cores)$rd.dat
      ref.rd.dat = rd.dat[ref.cells,,drop=F]
      idx = match(ref.cells, comb.dat$all.cells)
      index = NULL
      if(length(select.cells) >50000 & self.knn.method %in% c("Annoy.Euclidean")){
        require(BiocNeighbors)
        index = buildAnnoy(ref.rd.dat, distance ="Euclidean", transposed = FALSE)
      }      
      knn=get_knn_batch(dat=rd.dat, ref.dat = ref.rd.dat, k=k.tmp, method = self.knn.method, batch.size = batch.size, mc.cores=tmp.cores, index=index, transposed=FALSE)
      if(!is.null(index)){
        cleanAnnoyIndex(index)
      }
      knn = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))
      self.knn = knn
      index = NULL
      if(cross.knn.method  %in% c("Annoy.Euclidean", "Annoy.Cosine")){
        if(cross.knn.method=="Annoy.Cosine"){
          distance = "Cosine"
        }
        else{
          distance = "Euclidean"
        }
        if(length(select.cells)>50000){
          index = buildAnnoy(ref.dat, distance =distance, transposed = TRUE)
        }
      }
      knn =do.call("rbind", lapply(setdiff(select.sets,ref.set), function(set){
        cat("Set ", set, "\n")
        map.cells=  intersect(select.cells, colnames(dat.list[[set]]))
        if(length(map.cells)==0){
          return(NULL)
        }
        dat = dat.list[[set]][select.genes,map.cells,drop=F]      
        tmp.cores = mc.cores
        if(ncol(dat)< batch.size){
          tmp.cores = 1
        }         
        knn=get_knn_batch(dat=dat, ref.dat = ref.dat, k=k.tmp, method = cross.knn.method, batch.size = batch.size, mc.cores=tmp.cores, index=index, transposed=TRUE)
        #if(!is.null(comb.dat$cl.list)){
        #  test.knn = test_knn(knn, comb.dat$cl.list[[set]], colnames(ref.dat), comb.dat$cl.list[[ref.set]])          
        #  if(!is.null(test.knn)){
        #    cat("Knn", set, ref.set, cross.knn.method, "cl.score", test.knn$cl.score, "cell.score", test.knn$cell.score,"\n")
        #  }
        #}
        knn = matrix(idx[knn], nrow=nrow(knn), dimnames=list(row.names(knn), NULL))        
      }))
      if(!is.null(index)){
        cleanAnnoyIndex(index)
      }
      knn = rbind(self.knn, knn)
      knn.list[[ref.set]] = knn[select.cells,]      
    }
    knn.comb = do.call("cbind",knn.list)
    return(knn.comb)
  }

#' knn joint
#'
#' @param comb.dat 
#' @param ref.sets 
#' @param select.sets 
#' @param merge.sets 
#' @param select.cells 
#' @param select.genes 
#' @param cross.knn.method 
#' @param self.knn.method 
#' @param k 
#' @param sample.size 
#' @param cl.sample.size 
#' @param batch.size 
#' @param verbose 
#' @param mc.cores 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
knn_joint <- function(comb.dat, ref.sets=names(comb.dat$dat.list), select.sets= names(comb.dat$dat.list),merge.sets=ref.sets, select.cells=comb.dat$all.cells, select.genes=NULL, cross.knn.method="Annoy.Cosine", self.knn.method = "Annoy.Euclidean", method = "louvain", k=15,  sample.size = 5000, cl.sample.size = 100, batch.size = 10000, verbose=TRUE,mc.cores=1,...)
{
  #attach(comb.dat)
  with(comb.dat,{
  cat("Number of select cells", length(select.cells), "\n")
  cells.list = split(select.cells, meta.df[select.cells, "platform"])[select.sets]
  cells.list =  sample_sets_list(cells.list, cl.list[names(cl.list) %in% select.sets], sample.size=sample.size, cl.sample.size = cl.sample.size)
  ref.list = cells.list[ref.sets]
  
###Select genes for joint analysis
  if(is.null(select.genes)){
    select.genes = select_joint_genes(comb.dat, ref.list = ref.list,select.cells=select.cells, ...)
  }
  if(length(select.genes) < 5){
    return(NULL)
  }
  
  knn.comb= compute_knn(comb.dat, select.genes=select.genes, ref.list=ref.list, select.sets=select.sets, select.cells=select.cells, k=k, cross.knn.method =cross.knn.method, self.knn.method=self.knn.method, batch.size=batch.size, mc.cores=mc.cores)
  if(is.null(knn.comb)){
    return(NULL)
  }
  sampled.cells = unlist(cells.list)
  
  result = knn_jaccard_clust(knn.comb[sampled.cells,], method=method)
  result$knn = knn.comb
  ###preliminary clusters from louvain or leiden
  cl = result$cl
  if(length(cl) < nrow(result$knn)){
    diff.cells = setdiff(row.names(result$knn), names(cl))
    pred.df = predict_knn(result$knn[diff.cells,], all.cells, cl )$pred.df
    pred.cl= setNames(as.character(pred.df$pred.cl), row.names(pred.df))
    cl = c(cl, pred.cl[setdiff(names(pred.cl), names(cl))])     
  }
  cl.platform.counts = table(meta.df[names(cl), "platform"],cl)
  ###If a cluster is not present in reference sets, split the cells based on imputed cluster based on cells in reference set.

  ref.de.param.list = de.param.list[ref.sets]
  cl.min.cells = sapply(ref.de.param.list, function(x)x$min.cells)
  cl.big= cl.platform.counts[ref.sets,,drop=F] >= cl.min.cells
  bad.cl = colnames(cl.big)[colSums(cl.big) ==0]
  if(length(bad.cl) > 0){
    tmp.cells = names(cl)[cl %in% bad.cl]
    pred.df = predict_knn(result$knn[tmp.cells,,drop=F], all.cells, cl)$pred.df
    pred.cl= setNames(as.character(pred.df$pred.cl), row.names(pred.df))
    cl[names(pred.cl)]= pred.cl
  }
  merge.dat.list = comb.dat$dat.list[merge.sets]
  if(length(cl) < 5000  & length(cl) < length(comb.dat$all.cells)/2){
    for(x in names(merge.dat.list)){
      merge.dat.list[[x]] = merge.dat.list[[x]][, colnames(merge.dat.list[[x]]) %in% names(cl),drop=F]
    }
  }
  cl  = merge_cl_multiple(comb.dat, merge.dat.list=merge.dat.list, cl=cl, anchor.genes=select.genes)
  if(length(unique(cl))<=1){
    return(NULL)
  }
  if(verbose){
    print(table(cl))
  }
  result$ref.list = ref.list
  result$cl = cl
  result$markers = select.genes
  result$select.genes= select.genes
  result$ref.de.param.list = ref.de.param.list
  rm(merge.dat.list)
  gc()
  return(result)
})
}

#' Sim knn
#'
#' @param sim 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
sim_knn <- function(sim, k=15)
{
  require(matrixStats)
  th =  rowOrderStats(as.matrix(sim), which=ncol(sim)-k+1)
  select = sim >= th
  knn.index = t(apply(select, 1, function(x)head(which(x),k)))  
  knn.distance = do.call("rbind",lapply(1:nrow(sim), function(i) (1- sim[i,,drop=F])[knn.index[i,,drop=F]]))
  return(list(knn.index, knn.distance))
}

#' KNN cor
#'
#' @param ref.dat 
#' @param query.dat 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
knn_cor <- function(ref.dat, query.dat, k = 15)
{
  #sim = cor(as.matrix(query.dat), as.matrix(ref.dat), use="pairwise.complete.obs")
  sim = cor(as.matrix(query.dat), as.matrix(ref.dat))
  sim[is.na(sim)] = 0
  return(sim_knn(sim, k=k))
}

#' KNN cosine
#'
#' @param ref.dat 
#' @param query.dat 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
knn_cosine <- function(ref.dat, query.dat, k = 15)
  {
    library(qlcMatrix)
    sim=cosSparse(query.dat, ref.dat)
    sim[is.na(sim)] = 0
    knn.idx = sim_knn(sim, k=k)
    return(knn.idx)
  }



#' KNN Jaccard Louvain
#'
#' @param knn.index 
#'
#' @return
#' @export
#'
#' @examples
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param knn.index 
##' @param method 
##' @param prune 
##' @return 
##' @author Zizhen Yao
knn_jaccard_clust <- function(knn.index, method=c("louvain","leiden"),prune=0.05,...)
  {
    require(igraph)
    cat("Get jaccard\n")
    #sim=knn_jaccard(knn.index,...)
    sim = ComputeSNN(knn.index,prune=prune)
    rownames(sim) = colnames(sim) = row.names(knn.index)
    gr <- igraph::graph.adjacency(sim, mode = "undirected", 
                                  weighted = TRUE)
    if(method[1]=="louvain"){
      cat("Louvain clustering\n")
      result <- igraph::cluster_louvain(gr,...)
    }
    else{
      cat("Leiden clustering\n")
      library(leidenAlg)
      result <- leiden.community(gr,...)      
    }
    result$cl=membership(result)
    rm(sim)
    gc()
    return(result)
  }




#' Harmonize
#'
#' @param comb.dat 
#' @param prefix 
#' @param overwrite 
#' @param dir 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
harmonize <- function(comb.dat, prefix, overwrite=TRUE, dir=".",...)
  {
    if(!file.exists(dir)){
      dir.create(dir)
    }
    fn = file.path(dir, paste0(prefix, ".rda"))
    if(!overwrite){
      if(file.exists(fn)){
        load(fn)
        return(result)
      }
    }
    result = knn_joint(comb.dat, ...)
    save(result, file=fn)
    if(is.null(result)){
      return(NULL)
    }
    print("Cluster size")
    print(table(result$cl))
    #g = plot_cl_meta_barplot(result$cl, meta.df[names(result$cl), "platform"])
    #g = g + theme(axis.text.x = element_text(angle=45,hjust=1, vjust=1))
    #ggsave(paste0(prefix, ".platform.barplot.pdf"),g,height=5, width=12)
    #plot_confusion(result$cl, prefix,comb.dat)
    return(result)
  }




#' i harmonize
#'
#' @param comb.dat 
#' @param select.cells 
#' @param ref.sets 
#' @param prefix 
#' @param result 
#' @param overwrite 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
i_harmonize <- function(comb.dat, select.cells=comb.dat$all.cells, ref.sets=names(comb.dat$dat.list), prefix="", result=NULL, overwrite=TRUE, split.size = 50,...)
  {
    
    #attach(comb.dat)
    if(is.null(result)){
      result = harmonize(comb.dat=comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=prefix, overwrite=overwrite,...)
    }
    if(is.null(result)){
      return(NULL)
    }
    all.results= list(result)
    names(all.results) = prefix
    cl = result$cl
    for(i in as.character(sort(unique(result$cl)))){
      tmp.result = with(comb.dat, {
        tmp.prefix=paste(prefix, i,sep=".")
        print(tmp.prefix)
        select.cells= names(cl)[cl == i]
        platform.size = table(meta.df[select.cells, "platform"])        
        print(platform.size)
     
        sets = names(platform.size)
        
        pass.th = sapply(ref.sets, function(set)platform.size[[set]] >= de.param.list[[set]]$min.cells)
        pass.th2 = sapply(sets, function(set)platform.size[[set]] >= de.param.list[[set]]$min.cells*2)
        if(sum(pass.th) >= length(ref.sets) & sum(pass.th2) >= 1){
          tmp.result = i_harmonize(comb.dat, select.cells=select.cells, ref.sets=ref.sets, prefix=tmp.prefix, overwrite=overwrite, ...)
          }
        else{
          tmp.result = NULL
        }
      })
      if(!is.null(tmp.result)){
        all.results[names(tmp.result)] = tmp.result
      }           
    }
    return(all.results) 
  }



#' combine cl
#'
#' @param all.results 
#'
#' @return
#' @export
#'
#' @examples
combine_cl <- function(all.results)
  {
    cl = all.results[[1]]$cl
    cl = setNames(as.integer(cl),names(cl))
    markers=all.results[[1]]$markers
    n.cl = max(cl)
    for(i in 2:length(all.results)){
      #if(is.null(all.results[[i]]$cl) | length(unique(all.results[[i]]$cl)) < 2) next
      if(is.null(all.results[[i]]$cl)) next
      new.cl = all.results[[i]]$cl
      new.cl = setNames(as.integer(new.cl)+ n.cl,names(new.cl))
      cl[names(new.cl)] = new.cl
      n.cl = max(cl)
      cat(names(all.results)[i], n.cl, "\n")
    }
    return(cl)
  }




#' merge knn result
#'
#' @param split.results 
#'
#' @return
#' @export
#'
#' @examples
merge_knn_result <- function(split.results)
  {
    ref.cells = unlist(lapply(split.results, function(x)x$ref.cells))
    ref.cells = ref.cells[!duplicated(ref.cells)]
    markers =  unique(unlist(lapply(split.results, function(x)x$markers)))
    n.cl = 0
    cl = NULL
    cl.df = NULL
    knn = NULL
    knn.merge= NULL
    for(result in split.results){
      tmp.cl = setNames(as.integer(as.character(result$cl)) + n.cl, names(result$cl))
      tmp.cl.df = result$cl.df
      row.names(tmp.cl.df) = as.integer(row.names(tmp.cl.df)) + n.cl 
      cl = c(cl, tmp.cl)
      cl.df = rbind(cl.df, tmp.cl.df)
      n.cl = max(as.integer(as.character(cl)))
      orig.index = match(result$ref.cells, ref.cells)
      tmp.knn = result$knn[names(tmp.cl),]
      tmp.knn = matrix(orig.index[tmp.knn], nrow=nrow(tmp.knn))
      knn = rbind(knn, tmp.knn)
      tmp.knn = result$knn.merge[names(tmp.cl),]
      tmp.knn = matrix(orig.index[tmp.knn], nrow=nrow(tmp.knn))
      knn.merge = rbind(knn.merge, tmp.knn)
    }
    new.result = list(cl = as.factor(cl), cl.df = cl.df, markers=markers, knn=knn, ref.cells =ref.cells, knn.merge = knn.merge)
    return(new.result)
  }





#' Plot tsne
#'
#' @param cl 
#' @param cl.df 
#' @param comb.dat 
#' @param prefix 
#' @param tsne.df 
#' @param cex 
#' @param fn.size 
#' @param height 
#' @param width 
#'
#' @return
#' @export
#'
#' @examples
plot_tsne <- function(cl, cl.df, comb.dat, prefix, tsne.df, cex=0.3, fn.size=2, height=8, width=10)
  {
    library(ggplot2)
    library(gridExtra)
    with(comb.dat,{
    #cl.df$cluster_color = adjust_color(cl.df$cluster_color, adj.col)
    tmp = plot_tsne_cl(cl=cl, cl.df=cl.df,  tsne.df = tsne.df, cex=cex, fn.size = fn.size)
    tsne.df = tmp$tsne.df
    ggsave(paste("tsne.cl",prefix,"pdf", sep="."), tmp$g, height=height,width=width)
    
    tmp.df = meta.df[row.names(tsne.df),]
    tmp= setNames(as.factor(tmp.df$platform), row.names(tsne.df))
    meta.col=NULL
    if(length(levels(tmp))==2){
      meta.col = setNames(c("blue", "orange"), levels(tmp))
    }
    g= plot_tsne_meta(tsne.df, tmp, meta.col=meta.col, cex=cex)
    ggsave(paste("tsne",prefix, "platform.pdf", sep="."), g, height=height, width=width)

    plots = lapply(names(cl.list),function(x){
      tmp.cl = cl.list[[x]]
      tmp.cl = droplevels(tmp.cl[names(tmp.cl) %in% names(cl)])
      if(length(tmp.cl)==0){
        return(NULL)
      }
      g = plot_tsne_cl(cl=cl.list[[x]], cl.df=cl.df.list[[x]],  tsne.df = tsne.df[names(cl.list[[x]]),], cex=cex, fn.size = fn.size)$g
      g = g + ggtitle(x)
    })
    plots = plots[!sapply(plots, is.null)]
    ggsave(paste("tsne",prefix, "by.platform.cl.pdf", sep="."), marrangeGrob(grobs=plots, nrow=1, ncol=1), height=height, width=width)

    plots = lapply(names(cl.list),function(x){
      tmp.cells = names(cl)[meta.df[names(cl),"platform"]==x]
      if(length(tmp.cells)==0){
        return(NULL)
      }
      g = plot_tsne_cl(cl=droplevels(cl[tmp.cells]), cl.df=cl.df,  tsne.df = tsne.df[names(cl.list[[x]]),], cex=cex, fn.size = fn.size)$g
      g = g + ggtitle(x)
    })
    plots = plots[!sapply(plots, is.null)]
    ggsave(paste("tsne",prefix, "by.platform.pdf", sep="."), marrangeGrob(grobs=plots, nrow=1, ncol=1), height=height, width=width)
  })
  }



#' Plot markers cl means
#'
#' @param select.genes 
#' @param gene.ordered 
#' @param cl.means.list 
#' @param comb.dat 
#' @param cl 
#' @param cl.col 
#' @param prefix 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_markers_cl_means <- function(select.genes, gene.ordered=FALSE, cl.means.list = NULL, comb.dat=NULL, cl=NULL, cl.col=NULL, prefix="",...)
  {
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    blue.red <-colorRampPalette(c("blue", "white", "red"))
    if(is.null(cl.means.list)){
      cl.means.list=get_cl_means_list(comb.dat$dat.list,  comb.dat$de.param.list,select.genes=select.genes, cl=cl)
    }
    else{
      cl.means.list = sapply(cl.means.list, function(x)x[select.genes,],simplify=F)
    }
    if(!gene.ordered){
      gene.hc = hclust(dist(cl.means.list[[1]]), method="ward.D")
      select.genes = select.genes[gene.hc$order]
    }
    if(is.null(cl.col)){
      cl.col = jet.colors(length(unique(cl)))
    }
    cl.col = matrix(cl.col, nrow=1)
    colnames(cl.col) = levels(cl)
    pdf(paste0(prefix, ".cl.heatmap.pdf"),...)
    for(set in names(cl.means.list)){
      dat = cl.means.list[[set]][select.genes, ]
      cexCol = min(70/ncol(dat),1)
      cexRow = min(60/nrow(dat),1)
      heatmap.3(dat, Rowv=NULL, Colv=NULL, col=blue.red(100), trace="none",dendrogram="none", cexCol=cexCol, cexRow=cexRow, ColSideColors = cl.col, main=set)
    }
    dev.off()
  }


#' Get gene cl correlation
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
get_gene_cl_correlation <- function(cl.means.list)
  {
    sets=names(cl.means.list)
    gene.cl.cor = list()
    for(i in 1:(length(cl.means.list)-1)){
      for(j in (i+1):length(cl.means.list)){
        pair= paste(sets[i], sets[j], sep=":")
        common.cl = intersect(colnames(cl.means.list[[i]]), colnames(cl.means.list[[j]]))
        gene.cor =  pair_cor(cl.means[[i]][common.genes,common.cl],cl.means[[j]][common.genes,common.cl])
        gene.cl.cor[[pair]] = gene.cor
      }
    }
    return(gene.cl.cor)
  }


#' Simple dendrogram
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
simple_dend <- function(cl.means.list)
{
  levels = unique(unlist(lapply(cl.means.list, colnames)))
  n.counts = cl.cor=matrix(0, nrow=length(levels), ncol=length(levels))
  row.names(n.counts) = row.names(cl.cor)=levels
  colnames(n.counts)=colnames(cl.cor)=levels
  for(x in cl.means.list){
    tmp.cor= cor(x)
    tmp.cor[is.na(tmp.cor)] = 0
    cl.cor[colnames(x),colnames(x)] = cl.cor[colnames(x),colnames(x)] + tmp.cor
    n.counts[colnames(x),colnames(x)] =   n.counts[colnames(x),colnames(x)] +1
  }
  cl.cor = cl.cor/n.counts
  cl.cor[is.na(cl.cor)] = 0
  dend=as.dendrogram(hclust(as.dist(1-cl.cor)))
  dend = dend %>% set("labels_cex", 0.7)
  if (!is.null(l.color)) {
    dend = dend %>% set("labels_col", l.color[labels(dend)])
  }
  dend = dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 
    0.5)
  if (!is.null(l.color)) {
    dend = dend %>% set("leaves_col", l.color[labels(dend)])
  }
  if (!is.null(l.rank)) {
    dend = reorder_dend(dend, l.rank)
  }
  return(list(dend=dend, cl.cor=cl.cor))
}

#' Impute val cor
#'
#' @param dat 
#' @param impute.dat 
#'
#' @return
#' @export
#'
#' @examples
impute_val_cor <- function(dat, impute.dat)
  {
    gene.cor = pair_cor(dat, impute.dat)
    gene.cor[is.na(gene.cor)] = 0
    return(gene.cor)
  }


#' Build dend harmonize
#'
#' @param impute.dat.list 
#' @param cl 
#' @param cl.df 
#' @param ncores 
#'
#' @return
#' @export
#'
#' @examples
build_dend_harmonize <- function(impute.dat.list, cl, cl.df, ncores=1)
  {
    cl.means = do.call("rbind",sapply(impute.dat.list, function(x)get_cl_means(x,cl),simplify=F))
    l.rank = setNames(1:nrow(cl.df), row.names(cl.df))
    l.color = setNames(as.character(cl.df$cluster_color), row.names(cl.df))
    dend.result = build_dend(cl.means, l.rank = l.rank, l.color = l.color, nboot=100,ncores=ncores)
  }

