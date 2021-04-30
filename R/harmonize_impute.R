
#' Title
#'
#' @param knn.dist 
#' @param scale 
#' @param exclude.th 
#'
#' @return
#' @export
#'
#' @examples
get_knn_weight <- function(knn.dist, scale=0.2, exclude.th = 0.0001)
  {
    w = exp(-knn.dist*scale)
    if(exclude.th >= 0){
      w[knn.dist < exclude.th] = 0
    }
    return(w)
  }

#' Title
#'
#' @param knn.idx 
#' @param reference 
#' @param cl 
#'
#' @return
#' @export
#'
#' @examples
predict_knn_small <- function(knn.idx, reference, cl)
  {
    library(matrixStats)
    library(dplyr)
    query = row.names(knn.idx)
    df = data.frame(nn=as.vector(knn.idx), query=rep(row.names(knn.idx), ncol(knn.idx)))
    df = df[!is.na(df$nn),]
    df$nn.cl = cl[reference[df$nn]]
    tb=with(df, table(query, nn.cl))
    nn.size = table(df$query)
    tb = tb/as.vector(nn.size)
    pred.cl = setNames(colnames(tb)[apply(tb, 1, which.max)], row.names(tb))
    pred.score = setNames(rowMaxs(tb), row.names(tb))
    pred.df = data.frame(pred.cl, pred.score)
    return(list(pred.df=pred.df, pred.prob = tb))
  }

#' Title
#'
#' @param knn.idx 
#' @param reference 
#' @param cl 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
predict_knn <- function(knn.idx, reference, target, k=15, train.dat=NULL, test.dat=NULL, ...)
  {
    library(matrixStats)
    library(dplyr)
    if(is.null(knn.idx)){
      knn.idx = RANN::nn2(data=train.dat[reference,],query=test.dat,k=k)[[1]]        
      row.names(knn.idx) = row.names(test.dat)
    }
    dat = as.matrix(get_cl_mat(target))
    result = impute_knn(knn.idx, reference, dat, transpose_input=TRUE, transpose_output=TRUE,...)
    pred.target = setNames(colnames(result)[apply(result, 1, which.max)], row.names(result))
    pred.score = setNames(rowMaxs(result), row.names(result))
    pred.df = data.frame(pred.target, pred.score)
    return(list(pred.df=pred.df, pred.prob=result))
  }



#' Title
#'
#' @param knn.idx 
#' @param reference 
#' @param dat 
#' @param knn.dist 
#' @param w 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
impute_knn <- function(knn.idx, reference, dat, knn.dist=NULL, w=NULL, transpose_input=FALSE, transpose_output= FALSE)
  {
    cell.id = 1:nrow(knn.idx)
    if(transpose_input){
      reference.id = match(reference, row.names(dat))
      genes = colnames(dat)
    }
    else{
      reference.id = match(reference, colnames(dat))
      genes=row.names(dat)
    }
    
    if(transpose_output){
      impute.dat = matrix(0,length(cell.id), length(genes))
      row.names(impute.dat) = row.names(knn.idx)
      colnames(impute.dat) = genes
    }
    else{
      impute.dat = matrix(0, length(genes),length(cell.id))    
      row.names(impute.dat) = genes
      colnames(impute.dat) =  row.names(knn.idx)
    }
    ImputeKnn(knn.idx,cell.id, reference.id, dat, gene_idx_=NULL, w_mat_= w, impute.dat, transpose_input=transpose_input, transpose_output=transpose_output)
    impute.dat
}






impute_knn_global <- function(comb.dat, split.results, select.genes, select.cells, ref.list, sets=names(comb.dat$dat.list), max.dim=100, k=15, th=0.5, rm.eigen=NULL,rm.th=0.65,method="zscore",mc.cores=1,verbose=FALSE)
  {
    library(matrixStats)
    org.rd.dat.list <- list()
    knn.list <- list()
    impute.dat.list <- list()
    ###Impute the reference dataset in the original space globally
    for(x in names(ref.list))
      {
        print(x)
        impute.dat = matrix(0, nrow=length(select.genes), ncol=length(select.cells))
        dimnames(impute.dat) = list(select.genes, select.cells)
        impute.dat.list[[x]] = impute.dat        
        tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
        ref.cells = ref.list[[x]]
        rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th, method=method,mc.cores=mc.cores,verbose=verbose)
        org.rd.dat.list[[x]] = rd.result
        rd.result = org.rd.dat.list[[x]]
        if(!is.null(rm.eigen)){
          rd.dat  = filter_RD(rd.result$rd.dat, rm.eigen, rm.th,verbose=verbose)
        }
        #print(ncol(rd.dat))        
        knn = RANN::nn2(data=rd.dat[ref.cells,],query=rd.dat,k=k)[[1]]        
        row.names(knn) = row.names(rd.dat)
        cell.id = match(row.names(rd.dat), select.cells)
        dat = as.matrix(comb.dat$dat.list[[x]][select.genes,ref.cells])
        reference.id = 1:length(ref.cells)
        gene.id = 1:length(select.genes)
        ImputeKnn(knn, cell.id, reference.id, dat=dat, gene_idx_=gene.id, w_mat_= NULL, impute_dat=impute.dat.list[[x]], transpose_input=FALSE, transpose_output=FALSE)
      }

    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed. 
    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      if(x == names(split.results)[1]){
        impute.genes = select.genes
      }
      else{
        impute.genes=intersect(select.genes,c(result$markers, result$select.genes))
      }
      cl = result$cl
      knn = result$knn
      for(ref.set in intersect(names(result$ref.list),names(ref.list))){
        print(ref.set)
        tmp.cells = row.names(knn)
        query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
        select.cols = comb.dat$meta.df[comb.dat$all.cells[knn[1,]],"platform"] == ref.set
        cat("Impute genes", length(impute.genes),"\n")
        if(sum(select.cols)==0){
          next
        }
        if(length(query.cells)==0){
          next
        }

        select.knn = knn[query.cells,select.cols,drop=F]
        dat = impute.dat.list[[ref.set]]
        gene.id = match(impute.genes, row.names(dat))
        cell.id = match(query.cells, colnames(dat))
        reference.id = match(comb.dat$all.cells, colnames(dat))
        ImputeKnn(select.knn,cell.id, reference.id, dat=dat, gene_idx_=gene.id, w_mat_= NULL, impute_dat=dat, transpose_input=FALSE, transpose_output=FALSE)
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }


fast_knn <- function(query.dat, ref.dat=query.dat, distance="euclidean", k=15, M=16, ef=200,method="euclidean")
  {
    library("RcppHNSW")
    if(method=="euclidean"){
      p = new(HnswL2,ncol(ref.dat), nrow(ref.dat),M, ef)
    }
    else if(method=="cosine"){
      p = new(HnswCosine,ncol(ref.dat), nrow(ref.dat),M, ef)
    }
    else if(method=="ip"){
      p = new(HnswIP,ncol(ref.dat), nrow(ref.dat),M, ef)
    }
    p$addItems(ref.dat)
    knn.result = p$getAllNNsList(query.dat, k=15, include_distance=TRUE)
    row.names(knn.result[[1]]) = row.names(knn.result[[2]]) = row.names(query.dat)
    return(knn.result)
  }



