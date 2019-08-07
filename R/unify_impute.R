
get_knn_weight <- function(knn.dist, scale=0.2, exclude.th = 0.0001)
  {
    w = exp(-knn.dist*scale)
    if(exclude.th >= 0){
      w[knn.dist < exclude.th] = 0
    }
    return(w)
  }

predict_knn <- function(knn.idx, reference, cl)
  {
    library(matrixStats)
    library(dplyr)
    query = row.names(knn.idx)
    df = data.frame(nn=as.vector(knn.idx), query=rep(row.names(knn.idx), ncol(knn.idx)))
    df = df %>% filter(!is.na(nn))
    tmp.df = data.frame(nn=1:length(reference), nn.cl=cl[reference])
    df = df %>% left_join(tmp.df)
    tb=with(df, table(query, nn.cl))
    tb = tb/rowSums(!is.na(knn.idx))[row.names(tb)]
    pred.cl = setNames(colnames(tb)[apply(tb, 1, which.max)], row.names(tb))
    pred.score = setNames(rowMaxs(tb), row.names(tb))
    pred.df = data.frame(pred.cl, pred.score)
    return(list(pred.df=pred.df, pred.prob = tb))
  }

predict_knn_new <- function(knn.idx, reference, cl, ...)
  {
    library(matrixStats)
    library(dplyr)
    dat = as.matrix(get_cl_mat(cl))
    result = impute_knn(knn.idx, reference, dat, ...)
    pred.cl = setNames(colnames(result)[apply(result, 1, which.max)], row.names(result))
    pred.score = setNames(rowMaxs(result), row.names(result))
    pred.df = data.frame(pred.cl, pred.score)
    return(list(pred.df=pred.df, pred.prob=result))
  }






impute_knn <- function(knn.idx, reference, dat, knn.dist=NULL, ...)
  {
    query = row.names(knn.idx)
    impute.dat= matrix(0, nrow=nrow(knn.idx),ncol=ncol(dat))
    if(!is.null(knn.dist)){
      w = get_knn_weight(knn.dist,...)
    }
    else{
      w = matrix(1, nrow=nrow(knn.idx),ncol=ncol(dat))
    }
    total.w = rep(0, nrow(knn.idx))
    for(i in 1:ncol(knn.idx)){
      print(i)
      nn = reference[knn.idx[,i]]
      ###Ignore the neighbors not present in imputation reference
      select = nn %in% row.names(dat)   
      impute.dat[select,]= impute.dat[select,] +  dat[nn[select],] * w[select, i]
      total.w[select] = total.w[select]+ w[select,i]
    }
    impute.dat = impute.dat / total.w
    row.names(impute.dat) = row.names(knn.idx)
    colnames(impute.dat) = colnames(dat)
    return(impute.dat)
  }


iter_impute_knn <- function(knn, ref, dat, tol=10^-3,max.iter=100,...)
  {
    old.dat = NULL
    iter=0
    while(TRUE){
      iter= iter+1
      new.dat = impute_knn(knn, ref, dat,...)
      if(!is.null(old.dat)){
        diff = new.dat - old.dat
        diff.scaled = sum(abs(diff))/sum(abs(old.dat))
        print(iter)
        print(diff.scaled)
        if(diff.scaled < tol | iter > max.iter ){
          break
        }
      }
      dat = old.dat = new.dat
    }
    return(new.dat)
  }







#### assume within data modality have been performed
####
impute_knn_global <- function(comb.dat, split.results, select.genes, select.cells, ref.list, sets=comb.dat$sets, max.dim=80, th=0.5, rm.eigen=NULL,rm.th=0.65)
  {
    org.rd.dat.list <- list()
    knn.list <- list()
    impute.dat.list <- list()
    ###Impute the reference dataset in the original space globally
    for(x in names(ref.list))
      {
        print(x)
        tmp.cells= select.cells[comb.dat$meta.df[select.cells,"platform"]==x]
        ref.cells = intersect(ref.list[[x]],tmp.cells)
        rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.pca =max.dim, th=th)
        if(!is.null(rm.eigen)){
          rd.dat  = filter_RD(rd.result$rd.dat, rm.eigen, rm.th)
        }
        print(ncol(rd.dat))
        knn.result <- RANN::nn2(data=rd.dat[ref.cells,], query=rd.dat, k=15)
        knn <- knn.result[[1]]
        row.names(knn) = row.names(rd.dat)    
        org.rd.dat.list[[x]] = rd.result
        knn.list[[x]]=knn
        knn = knn.list[[x]]
        impute.dat.list[[x]] <- impute_knn(knn, ref.cells, as.matrix(t(comb.dat$dat.list[[x]][select.genes,ref.cells])))
      }
    ###cross-modality Imputation based on nearest neighbors in each iteraction of clustering using anchoring genes or genes shown to be differentiall expressed. 
    for(x in names(split.results)){
      print(x)
      result = split.results[[x]]
      cl = result$cl
      knn = result$knn
      for(ref.set in names(result$ref.list)){
        tmp.cells = row.names(knn)
        add.cells=FALSE
        query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
        if(any(!query.cells %in% row.names(impute.dat.list[[ref.set]]))){
          add.cells=TRUE
          impute.genes = select.genes
        }
        else{
          impute.genes=intersect(select.genes,c(result$select.markers, result$select.genes))
        }
        select.cols = comb.dat$meta.df[comb.dat$all.cells[knn[1,]],"platform"] == ref.set
        if(sum(select.cols)==0){
          next
        }
        else{
          ref.cells = intersect(comb.dat$all.cells[unique(as.vector(knn[, select.cols]))],select.cells)            
          knn = knn[query.cells,select.cols]
          impute.dat = impute_knn(knn, comb.dat$all.cells, impute.dat.list[[ref.set]][ref.cells,impute.genes])
        }
        if(!add.cells){
          impute.dat.list[[ref.set]][query.cells, impute.genes] <- impute.dat
        }
        else{
          impute.dat.list[[ref.set]] <- rbind(impute.dat.list[[ref.set]],impute.dat)
        }
        print("Impute dimension")
        print(dim(impute.dat.list[[ref.set]]))
        rm(impute.dat)
        gc()
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }

