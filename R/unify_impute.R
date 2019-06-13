impute_knn <- function(knn.idx, reference, dat)
  {
    query = row.names(knn.idx)
    impute.dat= matrix(0, nrow=nrow(knn.idx),ncol=ncol(dat))    
    k = rep(0, nrow(knn.idx))
    for(i in 1:ncol(knn.idx)){
      print(i)
      nn = reference[knn.idx[,i]]
      ###Ignore the neighbors not present in imputation reference
      select = nn %in% row.names(dat)
      impute.dat[select,]= impute.dat[select,] +  dat[nn[select],]
      k[select] = k[select]+1
    }
    impute.dat = impute.dat / k
    row.names(impute.dat) = row.names(knn.idx)
    colnames(impute.dat) = colnames(dat)
    return(impute.dat)
  }




#### assume within data modality have been performed
####
impute_knn_global <- function(comb.dat, split.results, select.genes, select.cells, ref.list, sets=comb.dat$sets, rm.eigen=NULL,max.dim=80, th=0.5)
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
        rd.result <- rd_PCA(comb.dat$dat.list[[x]], select.genes, select.cells=tmp.cells, sampled.cells = ref.cells, max.dim=max.dim, th=th, rm.eigen=rm.eigen, ncores=10)
        rd.dat  = rd.result$rd.dat
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
      for(ref.set in ref.sets){
        if(ref.set %in% names(result$ref.list)){
          tmp.cells = row.names(result$knn)
          add.cells=FALSE
          query.cells = intersect(tmp.cells[comb.dat$meta.df[tmp.cells,"platform"] != ref.set], select.cells)
          if(any(!query.cells %in% row.names(impute.dat.list[[ref.set]]))){
            add.cells=TRUE
            impute.genes = select.genes
          }
          else{
            impute.genes=intersect(select.genes,c(result$select.markers, result$select.genes))
          }
          select.cols = comb.dat$meta.df[comb.dat$all.cells[result$knn[1,]],"platform"] == ref.set
          if(sum(select.cols)==0){
            next
          }
          else{
            ref.cells = intersect(comb.dat$all.cells[unique(as.vector(knn[, select.cols]))],select.cells)            
            knn = result$knn[query.cells,select.cols]
            impute.dat = impute_knn(knn, comb.dat$all.cells, impute.dat.list[[ref.set]][ref.cells,impute.genes])
          }
          if(!add.cells){
            impute.dat.list[[ref.set]][query.cells, impute.genes] <- impute.dat
          }
          else{
            impute.dat.list[[ref.set]] <- rbind(impute.dat.list[[ref.set]],impute.dat)
          }
          rm(impute.dat)
          gc()
        }
      }
    }
    return(list(knn.list =knn.list, org.rd.dat.list = org.rd.dat.list,impute.dat.list=impute.dat.list, ref.list=ref.list))
  }

