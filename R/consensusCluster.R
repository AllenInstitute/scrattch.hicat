
collect_co_matrix <- function(cl.files,all.cells)
  {
    subsample.cl=list()
    co.matrix= matrix(0, nrow=length(all.cells),ncol=length(all.cells))
    pr.matrix = matrix(0, nrow=length(all.cells),ncol=length(all.cells))
    row.names(co.matrix)=row.names(pr.matrix)=colnames(co.matrix)=colnames(pr.matrix)=all.cells
    for(f in cl.files){
      tmp=load(f)
      cl = cl[intersect(names(cl),all.cells)]
      pr.matrix[names(cl),names(cl)]=   pr.matrix[names(cl),names(cl)] + 1
      for(x in unique(cl)){
        y=names(cl)[cl==x]
        co.matrix[y,y]= co.matrix[y,y]+1
      }
      subsample.cl[[f]]= cl
    }
    co.ratio = co.matrix/pr.matrix
    return(list(co.ratio,subsample.cl))
  }



collect_co_matrix_sparseM <- function(norm.dat,result.files,all.cells,max.cl.size=1000)
  {
    select.cells=c()
    cl.list=list()
    for(f in result.files){
      print(f)
      tmp=load(f)
      cl= result$cl
      test.cells = setdiff(all.cells, names(cl))
      markers=unique(result$markers)
      map.df = map.by.cor(norm.dat[markers,names(cl)],cl, norm.dat[markers,test.cells],method="means")$pred.df
      test.cl = setNames(map.df$pred.cl, row.names(map.df))
      all.cl = c(cl, test.cl)[all.cells]
      cl.size=table(all.cl)
      cl.size[cl.size > max.cl.size]=max.cl.size
      sampled.size = table(all.cl[select.cells])
      sampled.size[setdiff(names(cl.size), names(sampled.size))]= 0
      more.cells= cl.size - sampled.size[names(cl.size)]
      more.cells=more.cells[more.cells >0]
      more.cl = all.cl[setdiff(names(all.cl), select.cells)]
      add.cells= unlist(lapply(names(more.cells), function(cl){
        sample(names(more.cl)[more.cl==cl], more.cells[cl])
      }))
      select.cells= c(select.cells, add.cells)
      cl.list[[f]] = all.cl
    }
    cl.comb = do.call("rbind", sapply(cl.list, function(cl){
      tmp.df= data.frame(cell=select.cells, cl=all.cl[select.cells])
      tb=xtabs(~cl+cell, data=tmp.df)
      tb = Matrix(tb, sparse=TRUE)
    },simplify=F))
    cl.comb = t(cl.comb)
    save(cl.comb, file="cl.comb.rda")
    co.ratio = crossprod(t(cl.comb))
    co.ratio@x = co.ratio@x/length(result.files)
    return(list(co.ratio=co.ratio, cl.comb=cl.comb, cl.list=cl.list))
  }
    


merge_co_matrix <- function(co.ratio1, co.ratio2)
  {
    all.cells=c(row.names(co.ratio1),row.names(co.ratio2))
    tmp.co = Matrix(0,nrow= ncol(co.ratio1), ncol=ncol(co.ratio2))
    colnames(tmp.co)=colnames(co.ratio2)
    row.names(tmp.co)=row.names(co.ratio1)
    tmp.co1 = rbind(co.ratio1, t(tmp.co))
    tmp.co2 = rbind(tmp.co, co.ratio2)
    co.ratio = cbind(tmp.co1, tmp.co2)
    return(co.ratio)
  }


split_co_matrix <- function(cl, select.cl,...)
  {
    cl = setNames(as.integer(as.character(cl)),names(cl))
    n.cl=max(cl)
    new.cl=cl
    for(i in select.cl){
      tmp.cells=names(cl)[cl==i]
      tmp.cl= recursive_split(select.cells=tmp.cells,prefix=paste0("cl",".",i), ...)
      if(!is.null(tmp.cl)){
        new.cl[names(tmp.cl)] = tmp.cl + n.cl
        n.cl = max(new.cl)
      }
    }
    cl=new.cl
    cl = setNames(as.factor(cl), names(cl))
    cl = setNames(as.integer(cl),names(cl))
    return(cl)
  }


recursive_split <- function(co.ratio, norm.dat, select.cells=colnames(co.ratio), all.col=NULL, diff.th=0.25, prefix=NULL, min.cells=4,method="auto", de.param=de_param(),verbose=FALSE,cl=NULL)
  {
    print(prefix)
    if(is.null(cl)){
      if(length(select.cells)  < 2 * min.cells){
        return(NULL)
      }
      if(method=="auto"){
        if (length(select.cells)> 3000){
          select.method = "louvain"
        }
        else{
          if(!is.matrix(co.ratio)){
            co.ratio = as.matrix(co.ratio[select.cells, select.cells])
            norm.dat = norm.dat[,select.cells]
          }
          select.method="ward"
          print("Use ward")
        }
      }
      else{
        select.method = method
      }
      if(select.method=="ward"){
        tmp = init_co(co.ratio, select.cells, min.cells=min.cells, th = diff.th,method=select.method)
        if(is.null(tmp)){
          return(NULL)
        }
        tmp.cl = tmp$cl
      }
      else{###louvain
        adj.mat=co.ratio[select.cells, select.cells]
        gr = graph.adjacency(adj.mat, mode="undirected",weighted=TRUE)
        comm= cluster_louvain(gr)
        rm(gr)
        gc()
        if(pass_louvain(modularity(comm), adj.mat)){
          tmp.cl = setNames(comm$membership,select.cells)
        }
        else{
          return(NULL)
        }
      }
      if(length(unique(tmp.cl))==1){
        return(NULL)
      }
      print(table(tmp.cl))
      cell.cl.co.ratio =get_cell.cl.co.ratio(co.ratio, tmp.cl)
      tmp= merge_cl(norm.dat=norm.dat, cl=tmp.cl, rd.dat=cell.cl.co.ratio, min.cells=min.cells, max.cl.size=300,de.param= de.param)
      cat("Merge result:", names(tmp),"\n")
      if(is.null(tmp) | !is.list(tmp)) return(NULL)
      if (length(unique(tmp$cl))==1) return(NULL)
      tmp.cl= tmp$cl
      tmp.cl=merge_cl_by_co(tmp.cl, co.ratio, diff.th)
      tmp.cl = setNames(as.integer(tmp.cl),names(tmp.cl))
      if(length(unique(tmp.cl))==1) {
        return(NULL)
      }
      cl = tmp.cl[select.cells]
      print(table(cl))
    }
    n.cl=max(cl)
    new.cl=cl
    for(i in sort(unique(cl))){
      tmp.cells=names(cl)[cl==i]            
      tmp.cl= recursive_split(co.ratio=co.ratio, norm.dat=norm.dat, select.cells=tmp.cells,prefix=paste0(prefix,".",i),all.col=all.col, diff.th =diff.th, min.cells=min.cells, method=method, counts=counts, use.voom=use.voom, verbose=verbose,...)
      if(!is.null(tmp.cl)){
        new.cl[names(tmp.cl)] = tmp.cl + n.cl
        n.cl = max(new.cl)
      }
    }
    cl=new.cl
    cl = setNames(as.integer(as.factor(cl)), names(cl))
    if(verbose){
      tmp = display_cl(cl, norm.dat, prefix, col=all.col[,select.cells], max.cl.size=200,min.cells=min.cells, de.param= de.param())
      tmp.cl = tmp$cl
      tmp.cells = names(tmp.cl)
      tmp.cells = tmp.cells[order(tmp.cl)]
      sep = tmp.cl[tmp.cells]
      sep = which(sep[-1]!=sep[-length(sep)])
      pdf(paste0(prefix, ".co.pdf"))
      heatmap.3(as.matrix(co.ratio[tmp.cells, tmp.cells]), col = blue.red(100), trace="none", ColSideColors=all.col[,tmp.cells], Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
      dev.off()
    }
    print("Finish")
    return(cl)
  }

merge_cl_by_co <- function(cl, co.ratio, diff.th=0.25){
  cell.cl.co.ratio = get_cell.cl.co.ratio(co.ratio,cl)
  cl.co.ratio <- do.call("rbind",tapply(names(cl),cl, function(x)colMeans(cell.cl.co.ratio[x,])))
  co.within= diag(cl.co.ratio)
  co.df <- as.data.frame(as.table(cl.co.ratio),stringsAsFactors=FALSE)
  co.df = co.df[co.df[,1]<co.df[,2]& co.df[,3]>0.1,]
  co.df$within1 = co.within[co.df[,1]] 
  co.df$within2 = co.within[co.df[,2]]
  co.df$diff = pmax(co.df$within1, co.df$within2) - co.df[,3]
  co.df = co.df[co.df$diff  < diff.th,]
  co.df = co.df[order(co.df[,1],decreasing=T),]
  print(co.df)
  for(i in 1:nrow(co.df)){
    cl[cl==co.df[i,2]]=co.df[i,1]
  }
  return(cl)
}

get_cell.cl.co.ratio <- function(co.ratio,cl)
  {
    return(get_cl_means(co.ratio, cl))
  }

    
get_cl_co_stats <- function(co.ratio, cl)
  {
    cell.cl.co.ratio= get_cell.cl.co.ratio(co.ratio,cl)
    cl.co.ratio <- get_cl_means(t(cell.cl.co.ratio), cl)

    cell.cl.confusion <- unlist(sapply(1:ncol(cell.cl.co.ratio),function(i){
      select.cells=names(cl)[cl==colnames(cell.cl.co.ratio)[i]]
      cell.cl.co.ratio=setNames(rowMaxs(cell.cl.co.ratio[select.cells, -i])/cell.cl.co.ratio[select.cells, i,drop=F], select.cells)
    },simplify=F))
    cell.cl.confusion = cell.cl.confusion[names(cl)]
    cl.confusion = setNames(sapply(1:nrow(cl.co.ratio),function(i){
      max(cl.co.ratio[i, -i])/cl.co.ratio[i,i]
    }),row.names(cl.co.ratio))
    return(list(cell.cl.co.ratio=cell.cl.co.ratio, cell.cl.confusion=cell.cl.confusion, cl.co.ratio=cl.co.ratio, cl.confusion=cl.confusion))
  }


init_co <- function(co.ratio, select.cells, min.cells=4, th = 0.3,method="ward",verbose=FALSE)
  {
    tmp.dat = as.matrix(co.ratio[select.cells, select.cells])
    hc = hclust(as.dist(1-tmp.dat),method=method)
    tmp.cl=cut_co_matrix(tmp.dat, hc$order,w=min.cells-1, th = th)
    ord = colnames(tmp.dat)[hc$order]
    sep = which(tmp.cl[-1]!=tmp.cl[-length(tmp.cl)])
    if(verbose){
      pdf("co.pdf")
      heatmap.2(tmp.dat[ord,ord], Colv=NULL,Rowv=NULL,trace="none",col=blue.red(100),colsep=sep,sepcolor="black")
      dev.off()
    }    
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    tmp.cl=merge_cl_by_co(tmp.cl, tmp.dat, diff.th=th)
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    ###reassign each cells to the clusters, if it prefer to co-cluster with another cluster much better. 
    cell.cl.co.ratio= get_cell.cl.co.ratio(tmp.dat, tmp.cl)
    tmp.cl2=setNames(apply(cell.cl.co.ratio, 1, which.max),row.names(cell.cl.co.ratio))
    org.co.score = get_pair_matrix(cell.cl.co.ratio, names(tmp.cl), as.character(tmp.cl))
    select = tmp.cl!=tmp.cl2[names(tmp.cl)]  & rowMaxs(cell.cl.co.ratio) - org.co.score  > 0.1
    tmp.cl[select] = tmp.cl2[select]
    tmp.cl = setNames(as.integer(as.factor(tmp.cl)),names(tmp.cl))
    if(length(unique(tmp.cl))==1){
      return(NULL)    
    }
    if(length(unique(tmp.cl))>1){
      tmp <- do.call("cbind",tapply(names(tmp.cl), tmp.cl, function(x){
        rowMeans(tmp.dat[,x,drop=F])
      }))
      cl.hc = hclust(dist(t(tmp)),method="average")
      return(list(cl=tmp.cl, cl.hc=cl.hc, hc=hc))
    }
    return(NULL)    
  }
 
cut_co_matrix <- function(co.ratio, ord, w=3,th=0.25)
  {
    tmp=co.ratio[ord,ord]
    l = ncol(tmp)
    ##Compute the at every position the moving window
    co.w=sapply(1:(l-w),function(x)rowMeans(tmp[,x:(x+w)]))
    diff = colSums(abs(co.w[,1:(l-w*2-1),drop=F] - co.w[,(w+2):(l-w),drop=F]))
    diff = diff/ colSums(co.w[,1:(l-w*2-1),drop=F] + co.w[,(w+2):(l-w),drop=F])
    diff = c(rep(0,w),diff,rep(0,w))

    peak=slice(Rle(diff), th)
    sep=which.max(peak)
    diff.bin=cut(diff,breaks= c(seq(-0.05,1,length.out=100),2))
    diff.col= jet.colors(100)[diff.bin]
    cl = setNames(as.integer(cut(1:ncol(co.ratio), c(0,sep,ncol(co.ratio)+1))), colnames(co.ratio)[ord])
    return(cl)
  }

reassign_co_ratio <- function(co.ratio, cl, co.stats=NULL,confusion.th=0.4,min.cells=4)
  {
    if(is.null(co.stats)){
      co.stats = get_cl_co_stats(co.ratio, cl)
    }
    correct = 0
    while(1){
      cell.cl.co.ratio = co.stats$cell.cl.co.ratio
      cell.confusion = co.stats$cell.cl.confusion
      cl.confusion = co.stats$cl.confusion
      tmp.dat = cell.cl.co.ratio[names(cl),as.character(sort(unique(cl)))]
      pred.cl <- setNames(colnames(tmp.dat)[apply(tmp.dat, 1, which.max)], row.names(tmp.dat))      
      if(sum(cl==pred.cl) <= correct){
        break
      }
      correct = sum(cl==pred.cl)
      print(correct/length(cl))
      tmp.cells = names(pred.cl)[pred.cl!=cl[names(pred.cl)]]
      cl[tmp.cells]=pred.cl[tmp.cells]
      cl.size = table(cl)
      ###Remove small clusters with high average confusion score, assign cells to other cluster
      cl.small = names(cl.size)[cl.size <=10]
      cl.confusion[cl.small]
      rm.cl = cl.small[cl.confusion[cl.small] > confusion.th| cl.size[cl.small]< min.cells]
      tmp.cells = names(cl)[cl %in% rm.cl]
      tmp.dat = cell.cl.co.ratio[tmp.cells,setdiff(unique(cl),rm.cl),drop=F]
      pred.cl <- setNames(colnames(tmp.dat)[apply(tmp.dat, 1, which.max)], row.names(tmp.dat))
      cl[tmp.cells] = pred.cl[tmp.cells]
      ###Recompute co stats
      co.stats = get_cl_co_stats(co.ratio, cl)
    }
    return(list(cl=cl, co.stats=co.stats))
  }



