#' Title
#'
#' @param dat.list 
#' @param de.param.list 
#' @param cl 
#' @param pairs 
#' @param cl.means.list 
#' @param cl.present.list 
#' @param lfc.conservation.th 
#' @param de.genes.list 
#' @param max.cl.size 
#'
#' @return
#' @export
#'
#' @examples
de_genes_pairs_multiple <- function(dat.list, de.param.list, cl, pairs, cl.means.list=NULL, cl.present.list=NULL, cl.sqr.means.list=NULL, lfc.conservation.th=0.6, de.genes.list=NULL, max.cl.size=200, method="fast_limma")
  {
      cl.size = table(cl)
      if(is.null(de.genes.list)){
        de.genes.list = sapply(names(dat.list), function(x)list())
      }
      cl.size.platform=list()
      for(x in names(dat.list)){
        norm.dat = dat.list[[x]]
        tmp.cl = cl[names(cl) %in% colnames(norm.dat)]
        cl.size = table(tmp.cl)
        cl.size.platform[[x]]= cl.size
        select.cl = names(cl.size)[cl.size >= de.param.list[[x]]$min.cells]
        if(length(select.cl) < 2){
          next
        }
        tmp.cl = tmp.cl[tmp.cl %in% select.cl]
        if(is.factor(tmp.cl)){
          tmp.cl = droplevels(tmp.cl)
        }
        if(!is.null(max.cl.size)){
          tmp.cells = sample_cells(tmp.cl, max.cl.size)
          tmp.cl = tmp.cl[tmp.cells]
        }
        tmp.pairs= pairs[pairs[,1] %in% select.cl & pairs[,2] %in% select.cl,]
        if(nrow(tmp.pairs)==0){
          next
        }
        tmp = de_selected_pairs(norm.dat, cl= tmp.cl, pairs = tmp.pairs, de.param=de.param.list[[x]], method=method,cl.means=cl.means.list[[x]], cl.present = cl.present.list[[x]], cl.sqr.means=cl.sqr.means.list[[x]])
        de.genes.list[[x]] = c(de.genes.list[[x]], tmp)
      }               
      for(p in row.names(pairs)){
        lfc = sapply(names(cl.means.list), function(x){
          if(pairs[p,1] %in% colnames(cl.means.list[[x]]) & pairs[p,2] %in% colnames(cl.means.list[[x]])){
            cl.means.list[[x]][comb.dat$common.genes,pairs[p,1]] - cl.means.list[[x]][comb.dat$common.genes,pairs[p,2]]
          }
          else{
            NULL
          }
        })
        if(is.list(lfc)){
          lfc = lfc[!sapply(lfc,is.null)]
          lfc = do.call("cbind",lfc)
        }
        row.names(lfc) = comb.dat$common.genes
        sign1 = rowSums(lfc > 1)
        sign2 = rowSums(lfc < -1)
        frac = pmax(sign1, sign2)/ncol(lfc)
        select.genes = names(frac)[frac >= lfc.conservation.th]
        for(x in names(de.genes.list)){
          if(is.null(de.genes.list[[x]][[p]])){
            next
          }
          up.genes = de.genes.list[[x]][[p]]$up.genes
          down.genes = de.genes.list[[x]][[p]]$down.genes
          up.genes = up.genes[names(up.genes) %in% select.genes]
          down.genes = down.genes[names(down.genes) %in% select.genes]
          tmp = up.genes
          tmp[tmp > 20] = 20
          up.score <- sum(tmp)
          tmp = down.genes
          tmp[tmp > 20] = 20
          down.score <- sum(tmp)    
   
          de.genes.list[[x]][[p]]=list(
                              up.score = up.score,
                              down.score = down.score,
                              score = up.score + down.score,
                              up.num = length(up.genes),
                              down.num = length(down.genes),
                              num = length(up.genes) + length(down.genes)
                              )
          
        }
      }
      return(de.genes.list)
    }


#' Title
#'
#' @param cl.rd 
#'
#' @return
#' @export
#'
#' @examples
get_cl_sim <- function(cl.rd)
{
  if(ncol(cl.rd)>2 & nrow(cl.rd) > 2){
    sim=cor(cl.rd)
  }
  else{
    cl.diff=as.matrix(dist(t(cl.rd)))
    sim = 1 - cl.diff/max(cl.diff)
  }
  return(sim)
}

#' Title
#'
#' @param cl.rd.list 
#' @param FUN 
#'
#' @return
#' @export
#'
#' @examples
get_cl_sim_multiple <- function(cl.rd.list, FUN =pmax)
  {
    all.cl = unique(unlist(lapply(cl.rd.list, colnames)))
    cl.sim = matrix(-1, nrow=length(all.cl),ncol=length(all.cl))
    colnames(cl.sim) = row.names(cl.sim) = all.cl
    cl.rd.list = cl.rd.list[!sapply(cl.rd.list, is.null)]
    for(cl.rd in cl.rd.list){
      if(nrow(cl.rd) <= 1){
        next
      }
      tmp.sim = get_cl_sim(cl.rd)
      cl.sim[row.names(tmp.sim),colnames(tmp.sim)] = FUN(tmp.sim, cl.sim[row.names(tmp.sim),colnames(tmp.sim)])
    }
    diag(cl.sim)=1
    return(cl.sim)
  }


    


####Change criteria. If one of the platform shows significant DE genes, and the other platform show consistent fold change, keep the clusters seperate. 

#' Title
#'
#' @param comb.dat 
#' @param merge.dat.list 
#' @param cl 
#' @param anchor.genes 
#' @param verbose 
#' @param pairBatch 
#' @param de.genes.list 
#' @param lfc.conservation.th 
#' @param merge.type 
#'
#' @return
#' @export
#'
#' @examples
merge_cl_multiple <- function(comb.dat, merge.dat.list,  cl, anchor.genes, verbose=TRUE, pairBatch=40, de.genes.list=NULL, lfc.conservation.th=0.7, merge.type="undirectional", de.method="fast_limma")
{
  print("merge_cl_multiple")
  cl = setNames(as.character(cl),names(cl))
  merge_x_y <- function(x, y)
  {
    cl[cl==x]= y
    if(length(unique(cl))==1){
      return(NULL)
    }
    tmp.cells = names(cl)[cl==y]
    tmp = colnames(cl.sim)!=x
    cl.sim = cl.sim[tmp,tmp,drop=F]
    update.sim = setNames(rep(-1, ncol(cl.sim)), colnames(cl.sim))        
    for(set in names(merge.dat.list)){
      cl.rd = cl.rd.list[[set]]
      if(!is.null(cl.rd)){
        tmp = colnames(cl.rd)!=x   
        cl.rd = cl.rd[,tmp,drop=F]
      }    
      tmp.cells2 = intersect(tmp.cells, colnames(merge.dat.list[[set]]))
      if(length(tmp.cells2)==0){
        next
      }
      include.y = length(tmp.cells2) >= merge.de.param.list[[set]]$min.cells
      if(!is.null(cl.means.list)){
        tmp = colnames(cl.means.list[[set]])!=x
        cl.means.list[[set]] = cl.means.list[[set]][,tmp,drop=F]
        tmp.means = Matrix::rowMeans(merge.dat.list[[set]][,tmp.cells2,drop=F])
        if(include.y){
          if(!is.null(cl.means.list[[set]]) & nrow(cl.means.list[[set]])>0){
            cl.means.list[[set]][[y]] = tmp.means[row.names(cl.means.list[[set]])]
          }
          else{
            cl.means.list[[set]] = data.frame(tmp.means)
            colnames(cl.means.list[[set]])=y
          }
        }  
      }
      if(!is.null(cl.sqr.means.list) & de.method=="fast_limma"){
        tmp = colnames(cl.sqr.means.list[[set]])!=x
        cl.sqr.means.list[[set]] = cl.sqr.means.list[[set]][,tmp,drop=F]
        tmp.sqr.means = Matrix::rowMeans(merge.dat.list[[set]][,tmp.cells2,drop=F]^2)        
        if(include.y){
          if(!is.null(cl.sqr.means.list[[set]]) & nrow(cl.sqr.means.list[[set]])>0){
            cl.sqr.means.list[[set]][[y]] = tmp.sqr.means[row.names(cl.sqr.means.list[[set]])]
          }
          else{
            cl.sqr.means.list[[set]] = data.frame(tmp.means)
            colnames(cl.sqr.means.list[[set]])=y
          }
        }        
      }

      if(!is.null(cl.present.list)){
        
        tmp = colnames(cl.present.list[[set]])!=x
        cl.present.list[[set]] = cl.present.list[[set]][,tmp,drop=F]
        tmp.means = Matrix::rowMeans(merge.dat.list[[set]][,tmp.cells2,drop=F] >= merge.de.param.list[[set]]$low.th)
        if(include.y){
          if(!is.null(cl.present.list[[set]])&nrow(cl.present.list[[set]])>0){
            cl.present.list[[set]][[y]] =  tmp.means[row.names(cl.present.list[[set]])]
          }
          else{
            cl.present.list[[set]] = data.frame(tmp.means)
            colnames(cl.present.list[[set]])=y
          }   
        }
      }
      if(include.y){         
        tmp= Matrix::rowMeans(merge.dat.list[[set]][anchor.genes ,tmp.cells2,drop=FALSE])
        if(y %in% colnames(cl.rd)){
          cl.rd[,y]= tmp
        }
        else{
          cl.rd = cbind(cl.rd, tmp)
          colnames(cl.rd)[ncol(cl.rd)] = y
        }
        if(ncol(cl.rd)> 2){
          tmp.sim = cor(cl.rd[,y], cl.rd)
        }
        else{
          cl.diff=as.matrix(dist(t(cl.rd)))
          tmp.sim = (1 - cl.diff/max(cl.diff))[y,]          
        }
        update.sim[colnames(tmp.sim)] = pmax(update.sim[colnames(tmp.sim)],tmp.sim)
      }
      cl.rd.list[[set]]=cl.rd
    }
    cl.sim[y,] = update.sim
    return(list(cl=cl, cl.rd.list=cl.rd.list, cl.sim = cl.sim, cl.means.list=cl.means.list, cl.present.list = cl.present.list))
  }

  add_pairs_de_genes <- function(de.genes.list, cl, new.pairs)
    {
      if(verbose){
        print("Add de genes")
      }
      de.genes.list <- de_genes_pairs_multiple(merge.dat.list, merge.de.param.list, cl, pairs=new.pairs, cl.means.list=cl.means.list, cl.present.list=cl.present.list, cl.sqr.means.list= cl.sqr.means.list,lfc.conservation.th=lfc.conservation.th, de.genes.list=de.genes.list,method=de.method)
      if(verbose){
        print("Finish adding de genes")
      }
      return(de.genes.list)
    }
  
  rm_pairs_de_genes <- function(de.genes.list, rm.pairs)
    {
      for(x in names(de.genes.list)){
        de.genes.list[[x]] = de.genes.list[[x]][setdiff(names(de.genes.list[[x]]), rm.pairs)]
      }
      return(de.genes.list)
    }
  
  test_merge_multiple <-  function(de.genes.list, merge.type="undirectional")
    {
      merge.pairs <- unique(unlist(lapply(de.genes.list, names)))
      to.merge.df <- do.call("rbind", lapply(names(de.genes.list), function(x){
        to.merge <- sapply(merge.pairs, function(p){
          test_merge(de.genes.list[[x]][[p]], merge.de.param.list[[x]], merge.type=merge.type)
        })
        sc <- sapply(merge.pairs, function(p){de.genes.list[[x]][[p]]$score})
        sc[sapply(sc, is.null)] = 0
        sc = unlist(sc)
        df=data.frame(to.merge, sc, pair = merge.pairs,stringsAsFactors=FALSE)                     
      }))
      not.merged = with(to.merge.df, tapply(!to.merge, pair, sum))
      not.merged.pair = names(not.merged[not.merged > 0])
      to.merge.df = to.merge.df[!to.merge.df$pair %in% not.merged.pair,,drop=F]
      to.merge.sc =with(to.merge.df, tapply(sc, pair, max))
      return(sort(to.merge.sc))
    }
  
  
  merge.sets=names(merge.dat.list)
  merge.de.param.list = comb.dat$de.param.list[merge.sets]
  de.score.th = mean(sapply(merge.de.param.list, function(x)x$de.score.th))
  cl.platform.counts = table(comb.dat$meta.df[names(cl), "platform"],cl)[merge.sets,,drop=F]
  tmp = table(comb.dat$meta.df$platform)[merge.sets]
  cl.platform = cl.platform.counts / as.vector(tmp) 
  cl.platform = t(t(cl.platform) / colSums(cl.platform))  
  cl.min.cells = sapply(merge.de.param.list, function(x)x$min.cells)
  cl.big= cl.platform.counts >= cl.min.cells[rownames(cl.platform.counts)]
  cl.small = colnames(cl.big)[colSums(cl.big) == 0]
  cl.big =  colnames(cl.big)[colSums(cl.big) > 0]
  
  if(length(cl.big)==0){
    return(NULL)
  }
  #cl.rd.list = get_cl_means_list(merge.dat.list, merge.de.param.list, select.genes=anchor.genes, cl=cl)
  cl.rd.list = get_cl_means_list(merge.dat.list, cl=cl, select.genes=anchor.genes, de.param.list = merge.de.param.list)
  
  pairs=NULL
  ###Merge small clusters first
  cl.sim = get_cl_sim_multiple(cl.rd.list)
  while(length(cl.small)>0){
    knn = data.frame(cl=cl.small, nn=cl.big[sim_knn(cl.sim[cl.small, cl.big,drop=F],k=1)],stringsAsFactors=FALSE)
    knn$sim = get_pair_matrix(cl.sim, knn$cl, knn$nn)
    closest.pair = which.max(knn$sim)
    x = knn[closest.pair,1]
    y=  knn[closest.pair,2]
    if(verbose > 0){
        cat("Merge: ", x,y, "sim:", knn[closest.pair,3],"\n")
      }
    update.result=merge_x_y(x, y)
    if(is.null(update.result)){
      return(NULL)
    }
    cl = update.result$cl
    cl.rd.list = update.result$cl.rd.list
    cl.sim = update.result$cl.sim
    cl.means.list = update.result$cl.means.list
    cl.present.list = update.result$cl.present.list
    cl.sqr.means.list = update.result$cl.sqr.means.list
    cl.small = cl.small[cl.small!=x]
  }
  merge.de.param.list = comb.dat$de.param.list[merge.sets]

  cl.means.list = get_cl_means_list(merge.dat.list,  cl=cl, de.param.list=merge.de.param.list)
  cl.means.list = sapply(cl.means.list, as.data.frame, simplify=F)

  cl.sqr.means.list = get_cl_sqr_means_list(merge.dat.list, cl=cl, de.param.list=merge.de.param.list)
  cl.sqr.means.list = sapply(cl.sqr.means.list, as.data.frame, simplify=F)
  
  cl.present.list = get_cl_present_list(merge.dat.list, cl=cl, de.param.list=merge.de.param.list)
  cl.present.list = sapply(cl.present.list, as.data.frame, simplify=F)

  
  de.pairs = NULL
  de.genes.list = sapply(names(merge.dat.list), function(x)list(),simplify=F)
  while (length(unique(cl)) > 1) {
###Find pairs of nearest neighbrs as candidates for merging.
    k.tmp = pmin(4,ncol(cl.sim))
    nn=colnames(cl.sim)[sim_knn(cl.sim, k= k.tmp)]
    merge.pairs = data.frame(cl=rep(row.names(cl.sim), length(k.tmp)), nn=nn,stringsAsFactors=FALSE)
    merge.pairs = merge.pairs[merge.pairs[,1]!=merge.pairs[,2],]
    merge.pairs$sim = get_pair_matrix(cl.sim, merge.pairs$cl, merge.pairs$nn)
    
    tmp1 = pmin(merge.pairs[, 1], merge.pairs[, 2])
    tmp2 = pmax(merge.pairs[, 1], merge.pairs[, 2])
    merge.pairs[, 1:2] = cbind(tmp1, tmp2)
    
    p = paste(merge.pairs[, 1], merge.pairs[, 2], sep = "_")
    merge.pairs = merge.pairs[!duplicated(p), , drop = F]
    row.names(merge.pairs) = p[!duplicated(p)]
    merge.pairs = merge.pairs[order(merge.pairs$sim, decreasing = T), ,drop=F]
    merge.pairs = merge.pairs[!row.names(merge.pairs) %in% row.names(de.pairs),,drop=F]
    while(nrow(merge.pairs) > 0){
      new.pairs = head(row.names(merge.pairs),pairBatch)
      if(is.null(de.pairs)){
        de.pairs = merge.pairs[new.pairs,,drop=F]
      }else{       
        de.pairs = rbind(merge.pairs[new.pairs,],de.pairs)
      }
      de.genes.list = add_pairs_de_genes(de.genes.list, cl, merge.pairs[new.pairs,])
      if(sum(sapply(de.genes.list, length))==0){
        return(NULL)
      }
      merge.sc = test_merge_multiple(de.genes.list, merge.type=merge.type)
      if(length(merge.sc)>0){
        break
      }
      merge.pairs = merge.pairs[!row.names(merge.pairs)%in% new.pairs,]
    }
    if(length(merge.sc)==0){
      break
    }
    merged = c()   
    for(i in 1:length(merge.sc)){
      tmp = names(merge.sc[i])
      if(!tmp %in% row.names(de.pairs)){
        next
      }       
      p = de.pairs[tmp,,]
      x = p[1,1]
      y = p[1,2]
      sim = p[,3]
      if (i == 1 | merge.sc[i] < de.score.th/2 & sum(c(x,y) %in% merged) == 0){              
        if (verbose > 0) {
          cat("Merge ", x, y, merge.sc[i], sim, sum(cl == x),  "cells", sum(cl == y), "cells", "\n")
        }         
        update.result=merge_x_y(x=x, y=y)
        if(is.null(update.result)){
          return(NULL)
        }
        cl = update.result$cl
        cl.rd.list = update.result$cl.rd.list
        cl.sim = update.result$cl.sim
        cl.means.list = update.result$cl.means.list
        cl.present.list = update.result$cl.present.list
        cl.sqr.means.list = update.result$cl.sqr.means.list
        rm.pairs = de.pairs[, 1] %in% c(x,y) | de.pairs[, 2] %in% c(x,y)
        de.genes.list = rm_pairs_de_genes(de.genes.list, row.names(de.pairs)[rm.pairs])
        de.pairs = de.pairs[!rm.pairs,,drop=F]        
        merged = c(merged, c(x,y))
      }
    }
  }
  if (length(unique(cl)) < 2) {
    return(NULL)
  }
  if (verbose > 0) {
    print(table(cl))
  }
  return(cl)
}

