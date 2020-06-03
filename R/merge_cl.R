
#' Title
#'
#' @param de.pair 
#' @param de.param 
#' @param merge.type 
#'
#' @return
#' @export
#'
#' @examples
test_merge <- function(de.pair, de.param, merge.type="undirectional")
  {
    if(length(de.pair)==0){
      return(TRUE)
    }
    to.merge = FALSE
    if(merge.type=="undirectional"){
      if(!is.null(de.param$de.score.th)){
        to.merge=de.pair$score < de.param$de.score.th
      }
      if(!to.merge & !is.null(de.param$min.genes)){
        to.merge=de.pair$num < de.param$min.genes
      }
    }
    else{
      if(!is.null(de.param$de.score.th)){
        to.merge=de.pair$up.score < de.param$de.score.th | de.pair$down.score < de.param$de.score.th
      }
      if(!to.merge & !is.null(de.param$min.genes)){
        to.merge=de.pair$up.num < de.param$min.genes | de.pair$down.num < de.param$min.genes
      }
    }
    return(to.merge)
  }


#' Merge clusters based on pairwise differential expressed genes. 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param cl A vector of cluster membership with cell index as names, and cluster id as values. 
#' @param rd.dat Reduced dimensions for cells. Used to determine which clusters are close to each other. Clusters are merged among nearest neighbors first. 
#' @param de.param The DE gene criteria. See de_param for details. 
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately. 
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.   
#' @param de.method Use limma by default. We are still testing "chisq" mode.
#' @param de.genes If not null, use DE genes computated prevoiusly by DE_genes_pw or DE_genes_pairs to avoid recomputing.
#' @param return.markers If TRUE, compute the DE genes between very pairs of clusters as markers
#' @param pairBatch The number of pairs to be tested for merging in one batch. Avoid compairing many pairs at the same time to reduce memory comsumption. Default 40
#' @param sampled For big dataset, norm.dat may not include all cells from cl. If TRUE, norm.dat is the data matrix for downsampled cells, and no need for further down sampling. 
#'
#' @return A list with cl (cluster membership), de.genes (differentially expressed genes), sc (cluster pairwise de.score), markers (top cluster pairwise markers)
#' 
#' @export
#'
merge_cl<- function(norm.dat,
                    cl, 
                    rd.dat=NULL,
                    rd.dat.t = NULL,
                    de.param = de_param(), 
                    merge.type = c("undirectional","directional"), 
                    max.cl.size = 300,
                    de.method = "fast_limma",
                    de.genes = NULL, 
                    return.markers = FALSE,
                    pairBatch =40,
                    verbose = 0)
  {
    if(!is.integer(cl)){
      cl = setNames(as.integer(as.character(cl)), names(cl))
    }
    merge.type=merge.type[1]
    de.df=list()
    pairs=NULL
    if(!is.null(de.genes)){
      pairs=do.call("rbind",strsplit(names(de.genes), "_"))
      row.names(pairs)=names(de.genes)
    }
     ###Merge small clusters with the closest neighbors first.
    if(!is.null(rd.dat)){
      cl.rd = as.data.frame(get_cl_means(rd.dat,cl[names(cl) %in% row.names(rd.dat)]))
    }
    else{
      cl.rd = as.data.frame(get_cl_means(rd.dat.t,cl[names(cl) %in% colnames(rd.dat)]))
    }
    cl.size = table(cl)
    while(TRUE){
      cl.size = table(cl)
      if(length(cl.size)==1){
        break
      }
      cl.small =  names(cl.size)[cl.size < de.param$min.cells]
      if(length(cl.small)==0){
        break
      }
      ##Compute cluster similary on reduced dimension
      if(ncol(cl.rd)>2 & nrow(cl.rd) > 2){
        cl.sim = cor(as.matrix(cl.rd))
      }
      else{
        cl.diff=as.matrix(dist(t(as.matrix(cl.rd))))
        cl.sim = 1 - cl.diff/max(cl.diff)
      }
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
      cl[cl==y]= x
      p = as.character(c(x,y))
      cl.rd[[p[1]]] = get_weighted_means(as.matrix(cl.rd[,p]), as.vector(cl.size[p]))
      cl.rd[[p[2]]] = NULL
      cl.size[p[1]] = sum(cl.size[p])
      cl.size = cl.size[names(cl.size)!=p[2]]
    }		
    cl.means = as.data.frame(get_cl_means(norm.dat, cl))
    cl.present = as.data.frame(get_cl_present(norm.dat, cl,low.th=de.param$low.th))
    cl.sqr.means = as.data.frame(get_cl_sqr_means(norm.dat,cl))

    while(length(unique(cl)) > 1){
      if(length(unique(cl)) == 2){
        merge.pairs = as.data.frame(matrix(as.integer(colnames(cl.rd)), nrow=1))
        merge.pairs$sim = cor(cl.rd[,1], cl.rd[,2])
        row.names(merge.pairs) = paste(merge.pairs[,1], merge.pairs[,2],sep="_")
      }
      else{
        ##Compute cluster similary on reduced dimension
        if(ncol(cl.rd)>2 & nrow(cl.rd) > 2){
          cl.sim = cor(cl.rd)
        }
        else{
          cl.diff=as.matrix(dist(t(cl.rd)))
          cl.sim = 1 - cl.diff/max(cl.diff)
        }
        
        knn.matrix= sim_knn(cl.sim, k=pmin(4, ncol(cl.sim)))
        merge.pairs = do.call("rbind",apply(knn.matrix, 2,function(x)data.frame(c1=row.names(knn.matrix),c2=colnames(cl.sim)[x] , stringsAsFactors=FALSE)))
        merge.pairs = merge.pairs[merge.pairs[,1]!=merge.pairs[,2],]
        merge.pairs$sim = get_pair_matrix(cl.sim, merge.pairs[,1], merge.pairs[,2])
        merge.pairs[,1] = as.integer(merge.pairs[,1])
        merge.pairs[,2] = as.integer(merge.pairs[,2])
        tmp1 =pmin(merge.pairs[,1],merge.pairs[,2])
        tmp2 =pmax(merge.pairs[,1],merge.pairs[,2])
        merge.pairs[,1:2] = cbind(tmp1,tmp2)
        p = paste(merge.pairs[,1],merge.pairs[,2],sep="_")
        merge.pairs= merge.pairs[!duplicated(p),,drop=F]
        row.names(merge.pairs) = p[!duplicated(p)]
        merge.pairs = merge.pairs[order(merge.pairs$sim,decreasing=T),]
      }
     
      ###Determine the de score for these pairs
      if(nrow(merge.pairs)==0){
          break
        }
      
      #####get DE genes for new pairs
      new.pairs = setdiff(row.names(merge.pairs),names(de.genes))
      if(verbose > 0){
        cat("Compute DE genes\n")
      }
      tmp.de.genes =de_selected_pairs(norm.dat, cl=cl, pairs=merge.pairs[new.pairs,], de.param= de.param, method=de.method, cl.means=cl.means, cl.present=cl.present, cl.sqr.means=cl.sqr.means)      
      de.genes[names(tmp.de.genes)] = tmp.de.genes
      pairs = get_pairs(names(de.genes))
      
      tmp.pairs= intersect(names(de.genes), row.names(merge.pairs))
      sc = sapply(de.genes[tmp.pairs], function(x){
        if(length(x)>0){x$score}
        else{0}
      })
      sc = sort(sc)
                                        #print(head(sc,10))      
      to.merge = sapply(names(sc), function(p){
        to.merge = test_merge(de.genes[[p]], de.param, merge.type=merge.type)
      })
      if(sum(to.merge)==0){
        break
      }
      sc = sc[to.merge]
      to.merge= merge.pairs[names(sc),,drop=FALSE]
      to.merge$sc = sc          

        
      merged =c()
      ###The first pair in to.merge always merge. For the remaining pairs, if both clusters have already enough cells,
      ###or independent of previus merging, then they can be directly merged as well, without re-assessing DE genes. 
      for(i in 1:nrow(to.merge)){
        p = c(to.merge[i,1], to.merge[i,2])
        if(i == 1 | sc[i] < de.param$de.score.th /2  & length(intersect(p, merged))==0){
          cl[cl==p[2]] = p[1]

          p = as.character(p)
          if(verbose > 0){
            cat("Merge ",p[1], p[2], to.merge[i,"sc"], to.merge[i, "sim"], cl.size[p[1]],"cells", cl.size[p[2]],"cells", "\n")
          }
          
          cl.rd[[p[1]]] = get_weighted_means(as.matrix(cl.rd[,p]), cl.size[p])
          cl.rd[[p[2]]] = NULL          
          cl.means[[p[1]]] = get_weighted_means(as.matrix(cl.means[, p]), cl.size[p])
          cl.means[[p[2]]] = NULL
          cl.present[[p[1]]] = get_weighted_means(as.matrix(cl.present[, p]), cl.size[p])
          cl.present[[p[2]]] = NULL
          cl.sqr.means[[p[1]]] = get_weighted_means(as.matrix(cl.sqr.means[, p]), cl.size[p])
          cl.sqr.means[[p[2]]] = NULL
          cl.size[p[1]] = sum(cl.size[p])
          cl.size = cl.size[names(cl.size)!=p[2]]
          rm.pairs = row.names(pairs)[pairs[,1]%in% p | pairs[,2]%in% p]
          de.genes = de.genes[setdiff(names(de.genes),rm.pairs)]
          merged = c(merged,p)
        }
      }
    } 
    if(length(unique(cl))<2){
      return(NULL)
    }
    if(verbose > 0){
      print(table(cl))
    }
    markers = NULL
    if(return.markers){
      if(!is.null(max.cl.size)){
        sampled.cells = sample_cells(cl[names(cl) %in% colnames(norm.dat)],  max.cl.size)
        tmp.cl= cl[sampled.cells]
      }
      else{
        tmp.cl= cl
      }
      de.genes = de_all_pairs(norm.dat, cl=tmp.cl, de.genes=de.genes, de.param=de.param, cl.means=cl.means, cl.present=cl.present, cl.sqr.mans=cl.sqr.means)
    }
    markers = select_markers(norm.dat, cl, de.genes=de.genes, n.markers=50)$markers
    sc = sapply(de.genes, function(x){
      if(length(x)>0){x$score}
      else{0}
    })
    return(list(cl=cl, de.genes=de.genes,sc=sc, markers=markers))
  }

