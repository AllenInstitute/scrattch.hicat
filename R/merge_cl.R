
#' Merge clusters based on pairwise differential expressed genes. 
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1)
#' @param cl A vector of cluster membership with cells as names, and cluster id as values. 
#' @param rd.dat Reduced dimensions for cells. Used to determine which clusters are close to each other. Clusters are merged among nearest neighbors first. 
#' @param de.param The DE gene criteria. See de_param for details. 
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately. 
#' @param max.cl.size Sampled cluster size. This is to speed up limma DE gene calculation. Instead of using all cells, we randomly sampled max.cl.size number of cells for testing DE genes.   
#' @param de.method Use limma by default. We are still testing "chisq" mode. 
#' @return A list with cl (cluster membership), de.genes (differentially expressed genes), sc (cluster pairwise de.score), markers (top cluster pairwise markers)
#' @export
#'

merge_cl <- function(norm.dat, 
                    cl, 
                    rd.dat, 
                    de.param = de_param(), 
                    merge.type = c("undirectional","directional"), 
                    max.cl.size = 300,
                    de.method = "limma",
                    de.genes = NULL, 
                    return.markers = TRUE, 
                    verbose = 0)
  {
    merge.type=merge.type[1]
    print(merge.type)
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
      tmp.cells = sample_cells(cl,  max.cl.size,weights=cell.weights)
      tmp.dat = as.matrix(norm.dat[,tmp.cells])
      ###Merge small clusters with the closest neighbors first.
      cl.rd = Matrix::t(get_cl_means(Matrix::t(rd.dat),cl))
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
        else if(merge.type=="undirectional"){
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

