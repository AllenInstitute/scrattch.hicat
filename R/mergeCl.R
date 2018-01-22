###Merge multiple pairs are the same type if their merging are independent. 
mergeCl<- function(norm.dat, cl, rd.dat, min.padj=100, type="undirectional", min.cells=3, max.cl.size=300, rm.gene.mod=NULL,rm.eigen=NULL,rm.th=0.7, ...)
  {
    cl = setNames(as.integer(as.character(cl)), names(cl))
    de.genes=list()
    de.df=list()
    select.cells = names(cl)
    ###Remove genes correlate with rm.eigen
    if(!is.null(rm.eigen)){
      tmp.cells=sample(select.cells, min(length(select.cells), 4000))
      rm.kME=cor(t(as.matrix(norm.dat[,tmp.cells])), rm.eigen[tmp.cells,])
      rm.kME[is.na(rm.kME)] =0
      rm.score=setNames(rowMaxs(abs(rm.kME)),row.names(rm.kME))
      select = rm.score < rm.th 
      norm.dat = norm.dat[select,]
    }
    pairs=NULL
    while(length(unique(cl)) > 1){
      ###Merge small clusters with the closest neighbors first.
      ###Compute cluster-sums in reduced dimension.
      ###cluster-sum is used to compute cluster-means, save some time for re-computation after clusters are merged. 
      cl.rd <- do.call("rbind", tapply(names(cl),cl, function(x){
        colSums(rd.dat[x,,drop=F])
      },simplify=F))
      while(TRUE){
        cl.size = table(cl)
        cl.small = names(cl.size)[cl.size < min.cells]
        ###Merge small clusters with its closest neighbors.
        if(length(cl.small)>0){      
          ##Compute cluster similary on reduced dimension
          if(ncol(cl.rd)>2){
            cl.sim = cor(t(cl.rd))
          }
          else{
            cl.diff=as.matrix(dist(cl.rd/as.vector(cl.size[row.names(cl.rd)])))
            cl.sim = 1 - cl.diff/max(cl.diff)
          }
          tmp=as.data.frame(as.table(cl.sim[cl.small,,drop=F]))
          tmp[,1]=as.integer(as.character(tmp[,1]))
          tmp[,2]=as.integer(as.character(tmp[,2]))
          tmp = tmp[tmp[,1]!=tmp[,2],,drop=F]
          closest.pair = which.max(tmp$Freq)
          x = tmp[closest.pair,1]
          y=  tmp[closest.pair,2]
          cat("Merge: ", x,y, "sim:", tmp[closest.pair,3],"\n")
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
      ##Down sample cells for efficiency   
      if(!is.null(max.cl.size)){
        tmp.cells = unlist(tapply(names(cl),cl, function(x){
          if(length(x)>max.cl.size){
            x= sample(x, max.cl.size)
          }
          x
        },simplify=FALSE))
      }
      else{
        tmp.cells= names(cl)
      }
      tmp.dat = as.matrix(norm.dat[,tmp.cells])
      
      ####if the number of clusters is fewer than 10, compute all pairs. 
      if(length(cl.size)<10 & length(de.genes)==0){
        tmp.result=deScore(tmp.dat, cl=cl[tmp.cells], min.cells=min.cells, ...)
        de.genes=tmp.result$de.genes
        de.df = tmp.result$de.df
        sc = sapply(de.genes, function(x){
          if(length(x)>0){x$score}
          else{0}
        })
        merge.pairs=do.call("rbind",strsplit(names(de.genes), "_"))
        merge.pairs = gsub("cl","", merge.pairs)
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
        tmp.result=deScore.pairs(tmp.dat, cl=cl[tmp.cells], pairs=merge.pairs[new.pairs,,drop=F], min.cells=min.cells, ...)
        tmp.de.df = tmp.result[[1]]
        tmp.de.genes = tmp.result[[2]]
        de.genes[names(tmp.de.genes)] = tmp.de.genes
        de.df[names(tmp.de.df)] = tmp.de.df
        sc = sapply(de.genes[row.names(merge.pairs)], function(x){
          if(length(x)>0){x$score}
          else{0}
        })
      }
      sc = sort(sc)
      print(sc)

      to.merge = sapply(names(sc), function(p){
        x = de.genes[[p]]
        if(length(x)==0){
          to.merge = TRUE
        }
        else if(type=="undirectional"){
          to.merge=x$score < min.padj 
        }
        else{
          to.merge=x$up.score < min.padj | x$down.score < min.padj 
        }
        to.merge
      })
      if(sum(to.merge)==0){
        break
      }
      sc = sc[to.merge]
      to.merge= merge.pairs[names(sc),,drop=FALSE]
      merged =c()
      ###The first pair in to.merge always merge. 
      ###For the remaining pairs that are have score < min.adpj/2, and independent of previus merging, they can also be merged directly, without re-assessing DE genes. 
      for(i in 1:nrow(to.merge)){
        p = to.merge[i,]
        if(i == 1 | sc[i] < min.padj/2  & length(intersect(p, merged))==0){
          cat("Merge ",p[1], p[2], sc[i],"\n")
          
          cl[cl==p[2]] = p[1]
          rm.pairs = row.names(pairs)[pairs[,1]%in% p | pairs[,2]%in% p]
          de.genes = de.genes[setdiff(names(de.genes),rm.pairs)]
        }
        merged = c(merged, p)
      }
      pairs = pairs[names(de.genes),,drop=F]
    }
    return(list(cl=cl, de.genes=de.genes,sc=sc))
  }



