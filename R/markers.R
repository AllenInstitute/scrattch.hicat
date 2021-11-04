#' Title
#'
#' @param norm.dat 
#' @param cl 
#' @param n.markers 
#' @param de.genes 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
select_markers <- function(norm.dat, cl, n.markers=20,de.genes=NULL, mc.cores=1,...)
  {
    require(parallel)
    if(is.null(de.genes)){
      de.genes=de_all_pairs(norm.dat, cl, mc.cores=mc.cores,...)
    }
    pairs.df = get_pairs(names(de.genes))
    select.pairs = row.names(pairs.df)[pairs.df[,1] %in% cl & pairs.df[,2]%in% cl]
             
    de.markers = parallel::pvec(select.pairs, function(s){
      sapply(s, function(x){
        tmp = de.genes[[x]]
        select.genes = c(head(names(tmp$up.genes),n.markers), head(names(tmp$down.genes),n.markers))
      },simplify=F)
    },mc.cores=mc.cores)
    if(!is.null(norm.dat)){
      markers = intersect(unlist(de.markers),row.names(norm.dat))
    }
    else{
      markers = unique(unlist(de.markers))
    }
    return(list(markers=markers, de.genes=de.genes))
  }

#' Title
#'
#' @param de.genes 
#' @param top.n 
#' @param max.num 
#' @param bin.th 
#'
#' @return
#' @export
#'
#' @examples
get_gene_score <- function(de.genes,cl.means=NULL, all.genes=NULL, top.n=50, max.num=1000,bin.th=4, mc.cores=1)
  {
    require(Matrix)    
    require(parallel)
    
    if(is.null(all.genes)){
      all.genes <- parallel::pvec(names(de.genes), function(x){
        genes=lapply(x, function(p){
          de = de.genes[[p]]
          pair = strsplit(p, "_")[[1]]
          x = pair[[1]]
          y = pair[[2]]
          up.genes = names(de$up.genes)
          down.genes = names(de$down.genes)
                                        #Deal with pairs with too many DEX genes, include all binary markers
          up.binary.genes = down.binary.genes=NULL
          if(!is.null(cl.means)){
            lfc = cl.means[,x] - cl.means[,y]
            up.binary.genes = up.genes[lfc[up.genes] > bin.th]
            down.binary.genes = down.genes[lfc[down.genes] < -bin.th]            
          }
          up.genes=up.genes[up.genes %in% c(head(up.genes, top.n), up.binary.genes)]            
          down.genes=down.genes[down.genes %in% c(head(down.genes, top.n), down.binary.genes)]
          list(up=up.genes, down=down.genes)
        })
      },mc.cores=mc.cores)
      all.genes=unique(unlist(all.genes))
    } 
    up.gene.score = parallel::pvec(names(de.genes),function(x){
      mat=sapply(x, function(p){
        de = de.genes[[p]]
        tmp= match(all.genes, head(names(de$up.genes),max.num))
        tmp[is.na(tmp)]=max.num
        tmp = max.num - tmp
        tmp
      })
      mat = Matrix(mat, sparse=TRUE)
      row.names(mat) = all.genes
      mat
    }, mc.cores=mc.cores)
    if(is.list(up.gene.score)){
      up.gene.score = do.call("cbind",up.gene.score)
    }
    
    down.gene.score = parallel::pvec(names(de.genes),function(x){
      mat=sapply(x, function(p){
        de = de.genes[[p]]
        tmp= match(all.genes, head(names(de$down.genes),max.num))
        tmp[is.na(tmp)]=max.num
        tmp = max.num - tmp
        tmp[tmp < 0] = 0
        tmp
      })
      mat = Matrix(mat, sparse=TRUE)
      row.names(mat) = all.genes
      mat
    }, mc.cores=mc.cores)
    if(is.list(down.gene.score)){
      down.gene.score = do.call("cbind",down.gene.score)
    }
      
    row.names(up.gene.score)=row.names(down.gene.score)= all.genes
    return(list(up.gene.score=up.gene.score, down.gene.score=down.gene.score))
  }



#' Title
#'
#' @param de.genes 
#' @param add.genes 
#' @param gene.score 
#' @param rm.genes 
#' @param top.n 
#' @param max.num 
#'
#' @return
#' @export
#'
#' @examples
select_markers_pair <- function(de.genes, add.genes, cl.means, gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    pairs =do.call("rbind",strsplit(gsub("cl","",names(add.genes)), "_"))
    row.names(pairs)= names(add.genes)
    de.genes.list=list()
    if(is.null(gene.score)){
      de.genes= de.genes[names(add.genes)]
      tmp=get_gene_score(de.genes,cl.means=cl.means, top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score      
      gene.score=pmax(up.gene.score, down.gene.score)
      row.names(gene.score)=row.names(up.gene.score)
    }
    select.genes = setdiff(row.names(gene.score),rm.genes)
    gene.score= gene.score[,names(add.genes),drop=F]    
    while(sum(add.genes)>0 & nrow(gene.score)>0){
      g = order(get_row_means(gene.score, select.genes, ),decreasing=T)[1]
      g.n = row.names(gene.score)[g]
      select.pair = colnames(gene.score)[gene.score[g, ]< top.n]
      add.genes[select.pair]= add.genes[select.pair]-1
      for(p in select.pair){
        de.genes.list[[p]]=union(de.genes.list[[p]],g.n)
      }
      add.genes = add.genes[add.genes > 0]
      
      gene.score=gene.score[,names(add.genes),drop=F]
      select.genes = setdiff(select.genes, g)
    }
    return(de.genes.list)
  }


#' Title
#'
#' @param de.genes 
#' @param add.up 
#' @param add.down 
#' @param up.gene.score 
#' @param down.gene.score 
#' @param rm.genes 
#' @param top.n 
#' @param max.num 
#'
#' @return
#' @export
#'
#' @examples
select_markers_pair_direction <- function(de.genes, add.up,add.down,cl.means, up.gene.score=NULL,down.gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    up.genes = down.genes=NULL
    final.genes=c()
    if(is.null(up.gene.score)){
      pairs.n=union(names(add.up), names(add.down))
      pairs =do.call("rbind",strsplit(gsub("cl","",pairs.n), "_"))
      row.names(pairs)= pairs.n
      de.genes= de.genes[pairs.n]
      tmp=get_gene_score(de.genes,cl.means, top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score      
    }
    select.genes = setdiff(row.names(up.gene.score),rm.genes)
    while((sum(add.up)+sum(add.down)>0) & length(select.genes)>0){
      sc =0
      if(length(add.up)>0){
        sc = get_row_sums(up.gene.score, select.genes, names(add.up))
      }
      if(length(add.down)>0){
        sc = sc+ get_row_sums(down.gene.score, select.genes, names(add.down))
      }
      sc = sc[sc > 0]
      ##No more genes to choose from
      if(length(sc)==0){
        break
      }
      g = names(sc)[ which.max(sc)]
      select.up.pair = names(add.up)[up.gene.score[g,names(add.up)] > 0]
      select.down.pair = names(add.down)[down.gene.score[g,names(add.down)] > 0]
      if(length(select.up.pair) + length(select.down.pair)==0){
        break
      }
      add.up[select.up.pair]= add.up[select.up.pair]-1
      add.down[select.down.pair]= add.down[select.down.pair]-1
      add.up = add.up[add.up >0, drop=F]
      add.down = add.down[add.down >0, drop=F]

      tmp = data.frame(gene=rep(g, length(select.up.pair)), pair = select.up.pair)
      up.genes = rbind(up.genes, tmp)
      tmp = data.frame(gene=rep(g, length(select.down.pair)), pair = select.down.pair)
      down.genes = rbind(down.genes, tmp)      
      final.genes=c(final.genes, g)
      select.genes = setdiff(select.genes, final.genes)
      cat(g, length(add.up), length(add.down),length(select.genes),"\n")
    }
    markers= final.genes
    return(list(markers=markers, up.genes=up.genes, down.genes=down.genes))
  }

#' Title
#'
#' @param cl 
#' @param g1 
#' @param g2 
#' @param de.genes 
#' @param top.n 
#' @param max.num 
#' @param n.markers 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
# Always directory
select_markers_pair_group_top<- function(cl, g1,g2,de.genes,cl.means, top.n=50,max.num=1000,n.markers=20,up.gene.score=NULL, down.gene.score=NULL)
{
  require(matrixStats)
  pairs = do.call("rbind",strsplit(names(de.genes), "_"))
  pairs = gsub("cl", "",pairs)
  row.names(pairs)= names(de.genes)
  up.pairs = row.names(pairs)[pairs[,1] %in% g1 & pairs[,2] %in% g2]
  down.pairs = row.names(pairs)[pairs[,1] %in% g2 & pairs[,2] %in% g1]
  select.pairs = c(up.pairs, down.pairs)
  if(is.null(up.gene.score)){
    tmp=get_gene_score(de.genes[select.pairs],cl.means=cl.means, top.n=top.n, max.num=max.num,bin.th=4)
    up.gene.score=tmp$up.gene.score
    down.gene.score=tmp$down.gene.score
  }
  all.genes = row.names(up.gene.score)
  tmp.up.gene.score = cbind(up.gene.score[,up.pairs,drop=F], down.gene.score[,down.pairs,drop=F])
  tmp.down.gene.score = cbind(down.gene.score[,up.pairs,drop=F], up.gene.score[,down.pairs,drop=F])
  tmp.up.gene.score.total = Matrix::rowSums(tmp.up.gene.score)
  tmp.down.gene.score.total = Matrix::rowSums(tmp.down.gene.score)
  
  tmp.up.gene.score.total = sort(tmp.up.gene.score.total[tmp.up.gene.score.total>0],decreasing=TRUE)
  tmp.down.gene.score.total = sort(tmp.down.gene.score.total[tmp.down.gene.score.total>0],decreasing=TRUE)

  up.genes = names(head(tmp.up.gene.score.total, n.markers))
  down.genes = names(head(tmp.down.gene.score.total,n.markers))

  up.num = Matrix::colSums(tmp.up.gene.score[up.genes, , drop = F] > 0)  
  down.num = Matrix::colSums(tmp.down.gene.score[down.genes, , drop = F] > 0)   
  genes = union(up.genes, down.genes)
  select = !row.names(tmp.up.gene.score) %in% genes
  tmp.up.gene.score = tmp.up.gene.score[select,,drop=F]
  tmp.down.gene.score = tmp.down.gene.score[select,,drop=F]

  return(list(up.genes=up.genes, down.genes=down.genes, up.num=up.num, down.num=down.num, up.gene.score=tmp.up.gene.score, down.gene.score=tmp.down.gene.score, up.gene.score.total=tmp.up.gene.score.total, down.gene.score.total = tmp.down.gene.score.total))
}

select_markers_pair_group<- function(cl, g1,g2,de.genes,cl.means, top.n=50,max.num=1000,n.markers=20,up.gene.score=NULL, down.gene.score=NULL)
  {
    result=select_markers_pair_group_top(cl, g1,g2,de.genes,cl.means, top.n=50,max.num=1000,n.markers=n.markers,up.gene.score, down.gene.score)
    up.gene.score= result$up.gene.score
    down.gene.score=result$down.gene.score
    up.num = result$up.num
    down.num  = result$down.num
    up.genes = result$up.genes
    down.genes = result$down.genes
    genes = c(up.genes, down.genes)
    
    add.up = setNames(rep(n.markers, length(up.num)),names(up.num)) - up.num
    add.down = setNames(rep(n.markers, length(down.num)),names(down.num)) - up.num

    if(sum(add.up) + sum(add.down) > 0){    
      tmp = select_markers_pair_direction(add.up,add.down,de.genes=de.genes, cl.means=cl.means,up.gene.score=up.gene.score,down.gene.score=down.gene.score)
      genes = c(genes, tmp$markers)      
    }
    return(list(genes=genes, up.genes=up.genes, down.genes=down.genes))
  }



#' Title
#'
#' @param de.genes 
#' @param up.gene.score 
#' @param down.gene.score 
#' @param default.markers 
#' @param pair.num 
#' @param add.up 
#' @param add.down 
#' @param rm.genes 
#' @param pairs 
#'
#' @return
#' @export
#'
#' @examples
select_N_markers <- function(de.genes, cl.means, up.gene.score=NULL, down.gene.score=NULL, default.markers=NULL, pair.num =1, add.up=pair.num, add.down=pair.num, rm.genes=NULL, pairs=names(de.genes), mc.cores=20)
  {
   add.up = setNames(rep(add.up, length(pairs)), pairs)
   add.down= setNames(rep(add.down, length(pairs)), pairs)
   up.default = down.default = c()
   if(!is.null(default.markers)){
     tmp= intersect(default.markers, row.names(up.gene.score))
     up.default=up.gene.score[tmp,]
     up.num = colSums(up.default)[pairs]     
     down.default = down.gene.score[tmp,]
     down.num = colSums(down.default)[pairs]
     add.up = pmax(add.up -  up.num,0)
     add.down = pmax(add.down -  down.num,0)
   }
   add.up = add.up[add.up>0, drop=F]
   add.down = add.down[add.down>0, drop=F]
   result = select_markers_pair_direction(add.up,add.down,de.genes=de.genes, cl.means=cl.means,up.gene.score=up.gene.score,down.gene.score=down.gene.score,rm.genes=c(rm.genes,default.markers),top.n=50,max.num=2000)
   up.genes = result$up.genes
   down.genes = result$down.genes
   if(length(default.markers)>0){
     up.default = as(up.default,"TsparseMatrix")
     up.default = data.frame(gene = row.names(up.default)[up.default@i+1], pair= colnames(up.default)[up.default@j+1])
     up.genes = rbind(up.default,up.genes)
     
     down.default = as(down.default,"TsparseMatrix")
     down.default = data.frame(gene = row.names(down.default)[down.default@i+1], pair= colnames(down.default)[down.default@j+1])
     down.genes = rbind(down.default,down.genes)
   }
   up.genes = split(up.genes$gene, up.genes$pair)
   down.genes = split(down.genes$gene, down.genes$pair)

   return(list(up.genes=up.genes, down.genes=down.genes, markers=c(default.markers, result$markers)))
 }



#' Title
#'
#' @param de.genes 
#' @param cl 
#' @param n.markers 
#' @param default.markers 
#' @param rm.genes 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
select_pos_markers <- function(de.genes, cl, select.cl=unique(cl), cl.means=NULL, n.markers=3, default.markers=NULL, rm.genes=NULL, up.gene.score=NULL, down.gene.score=NULL,mc.cores=1)
  {
    library(parallel)
    if(!is.null(de.genes)){
      pairs.df = get_pairs(names(de.genes))
    }
    else if(!is.null(up.gene.score)){
      pairs.df = get_pairs(colnames(up.gene.score))
    }else{
      stop("both de.genes and up.gene.score are null")
    }

    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types
    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      print(x)
      up.pairs = row.names(pairs.df)[pairs.df[,1] == x]
      down.pairs = row.names(pairs.df)[pairs.df[,2] == x]
      if(is.null(up.gene.score) | is.null(down.gene.score)){
        tmp=get_gene_score(de.genes[c(up.pairs,down.pairs)], cl.means=cl.means)
        up.gene.score=tmp$up.gene.score
        down.gene.score=tmp$down.gene.score
      }
      add.up = setNames(rep(n.markers, length(up.pairs)), up.pairs)
      add.down = setNames(rep(n.markers, length(down.pairs)), down.pairs)
      
      up.default = sapply(up.pairs, function(p){intersect(names(de.genes[[p]]$up.genes), default.markers)},simplify=F)
      down.default = sapply(down.pairs, function(p){intersect(names(de.genes[[p]]$down.genes), default.markers)},simplify=F)
      if(length(add.up) > 0 & length(up.default)>0){
        add.up = pmax(add.up -  sapply(up.default, length),0)
      }
      if(length(add.down) > 0 & length(down.default)>0){
        add.down = pmax(add.down -  sapply(down.default, length),0)
      }
      tmp.result = select_markers_pair_direction(add.up=add.up, add.down=add.down,de.genes=de.genes, cl.means=cl.means, up.gene.score=up.gene.score,down.gene.score=down.gene.score,rm.genes=c(rm.genes,default.markers))
      result=list(unique(c(tmp.result$markers, unlist(up.default), unlist(down.default))))
      names(result)=x
      result
    }
    return(cl.markers)
  }


select_top_pos_markers <- function(de.genes, cl, n.markers=3, up.gene.score, down.gene.score,mc.cores=10)
  {
    library(parallel)
    pairs.df = get_pairs(names(de.genes))    
    
    ###for each cluster, find markers that discriminate it from other types
    cl.markers <- parallel::pvec(levels(cl), function(x){
      sapply(x, function(tmp.cl){
        up.pairs = row.names(pairs.df)[pairs.df[,1] == tmp.cl]
        down.pairs = row.names(pairs.df)[pairs.df[,2] == tmp.cl]
        if(length(up.pairs)>0){
          sc = get_row_sums(up.gene.score, select.col=up.pairs)
        }
        else{
          sc = 0
        }
        if(length(down.pairs)>0){
          sc = sc+ get_row_sums(down.gene.score, select.col=down.pairs)
        }
        head(names(sort(sc,decreasing=T)), n.markers)        
      },simplify=F)
    },mc.cores=mc.cores)
    return(cl.markers)
  }


#' Title
#'
#' @param cl.dat 
#' @param th 
#' @param tau.th 
#' @param n.markers 
#'
#' @return
#' @export
#'
#' @examples
markers_max_tau <- function(cl.dat, th=0.5, tau.th=0.7,n.markers=3)
  {
    tau = calc_tau(cl.dat)
    tmp=split(row.names(cl.dat), colnames(cl.dat)[apply(cl.dat,1, which.max)])
    tau.df = do.call("rbind", sapply(names(tmp),function(x){
      g = tmp[[x]]
      g = g[cl.dat[g,x] > th]
      g = g[order(tau[g],decreasing=T)]
      if(length(g)>0){
        df=data.frame(g=g, cl=x, cl.dat=cl.dat[g, x],tau=tau[g])
        df = df[df$tau>tau.th,]
        head(df,n=n.markers)
      }
      else{
        NULL
      }
    },simplify=F))
  }

#' Title
#'
#' @param cl 
#' @param cl.present.counts 
#' @param present.th 
#' @param tau.th 
#' @param top.n 
#'
#' @return
#' @export
#'
#' @examples
markers_tau_one_vs_other <- function(cl, cl.present.counts, present.th=0.4, tau.th=0.8,top.n=10)
  {
    all.g = row.names(cl.present.counts)
    all.g = setdiff(all.g,c(grep("LOC", all.g,value=T),grep("Rik$",all.g, value=T)))
    cl.size=table(cl)[colnames(cl.present.counts)]
    cl.present.prob = t(t(cl.present.counts)/as.vector(cl.size))
    tau.genes = sapply(colnames(cl.present.counts),function(x){
      fg = cl.present.counts[,x]/sum(cl==x)
      bg = (rowSums(cl.present.counts) - cl.present.counts[,x])/ (length(cl) - sum(cl==x))
      tau = (fg - bg)/pmax(bg,fg)
      g= all.g[cl.present.prob[all.g,x]>0.5]
      g = g[order(tau[g],decreasing=T)]
      g = g[tau[g]>0.97 | (tau[g]>0.7 & g %in% head(g,top.n))]
      #paste(g,collapse=" ")
    },simplify=F)
  }




#' Title
#'
#' @param cl.g 
#' @param norm.dat 
#' @param cl 
#' @param de.param 
#' @param n.markers 
#' @param cl.present.counts 
#'
#' @return
#' @export
#'
#' @examples
group_specific_markers <- function(cl.g, 
                                   norm.dat, 
                                   cl, 
                                   de.param, 
                                   n.markers = 5, 
                                   cl.present.counts = NULL) {
  cl <- droplevels(cl)
  select.cells <- names(cl)[cl %in% cl.g]
  not.select.cells <- setdiff(names(cl), select.cells)
  
  if(is.null(cl.present.counts)){
    fg <- Matrix::rowSums(norm.dat[, select.cells] > de.param$low.th)
    bg <- Matrix::rowSums(norm.dat[, not.select.cells] > de.param$low.th)
  } else{
    fg <- Matrix::rowSums(cl.present.counts[, cl.g, drop = F])
    bg <- Matrix::rowSums(cl.present.counts[, levels(cl), drop = F]) - fg
  }
  
  bg.freq <- bg / length(not.select.cells)
  fg.freq <- fg / length(select.cells)
  
  tau <- (fg.freq - bg.freq) / pmax(bg.freq, fg.freq)
  ratio <- fg / (fg + bg)
  
  stats <- vec_chisq_test(fg, 
                          rep(length(select.cells), length(fg)), 
                          bg, 
                          rep(length(not.select.cells), length(bg))
  )
  
  g <- names(fg.freq)[fg.freq > de.param$q1.th & tau > de.param$q.diff.th]
  g <- g[order(tau[g] + ratio[g] / 4 + fg.freq[g] / 5, decreasing = T)]
  select.g <- c(g[tau[g] > 0.95], head(g, n.markers))
  g <- g[g %in% select.g]
  
  if(length(g) > 0){
    df <- data.frame(g = g,
                     specificity = round(tau[g], digits = 2), 
                     fg.freq = round(fg.freq[g], digits = 2), 
                     bg.freq = round(bg.freq[g], digits = 2), 
                     fg.counts = fg[g],
                     bg.counts = bg[g],
                     pval = stats[g, "pval"])
    return(df)
  } else {
    return(NULL)
  }
}


#' Title
#'
#' @param cl.g 
#' @param norm.dat 
#' @param cl 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
within_group_specific_markers <- function(cl.g, norm.dat, cl, ...)
  {
    cl = as.factor(cl)
    cl = droplevels(cl[cl %in% cl.g])
    df = do.call("rbind", sapply(cl.g, function(x){
      df=group_specific_markers(x, norm.dat, cl,...)
      if(!is.null(df)){
        df$cl = x
        df
      }
      else{
        NULL
      }
    },simplify=F))
    return(df)
  }


###Beta score from Trygve
#' Title
#'
#' @param propExpr 
#' @param spec.exp 
#' @param mcores 
#'
#' @return
#' @export
#'
#' @examples
get_beta_score <- function(propExpr, spec.exp = 2, mcores=1){
  library(doMC)
  if(mcores ==1){
    registerDoSEQ()
  }
  else{
    registerDoMC(cores=mc.cores)
    on.exit(parallel::stopCluster(), add = TRUE)
  }
  calc_beta <- function(y, spec.exp = 2, eps1 = 1e-10) {
    d1 <- as.matrix(dist(y))
    score1 <- sum(d1^spec.exp) / (sum(d1) + eps1)
    return(score1)
  }
  res <- foreach(i = 1:nrow(propExpr), .combine="c") %dopar% {
    print(i)
    calc_beta(propExpr[i,])
  }
  names(res) = row.names(propExpr)
  return(res)
}



###select_markers_groups 
#' Title
#'
#' @param de.genes 
#' @param cl.group Assignment of clusters to groups cluster as names, and group id as values.
#'
#' @return
#' @export
#'
#' @examples


select_markers_groups <- function(de.genes, cl.group, cl.means, n.markers=1,...)
  {
    pairs = create_pairs(as.character(names(cl.group)))
    pairs = pairs[cl.group[pairs[,1]]!= cl.group[pairs[,2]], ]
    markers=select_N_markers(de.genes, cl.means=cl.means, pairs=row.names(pairs), ...)
    return(markers)
  }

  
