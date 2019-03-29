select_markers <- function(norm.dat, cl, n.markers=20,de.genes=NULL, ...)                           
  {
    if(is.null(de.genes)){
      de.genes=de_score(norm.dat, cl, ...)
    }
    pairs = names(de.genes)
    pairs.df = gsub("cl","", do.call("rbind",strsplit(pairs, "_")))
    row.names(pairs.df)=pairs
    select.pairs = pairs[pairs.df[,1] %in% cl & pairs.df[,2]%in% cl]
    de.markers = sapply(select.pairs, function(s){
      tmp = de.genes[[s]]
      c(head(tmp$up.genes,n.markers), head(tmp$down.genes,n.markers))
    },simplify=F)
    markers = intersect(unlist(de.markers),row.names(norm.dat))
    return(list(markers=markers, de.genes=de.genes[select.pairs]))
  }

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




get_gene_score <- function(de.genes,top.n=50, max.num=1000,bin.th=4)
  {
    select.genes <- sapply(de.genes, function(x){
      up.genes = x$up.genes
      #Deal with pairs with too many DEX genes, include all binary markers
      up.binary.genes = up.genes[x$q.stats[up.genes, 3] < 1 & x$q.stats[up.genes,2]>bin.th]
      up.genes=up.genes[up.genes %in% c(head(up.genes, top.n), up.binary.genes)]
      down.genes = x$down.genes
      down.binary.genes = down.genes[x$q.stats[down.genes, 1] < 1 & x$q.stats[down.genes,4]>bin.th]
      down.genes=down.genes[down.genes %in% c(head(down.genes, top.n), down.binary.genes)]
      list(up=up.genes, down=down.genes)
    },simplify=F)
    all.genes=unique(unlist(select.genes))
    up.gene.score = sapply(select.genes, function(x){   
      tmp=match(all.genes, x$up)
      tmp[is.na(tmp)]=max.num
      tmp
    })
    down.gene.score = sapply(select.genes, function(x){   
      tmp=match(all.genes, x$down)
      tmp[is.na(tmp)]=max.num
      tmp
    })
    row.names(up.gene.score)=row.names(down.gene.score)= all.genes
    return(list(up.gene.score=up.gene.score, down.gene.score=down.gene.score))
  }



select_markers_pair <- function(norm.dat, de.genes, add.genes, gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    pairs =do.call("rbind",strsplit(gsub("cl","",names(add.genes)), "_"))
    row.names(pairs)= names(add.genes)
    de.genes.list=list()
    if(is.null(gene.score)){
      de.genes= de.genes[names(add.genes)]
      tmp=get_gene_score(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score      
      gene.score=pmin(up.gene.score, down.gene.score)
      row.names(gene.score)=row.names(up.gene.score)
    }
    select.genes = setdiff(row.names(gene.score),rm.genes)
    gene.score= gene.score[select.genes,names(add.genes),drop=F]
    final.genes=list()
    while(sum(add.genes)>0 & nrow(gene.score)>0){
      g = order(rowSums(gene.score))[1]
      g.n = row.names(gene.score)[g]
      select.pair = colnames(gene.score)[gene.score[g, ]< top.n]
      add.genes[select.pair]= add.genes[select.pair]-1
      for(p in select.pair){
        de.genes.list[[p]]=union(de.genes.list[[p]],g.n)
      }
      add.genes = add.genes[add.genes > 0]
      gene.score=gene.score[-g,names(add.genes),drop=F]
    }
    return(de.genes.list)
  }


select_markers_pair_direction <- function(de.genes, add.up,add.down,up.gene.score=NULL,down.gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    up.genes = down.genes=list()
    final.genes=c()
    if(is.null(up.gene.score)){
      pairs.n=union(names(add.up), names(add.down))
      pairs =do.call("rbind",strsplit(gsub("cl","",pairs.n), "_"))
      row.names(pairs)= pairs.n
      de.genes= de.genes[pairs.n]
      tmp=get_gene_score(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score      
    }
    select.genes = setdiff(row.names(up.gene.score),rm.genes)
    while((sum(add.up)+sum(add.down)>0) & length(select.genes)>0){
      sc =0
      if(length(add.up)>0){
        sc = rowSums(up.gene.score[select.genes, names(add.up),drop=F])
      }
      if(length(add.down)>0){
        sc = sc+rowSums(down.gene.score[select.genes, names(add.down),drop=F])
      }
      g = select.genes[ which.min(sc)]
      select.up.pair = names(add.up)[up.gene.score[g,names(add.up)]< top.n]
      select.down.pair = names(add.down)[down.gene.score[g,names(add.down)]< top.n]
      if(length(select.up.pair) + length(select.down.pair)==0){
        break
      }
      add.up[select.up.pair]= add.up[select.up.pair]-1
      add.down[select.down.pair]= add.down[select.down.pair]-1
      add.up = add.up[add.up >0, drop=F]
      add.down = add.down[add.down >0, drop=F]
      for(p in select.up.pair){
        up.genes[[p]]=c(up.genes[[p]],g)
      }
      for(p in select.down.pair){
        down.genes[[p]]=c(down.genes[[p]],g)
      }
      final.genes=c(final.genes, g)
      select.genes = setdiff(select.genes, final.genes)
      if(length(select.genes)>0){
        select = rep(TRUE,length(select.genes))
        if(length(add.up)>0){
          select = rowSums(up.gene.score[select.genes,names(add.up),drop=F] < top.n)>0
        }
        if(length(add.down)>0){
          select = select | rowSums(down.gene.score[select.genes,names(add.down),drop=F] < top.n)>0
        }
        select.genes= select.genes[select]
      }
      cat(g, length(add.up), length(add.down),"\n")
    }
    markers= final.genes
    return(list(markers=markers, up.genes=up.genes, down.genes=down.genes))
  }


select_N_markers <- function(de.genes, up.gene.score=NULL, down.gene.score=NULL, default.markers=NULL, pair.num = 1,rm.genes=NULL)
  {
   add.up = add.down=setNames(rep(pair.num, length(de.genes)), names(de.genes))
   if(!is.null(default.markers)){
     up.default = sapply(de.genes, function(x){intersect(x$up.genes, default.markers)},simplify=F)
     down.default = sapply(de.genes, function(x){intersect(x$down.genes, default.markers)},simplify=F)
     add.up = pmax(add.up -  sapply(up.default, length),0)
     add.down = pmax(add.down -  sapply(down.default, length),0)
   }
   add.up = add.up[add.up>0, drop=F]
   add.down = add.down[add.down>0, drop=F]
   result = select_markers_pair_direction(add.up,add.down,de.genes=de.genes, up.gene.score=up.gene.score,down.gene.score=down.gene.score,rm.genes=c(rm.genes,default.markers),top.n=50,max.num=2000)
   up.genes = up.default
   down.genes=down.default
   for(x in names(result$up.genes)){
     up.genes[[x]]=c(up.genes[[x]], result$up.genes[[x]])
   }
   for(x in names(result$down.genes)){
     down.genes[[x]]=c(down.genes[[x]], result$down.genes[[x]])
   }
   return(list(up.genes=up.genes, down.genes=down.genes, markers=c(default.markers, result$markers)))
 }


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



node_vs_sibling_markers <- function(dend.list, norm.dat, cl,cl.df,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      dend = dend.list[[x]]
      print(labels(dend)) 
      all.cl =  droplevels(cl[cl %in% labels(dend)])
      if(length(dend)>1){
        do.call("rbind",sapply(1:length(dend),function(i){          
          cl.g =  labels(dend[[i]])
          df=group_specific_markers(cl.g, norm.dat, all.cl,...)
          if(!is.null(df)){
            df$node = attr(dend[[i]],"label")
            df$parent = attr(dend,"label")
          }
          df
        },simplify=F))
      }
      else{
        NULL
      }
    },simplify=F))
  }


node_specific_markers <- function(dend.list, norm.dat, cl,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      print(x)
      df=group_specific_markers(labels(dend.list[[x]]), norm.dat, cl,...)
      if(!is.null(df)){
        df$cl = x
      }
      df
    },simplify=F))
  }


select_pos_markers <- function(de.genes, cl, n.markers=3, up.gene.score=NULL, down.gene.score=NULL)
  {
    pairs = names(de.genes)
    pairs.df = as.data.frame(do.call("rbind", strsplit(pairs, "_")))
    row.names(pairs.df) = pairs
    pairs.df = pairs.df[pairs.df[,1]%in% levels(cl) & pairs.df[,2]%in% levels(cl),]
    if(is.null(up.gene.score) | is.null(down.gene.score)){
      tmp=get_gene_score(de.genes)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score
    }
    cl.df$markers=""
###for each cluster, find markers that discriminate it from other types
    cl.markers <- sapply(levels(cl), function(tmp.cl){
      up.pairs = row.names(pairs.df)[pairs.df[,1] == tmp.cl]
      down.pairs = row.names(pairs.df)[pairs.df[,2] == tmp.cl]
      add.up = setNames(rep(n.markers, length(up.pairs)), up.pairs)
      add.down = setNames(rep(n.markers, length(down.pairs)), down.pairs)
      tmp.result = select_markers_pair_direction(de.genes, add.up=add.up, add.down = add.down,up.gene.score=up.gene.score, down.gene.score=down.gene.score )
      tmp.result$markers
    },simplify=F)  
  }



###Beta score from Trygve

get_beta_score <- function(propExpr, spec.exp = 2, mcores=1){
  if(mcores ==1){
    registerDoSEQ()
  }
  else{
    cl <- parallel::makeForkCluster(mcores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
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
