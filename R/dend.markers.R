selectMarkersPair <- function(norm.dat, add.genes,de.genes=NULL, gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    pairs =do.call("rbind",strsplit(gsub("cl","",names(add.genes)), "_"))
    row.names(pairs)= names(add.genes)
    de.genes.list=list()
    if(is.null(gene.score)){
      de.genes= de.genes[names(add.genes)]
      tmp=getGeneScore(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score      
      gene.score=pmin(up.gene.score, down.gene.score)
      row.names(gene.score)=all.genes
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

selectMarkersPairDirection <- function(add.up,add.down,de.genes=NULL, up.gene.score=NULL,down.gene.score=NULL,rm.genes=NULL,top.n=50,max.num=2000)
  {
    up.genes = down.genes=list()
    final.genes=c()
    if(is.null(up.gene.score)){
      pairs.n=union(names(add.up), names(add.down))
      pairs =do.call("rbind",strsplit(gsub("cl","",pairs.n), "_"))
      row.names(pairs)= pairs.n
      de.genes= de.genes[pairs.n]
      tmp=getGeneScore(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
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


selectNMarkers <- function(de.genes, up.gene.score=NULL, down.gene.score=NULL, default.markers=NULL, pair.num = 1,rm.genes=NULL)
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
   result = selectMarkersPairDirection(add.up,add.down,de.genes=de.genes, up.gene.score=up.gene.score,down.gene.score=down.gene.score,rm.genes=c(rm.genes,default.markers),top.n=50,max.num=2000)
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


selectMarkersPairGroup <- function(norm.dat,cl, g1,g2,de.genes,top.n=50,max.num=1000,n.markers=20,up.gene.score=NULL, down.gene.score=NULL)
  {
    pairs = do.call("rbind",strsplit(names(de.genes), "_"))
    pairs = gsub("cl", "",pairs)
    row.names(pairs)= names(de.genes)
    up.pairs = row.names(pairs)[pairs[,1] %in% g1 & pairs[,2] %in% g2]
    down.pairs = row.names(pairs)[pairs[,1] %in% g2 & pairs[,2] %in% g1]
    select.pairs = c(up.pairs, down.pairs)
    if(is.null(up.gene.score)){
      tmp=get_gene_score(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score
      row.names(down.gene.score)= all.genes
    }
    all.genes = row.names(up.gene.score)
    tmp.up.gene.score = cbind(up.gene.score[,up.pairs,drop=F], down.gene.score[,down.pairs,drop=F])
    tmp.down.gene.score = cbind(down.gene.score[,up.pairs,drop=F], up.gene.score[,down.pairs,drop=F])
    
    up.genes = row.names(tmp.up.gene.score)[head(order(rowSums(tmp.up.gene.score)), n.markers)]
    down.genes = row.names(tmp.down.gene.score)[head(order(rowSums(tmp.down.gene.score)), n.markers)]
        
    up.num = colSums(tmp.up.gene.score[up.genes,,drop=F] < max.num)
    down.num = colSums(tmp.down.gene.score[down.genes,,drop=F] < max.num)
    total.num = up.num + down.num
    add.genes = setNames(rep(n.markers, ncol(tmp.up.gene.score)), colnames(tmp.up.gene.score)) - total.num
    add.genes = add.genes[add.genes > 0]
    
    up.genes = up.genes[rowMins(tmp.up.gene.score[up.genes,,drop=F]) < max.num]
    down.genes = down.genes[rowMins(tmp.down.gene.score[down.genes,,drop=F]) < max.num]
    genes = union(up.genes, down.genes)
    if(length(add.genes)>0){
      tmp=selectMarkersPair(norm.dat, add.genes= add.genes,de.genes= de.genes, gene.score=pmin(tmp.up.gene.score, tmp.down.gene.score), rm.genes=c(up.genes, down.genes),top.n=top.n)
      genes=union(genes, unlist(tmp))
    }
    return(genes)
  }




findNodeSpecificMarkers <- function(dend,norm.dat, cl, cl.df, n.markers=10,de.genes=NULL,up.gene.score=NULL, down.gene.score=NULL,top.n=50,max.num=2000)
{
  m=list()
  g1= labels(dend)
  g2=  setdiff(levels(cl), g1)
  all.genes = row.names(up.gene.score)
  pairs = do.call("rbind",strsplit(colnames(up.gene.score), "_"))
  row.names(pairs)= colnames(up.gene.score)
  up.pairs = row.names(pairs)[pairs[,1] %in% g1 & pairs[,2] %in% g2]
  down.pairs = row.names(pairs)[pairs[,1] %in% g2 & pairs[,2] %in% g1]
  if(length(up.pairs)>0 & length(down.pairs)>0){
    tmp.up.gene.score = cbind(up.gene.score[,up.pairs,drop=F], down.gene.score[,down.pairs,drop=F])
    tmp.down.gene.score = cbind(down.gene.score[,up.pairs,drop=F], up.gene.score[,down.pairs,drop=F])
    up.genes = row.names(tmp.up.gene.score)[head(order(rowSums(tmp.up.gene.score)), n.markers)]
    up.genes = up.genes[rowSums(tmp.up.gene.score[up.genes,,drop=F] < max.num) > 0]
    if(length(up.genes)>0){
      up.num = colSums(tmp.up.gene.score[up.genes,,drop=F] < max.num)
      markers=up.genes
      m[[attr(dend,"label")]]=markers
    }
  }
  if(length(dend)>1){
    for(i in 1:length(dend)){
      m=c(m,findNodeSpecificMarkers(dend[[i]],norm.dat, cl, cl.df, de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,top.n=top.n, max.num=max.num, n.markers=n.markers))
    }
  }
  return(m)
}

findDendMarkers <- function(dend,norm.dat, cl, cl.df,de.genes,up.gene.score=NULL, down.gene.score=NULL,...)
  {
    require(dendextend)
    require(randomForest)
    print(dend)
    if(length(dend)>1){
      cl.g = sapply(1:length(dend),function(i){
        g = row.names(cl.df)[cl.df$cluster_label %in% (dend[[i]] %>% labels)]
      },simplify=F)
      
      markers = c()
      for(i in 1:(length(cl.g)-1)){
        for(j in (i+1):length(cl.g)){
          g = selectMarkersPairGroup(norm.dat,cl, cl.g[[i]],cl.g[[j]],de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)
          markers=union(markers, g)          
        }
      }
      
      tmp.cl = setNames(rep(1:length(cl.g),sapply(cl.g,length)),unlist(cl.g))
      select.cells = names(cl)[cl %in% names(tmp.cl)]
      select.cells = unlist(tapply(select.cells, cl[select.cells],function(x)sample(x,min(length(x), 50))))
      tmp.cl = setNames(tmp.cl[as.character(cl[select.cells])],select.cells)
      ###if too few genes selected, add more genes
      if(length(markers)< 4){
        de.df= DE.genes.pw(norm.dat[,select.cells], paste0("cl",tmp.cl))
        markers = selectMarkers(norm.dat, tmp.cl, de.df, n.markers=10, q=0.4, q.b=0.7)$markers
      }
      
      rf = randomForest(t(norm.dat[markers,select.cells]),factor(tmp.cl))
      w = importance(rf)
      attr(dend, "markers")=w[,1]
      for(i in 1:length(dend)){
        dend[[i]] = findDendMarkers(dend[[i]], norm.dat, cl, cl.df,de.genes=de.genes,up.gene.score=up.gene.score, down.gene.score=down.gene.score,...)
      }
    }
    return(dend)
  }

findDendMarkersSpecficity <- function(dend,norm.dat, cl, cl.label)
{
  if(length(dend)>1){
    cl.g = sapply(1:length(dend),function(i){
      g = names(cl.label)[cl.label %in% (dend[[i]] %>% labels)]
    },simplify=F)
    markers=names(sort(attr(dend, "markers"),decreasing=T))
    cl.g.mean = sapply(cl.g, function(x){
      tmp.cells = names(cl)[cl %in% x]
      rowMeans(norm.dat[markers,tmp.cells])
    })
    colnames(cl.g.mean)=1:ncol(cl.g.mean)
    cl.g.mean.diff1 = cl.g.mean - rowMaxs(cl.g.mean)
    cl.g.mean.diff2 = cl.g.mean - rowMins(cl.g.mean)
    cl.g.list = lapply(colnames(cl.g.mean), function(x){row.names(cl.g.mean)[cl.g.mean.diff1[,x] > -1 & cl.g.mean.diff2[,x] > 1]})
    names(cl.g.list) = sapply(1:length(dend),function(i){attr(dend[[i]],"label")})
    attr(dend, "markers.byCl")= cl.g.list
    for(i in 1:length(dend)){
      dend[[i]] = findDendMarkersSpecficity(dend[[i]],norm.dat, cl, cl.label)
    }
  }
  return(dend)
}
  


mapDendMarkers <- function(dend.list, map.dat,select.cells,th=0.5)
  {
    map.gene.num <- matrix(0, nrow=ncol(map.dat), ncol= length(dend.list))
    row.names(map.gene.num)= colnames(map.dat)
    colnames(map.gene.num)=names(dend.list)
    for(x in names(dend.list)){
      node = dend.list[[x]]
      if(length(node)>1){
        for(i in 1:length(node)){
          l = attr(node[[i]], "label")
          print(l)
          tmp.genes= intersect(attr(node, "markers.byCl")[[i]],row.names(map.dat))
          map.gene.num[select.cells, l] = colSums(map.dat[tmp.genes,select.cells]>th)
        }
      }
    }
    return(map.gene.num)
  }

getNodeSpecificMarkers <- function(dend.list, norm.dat, cl,cl.df,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      print(x)
      cl.g = row.names(cl.df)[cl.df$cluster_label %in% labels(dend.list[[x]])]      
      df=getGroupSpecificMarkers(cl.g, norm.dat, cl,...)
      if(!is.null(df)){
        df$cl = x
      }
      df
    },simplify=F))
  }


getNodeVsSiblingMarkers <- function(dend.list, norm.dat, cl,cl.df,...)
  {    
    do.call("rbind",sapply(names(dend.list), function(x){
      dend = dend.list[[x]]
      all.cl =  droplevels(cl[cl %in% row.names(cl.df)[cl.df$cluster_label %in% labels(dend)]])
      if(length(dend)>1){
        do.call("rbind",sapply(1:length(dend),function(i){          
          cl.g = row.names(cl.df)[cl.df$cluster_label %in% labels(dend[[i]])]
          df=getGroupSpecificMarkers(cl.g, norm.dat, all.cl,...)
          if(!is.null(df)){
            df$cl = attr(dend[[i]],"label")
          }
          df
        },simplify=F))
      }
      else{
        NULL
      }
    },simplify=F))
  }

getGroupSpecificMarkers <- function(cl.g, norm.dat, cl,cl.present.counts=NULL,low.th=1,q1.th=0.5, q.diff.th=0.7,n.markers=5)
  {
    cl= droplevels(cl)
    select.cells = names(cl)[cl %in% cl.g]
    not.select.cells = setdiff(names(cl), select.cells)
    if(is.null(cl.present.counts)){
      fg = rowSums(norm.dat[,select.cells]> low.th)/length(select.cells)
      bg = rowSums(norm.dat[,not.select.cells]> low.th)/length(not.select.cells)
    }
    else{
      fg = rowSums(cl.present.counts[,cl.g,drop=F])
      bg = rowSums(cl.present.counts[,levels(cl),drop=F]) - fg
      bg.freq= bg/length(not.select.cells)
      fg.freq = fg/length(select.cells)
    }
    tau = (fg.freq - bg.freq)/pmax(bg.freq,fg.freq)
    g = names(fg.freq)[fg.freq > q1.th & tau > q.diff.th]
    g = g[order(tau[g]+fg.freq[g],decreasing=T)]
    g = head(g, n.markers)
    if(length(g > 0)){
      df=data.frame(g=g,specifity=tau[g], fg.freq=fg.freq[g],bg.freq = bg.freq[g], fg.counts=fg[g],bg.counts=bg[g])
      return(df)
    }
    return(NULL)
  }


calc_tau <- function(m, byRow=TRUE)
{
  if(!byRow){
    m = t(m)
  }
  m = m/rowMaxs(m)
  tau = rowSums(1 - m)/(ncol(m) - 1)
  tau[is.na(tau)]=0
  return(tau)
}

cl_max_tau <- function(cl.dat, th=0.5, tau.th=0.7,n.markers=3)
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

tau_one_vs_other_genes <- function(cl, cl.present.counts, present.th=0.4, tau.th=0.8,top.n=10)
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

