#' Title
#'
#' @param all.results 
#'
#' @return
#' @export
#'
#' @examples
combine_cl <- function(all.results)
  {
    cl = all.results[[1]]$cl
    cl = setNames(as.integer(cl),names(cl))
    markers=all.results[[1]]$markers
    n.cl = max(cl)
    for(i in 2:length(all.results)){
      if(is.null(all.results[[i]]$cl) | length(unique(all.results[[i]]$cl)) < 2) next
      new.cl = all.results[[i]]$cl
      new.cl = setNames(as.integer(new.cl)+ n.cl,names(new.cl))
      cl[names(new.cl)] = new.cl
      n.cl = max(cl)
      cat(names(all.results)[i], n.cl, "\n")
    }
    return(cl)
  }

#' Title
#'
#' @param cl 
#' @param cl.df 
#' @param comb.dat 
#' @param prefix 
#' @param tsne.df 
#' @param cex 
#' @param fn.size 
#' @param height 
#' @param width 
#'
#' @return
#' @export
#'
#' @examples
plot_tsne <- function(cl, cl.df, comb.dat, prefix, tsne.df, cex=0.3, fn.size=2, height=8, width=10)
  {
    library(ggplot2)
    library(gridExtra)
    with(comb.dat,{
    #cl.df$cluster_color = adjust_color(cl.df$cluster_color, adj.col)
    tmp = plot_tsne_cl(cl=cl, cl.df=cl.df,  tsne.df = tsne.df, cex=cex, fn.size = fn.size)
    tsne.df = tmp$tsne.df
    ggsave(paste("tsne.cl",prefix,"pdf", sep="."), tmp$g, height=height,width=width)
    
    tmp.df = meta.df[row.names(tsne.df),]
    tmp= setNames(as.factor(tmp.df$platform), row.names(tsne.df))
    meta.col=NULL
    if(length(levels(tmp))==2){
      meta.col = setNames(c("blue", "orange"), levels(tmp))
    }
    g= plot_tsne_meta(tsne.df, tmp, meta.col=meta.col, cex=cex)
    ggsave(paste("tsne",prefix, "platform.pdf", sep="."), g, height=height, width=width)

    plots = lapply(names(cl.list),function(x){
      tmp.cl = cl.list[[x]]
      tmp.cl = droplevels(tmp.cl[names(tmp.cl) %in% names(cl)])
      if(length(tmp.cl)==0){
        return(NULL)
      }
      g = plot_tsne_cl(cl=cl.list[[x]], cl.df=cl.df.list[[x]],  tsne.df = tsne.df[names(cl.list[[x]]),], cex=cex, fn.size = fn.size)$g
      g = g + ggtitle(x)
    })
    plots = plots[!sapply(plots, is.null)]
    ggsave(paste("tsne",prefix, "by.platform.cl.pdf", sep="."), marrangeGrob(grobs=plots, nrow=1, ncol=1), height=height, width=width)

    plots = lapply(names(cl.list),function(x){
      tmp.cells = names(cl)[meta.df[names(cl),"platform"]==x]
      if(length(tmp.cells)==0){
        return(NULL)
      }
      g = plot_tsne_cl(cl=droplevels(cl[tmp.cells]), cl.df=cl.df,  tsne.df = tsne.df[names(cl.list[[x]]),], cex=cex, fn.size = fn.size)$g
      g = g + ggtitle(x)
    })
    plots = plots[!sapply(plots, is.null)]
    ggsave(paste("tsne",prefix, "by.platform.pdf", sep="."), marrangeGrob(grobs=plots, nrow=1, ncol=1), height=height, width=width)
  })
  }



#' Title
#'
#' @param consensus.cl 
#' @param prefix 
#' @param comb.dat 
#' @param consensus.cl.df 
#' @param do.droplevels 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_confusion <- function(consensus.cl, prefix, comb.dat,consensus.cl.df = NULL, do.droplevels = FALSE,cex.x=6, cex.y=6,...)
{
  g.list=list()
  for(x in names(comb.dat$cl.list)){
    if(sum(names(comb.dat$cl.list[[x]]) %in% names(consensus.cl)) > 0){
      if(is.null(consensus.cl.df)){
        g = compare_annotate(consensus.cl, comb.dat$cl.list[[x]], comb.dat$cl.df.list[[x]], rename = FALSE, do.droplevels=do.droplevels,cex.x=cex.x, cex.y=cex.y)$g
        g = g + xlab("consensus cluster") + ylab(x)
      }
      else{
        g = compare_annotate(comb.dat$cl.list[[x]], consensus.cl, consensus.cl.df, rename = FALSE, do.droplevels=do.droplevels)$g
        g = g + ylab("consensus cluster") + xlab(x)
      }
      g.list[[x]] <- g
      ggsave(paste(prefix, x, "pdf", sep="."), g,...)
    }
  }
  return(g.list)
}



#' Title
#'
#' @param dat.list 
#' @param cl 
#' @param de.param.list
#' @param min.cells
#' @param select.genes 
#' @param sets 
#'
#' @return
#' @export
#'
#' @examples

get_cl_means_list <- function(dat.list, cl, de.param.list=NULL, min.cells=NULL, select.genes=NULL, sets=names(dat.list))
  {
    if(is.null(min.cells)){
      if(!is.null(de.param.list)){
        min.cells = lapply(de.param.list,"[[", "min.cells")
      }
      else{
        min.cells = setNames(rep(1, length(dat.list)), names(dat.list))
      }
    }
    else{
      if(length(min.cells) ==1){
        min.cells = setNames(rep(min.cells, length(dat.list)), names(dat.list))
      }
    }
    if(is.null(min.cells)){
      min.cells = setNames(rep(1, length(dat.list)), names(dat.list))
    }
    cl.means.list = list()
    for(x in sets){
      tmp.cells = intersect(names(cl), colnames(dat.list[[x]]))
      tmp.cl = cl[tmp.cells]
      cl.size = table(tmp.cl)
      select.cl = names(cl.size)[cl.size >= min.cells[[x]]]
      if(length(select.cl)==0){
        next
      }
      tmp.cl = tmp.cl[tmp.cl %in% select.cl]
      if(is.factor(tmp.cl)){
        tmp.cl=droplevels(tmp.cl)
      }
      if(is.null(select.genes)){
        tmp=get_cl_means(dat.list[[x]], tmp.cl)
      }
      else{
        tmp=get_cl_means(dat.list[[x]], tmp.cl)[select.genes,,drop=F]
      }
      cl.means.list[[x]]= tmp
    }
    return(cl.means.list)
  }


get_cl_sqr_means_list <- function(dat.list, cl, de.param.list=NULL, min.cells=NULL, select.genes=NULL, sets=names(dat.list))
  {
    if(is.null(min.cells)){
      if(!is.null(de.param.list)){
        min.cells = lapply(de.param.list,"[[", "min.cells")
      }
      else{
        min.cells = setNames(rep(1, length(dat.list)), names(dat.list))
      }
    }
    else{
      if(length(min.cells) ==1){
        min.cells = setNames(rep(min.cells, length(dat.list)), names(dat.list))
      }
    }
    cl.sqr.means.list = list()
    for(x in sets){
      tmp.cells = intersect(names(cl), colnames(dat.list[[x]]))
      tmp.cl = cl[tmp.cells]
      cl.size = table(tmp.cl)
      select.cl = names(cl.size)[cl.size >= min.cells[[x]]]
      if(length(select.cl)==0){
        return(NULL)
      }
      tmp.cl = tmp.cl[tmp.cl %in% select.cl]
      if(is.factor(tmp.cl)){
        tmp.cl=droplevels(tmp.cl)
      }
      tmp=get_cl_sqr_means(dat.list[[x]], tmp.cl)
      if(!is.null(select.genes)){
        tmp=tmp[select.genes,,drop=F]
      }
      cl.sqr.means.list[[x]]= tmp
    }
    return(cl.sqr.means.list)
  }



#' Title
#'
#' @param dat.list 
#' @param de.param.list 
#' @param select.genes 
#' @param cl 
#' @param sets 
#'
#' @return
#' @export
#'
#' @examples
get_cl_present_list <- function(dat.list, de.param.list=NULL, min.cells=NULL, select.genes=NULL, cl, sets=names(dat.list))
  {
    if(is.null(min.cells)){
      if(!is.null(de.param.list)){
        min.cells = lapply(de.param.list,"[[", "min.cells")
      }
      else{
        min.cells = setNames(rep(1, length(dat.list)), names(dat.list))
      }
    }
    else{
      if(length(min.cells) ==1){
        min.cells = setNames(rep(min.cells, length(dat.list)), names(dat.list))
      }
    }
    cl.present.list = list()
    for(x in sets){
      tmp.cells = intersect(names(cl), colnames(dat.list[[x]]))
      tmp.cl = cl[tmp.cells]
      cl.size = table(tmp.cl)
      select.cl = names(cl.size)[cl.size >= min.cells[[x]]]
      if(length(select.cl)==0){
        return(NULL)
      }
      tmp.cl = tmp.cl[tmp.cl %in% select.cl]
      if(is.factor(tmp.cl)){
        tmp.cl=droplevels(tmp.cl)
      }
      tmp=get_cl_present(dat.list[[x]], tmp.cl, de.param.list[[x]]$low.th)
      if(!is.null(select.genes)){
        tmp=tmp[select.genes,,drop=F]
      }
      cl.present.list[[x]]= tmp
    }
    return(cl.present.list)
  }



#' Title
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
get_gene_cl_correlation <- function(cl.means.list)
  {
    sets=names(cl.means.list)
    gene.cl.cor = list()
    for(i in 1:(length(cl.means.list)-1)){
      for(j in (i+1):length(cl.means.list)){
        pair= paste(sets[i], sets[j], sep=":")
        common.cl = intersect(colnames(cl.means.list[[i]]), colnames(cl.means.list[[j]]))
        common.genes = intersect(row.names(cl.means.list[[i]]), row.names(cl.means.list[[j]]))
        gene.cor =  pair_cor(cl.means.list[[i]][common.genes,common.cl],cl.means.list[[j]][common.genes,common.cl])
        gene.cl.cor[[pair]] = gene.cor
      }
    }
    return(gene.cl.cor)
  }

#' Title
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
get_de_lfc_list <- function(cl.means.list)
  {
    sets=names(cl.means.list)
    de.gene.sign = NULL
    de.lfc.list = sapply(sets, function(set){
      cl.means = cl.means.list[[set]]
      cn = colnames(cl.means)
      cl.n = length(cn)
      pairs = cbind(rep(cn, rep(cl.n, cl.n)), rep(cn, cl.n))
      pairs = pairs[pairs[, 1] < pairs[, 2], , drop = F]
      row.names(pairs)= paste0(pairs[,1],"_", pairs[,2])
      lfc = cl.means[,pairs[,1]] - cl.means[,pairs[,2]]
      colnames(lfc) = row.names(pairs)
      lfc
    },simplify=F)
    return(de.lfc.list)
  }



#' Title
#'
#' @param dat.list 
#' @param de.param.list 
#' @param cl 
#' @param select.sets 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
get_de_result  <- function(dat.list, de.param.list,  cl, select.sets=names(de.param.list), max.cl.size=200,...)
  {
    de.genes.list <- sapply(select.sets, function(x){
      tmp.cl =cl[names(cl) %in% colnames(dat.list[[x]])]
      if(is.factor(tmp.cl)){
        tmp.cl = droplevels(tmp.cl)
      }
      tmp.cells = sample_cells(tmp.cl, max.cl.size)
      tmp.cl = tmp.cl[tmp.cells]

      if(length(unique(tmp.cl)) > 1){
        de.result = de_all_pairs(norm.dat= dat.list[[x]], cl=tmp.cl, de.param= de.param.list[[x]],...)
      }
      else{
        return(NULL)
      }
    },simplify=F)
  }

#' Title
#'
#' @param de.genes.list 
#' @param cl.means.list 
#' @param common.genes 
#' @param max.num 
#' @param pairs 
#' @param frac.th 
#'
#' @return
#' @export
#'
#' @examples
comb_de_result <- function(de.genes.list, cl.means.list, common.genes=NULL, max.num=10000, pairs=NULL, frac.th=0.7)
{
  sets = names(cl.means.list)
  if(is.null(common.genes)){
    common.genes = row.names(cl.means.list[[1]])
    for(x in 2:length(cl.means.list)){
      common.genes= intersect(common.genes, row.names(cl.means.list[[x]]))
    }
  }
  if(is.null(pairs)){
    pairs = unique(unlist(lapply(de.genes.list, names)))
  }
  de.genes = list()
  for(p in pairs){
    de.counts = table(unlist(lapply(names(de.genes.list), function(set){
      de = de.genes.list[[set]][[p]]
      c(names(de$up.genes),names(de$down.genes))})))      
    g = intersect(names(de.counts), common.genes)
    pair = unlist(strsplit(p,"_"))
    lfc = lapply(sets, function(set){
      cl.means = cl.means.list[[set]]
      if(all(pair %in% colnames(cl.means))){
        lfc = cl.means[g,pair[1]] - cl.means[g,pair[2]]
      }
      else{
        NULL
      }
    })
    lfc = do.call("cbind",lfc)
    if(is.null(lfc)){
      next
    }
    lfc = as.matrix(lfc)

    
    rank =  do.call("cbind",lapply(names(de.genes.list), function(set){
      de = de.genes.list[[set]][[p]]
      
      tmp1=match(g, names(de$up.genes))
      tmp2=match(g, names(de$down.genes))
      tmp1[is.na(tmp1)] = max.num
      tmp2[is.na(tmp2)] = max.num
      tmp = pmin(tmp1, tmp2)
    })) 
    rank[rank > max.num] = max.num
    row.names(rank)=g
    rank.mean = rowMeans(rank)
    sign = rowSums(lfc > 0)
    sign1 = rowSums(lfc > 1)
    sign2 = rowSums(lfc < -1)
    frac = pmax(sign1, sign2)/ncol(lfc)
    lfc = rowMeans(lfc)
    select.g = names(frac)[frac > frac.th]
    df = data.frame(lfc =lfc, frac=frac, counts = as.vector(de.counts[g]), rank.mean=rank.mean)
    
    df = df[select.g,,drop=FALSE]
    ord = order(with(df, rank.mean))
    df = df[ord,]
    
    up = row.names(df)[which(df$lfc > 0)]
    down = row.names(df)[which(df$lfc < 0)]
    select = c(up,down)
    de.genes[[p]] = list(up.genes =setNames(df[up,"counts"], up), down.genes=setNames(df[down,"counts"], down), up.num = length(up),down.num=length(down),num=length(select),genes=select, de.df = df)    
  }
  return(de.genes)
}


#' Title
#'
#' @param comb.dat 
#' @param cl 
#' @param de.genes.list 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
select_joint_markers <- function(comb.dat, cl, de.genes.list,...)
  {
    tmp <- sapply(names(de.genes.list), function(x){
      tmp.cl =cl[names(cl) %in% colnames(comb.dat$dat.list[[x]])]
      if(is.factor(tmp.cl)){
        tmp.cl = droplevels(tmp.cl)
      }
      if(length(levels(tmp.cl)) > 1){
        de.result = display_cl(tmp.cl, comb.dat$dat.list[[x]], max.cl.size = 200, de.genes=de.genes.list[[x]],...)$markers
      }
    },simplify=F)
    marker.counts <- table(unlist(tmp))
    return(marker.counts)
  }


#' Title
#'
#' @param cl.means.list 
#'
#' @return
#' @export
#'
#' @examples
build_dend_with_means <- function(cl.means.list)
{
  levels = unique(unlist(lapply(cl.means.list, colnames)))
  n.counts = tmp.cor=matrix(0, nrow=length(levels), ncol=length(levels))
  row.names(n.counts) = row.names(tmp.cor)=levels
  colnames(n.counts)=colnames(tmp.cor)=levels
  for(x in cl.means.list){
    tmp.cor[colnames(x),colnames(x)] = tmp.cor[colnames(x),colnames(x)] + cor(x)
    n.counts[colnames(x),colnames(x)] =   n.counts[colnames(x),colnames(x)] +1
  }
  tmp.cor = tmp.cor/n.counts
  tmp.cor[is.na(tmp.cor)]=0
  hclust(as.dist(1-tmp.cor))
}

#' Title
#'
#' @param dat 
#' @param impute.dat 
#'
#' @return
#' @export
#'
#' @examples
impute_val_cor <- function(dat, impute.dat)
  {
    gene.cor = pair_cor(dat, impute.dat)
    gene.cor[is.na(gene.cor)] = 0
    return(gene.cor)
  }


#' Title
#'
#' @param dat.list 
#' @param select.genes 
#' @param select.cells 
#' @param pairs 
#'
#' @return
#' @export
#'
#' @examples
gene_gene_cor_conservation <- function(dat.list, select.genes, select.cells,pairs=NULL)
  {
    sets = names(dat.list)
    gene.cor.list = sapply(sets, function(set){
      dat = dat.list[[set]]
      gene.cor = cor(t(as.matrix(dat[select.genes,intersect(colnames(dat),select.cells)])))
      gene.cor[is.na(gene.cor)] = 0
      gene.cor
    },simplify=F)
    if(is.null(pairs)){
      n.sets = length(sets)	
      pairs = cbind(rep(sets, rep(n.sets,n.sets)), rep(sets, n.sets))
      pairs = pairs[pairs[,1]<pairs[,2],,drop=F]
    }
    gene.cor.mat= sapply(1:nrow(pairs), function(i){
      p = pairs[i,]
      pair_cor(gene.cor.list[[p[1]]], gene.cor.list[[p[2]]])
    })
    colnames(gene.cor.mat) = paste0(pairs[,1],":",pairs[,2])
    return(gene.cor.mat)
  }



#' Title
#'
#' @param dat.list 
#' @param cl 
#' @param de.param.list 
#' @param prefix 
#' @param common.genes 
#' @param comb.de.genes 
#' @param cl.means.list 
#' @param col.list 
#' @param cl.col 
#' @param select.genes 
#' @param save.matrix 
#' @param n.markers 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_markers <- function(dat.list, cl,  de.param.list,prefix, common.genes, comb.de.genes=NULL, ref.set=names(dat.list)[[1]], cl.means.list=NULL, col.list=NULL, cl.col=NULL, select.genes=NULL, save.matrix=FALSE,n.markers = 20,...)
  {
    sets=names(dat.list)
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    blue.red <-colorRampPalette(c("blue", "white", "red"))
    if(is.null(cl.means.list)){
      cl.means.list = get_cl_means_list(dat.list, de.param.list = de.param.list, cl=cl)
    }
    if(is.null(select.genes)){
      if(is.null(comb.de.genes)){
        de.genes.list = sapply(sets, function(set){
          dat = dat.list[[set]]
          tmp.cells=  intersect(names(cl), colnames(dat))
          if(length(tmp.cells)==0){
            return(NULL)
          }
          tmp.cl = cl[tmp.cells]
          if(is.factor(tmp.cl)){
            tmp.cl = droplevels(tmp.cl)
          }
          tmp=display_cl(cl=tmp.cl, norm.dat=dat, max.cl.size = 200, de.param=de.param.list[[set]], n.markers=n.markers)$de.genes
        },simplify=F)
        comb.de.genes = comb_de_result(de.genes.list, cl.means.list, common.genes=common.genes)
      }
      select.genes = select_markers(dat.list[[1]], cl, n.markers=n.markers, de.genes=comb.de.genes)$markers
    }
    gene.hc = hclust(dist(cl.means.list[[ref.set]][select.genes,]), method="ward.D")
    if(is.null(cl.col)){
      cl.col = jet.colors(length(unique(cl)))
    }
    cl.col = cl.col[as.factor(cl)]
    cl.col =t(as.matrix(cl.col, ncol=1))
    colnames(cl.col)= names(cl)
    if(save.matrix){
      dat.matrix = list()
    }
    else{
      dat.matrix=NULL
    }
    for(set in sets){
      dat = dat.list[[set]]
      tmp.cl = cl[intersect(names(cl), colnames(dat))]
      if(is.factor(tmp.cl)){
        tmp.cl=droplevels(tmp.cl)
      }
      tmp.cells=  sample_cells(tmp.cl, 100)
      tmp.cl = tmp.cl[tmp.cells]
      tmp.col = t(as.matrix(cl.col[,tmp.cells]))
      if(!is.null(col.list) & !is.null(col.list[[set]])){
        tmp.col = rbind(tmp.col,col.list[[set]][,tmp.cells])
      }
      cells = plot_cl_heatmap(dat, cl=tmp.cl, markers= select.genes, gene.hc=gene.hc, prefix=paste(prefix, set, "markers.pdf", sep="."),ColSideColors=tmp.col,...)
      if(save.matrix){
        dat.matrix[[set]] = dat[gene.hc$labels[gene.hc$order], cells]
      }
    }
    return(list(select.genes=select.genes, dat.matrix = dat.matrix, comb.de.genes= comb.de.genes,gene.hc=gene.hc))
  }



get_incident_matrix <- function(cl1, cl2, consensus.cl,pseudo=0.01)
{
  tb1 = table(cl1, consensus.cl[names(cl1)]) + pseudo
  freq1 = as.matrix(tb1/rowSums(tb1))
  tb2 = table(cl2, consensus.cl[names(cl2)]) + pseudo
  freq2 = as.matrix(tb2/rowSums(tb2))
  kl = matrix(0, nrow=nrow(freq1), ncol=nrow(freq2))
  row.names(kl) = row.names(freq1)
  colnames(kl) = row.names(freq2)
  library(entropy)
  for(i in 1:nrow(freq1))
    for(j in 1:nrow(freq2)){
      kl[i,j] = KL.plugin(freq1[i,], freq2[j,]) 
    }

  kl.exp = exp(-kl)
  kl.sim = kl.exp/rowSums(kl.exp)
  
  pairs.df = melt(kl)
  colnames(pairs.df)=c("cl1","cl2","KL")
  pairs.df$sim = as.vector(kl.sim)
  return(list(kl, kl.sim,pairs.df))
}


