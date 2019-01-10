source("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat_big_no_names/R/merge_cl.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/transcript_analysis/scrattch.hicat/R/annotate.R")


jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))




###comb.dat include the following elements
###dat.list a list of data matrix
###ref.list a list of sampled cells used as reference for each dataset.
###ref.de.param.list the DE gene criteria for each reference dataset (optional)
###meta.df merged meta data for all datasets. 
###cl.list clusters for each dataset (optional)
###cl.df.list cluster annotations for each dataset (optional) 


prepare_joint  <- function(dat.list, ref.sets, ref.list=NULL,ref.de.param.list=NULL, meta.df=NULL, cl.list=NULL, cl.df.list = NULL, de.genes.list=NULL, sampled.size = 10000, cl.sample.size = 200)
  {
    all.cells = unlist(lapply(dat.list, colnames))
    common.genes = row.names(dat.list[[1]])
    for(x in 2:length(dat.list)){
      common.genes= intersect(common.genes, row.names(dat.list[[x]]))
    }
    for(x in names(dat.list)){
      colnames(dat.list[[x]]) = paste(x, colnames(dat.list[[x]]), sep=".")
    }
    if(!is.null(cl.list)){
        for(x in names(cl.list)){
          names(cl.list[[x]]) = paste(x, names(cl.list[[x]]), sep=".")
        }
    }
    if(is.null(ref.list)){
      ref.list = sapply(ref.sets, function(ref.set){
        all.cells = colnames(dat.list[[ref.set]])
        if(is.null(cl.list[[ref.set]])){
          sample(all.cells, pmin(sample.size, length(all.cells)))
        }
        else{
          sample_cells(cl.list[[ref.set]], cl.sample.size)
        }
      },simplify=F)
    }
    if(is.null(ref.de.param.list)){
      ref.de.param.list = sapply(ref.sets, function(ref.set) de_param(), simplify=F)
    }
    platform = do.call("c",lapply(names(dat.list), function(p){
      dat = dat.list[[p]]
      setNames(rep(p, ncol(dat)), colnames(dat))
    }))
    gene.counts <- do.call("c",lapply(names(dat.list), function(p){
      dat = dat.list[[p]]
      setNames(colSums(dat > 0), colnames(dat))
    }))
    meta.df = cbind(meta.df, data.frame(platform, gene.counts))
    comb.dat = list(dat.list=dat.list, ref.list = ref.list, meta.df = meta.df, ref.de.param.list = ref.de.param.list, cl.list=cl.list, cl.df.list = cl.df.list, de.genes.list = de.genes.list(), common.genes=common.genes, all.cells= all.cells)
  }

get_knn <- function(ref.dat, dat, k, method ="cor")
  {
    if(method=="cor"){
      knn.index = knn_cor(ref.dat, dat,k=k)  
    }
    else if(method=="RANN"){
      knn.index = RANN::nn2(t(ref.dat), t(dat), k=k)[[1]]
    }
    else{
      stop(paste(method, "method unknown"))
    }
    return(knn.index)
  }


select_joint_genes  <-  function(comb.dat, maxGenes=2000, vg.padj.th=0.5, max.dim=20)
  {
    attach(comb.dat)
    select.genes = lapply(names(ref.list), function(ref.set){
      ref.cells = ref.list[[ref.set]]
      tmp.cells=  intersect(select.cells, ref.cells)
      de.param = ref.de.param.list[[ref.set]]
      ###if cluster membership is available, use cluster DE genes
      if(!is.null(cl.list[[ref.set]])){
        cl = droplevels(cl.list[[ref.set]][tmp.cells])
        cl.size = table(cl)
        cl = droplevels(cl[cl %in% names(cl.size)[cl.size > de.param$min.cells]])
        if(length(levels(cl))==1){
          return(NULL)
        }
        de.genes = de.genes.list[[ref.set]]
        select.genes = display_cl(cl, norm.dat=ref.dat, max.cl.size = 200, de.param=de.param, n.markers=50, de.genes= de.genes)$markers
        select.genes = intersect(select.genes, common.genes)
      }
      ####if cluster membership is not available, use high variance genes and genes with top PCA loading. 
      else{
        norm.dat = dat.list[[ref.set]][common.genes,ref.cells]
        counts =  2^norm.dat -1
        vg = findVG(as.matrix(counts[select.genes,sampled.cells]),plot.fig=plot.fig)
        select.genes = row.names(vg)[which(vg$loess.padj < vg.padj.th | vg$dispersion >3)]
        select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
        rd = rd_PCA(norm.dat,select.genes, select.cells, max.pca = max.dim)
        gene.rank = apply(-abs(rd$pca$rotation[,i]), 2, rank) 
        select.genes = row.names(gene.rank)[rowMins(gene.rank) < 100]
      }
    })
    gene.score = table(unlist(select.genes))
    if(length(gene.score)==0){
      return(NULL)
    }
    select.genes= names(head(sort(gene.score, decreasing=T), maxGenes))
    return(select.genes)
  }

knn_joint <- function(comb.dat, select.sets= names(comb.dat$dat.list),select.cells=comb.dat$all.cells, select.genes=NULL, method="cor", k=15,  n.sampled.cells = 10000,...)
{	  
  attach(comb.dat)
  ###Select genes for joint analysis
  if(is.null(select.genes)){
    select.genes = select_joint_genes(comb.dat, ...)
  }
  cat("Number of select genes", length(select.genes), "\n")
  cat("Get knn\n")
  ###index is the index of knn from all the cells
  knn.comb = do.call("rbind",lapply(select.sets, function(set){
    dat = dat.list[[set]]
    tmp.cells=  intersect(select.cells, colnames(dat))
    print(length(tmp.cells))
    do.call("cbind", lapply(names(ref.list), function(ref.set){
      ref.dat = dat.list[[ref.set]][select.genes, ref.list[[ref.set]]]
      knn=get_knn(ref.dat, dat[select.genes, tmp.cells], k, method = method)
      idx = match(colnames(ref.dat), all.cells)
      knn = matrix(idx[knn], nrow=nrow(knn))
      row.names(knn) = tmp.cells
    }))
  }))
  sets = setNames(rep(names(dat.list), sapply(dat.list, ncol)), row.names(knn.comb))
  sampled.cells = union(sample_cells(sets, n.sampled.cells), ref.cells)
  result = knn_jaccard_louvain(knn.comb[sampled.cells,])
  result$cl.mat = t(result$louvain$memberships)
  row.names(result$cl.mat) = sampled.cells
  result$knn = knn.comb
  cl = setNames(result$cl.mat[,1], row.names(result$cl.mat))
  rm.cells = c()
  for(ref.set in names(ref.list)){
    tmp.cells=  intersect(names(cl), ref.list[[ref.set]])
    tmp.cl = cl[tmp.cells]
    ###these cells are present in clusters missing cells from the referece dataset. 
    rm.cells = names(cl)[!cl %in% tmp.cl]
    cl = cl[cl %in% tmp.cl]
    cat("Merge clusters\n")
    rd.dat.t = ref.dat[select.genes, names(tmp.cl)]
    de.param = ref.de.param.list[[ref.set]]
    merge.result = merge_cl(ref.dat[,names(tmp.cl)], cl, rd.dat.t = rd.dat.t, de.param = de.param, sampled=TRUE)
    if(is.null(merge.result)){
      return(NULL)
    }
    cl = merge.result$cl
  }
  if(length(cl) < nrow(result$knn)){
    pred.df = predict_knn(result$knn, all.cells, cl)
    pred.cl= setNames(pred.df$pred.cl, row.names(pred.df))
    cl = c(cl, pred.cl[setdiff(names(pred.cl), names(cl))])
  }
  result$cl = cl
  result$markers = select.genes
  result$select.genes= select.genes
  result$rm.cells = rm.cells
  return(result)
}

sim_knn <- function(sim, k=15)
{
  
  th =  rowOrderStats(as.matrix(sim), which=ncol(sim)-k+1)
  select = sim >= th
  knn.idx = t(apply(select, 1, function(x)head(which(x),k)))
  return(knn.idx)
}

knn_cor <- function(ref.dat, query.dat, k = 15)
{
  require(matrixStats)
  tmp.cor = cor(as.matrix(query.dat), as.matrix(ref.dat))
  tmp.cor[is.na(tmp.cor)] = 0
  knn.idx = sim_knn(tmp.cor, k=k)
  return(knn.idx)
}



jaccard2 <- function(m) {
  library(Matrix)
  
  ## common values:
  A <-  tcrossprod(m)
  B <- as(A, "dgTMatrix")
  
  ## counts for each row
  b <- Matrix::rowSums(m)  
  
   
  ## Jacard formula: #common / (#i + #j - #common)
  x = B@x / (b[B@i+1] + b[B@j+1] - B@x)
  B@x = x
  return(B)
}


knn_jaccard <- function(knn.index)
  {
    knn.df = data.frame(i = rep(1:nrow(knn.index), ncol(knn.index)), j=as.vector(knn.index))
    knn.mat = sparseMatrix(i = knn.df[[1]], j=knn.df[[2]], x=1)
    sim= jaccard2(knn.mat)
    row.names(sim) = colnames(sim) = row.names(knn.index)
    return(sim)
  }


knn_jaccard_louvain <- function(knn.index)
  {
    require(igraph)
    cat("Get jaccard\n")
    sim=knn_jaccard(knn.index)
    cat("Louvain clustering\n")
    gr <- igraph::graph.adjacency(sim, mode = "undirected", 
                                  weighted = TRUE)
    result <- igraph::cluster_louvain(gr)
    return(list(sim=sim,louvain=result))
  }


predict_knn <- function(knn.idx, reference, cl)
  {
    query = row.names(knn.idx)
    df = data.frame(nn=as.vector(knn.idx), query=rep(row.names(knn.idx), ncol(knn.idx)))
    df$nn.cl = cl[reference[df$nn]]
    tb=with(df, table(query, nn.cl))
    pred.cl = setNames(colnames(tb)[apply(tb, 1, which.max)], row.names(tb))
    pred.score = setNames(rowMaxs(tb)/ncol(knn.idx), row.names(tb))
    pred.df = data.frame(pred.cl, pred.score)
    return(pred.df)
  }



impute_knn <- function(knn.idx, reference, dat)
  {
    query = row.names(knn.idx)
    impute.dat= sapply(1:ncol(dat), function(x){
      print(x)
      tmp.dat = sapply(1:ncol(knn.idx), function(i){
        dat[reference[knn.idx[,i]],x]
      })
      rowMeans(tmp.dat)
    })
    row.names(impute.dat) = row.names(knn.idx)
    colnames(impute.dat) = colnames(dat)
    return(impute.dat)
  }



process <- function(comb.dat, select.cells, prefix, ...)
  {
    result = knn_joint(comb.dat,select.cells = select.cells, ...)
    if(is.null(result)){
      return(NULL)
    }
    print(table(result$cl))

    if(!is.null(cl.list)){
      select.ref = names(ref.list)[[1]]
      ref.cl = cl.list[[select.ref]]
      ref.cl.df = cl.df.list[[select.ref]]
      tmp = compare_annotate(result$cl, ref.cl, ref.cl.df)
      result$cl = tmp$cl
      result$cl.df = tmp$cl.df
    }
    else{
      cl.df = get_cl_df(result$cl)
      result$cl.df = cl.df
    }
    cl = result$cl
    tmp = factor(cl.df[as.character(cl),"cluster_label"], levels=cl.df$cluster_label)
    g = plot_cl_meta_barplot(setNames(meta.df[names(cl), "platform"], names(cl)), tmp)
    g = g + theme(axis.text.x = element_text(angle=45,hjust=1, vjust=1))
    ggsave(paste0(prefix, ".platform.barplot.pdf"),g)

    if(length(result$cl) < length(select.cells)){
      pred.df = predict_knn(result$knn, result$ref.cells, cl)
      pred.cl = setNames(factor(as.character(pred.df$pred.cl), levels = row.names(result$cl.df)),  row.names(pred.df))
      result$cl = pred.cl
    }
    #plot_confusion(result$cl, prefix)
    return(result)
  }

iter_process <- function(comb.dat, select.cells, prefix, result=NULL, ...)
  {
    if(is.null(result)){
      result = process(select.cells, prefix, ...)
      plot_confusion(result$cl, prefix=prefix)      
    }
    finer.results <- sapply(as.character(sort(unique(result$cl))), function(i){
      tmp.prefix=paste(prefix, i,sep=".")
      print(tmp.prefix)
      select.cells= names(result$cl)[result$cl == i]
      platform.size = table(meta.df[select.cells, "platform"])
      if(sum(platform.size >= 10) > 1){
        tmp.result = process(select.cells=select.cells,  prefix=tmp.prefix, ...)
      }
      else{
        return(NULL)
      }
    },simplify=F)
    new.result <- combine_knn_result(result, finer.result)
  }



combine_knn_result <- function(result, split.result)
  {
    ref.cells = result$ref.cells
    knn= result$knn
    cl = result$cl
    markers=result$markers
    n.cl = 0
    new.cl = setNames(rep(1, length(cl)), names(cl))
    knn.merge = knn
    for(x in as.character(sort(unique(cl)))){
      print(x)
      tmp.cells = names(cl)[cl==x]
      if(!x %in% names(split.result)){
        new.cl[tmp.cells]=n.cl+1
        n.cl = n.cl+1
      }
      else{
        print("finer split")       
        tmp.cl = split.result[[x]]$cl
        new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
        markers=union(markers, split.result[[x]]$markers)
        ####use finer split knn to replace the original knn
        orig.index = match(split.result[[x]]$ref.cells, result$ref.cells)
        tmp.knn = split.result[[x]]$knn[names(tmp.cl),]
        tmp.knn = matrix(orig.index[tmp.knn], nrow=nrow(tmp.knn))
        row.names(tmp.knn) = names(tmp.cl)
        knn.merge[row.names(tmp.knn),] = tmp.knn
        n.cl = max(new.cl)        
      }
      cat("New #", n.cl, "\n")
    }
    new.result = list(cl = new.cl, markers=markers, knn=knn)
    new.result$ref.cells = result$ref.cells
    new.result$cl.mat = cbind(result$cl, new.result$cl[names(result$cl)])
    new.result$knn.merge = cbind(result$knn, knn.merge)
    return(new.result)
  }




merge_knn_result <- function(split.results)
  {
    ref.cells = unlist(lapply(split.results, function(x)x$ref.cells))
    ref.cells = ref.cells[!duplicated(ref.cells)]
    markers =  unique(unlist(lapply(split.results, function(x)x$markers)))
    n.cl = 0
    cl = NULL
    cl.df = NULL
    knn = NULL
    knn.merge= NULL
    for(result in split.results){
      tmp.cl = setNames(as.integer(as.character(result$cl)) + n.cl, names(result$cl))
      tmp.cl.df = result$cl.df
      row.names(tmp.cl.df) = as.integer(row.names(tmp.cl.df)) + n.cl 
      cl = c(cl, tmp.cl)
      cl.df = rbind(cl.df, tmp.cl.df)
      n.cl = max(as.integer(as.character(cl)))
      orig.index = match(result$ref.cells, ref.cells)
      tmp.knn = result$knn[names(tmp.cl),]
      tmp.knn = matrix(orig.index[tmp.knn], nrow=nrow(tmp.knn))
      knn = rbind(knn, tmp.knn)
      tmp.knn = result$knn.merge[names(tmp.cl),]
      tmp.knn = matrix(orig.index[tmp.knn], nrow=nrow(tmp.knn))
      knn.merge = rbind(knn.merge, tmp.knn)
    }
    new.result = list(cl = as.factor(cl), cl.df = cl.df, markers=markers, knn=knn, ref.cells =ref.cells, knn.merge = knn.merge)
    return(new.result)
  }





plot_tsne <- function(sim, cl, cl.df, meta.df, prefix, tsne.df = NULL, cex=0.3, fn.size=2, adj.col=100)
  {
    require(Rtsne)
    if(is.null(tsne.df)){
      sim= as.matrix(sim)
      tsne.result = Rtsne(1-sim, is_distance=TRUE, theta=0.1)
      tsne.df = as.data.frame(tsne.result$Y)
      colnames(tsne.df) = c("Lim1", "Lim2")
      row.names(tsne.df) = row.names(sim)
    }
    #cl.df$cluster_color = adjust_color(cl.df$cluster_color, adj.col)
    tmp = plot_tsne_cl(cl=cl, cl.df=cl.df,  tsne.df = tsne.df, cex=cex, fn.size = fn.size)
    tsne.df = tmp$tsne.df
    ggsave(paste("tsne.cl",prefix,"pdf", sep="."), tmp$g)

    tmp.df = meta.df[row.names(tsne.df),]
    tmp= setNames(as.factor(tmp.df$platform), row.names(tsne.df))
    meta.col=NULL
    if(length(levels(tmp))==2){
      meta.col = setNames(c("blue", "orange"), levels(tmp))
    }
    g= plot_tsne_meta(tsne.df, tmp, meta.col=meta.col)
    ggsave(paste("tsne",prefix, "platform.pdf", sep="."), g, height=8, width=10)
    
    org.cl = cl[row.names(sim)]
    org.cl.df = cl.df
    
    for(x in intersect(sets,meta.df[row.names(tsne.df),"platform"])){
      load(file.path("../", x, "cl.final.rda"))
      names(cl) = paste(x, names(cl), sep=".")
      tmp.cl = droplevels(cl[intersect(names(cl), row.names(tsne.df))])
      cl.size = table(tmp.cl)
      select.cl = names(cl.size)[cl.size > 2 & cl.size > sum(cl.size) * 0.01]
      tmp.cl = droplevels(tmp.cl[tmp.cl %in% select.cl])
      tmp.cl.df = cl.df[levels(tmp.cl),]
      if(nrow(tmp.cl.df) > 1){
        tmp.cl.df$cluster_color = adjust_color(as.character(tmp.cl.df$cluster_color),adj.col)
      }
      g = plot_tsne_cl(cl=tmp.cl, cl.df=tmp.cl.df,  tsne.df = tsne.df[names(tmp.cl),], cex=cex, fn.size = fn.size)$g
      ggsave(paste("tsne.",prefix, x, "pdf", sep="."), g)
    }

    #save(tsne.df, file=paste("tsne.df", prefix, "rda",sep="."))
    
    return(tsne.df)
  }



plot_confusion <- function(consensus.cl, prefix, comb.dat)
  {
    attach(comb.dat)
    for(x in names(cl.list)){
      names(cl)  = paste(x, names(cl), sep=".")
      ggsave(paste(prefix, x, "pdf", sep="."), compare_annotate(consensus.cl, cl, cl.df, rename = FALSE)$g)
    }
  }

plot_markers <- function(comb.dat, cl, prefix, maxGenes=500, cl.col=NULL, select.genes=NULL, save.matrix=FALSE,...)
  {
    attach(comb.dat)
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    blue.red <-colorRampPalette(c("blue", "white", "red"))
    de.result=NULL
    if(is.null(select.genes)){
      de.result = lapply(names(ref.list), function(ref.set){
        print(x)
        dat = dat.list[[ref.set]]
        tmp.cells=  intersect(names(cl), colnames(dat))
        tmp.cl = droplevels(cl[tmp.cells])
        print(table(tmp.cl))
        display_cl(tmp.cl, norm.dat=dat, max.cl.size = 200, de.param=ref.de.param.list[[ref.set]], n.markers=20)
      })
      gene.score = table(unlist(lapply(de.result, function(x)intersect(x$markers,common.genes))))
      select.genes= names(head(sort(gene.score, decreasing=T), maxGenes))
      if(!is.null(all.genes)){
        select.genes = intersect(select.genes, all.genes)
      }
    }
    tmp.dat= dat.list[[names(ref.list)[[1]]]]
    tmp.cl.means = get_cl_means(tmp.dat[select.markers,], cl[intersect(names(cl), colnames(tmp.dat))])
    gene.hc = hclust(dist(tmp.cl.means), method="ward.D")
    if(is.null(cl.col)){
      cl.col = jet.colors(length(unique(cl)))
    }
    cl.col = cl.col[as.factor(cl)]
    tmp.col =t(as.matrix(cl.col, ncol=1))
    colnames(tmp.col)= names(cl)
    if(save.matrix){
      dat.matrix = list()
    }
    else{
      dat.matrix=NULL
    }
    for(x in names(dat.list)){
      dat = dat.list[[x]]
      tmp.cells=  sample_cells(droplevels(cl[intersect(names(cl), colnames(dat))]), 100)
      tmp.cl = droplevels(cl[tmp.cells])      
      cells = plot_cl_heatmap(dat, cl=tmp.cl, markers= select.genes, gene.hc=gene.hc, prefix=paste(prefix, x, "markers.pdf", sep="."),ColSideColors=tmp.col,...)
      if(save.matrix){
        dat.matrix[[x]] = dat[gene.hc$labels[gene.hc$order], cells]
      }
    }
    return(list(select.genes=select.genes, dat.matrix = dat.matrix, de.result= de.result))
  }

