fast_tsne <- function(dat,theta=0.05, nthreads=12, perplexity=20, fast.tsne.path="~/src/FIt-SNE", ...){
  source(file.path(fast.tsne.path, "fast_tsne.R"))
  fast.tsne.df <- fftRtsne(dat, theta=theta,nthreads = nthreads, perplexity=perplexity, fast_tsne_path=file.path(fast.tsne.path,"bin/fast_tsne"),...)
  row.names(fast.tsne.df)= row.names(dat)
  colnames(fast.tsne.df)=c("Lim1","Lim2")
  fast.tsne.df = as.data.frame(fast.tsne.df)
}


plot_tsne_cl <- function(norm.dat, select.genes, cl, cl.df, tsne.df = NULL, show.legend=FALSE, cex=0.15, fn.size=2, alpha.val=1, ...)
  {
    library(ggplot2)
    require(Rtsne)
    if(is.null(tsne.df)){
      tsne.result = Rtsne(t(as.matrix(norm.dat[select.genes,names(cl)])),...)$Y
      row.names(tsne.result)=names(cl)
      tsne.df = as.data.frame(tsne.result[names(cl),])
      colnames(tsne.df)=c("Lim1","Lim2")
    }
    cl.color = setNames(cl.df$cluster_color, cl.df$cluster_label)
    cl.label = setNames(cl.df$cluster_label, row.names(cl.df))
    
    g <- plot_RD_cl(tsne.df, cl, cl.color, cl.label, cex=cex, fn.size =fn.size, alpha.val=alpha.val,show.legend=show.legend)
    return(list(tsne.df=tsne.df, g=g))    
  }

###Backward compatible to with the original function definition.
plot_tSNE_meta <- function(tsne.df, ...)
  {
    plot_RD_meta(tsne.df, ...)
  }

###Backward compatible to with the original function definition.
plot_tSNE_gene <- function(tsne.df, ...)
  {
    plot_RD_gene(tsne.df, ...)
  }


