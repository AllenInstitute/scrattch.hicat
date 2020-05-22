#' Title
#'
#' @param dat 
#' @param theta 
#' @param nthreads 
#' @param perplexity 
#' @param fast.tsne.path 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fast_tsne <- function(dat,theta=0.05, nthreads=12, perplexity=20, fast.tsne.path="~/src/FIt-SNE", ...){
  source(file.path(fast.tsne.path, "fast_tsne.R"))
  fast.tsne.df <- fftRtsne(dat, theta=theta,nthreads = nthreads, perplexity=perplexity, fast_tsne_path=file.path(fast.tsne.path,"bin/fast_tsne"),...)
  row.names(fast.tsne.df)= row.names(dat)
  colnames(fast.tsne.df)=c("Lim1","Lim2")
  fast.tsne.df = as.data.frame(fast.tsne.df)
}


#' Title
#'
#' @param norm.dat 
#' @param select.genes 
#' @param cl 
#' @param cl.df 
#' @param tsne.df 
#' @param show.legend 
#' @param cex 
#' @param fn.size 
#' @param alpha.val 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param tsne.df 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_tSNE_meta <- function(tsne.df, ...)
  {
    plot_RD_meta(tsne.df, ...)
  }

###Backward compatible to with the original function definition.
#' Title
#'
#' @param tsne.df 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_tSNE_gene <- function(tsne.df, ...)
  {
    plot_RD_gene(tsne.df, ...)
  }


#' Title
#'
#' @param dat 
#' @param select.genes 
#' @param select.samples 
#' @param dims 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
build_tsne <- function(dat, select.genes = rownames(dat), select.samples = colnames(dat), dims=2, ...)
{
  library(Rtsne)
  if(length(select.genes)<2)   select.genes   <- 1:dim(dat)[1]
  if(length(select.samples)<2) select.samples <- 1:dim(dat)[2]
  tsne.result = Rtsne(t(as.matrix(dat[select.genes,select.samples])), dims=dims, ...)$Y
  row.names(tsne.result) = select.samples
  tsne.df = as.data.frame(tsne.result[select.samples, ])
  colnames(tsne.df) = letters[c(24:26,1:23)][1:dims]
  tsne.df$sample_name = rownames(tsne.df)
  tsne.df
}

