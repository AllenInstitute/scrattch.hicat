fast_tsne <- function(dat,theta=0.05, nthreads=12, perplexity=20, fast.tsne.path="~/src/FIt-SNE", ...){
  source(file.path(fast.tsne.path, "fast_tsne.R"))
  fast.tsne.df <- fftRtsne(dat, theta=theta,nthreads = nthreads, perplexity=perplexity, fast_tsne_path=file.path(fast.tsne.path,"bin/fast_tsne"),...)
  row.names(fast.tsne.df)= row.names(dat)
  colnames(fast.tsne.df)=c("Lim1","Lim2")
  fast.tsne.df = as.data.frame(fast.tsne.df)
}

