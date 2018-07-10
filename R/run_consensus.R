#' Wrapper function to repeatively run clustering on subsampled cells and infer consensus clusters
#'
#' @param norm.dat normalized expression data matrix in log transform, using genes as rows, and cells and columns. Users can use log2(FPKM+1) or log2(CPM+1).
#' @param select.cells The cells to be clustered. Default: columns of norm.dat
#' @param niter The number of iteractions to run. Default 100.
#' @param sample.frac The fraction of of cells sampled per run. Default: 0.8. 
#' @param de.param The differential gene expression threshold. See de_param() function for details. 
#' @param output_dir The output directory to store clutering results for each iteraction. 
#' @param mc.cores The number of cores to be used for parallel processing. 
#' @param override binary variable determine if the clustering results already stored in output_dir should be overriden.  
#' @param init.result The pre-set high level clusters. If set, the function will only find finer splits of the current clusters.  
#' @param ... Other parameters passed to iter_clust
#'
#' @return
#' @export
#'
#' @examples
run_consensus_clust <- function(norm.dat, select.cells=colnames(norm.dat), niter=100, sample.frac=0.8, de.param=de_param(), output_dir="subsample_result",mc.cores=1, override=FALSE, init.result=NULL, ...)
{
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  all.cells=select.cells
  if(!is.null(init.result)){
    all.cells= intersect(all.cells, names(init.result$cl))
  }
  run <- function(i,...){
    prefix = paste("iter",i,sep=".")
    print(prefix)
    outfile= file.path(output_dir, paste0("result.",i,".rda"))
    if(file.exists(outfile)& !override){
      return(NULL)
    }
    select.cells=sample(all.cells, round(length(all.cells)*sample.frac))
    save(select.cells, file=file.path(output_dir, paste0("cells.",i,".rda")))
    result <- iter_clust(norm.dat=norm.dat, select.cells=select.cells,prefix=prefix, de.param = de.param, result=init.result, ...)
    save(result, file=outfile)
  }
  
  if (mc.cores==1){
    sapply(1:niter, function(i){run(i,...)})
  }
  else{
    require(foreach)
    require(doParallel)
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    foreach(i=1:niter, .combine='c') %dopar% run(i)
    stopCluster(cl)
  }
  result.files=file.path(output_dir, dir(output_dir, "result.*.rda"))
  if(length(all.cells) < 100000){
    co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
    co.ratio = co.result$co.ratio
    consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, norm.dat, select.cells=all.cells, de.param = de.param)
    refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio, tol.th=0.01, confusion.th=0.6, min.cells= de.param$min.cells)
    markers = consensus.result$markers
  }
  else{
    result <- iter_clust(norm.dat=norm.dat, select.cells=all.cells, de.param = de.param, ...)
    co.result <- collect_subsample_cl_matrix(norm.dat,result.files,all.cells)
    cl=merge_cl_by_co(result$cl, co.ratio=co.result$co.ratio, cl.mat=co.result$cl.mat, diff.th=0.25)
    refine.result = refine_cl(cl, cl.mat = co.result$cl.mat, tol.th=0.01, confusion.th=0.6, min.cells=de.param$min.cells)
    markers=result$markers      
  }
  cl = refine.result$cl
  merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat=Matrix::t(norm.dat[markers,]), de.param = de.param,return.markers=FALSE)
  return(list(co.result=co.result, cl.result=merge.result))
}
