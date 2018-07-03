#' Run bootstrapped, iterative consensus clustering
#'
#' @param norm.dat A matrix of normalized data
#' @param select.cells A selected subset of sample IDs. Default is colnames(norm.dat).
#' @param niter Number of bootstrapping rounds to run. Default is 100.
#' @param sample.frac Fraction of samples to use per bootstrapping round. Default is 0.8
#' @param de.param DE parameters set using `de_param()`.
#' @param output_dir Output directory for intermediate clustering results. Set to NULL to run in memory-only mode. This may take a significant quantity of RAM.
#' @param overwrite Overwrite flag for file output. Default is FALSE.
#' @param mc.cores Number of cores to use for parallel processing. If set to 1 (default), will use single-core mode. If "auto", will use detected cores - 1.
#' @param init.result Initial results.
#' @param ... Additional parameters passed to `iter_clust()`
#'
#'
run_consensus_clust <- function(norm.dat, 
                                select.cells = colnames(norm.dat), 
                                niter = 100, 
                                sample.frac = 0.8, 
                                de.param = de_param(), 
                                output_dir = "subsample_result",
                                overwrite = FALSE, 
                                mc.cores = 1, 
                                init.result = NULL,
                                ...)
  {
    
    if(mc.cores == "auto") {
      library(foreach)
      library(doParallel)
      
      mc.cores <- detectCores() - 1
    } else if(mc.cores > 1) {
      library(foreach)
      library(doParallel)
    }
  
    if(!dir.exists(output_dir) & !is.null(output_dir)){
      dir.create(output_dir)
    }
  
    all.cells <- select.cells
    
    if(!is.null(init.result)){
      all.cells <- intersect(all.cells, names(init.result$cl))
    }
    
    run <- function(i,...){
      prefix <- paste("iter",i,sep=".")
      print(prefix)
      
      if(!is.null(output_dir)) {
        outfile <- file.path(output_dir, paste0("result.",i,".rda"))
        if(file.exists(outfile) & overwrite == FALSE){
          stop(paste("Output file",outfile,"exists and overwrite = FALSE. Stopping."))
        }
      }
      
      select.cells <- sample(all.cells, round(length(all.cells) * sample.frac))
      
      if(!is.null(output_dir)) {
        save(select.cells, file = file.path(output_dir, paste0("cells.",i,".rda")))
      }
      
      result <- iter_clust(norm.dat = norm.dat, 
                           select.cells = select.cells,
                           prefix = prefix, 
                           de.param = de.param, 
                           result = init.result, 
                           ...)
      
      if(!is.null(output_dir)) {
        save(result, file = outfile)
      }
      
      result
    }
    
    # If running in single-core mode, run sapply
    if(mc.cores == 1){
      results <- sapply(1:niter, function(i){run(i,...)})
    } else {
      # If running in multi-core mode, run doParallel
      cl <- makeCluster(mc.cores)
      registerDoParallel(cl)
      results <- foreach(i = 1:niter, .combine = 'c') %dopar% run(i)
      stopCluster(cl)
    }
    
    if(!is.null(output_dir)) {
      result.files <- file.path(output_dir, dir(output_dir, "result.*.rda"))
    }
    
    if(length(all.cells) < 100000) {
      co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
      co.ratio <- co.result$co.ratio
      consensus.result <- iter_consensus_clust(co.ratio, 
                                               co.result$cl.list, 
                                               norm.dat, 
                                               select.cells = all.cells, 
                                               de.param = de.param)
      refine.result <- refine_cl(consensus.result$cl, 
                                 co.ratio = co.ratio, 
                                 tol.th = 0.01, 
                                 confusion.th = 0.6)
      markers <- consensus.result$markers
    } else {
      result <- iter_clust(norm.dat = norm.dat, 
                           select.cells = all.cells, 
                           de.param = de.param, 
                           ...)
      co.result <- collect_subsample_cl_matrix(norm.dat, 
                                               result.files, 
                                               all.cells)
      cl <- merge_cl_by_co(result$cl, 
                           co.ratio = co.result$co.ratio, 
                           cl.mat = co.result$cl.mat, 
                           diff.th = 0.25)
      refine.result <- refine_cl(cl, 
                                cl.mat = co.result$cl.mat, 
                                tol.th = 0.01, 
                                confusion.th = 0.6)
      markers <- result$markers      
    }
    cl <- refine.result$cl
    merge.result <- merge_cl(norm.dat = norm.dat, 
                             cl = cl, 
                             rd.dat = t(norm.dat[markers,]), 
                             de.param = de.param,
                             return.markers = FALSE)
    
    return(list(co.result = co.result, 
                cl.result = merge.result))
  }



