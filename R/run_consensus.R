run_consensus_clust <- function(norm.dat, niter=100, sample.frac=0.8, de.param=de_param(), output_dir="subsample_result",mc.cores=6, override=FALSE, ...)
  {
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    all.cells= colnames(norm.dat)
    run <- function(i,...){
      require(iterclust)    
      prefix = paste("iter",i,sep=".")
      outfile= file.path(output_dir, paste0("result.",i,".rda"))
      if(file.exists(outfile)& !override){
        return(NULL)
      }
      select.cells=sample(all.cells, round(length(all.cells)*sample.frac))
      save(select.cells, file=file.path(output_dir, paste0("cells.",i,".rda")))
      result <- iter_clust(norm.dat=norm.dat, select.cells=select.cells,prefix=prefix, de.param = de.param, ...)
      save(result, file=outfile)
    }
    if (mc.cores==1){
      sapply(1:n.iter, function(i){run(i,...)})
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
    co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
    co.ratio = co.result$co.ratio
    consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, norm.dat, select.cells=all.cells, de.param = de.param)
    refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio, tol.th=0.01, confusion.th=0.6)
    merge.result= merge_cl(norm.dat=norm.dat, cl=refine.result$cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE)
    return(list(co.result=co.result, cl.result=merge.result))
  }



