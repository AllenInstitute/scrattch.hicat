# rd_Rumap (dat, param)
#
# dat (nsample x nfeat) : full dim. data
# umap.param : umap parameter
#
# contact : changkyul@alleninstitute.org
#
#
# Default umap configuration parameters
#            n_neighbors: 15
#           n_components: 2
#                 metric: euclidean
#               n_epochs: 200
#                  input: data
#                   init: spectral
#               min_dist: 0.1
#       set_op_mix_ratio: 1
#     local_connectivity: 1
#              bandwidth: 1
#                  alpha: 1
#                  gamma: 1
#   negative_sample_rate: 5
#                      a: NA
#                      b: NA
#                 spread: 1
#           random_state: NA
#        transform_state: NA
#                    knn: NA
#            knn_repeats: 1
#                verbose: FALSE
#        umap_learn_args: NA
#
rd_Rumap <- function (dat, umap.param,method="umap") {
  if(method=="umap"){
    library(umap)
    # set up parameters                                        
    custom.config = umap.defaults
    if (length(umap.param) > 0) {
      inparam = names(umap.param)
      for (key in inparam) custom.config[[key]] = umap.param[[key]]
    }
    
    # run ummap
    if (custom.config$metric == "correlation") {
      dat.dist = (1 - cor(t(dat)))/2
      custom.config$metric = "euclidean"
      tmp = umap(dat.dist, config=custom.config, random_state=umap.param$random_state, input="dist")
    } else {
      tmp = umap(dat, config=custom.config, random_state=umap.param$random_state)
    }
    rd.umap = tmp$layout
  }

  return(rd.umap)
}

umap_param <- function(n_neighbors=25, metric="correlation", min_dist=0.4, random_seed=123)
  {
    umap.param <- list()
    umap.param$n_neighbors  = n_neighbors
    umap.param$metric       = metric
    umap.param$min_dist     = min_dist
    umap.param$random_state = random_seed
    return(umap.param)
  }

#
# PCA_umap (rd.dat.list, cl, Comb.dat, ref.set, N.sampled.cells, umap.param)
#
# rd.dat.list[[ref.set]] : dimension reduction data for each platform
# cl                     : cluster result
# Comb.dat               : combined data set 
# ref.set                : platform to be used for reference
# N.sampled.cells        : max number of cells selected in each cluster
# umap.param             : umap parameter
#
# contact : changkyul@alleninstitute.org
#
PCA_umap <- function(dat,
                     cl, 
                     rm.eigen=NULL,
                     rm.th=0.6,
                     umap.param = NULL,... )
{
  rd.dat = rd_PCA(dat, select.cells = names(cl),...)
  if(!is.null(rm.eigen)){
    rd.dat = filter_RD(rd.dat, rm.eigen, rm.th)
  }  
  ######################################
  # UMAP
  if (is.null(umap.param)) {
    umap.param= umap_param()    
  } 
  umap = rd_Rumap(rd.dat[sampled.cells,], umap.param=umap.param)
  return(umap)
}


