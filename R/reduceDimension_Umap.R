#' @importFrom reticulate py_module_available py_set_seed import
#'
#' @rdname RunUMAP
#' @method RunUMAP default
#' @export
#'
rd_Umap <- function(
                    mat,                             
                    n.neighbors = 30L,
                    n.components = 2L,
                    metric = "correlation",
                    n.epochs = NULL,
                    learning.rate = 1.0,
                    min.dist = 0.3,
                    spread = 1.0,
                    set.op.mix.ratio = 1.0,
                    local.connectivity = 1L,
                    repulsion.strength = 1,
                    negative.sample.rate = 5,
                    a = NULL,
                    b = NULL,
                    seed.use = 42,
                    metric.kwds = NULL,
                    angular.rp.forest = FALSE,
                    reduction.key = 'UMAP_',
                    verbose = TRUE)
{
  library(reticulate)
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs_as.integer(x = n.epochs)
  }
  umap_import_import(module = "umap", delay_load = TRUE)
  umap_umap_import$UMAP(
                        n_neighbors = as.integer(x = n.neighbors),
                        n_components = as.integer(x = n.components),
                        metric = metric,
                        n_epochs = n.epochs,
                        learning_rate = learning.rate,
                        min_dist = min.dist,
                        spread = spread,
                        set_op_mix_ratio = set.op.mix.ratio,
                        local_connectivity = local.connectivity,
                        repulsion_strength = repulsion.strength,
                        negative_sample_rate = negative.sample.rate,
                        a = a,
                        b = b,
                        metric_kwds = metric.kwds,
                        angular_rp_forest = angular.rp.forest,
                        verbose = verbose)
  umap_output= umap$fit_transform(mat)
  colnames(umap_output)=paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(umap_output)=rownames(object)
  return(umap_output)
}
