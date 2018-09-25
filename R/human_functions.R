# 1. check_neun

#' Check whether clusters are neun/non-neun
#'
#' @param anno anno dataframe with the following required colnames - "facs_population_plan", "cluster"
#' @param cluster "cluster" column from anno dataframe
#' @param neun.thresh default is 0.5
#'
#' @return returns all the clusters and true/false of whether they are neun positive
#' @export
#'
#' @examples check_neun(anno, anno$cluster, neun.thresh = 0.5)
check_neun <- function(anno, cluster, neun.thresh = 0.5) 
{
  neun.pos.cnt <- tapply(anno$facs_population_plan, cluster, 
                         function(x) {
                           sum(grepl("NeuN-pos", x, ignore.case = TRUE), na.rm = TRUE)/length(x[!is.na(x)])
                         })
  neun.check <- neun.pos.cnt > neun.thresh
  return(neun.check)
}

#---------------------------------------------------------------------------------------------------------------------------
# 2. check_qc

#' Check qc of some numeric metrics
#'
#' @param x numeric metric to check qc on
#' @param qc.iqr.mult default is 3
#'
#' @return returns binary result of whether qc failed or passed
#' @export
#'
#' @examples check_qc(x,  qc.iqr.mult = 3)
check_qc <- function(x, qc.iqr.mult = 3) 
{
  med.diff <- abs(x - median(x, na.rm = TRUE))
  iqr <- abs(quantile(x, 0.75, na.rm = TRUE) - quantile(x, 
                                                        0.25, na.rm = TRUE))
  outlier.check <- ifelse(med.diff > qc.iqr.mult * iqr, 1, 
                          0)
  return(outlier.check)
}

#---------------------------------------------------------------------------------------------------------------------------
# 3. check_outlier

#' check for outlier clusters
#'
#' @param anno anno dataframe with the following required columns - facs_population_plan", "cluster"
#' @param cluster clusters labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names and cpm normalized
#' @param select.cells column nmaes of norm.dat
#' @param keep.cl default NULL
#' @param neun.thresh default 0.5
#' @param qc.metrics required columns from anno dataframe. default is Genes.Detected.CPM", "percent_reads_aligned_total", "complexity_cg"
#' @param test.genes default is "SNAP25", "GAD1", "GAD2", "SLC17A7", "SLC17A6"
#' @param plot default is TRUE
#' @param plot.path path of plot, default is ./output/
#'
#' @return gives outlier clusters and exploratory plots
#' @export
#'
#' @examples check_outlier()
                                                                                                                 
check_outlier <- function(anno, cluster, norm.dat, 
                          select.cells = colnames(norm.dat), 
                          keep.cl = NULL, 
                          neun.thresh = 0.5, 
                          qc.metrics = c("Genes.Detected.CPM", "percent_reads_aligned_total", "complexity_cg"), 
                          test.genes = c("SNAP25",  "GAD1", "GAD2", "SLC17A7", "SLC17A6"), 
                          plot = TRUE, 
                          plot.path = "output/") 
{
  select.cells <- intersect(names(cluster), select.cells)
  select.id <- match(select.cells, colnames(norm.dat))
  norm.dat <- norm.dat[, select.id]
  anno <- droplevels(anno[select.id, ])
  neuronal.cl <- check_neun(anno, cluster, neun.thresh = neun.thresh)
  qc.metrics <- intersect(qc.metrics, colnames(anno))
  qc.median <- apply(anno[, qc.metrics], 2, function(x) {
    tapply(x, cluster, function(y) median(as.numeric(y), na.rm = TRUE))
  })
  qc.check <- apply(qc.median, 2, check_qc)
  qc.outlier <- neuronal.cl & apply(qc.check, 1, sum) > 0
  if (plot == TRUE & sd(qc.check) > 0) {
    pdf(file = paste0(plot.path, "/cluster_qc_check.pdf"), 
        height = nrow(qc.check)/4 + 1, width = 2, onefile = FALSE)
    pheatmap::pheatmap(qc.check, cluster_rows = FALSE, color = c("grey90",  "black"), legend = FALSE)
    dev.off()
  }
  expr.cnt <- apply(norm.dat[test.genes, ], 1, function(x) {
    expr.thresh <- 3
    tapply(x, cluster, function(x) sum(x > expr.thresh)/length(x))
  })
  expr.outlier <- neuronal.cl & (expr.cnt[, "SNAP25"] < 0.5 | 
                                (expr.cnt[, "GAD1"] > 0.5 | expr.cnt[, "GAD2"] > 0.5) & 
                                (expr.cnt[, "SLC17A7"] > 0.5 | expr.cnt[, "SLC17A6"] > 0.5) | 
                                (expr.cnt[, "GAD1"] < 0.8 & expr.cnt[, "GAD2"] < 0.8 & expr.cnt[, "SLC17A7"] < 0.8 & expr.cnt[,"SLC17A6"] < 0.8))
  if (plot == TRUE) {
    pdf(file = paste0(plot.path, "/cluster_marker_check.pdf"), 
        width = 8, height = 8)
    gene.pairs <- list(c(1, 4), c(2, 4), c(3, 5))
    for (pair1 in gene.pairs) {
      plot(expr.cnt[, pair1], type = "n")
      text(expr.cnt[, pair1], labels = rownames(expr.cnt))
    }
    dev.off()
  }
  outlier.cl <- unique(c(names(qc.outlier)[qc.outlier == TRUE], 
                         names(expr.outlier)[expr.outlier == TRUE]))
  outlier.cl.final <- setdiff(outlier.cl, keep.cl)
  return(outlier.cl.final)
}

#---------------------------------------------------------------------------------------------------------------------------
# 4. group_cl

#' Title
#'
#' @param anno anno dataframe with the following required columns - facs_population_plan", "cluster"
#' @param cluster clusters labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names and cpm normalized
#' @param select.cells column nmaes of norm.dat
#' @param neun.thresh default 0.5
#' @param outlier.cl default is NULL
#' @param donor.cl default is NULL
#'
#' @return gives the clusters names in each class (exc, inh, etc. ) as a list
#' @export
#'
#' @examples group_cl()
group_cl <- function(anno, cluster, norm.dat, select.cells = colnames(norm.dat), 
                      neun.thresh = 0.5, outlier.cl = NULL, donor.cl = NULL) 
{
  select.cells <- intersect(names(cluster), select.cells)
  select.id <- match(select.cells, colnames(norm.dat))
  norm.dat <- norm.dat[, select.id]
  anno <- droplevels(anno[select.id, ])
  neuronal.cl <- check_neun(anno, cluster, neun.thresh = neun.thresh)
  test.genes = c("GAD1", "GAD2", "SLC17A7", "SLC17A6")
  expr.cnt <- apply(norm.dat[test.genes, ], 1, function(x) {
    expr.thresh <- 3
    tapply(x, cluster, function(x) sum(x > expr.thresh)/length(x))
  })
  cl.class <- list()
  cl.class[["inh"]] <- setdiff(names(neuronal.cl)[neuronal.cl & 
                                                    (expr.cnt[, "GAD1"] > 0.5 | expr.cnt[, "GAD2"] > 0.5)], 
                               c(outlier.cl, donor.cl))
  cl.class[["exc"]] <- setdiff(names(neuronal.cl)[neuronal.cl & 
                                                    (expr.cnt[, "SLC17A7"] > 0.2 | expr.cnt[, "SLC17A6"] > 
                                                       0.2)], c(outlier.cl, donor.cl, cl.class[["inh"]]))
  cl.class[["glia"]] <- setdiff(names(neuronal.cl)[!neuronal.cl], 
                                c(outlier.cl, donor.cl))
  cl.class[["donor"]] <- setdiff(donor.cl, outlier.cl)
  cl.class[["outlier"]] <- setdiff(unique(cluster), unlist(cl.class[c("inh", 
                                                                      "exc", "glia", "donor")]))
  return(cl.class)
}


#---------------------------------------------------------------------------------------------------------------------------
# 5. check_donor

#' Title
#'
#' @param anno anno dataframe with the following required columns - facs_population_plan", "cluster"
#' @param cluster clusters labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names and cpm normalized
#' @param select.cells column nmaes of norm.dat
#' @param keep.cl default is NULL
#' @param meta1 default is "roi" (colname should be present in anno dataframe)
#' @param meta2 default is "external_donor_name" (colname should be present in anno dataframe)
#' @param plot default is TRUE
#' @param plot.path path of plot, default is "./output/"
#'
#' @return Gives cluster ids vs donor expression/region heatmap
#' @export
#'
#' @examples check_donor()
check_donor <- function(anno, cluster, norm.dat, select.cells = names(cluster), 
                         keep.cl = NULL, meta1 = "roi", meta2 = "external_donor_name", 
                         plot = TRUE, plot.path = "output/") 
{
  select.cells <- intersect(names(cluster), select.cells)
  select.id <- match(select.cells, colnames(norm.dat))
  norm.dat <- norm.dat[, select.id]
  anno <- droplevels(anno[select.id, ])
  layer.donor <- as.matrix(table(anno[, meta1], anno[, meta2]))
  layer.donor.prop <- sweep(layer.donor, 1, rowSums(layer.donor), 
                            "/")
  if (plot == TRUE & sd(layer.donor.prop) > 0) {
    hm.colors <- colorRampPalette(c("white", RColorBrewer::brewer.pal(9, 
                                                                      "YlOrRd")))(100)
    pdf(file = paste0(plot.path, "/layer_donor_prop.pdf"), 
        width = nrow(layer.donor.prop)/6 + 1, height = ncol(layer.donor.prop)/6 + 
          1, onefile = FALSE)
    pheatmap::pheatmap(t(layer.donor.prop), cluster_rows = FALSE, 
                       cluster_cols = FALSE, color = hm.colors)
    dev.off()
  }
  layer.cl <- as.matrix(table(anno[, meta1], cluster))
  layer.cl.prop <- sweep(layer.cl, 2, colSums(layer.cl), "/")
  cl.donor.exp.cnt <- round(t(t(layer.cl.prop) %*% layer.donor))
  cl.donor.exp.prop <- sweep(cl.donor.exp.cnt, 2, colSums(cl.donor.exp.cnt), 
                             "/")
  cl.donor.obs.cnt <- t(as.matrix(table(cluster, anno[, meta2])))
  cl.donor.obs.prop <- sweep(cl.donor.obs.cnt, 2, colSums(cl.donor.obs.cnt), 
                             "/")
  if (plot == TRUE & sd(cl.donor.obs.prop) > 0) {
    hm.colors <- colorRampPalette(c("white", RColorBrewer::brewer.pal(9, 
                                                                      "YlOrRd")))(100)
    pdf(file = paste0(plot.path, "/cluster_donor_prop.pdf"), 
        width = nrow(cl.donor.obs.prop)/6 + 1, height = ncol(cl.donor.obs.prop)/6 + 
          1, onefile = FALSE)
    pheatmap::pheatmap(t(cl.donor.obs.prop), cluster_rows = FALSE, 
                       cluster_cols = FALSE, color = hm.colors)
    dev.off()
  }
  cl.donor.prop.dev <- cl.donor.obs.prop - cl.donor.exp.prop
  cl.donor.prop.dev.bin <- ifelse(cl.donor.prop.dev > 0.5, 
                                  1, ifelse(cl.donor.prop.dev < -0.5, -1, ifelse(cl.donor.obs.prop > 
                                                                                   0.9, 1, 0)))
  donor.cl <- names(which(apply(cl.donor.prop.dev.bin, 2, function(x) any(x == 
                                                                            1))))
  donor.cl.final <- setdiff(donor.cl, keep.cl)
  return(donor.cl.final)
}

#---------------------------------------------------------------------------------------------------------------------------