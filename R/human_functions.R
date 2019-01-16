# 1. check_neun

#' Check whether clusters are Neun+ or Neun- 
#' 
#' This function is very specific to brain cell/nuclear data where Neun is used as a
#'   stain for neuronal vs. non-neuronal identity.  This function determine's whether 
#'   a cluster is predicted to be neuronal or non-neuronal based on expression of Neun.
#'
#' @param anno anno dataframe which must include column name listed in `neun.colname`
#' @param cluster cluster labels for all cells along with sample_id as their names
#' @param neun.thresh fraction of cells expressing NeuN to be considered 
#'   NeuN positive (default is 0.5)
#' @param neun.colname column name in `anno` with the Neun information
#' @param neun.val value corresponding to non-neuronal marker in neun.colname in anno
#'
#' @return returns all the clusters and true/false of whether they are Neun positive
#' @export
#'
#' @examples check_neun(anno, anno$cluster, neun.thresh = 0.5)
check_neun <- function(anno, cluster,
                       neun.thresh = 0.5,
                       neun.colname = "facs_population_plan",
                       neun.val = "NeuN-pos") {
  neun.pos.cnt <- tapply(
    anno[, neun.colname], cluster,
    function(x) {
      sum(grepl(neun.val, x, ignore.case = TRUE), na.rm = TRUE) / length(x[!is.na(x)])
    }
  )
  neun.check <- neun.pos.cnt > neun.thresh
  return(neun.check)
}

#---------------------------------------------------------------------------------------------------------------------------
# 2. check_qc

#' Check qc of numeric metric
#' 
#' Identifies values in a numeric vector that are sufficiently higher than expected.
#'
#' @param x numeric vector corresponding to metric to check qc on
#' @param qc.iqr.mult How many interquartile ranges about the median must a value 
#'   be to be considered an outlier? (default is 3)
#'
#' @return returns binary result of whether qc failed or passed
#' @export
#'
#' @examples check_qc(x,  qc.iqr.mult = 3)
check_qc <- function(x, qc.iqr.mult = 3) {
  med.diff <- abs(x - median(x, na.rm = TRUE))
  iqr      <- abs(quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))
  outlier.check <- ifelse(med.diff > qc.iqr.mult * iqr, 1, 0)
  return(outlier.check)
}

#---------------------------------------------------------------------------------------------------------------------------
# 3. check_outlier

#' Check for outlier clusters
#' 
#' This function checks for outliers looks for unexpected combinations of marker gene 
#'   expression (e.g., GAD1 + SLC17A7) and for particularly high or low expression of
#'   indicated QC metrics, and flags any of the clusters meeting those criteria as 
#'   potential outliers.  This should (in theory) find things like poor quality 
#'   clusters and clusters of doublets.  Specific genes and thresholds currnetly
#'   hard-coded in, but might be updated in later iterations.
#'
#' @param anno anno dataframe which must include column names listed in `neun.colname`
#'   and `qc.metrics`.  "cluster" is added from `cluster` parameter below.
#' @param cluster cluster labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names 
#'   and cpm normalized
#' @param select.cells column nmaes of norm.dat
#' @param keep.cl clusters to definitely keep in analysis (e.g., to exclude from 
#'   consideration as an outlier cluster) default is NULL
#' @param neun.thresh fraction of cells expressing NeuN to be considered 
#'   NeuN positive (default is 0.5)
#' @param neun.colname column name in anno with the Nuen information
#' @param neun.val value corresponding to non-neuronal marker in neun.colname in anno
#' @param qc.metrics required columns from anno dataframe. default is Genes.Detected.CPM", 
#'   "percent_reads_aligned_total", "complexity_cg"
#' @param test.genes CURRENTLY NOT USED.  This function will eventually allow for a 
#'   pre-defined set of genes to be entered. default is "SNAP25", "GAD1", "GAD2", 
#'   "SLC17A7", "SLC17A6"
#' @param expr.th expression threshold for detecting test genes
#' @param prop.th proportion threshold of detected genes by cluster (default is 0.5, 0.5, 0.2, 0.2)
#' @param min.prop.th one of last 4 test genes should have detection at least at this amount
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
                          neun.colname = "facs_population_plan",
                          neun.val = "NeuN-pos", 
                          qc.metrics = c("Genes.Detected.CPM", "percent_reads_aligned_total", "complexity_cg"),
                          test.genes = c("SNAP25", "GAD1", "GAD2", "SLC17A7", "SLC17A6"),
                          expr.th = 3,
                          prop.th = c(0.5, 0.5, 0.5, 0.5, 0.5), 
                          min.prop.th = 0.8,
                          plot = TRUE,
                          plot.path = "output/") {
  select.cells <- intersect(names(cluster), select.cells)
  select.id    <- match(select.cells, colnames(norm.dat))
  norm.dat     <- norm.dat[, select.id]
  anno         <- droplevels(anno[select.id, ])
  anno$cluster <- cluster[select.cells]
  neuronal.cl  <- check_neun(anno, cluster, neun.thresh = neun.thresh, 
                             neun.colname = neun.colname, neun.val = neun.val)
  qc.metrics   <- intersect(qc.metrics, colnames(anno))
  qc.median    <- apply(anno[, qc.metrics], 2, function(x) {
    tapply(x, cluster, function(y) median(as.numeric(y), na.rm = TRUE))
  })
  qc.check     <- apply(qc.median, 2, check_qc)
  qc.outlier   <- neuronal.cl & apply(qc.check, 1, sum) > 0
  if (plot == TRUE & sd(qc.check) > 0) {
    pdf(
      file = paste0(plot.path, "/cluster_qc_check.pdf"),
      height = nrow(qc.check) / 4 + 1, width = 2, onefile = FALSE
    )
    pheatmap::pheatmap(qc.check, cluster_rows = FALSE, 
                       color = c("grey90", "black"), legend = FALSE)
    dev.off()
  }
  expr.cnt <- apply(norm.dat[test.genes, ], 1, function(x) {
    tapply(x, cluster, function(x) sum(x > expr.th) / length(x))
  })
  expr.outlier <- neuronal.cl & (expr.cnt[, test.genes[1]] < prop.th[1] |
    (expr.cnt[, test.genes[2]] > prop.th[2] | expr.cnt[, test.genes[3]] > prop.th[3]) &
      (expr.cnt[, test.genes[4]] > prop.th[4] | expr.cnt[, test.genes[5]] > prop.th[5]) |
    (expr.cnt[, test.genes[2]] < min.prop.th & expr.cnt[, test.genes[3]] < min.prop.th & 
       expr.cnt[, test.genes[4]] < min.prop.th & expr.cnt[, test.genes[5]] < min.prop.th))
  if (plot == TRUE) {
    pdf(
      file = paste0(plot.path, "/cluster_marker_check.pdf"),
      width = 8, height = 8
    )
    gene.pairs <- list(c(1, 4), c(2, 4), c(3, 5))
    for (pair1 in gene.pairs) {
      plot(expr.cnt[, pair1], type = "n")
      text(expr.cnt[, pair1], labels = rownames(expr.cnt))
    }
    dev.off()
  }
  outlier.cl <- unique(c(
    names(qc.outlier)[qc.outlier == TRUE],
    names(expr.outlier)[expr.outlier == TRUE]
  ))
  outlier.cl.final <- setdiff(outlier.cl, keep.cl)
  return(outlier.cl.final)
}

#---------------------------------------------------------------------------------------------------------------------------
# 4. group_cl

#' Assign clusters to a group
#' 
#' This function assigns each cluster to a predefined set of classes (exc, inh, glia, 
#'   donor, outlier) using test.genes (default GAD1, GAD2, SLC17A7, SLC17A6).
#'   Later iterations of the function will allow for user-defined classes and genes.
#'
#' @param anno anno dataframe which must include column names listed in `neun.colname`
#'   "cluster" is added from `cluster` parameter below.
#' @param cluster clusters labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names 
#'   and cpm normalized
#' @param select.cells column names of norm.dat
#' @param test.genes genes used to define cell classes (default is GAD1, GAD2, SLC17A7, SLC17A6)
#' @param expr.th expression threshold for detecting test genes
#' @param prop.th proportion threshold of detected genes by cluster (default is 0.5, 0.5, 0.2, 0.2)
#' @param neun.thresh fraction of cells expressing NeuN to be considered 
#'   NeuN positive (default is 0.5)
#' @param neun.colname column name in anno with the Neun information
#' @param neun.val value corresponding to non-neuronal marker in neun.colname in anno
#' @param outlier.cl vector of clusters previously defined as "outlier" (default is NULL)
#' @param donor.cl vector of clusters previously defined as "donor" (default is NULL)
#'
#' @return gives the clusters names in each class (exc, inh, glia, donor, outlier)
#'   as a list
#' @export
#'
#' @examples group_cl()
group_cl <- function(anno, cluster, norm.dat,
                     select.cells = colnames(norm.dat),
                     neun.thresh = 0.5,
                     neun.colname = "facs_population_plan",
                     neun.val = "NeuN-pos"
                     test.genes = c("GAD1", "GAD2", "SLC17A7", "SLC17A6"), 
                     expr.th = 3,
                     prop.th = c(0.5, 0.5, 0.2, 0.2), 
                     outlier.cl = NULL,
                     donor.cl = NULL) {
  select.cells <- intersect(names(cluster), select.cells)
  select.id    <- match(select.cells, colnames(norm.dat))
  norm.dat     <- norm.dat[, select.id]
  anno         <- droplevels(anno[select.id, ])
  anno$cluster <- cluster[select.cells]
  neuronal.cl  <- check_neun(anno, cluster, neun.thresh = neun.thresh, 
                             neun.colname = neun.colname, neun.val = neun.val)
  expr.cnt     <- apply(norm.dat[test.genes, ], 1, function(x) {
    tapply(x, cluster, function(x) sum(x > expr.th) / length(x))
  })
  cl.class <- list()
  cl.class[["inh"]] <- setdiff(
    names(neuronal.cl)[
      neuronal.cl & (expr.cnt[, test.genes[1]] > prop.th[1] | 
                       expr.cnt[, test.genes[2]] > prop.th[2])
    ],
    c(outlier.cl, donor.cl)
  )
  cl.class[["exc"]] <- setdiff(
    names(neuronal.cl)[
      neuronal.cl & (expr.cnt[, test.genes[3]] > prop.th[3] | 
                       expr.cnt[, test.genes[4]] > prop.th[4])
    ],
    c(outlier.cl, donor.cl, cl.class[["inh"]])
  )
  cl.class[["glia"]]    <- setdiff(names(neuronal.cl)[!neuronal.cl], 
                                   c(outlier.cl, donor.cl))
  cl.class[["donor"]]   <- setdiff(donor.cl, outlier.cl)
  cl.class[["outlier"]] <- setdiff(unique(cluster), 
                                   unlist(cl.class[c("inh", "exc", "glia", "donor")]))
  return(cl.class)
}


#---------------------------------------------------------------------------------------------------------------------------
# 5. check_donor

#' Find donor clusters
#' 
#' This function identifies clusters that are nearly exclusively expressed in one donor.
#'   Currently it is hard-coded to define clusters with >90% enrichment in a single donor
#'   or >50% more than expected based on the total number of cells per donor.  Later
#'   function iterations will include these options as parameters.
#'
#' @param anno anno dataframe which must include column names listed in `neun.colname`
#'   `meta1_area` and `meta2_donor` below. "cluster" is added from `cluster` parameter below.
#' @param cluster clusters labels for all cells along with sample_id as their names
#' @param norm.dat expression dataframe with columns as cells and rows as gene names 
#'   and cpm normalized
#' @param select.cells column nmaes of norm.dat
#' @param keep.cl clusters to definitely keep in analysis (e.g., to exclude from 
#'   consideration as a donor cluster) default is NULL
#' @param meta1_area The metadata column that contains information about the area 
#'   (e.g. cortical region and layer).  default is "roi" 
#' @param meta2_donor The metadata column that contains information about the donor. 
#'   default is "external_donor_name" 
#' @param plot default is TRUE
#' @param plot.path path of plot, default is "./output/"
#'
#' @return Gives cluster ids vs donor expression/region heatmap
#' @export
#'
#' @examples check_donor()
check_donor <- function(anno, cluster, norm.dat, select.cells = names(cluster),
                        keep.cl = NULL, meta1_area = "roi", meta2_donor = "external_donor_name",
                        plot = TRUE, plot.path = "output/") {
  select.cells <- intersect(names(cluster), select.cells)
  select.id    <- match(select.cells, colnames(norm.dat))
  norm.dat     <- norm.dat[, select.id]
  anno         <- droplevels(anno[select.id, ])
  layer.donor  <- as.matrix(table(anno[, meta1_area], anno[, meta2_donor]))
  layer.donor.prop <- sweep(layer.donor, 1, rowSums(layer.donor),"/")
  if (plot == TRUE & sd(layer.donor.prop) > 0) {
    hm.colors <- colorRampPalette(c("white", RColorBrewer::brewer.pal(9,"YlOrRd")))(100)
    pdf(
      file = paste0(plot.path, "/layer_donor_prop.pdf"),
      width = nrow(layer.donor.prop) / 6 + 1, height = ncol(layer.donor.prop) / 6 +
        1, onefile = FALSE
    )
    pheatmap::pheatmap(t(layer.donor.prop),
      cluster_rows = FALSE,
      cluster_cols = FALSE, color = hm.colors
    )
    dev.off()
  }
  layer.cl          <- as.matrix(table(anno[, meta1_area], cluster))
  layer.cl.prop     <- sweep(layer.cl, 2, colSums(layer.cl), "/")
  cl.donor.exp.cnt  <- round(t(t(layer.cl.prop) %*% layer.donor))
  cl.donor.exp.prop <- sweep(cl.donor.exp.cnt, 2, colSums(cl.donor.exp.cnt), "/")
  cl.donor.obs.cnt  <- t(as.matrix(table(cluster, anno[, meta2_donor])))
  cl.donor.obs.prop <- sweep(cl.donor.obs.cnt, 2, colSums(cl.donor.obs.cnt), "/")
  if (plot == TRUE & sd(cl.donor.obs.prop) > 0) {
    hm.colors <- colorRampPalette(c("white", RColorBrewer::brewer.pal(9,"YlOrRd")))(100)
    pdf(
      file = paste0(plot.path, "/cluster_donor_prop.pdf"),
      width = nrow(cl.donor.obs.prop) / 6 + 1, 
      height = ncol(cl.donor.obs.prop) / 6 + 1, 
      onefile = FALSE
    )
    pheatmap::pheatmap(t(cl.donor.obs.prop),
      cluster_rows = FALSE, cluster_cols = FALSE, color = hm.colors
    )
    dev.off()
  }
  cl.donor.prop.dev <- cl.donor.obs.prop - cl.donor.exp.prop
  cl.donor.prop.dev.bin <- ifelse(cl.donor.prop.dev > 0.5,
    1, ifelse(cl.donor.prop.dev < -0.5, -1, ifelse(cl.donor.obs.prop > 0.9, 1, 0))
  )
  donor.cl <- names(which(apply(cl.donor.prop.dev.bin, 2, function(x) any(x == 1))))
  donor.cl.final <- setdiff(donor.cl, keep.cl)
  return(donor.cl.final)
}

#---------------------------------------------------------------------------------------------------------------------------



#' Compare and plot two sets of cluster assignments
#' 
#' This is the subset of the `compare_annotate` function that does the plotting.
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param ref.cl A cluster factor object for the reference clusters
#'
#' @return g A ggplot2 dot plot object for the comparison.
#'
#' @export
compare_plot <- function(cl,ref.cl){
  
  library(ggplot2)
  
  common.cells <- intersect(names(cl),names(ref.cl))
  # compare predicted cluster member with the new clustering result 
  tb <- table(cl[common.cells], ref.cl[common.cells])
  
  # Plot the mapping
  tb.df <- as.data.frame(tb)
  tb.df <- tb.df[tb.df$Freq > 0,]
  
  select.cells <- names(cl)
  
  # Compute Jaccard statistics for each pair of clusters
  tb.df$jaccard <- 0
  for(i in 1:nrow(tb.df)){
    n_ol <- length(union(common.cells[cl[common.cells] == as.character(tb.df[i,1])],
                         common.cells[ref.cl[common.cells] == as.character(tb.df[i,2])]))
    
    tb.df$jaccard[i] <- tb.df$Freq[i] / n_ol
  }
  
  colnames(tb.df) <- c("cl","ref.cl","Freq","jaccard")
  
  g <- ggplot(tb.df, 
              aes(x = cl, 
                  y = ref.cl)) + 
    geom_point(aes(size = sqrt(Freq),
                   color = jaccard)) + 
    theme(axis.text.x = element_text(vjust = 0.1,
                                     hjust = 0.2, 
                                     angle = 90,
                                     size = 7),
          axis.text.y = element_text(size = 6)) + 
    scale_color_gradient(low = "yellow", high = "darkblue") + 
    scale_size(range=c(0,3))
  
  g
}


#---------------------------------------------------------------------------------------------------------------------------


#' Reorder factors of one annotation to match another
#' 
#' This function takes two identical annotations (e.g, two cluster calls for the 
#'   same cells) and updates the levels of the match cluster to correspond to 
#'   the levels of the reference cluster.  This is useful to run prior to
#'   `compare_plot` if clusters are not already ordered.
#'   
#' `clRef` and `clMatch` should be factors, but can be coerced from characters
#'   if needed.
#' 
#' @param clMatch A cluster factor object to compare to a reference
#' @param clRef A cluster factor object for the reference clusters
#'
#' @return The `clMatch` vector with levels rearranged to match `clRef` in a 
#'   reasonable way.
#'
#' @export
reorder_factor <- function(clMatch,clRef){
  out <- table(clMatch,clRef)
  out <- t(out)/colSums(out)
  fac <- apply(out,2,function(x) which.max(x)+0.01*sum(cumsum(x)))
  fac <- names(sort(fac))
  factor(clMatch,levels=fac)
}
