# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# de_param()
#
# vec_chisq_test()
#
# score_pair_limma() Tests for one pair of clusters
#
# score_pair_chisq() Tests for one pair of clusters
#   vec_chisq_test() de.genes.R
#
# de_selected_pairs() Tests for multiple pairs of clusters. Was DE_genes_pairs().
#   get_cl_means() util.R
#   score_pair_limma() de.genes.R
#   score_pair_chisq() de.genes.R
#
# de_all_pairs() Tests for every pair of clusters. Was DE_genes_pw().
#   DE_genes_pairs() de.genes.R
#
# compute_pair_deScore() Compute deScores based on score_pair_X() results and de_param(). Was de_pairs().
# 
# de_score()
#   de_score_pairs() de.genes.R
#
# de_score_pairs()
#   DE_genes_pairs de.genes.R
#
# get_de_matrix()
#   get_pairs() de.genes.R
#   convert_pair_matrix() util.R
#
# plot_de_num()
#   get_de_matrix() de.genes.R
#   heatmap.3() heatmap.R
#
#  DE_genes_cat_by_cl
#    deScore.pairs() ??
#
# plot_de_lfc_num()

#' Set differential expression (DE) parameters for genes and clusters.
#'
#' This function provides a convenient way to manage settings for differential expression tests in scrattch.hicat.
#'
#' Calling \code{de.param()} without additional parameters provides reasonable defaults for high depth (e.g. SMART-seq) datasets.
#'
#' @param low.th  Lower boundary for normalized gene expression. Default = 1. See details.
#' @param padj.th Upper boundary threshold for adjusted p-values for differential expression tests. Default = 0.01. See details. 
#' @param lfc.th  Lower boundary threshold for log2(fold change) for differential expression tests. Default = 1 (i.e. 2-fold). See details
#' @param q1.th Lower boundary threshold for foreground detection in proportional comparisons. Default = 0.5. See details.
#' @param q2.th Upper boundary threshold for background detection in proportional comparisons. Default = NULL. See details.  
#' @param q.diff.th Threshold for scaled difference in proportions. Default = 0.7. See details.
#' @param de.score.th Lower boundary of total differential expression scores for cluster comparisons. Default = 150. See details.
#' @param min.cells The minimum number of cells allowed in each cluster. Default = 4. See details.
#' @param min.genes The minimum number of differentially expressed genes required to separate clusters. Default = 5. See details.
#' 
#' @details
#' \strong{Gene detection threshold:}
#' 
#' \code{low.th} sets a lower bound for normalized gene expression to determine whether or not a gene is considered to be detected.
#' This is used to filter genes that are too low in expression to be reliably detected.\cr
#' This parameter can be set globally by providing a single value, or per-gene by providing a named vector.
#' 
#' \strong{Differential expression test thresholds:}
#' 
#' \code{scrattch.hicat} utilizes \code{limma}'s \link[limma]{eBayes} or Chi-Square tests for differential gene expression. These parameters are used to 
#' determine which genes are considered differentially expressed:\cr
#' \code{padj.th} is the threshold for adjusted p-values. Adjusted p-values must be below this threshold to be considered significant.\cr 
#' \code{lfc.th} is the threshold for abs(log2(Fold Change)).
#' 
#' \strong{Cluster proportion thresholds:}
#' 
#' We use \code{q1.th}, \code{q2.th} and \code{q.diff.th} for additional tests based on the proportion of cells in each cluster that express each gene.
#' For every pair of clusters, we define \emph{q1} and \emph{q2} as the proportion of cells with expression greater than \code{low.th} 
#' (above) in the foregound and background cluster, respectively. 
#' We use \code{q1.th} to select genes in a high proportion of foreground clusters, and \code{q2.th} to select genes in a low proportion of background clusters. 
#' Finally, we use \code{q.diff.th} to test for the difference between the foreground and background proportions.
#' 
#' \code{q1.th}: The minimum proportion of cells in the foreground cluster with expression greater than \code{low.th}.\cr
#' \code{q2.th}: The maximum proportion of cells in the background cluster with expression greater than \code{low.th}.\cr
#' \code{q.diff.th}: The scaled proportional difference between \emph{q1} and \emph{q2}, defined as \eqn{abs(q1 - q2) / max(q1, q2)} .
#' 
#' \strong{Cluster-wise p-value threshold:}
#' 
#' After performing differential expression tests between a pair of clusters, we use \emph{de.score} as a way to determine if enough overall 
#' differential expression is observed to consider the two clusters distinct from each other.
#' 
#' We define \emph{de.score} for each gene as \eqn{min(-log10(p.adj), 20)}. This sets a cap on the contribution of each gene to the cluster-wise \emph{de.score} value at 20.\cr
#' The \emph{de.score} for a pair of clusters is the sum of the gene-wise \emph{de.score} values. 
#' 
#' Only genes passing the \code{padj.th} and \code{lfc.th} thresholds (above) contribute to the \emph{de.score}.
#' 
#' \code{de.score.th} is used as a minimum value for the cluster-wise \emph{de.score} in a pairwise comparison between clusters.
#' 
#' \strong{Cell and gene count thresholds:}
#' 
#' \code{min.cells} is the minimum size allowed for a cluster. If a cluster size is below \code{min.cells}, it will be merged with the nearest cluster.\cr
#' 
#' \code{min.genes} is the minimum number of differentially expressed genes (passing the \code{padj.th} and \code{lfc.th} thresholds, above) 
#' required to consider two clusters separate.
#' 
#' @return returns a list of parameters for reuse
#' @export
#' 
#' @examples 
#' 
#' # Recommended initial parameters for SMART-Seq (> 8,000 genes per sample):
#' 
#' sm_param <- de_param(low.th = 1,
#'                      padj.th = 0.01,
#'                      lfc.th = 1,
#'                      q1.th = 0.5,
#'                      q2.th = NULL,
#'                      q.diff.th = 0.7,
#'                      de.score.th = 150,
#'                      min.cells = 4,
#'                      min.genes = 5)
#' 
#' # Recommended initial parameters for 10x Cells (> 3,000 genes per sample):
#' 
#' tx_param <- de_param(low.th = 1,
#'                      padj.th = 0.01,
#'                      lfc.th = 1,
#'                      q1.th = 0.4, # Reduced due to dropout
#'                      q2.th = NULL,
#'                      q.diff.th = 0.7,
#'                      de.score.th = 150,
#'                      min.cells = 10, # Increased due to higher number of cells
#'                      min.genes = 5)
#' 
#' # Recommended initial parameters for 10x Nuclei (> 1,000 genes per sample):
#'
#' tx_param <- de_param(low.th = 1,
#'                      padj.th = 0.01,
#'                      lfc.th = 1,
#'                      q1.th = 0.3, # Reduced due to dropout
#'                      q2.th = NULL,
#'                      q.diff.th = 0.7,
#'                      de.score.th = 100, # Reduced due to decreased detection
#'                      min.cells = 10, # Increased due to higher number of cells
#'                      min.genes = 5) 
#'
de_param <- function(low.th = 1,
                     padj.th = 0.01, 
                     lfc.th = 1, 
                     q1.th = 0.5, 
                     q2.th = NULL,
                     q.diff.th = 0.7, 
                     de.score.th = 150, 
                     min.cells = 4, 
                     min.genes = 5) {
  
  if(padj.th > 1) {
    stop("padj.th must be a value <= 1.")
  }
  if(padj.th <= 0) {
    stop("padj.th must be a value > 0.")
  }
  if(lfc.th < 0) {
    stop("lfc.th must be a non-negative value.")
  }
  if(low.th < 0) {
    stop("low.th must be a non-negative value.")
  }
  if(q1.th < 0 | q1.th > 1) {
    stop("q1.th must be a non-negative value between 0 and 1.")
  }
  if(!is.null(q2.th)) {
    if(q2.th < 0 | q2.th > 1) {
      stop("q2.th must be a non-negative value between 0 and 1.")
    }
  }
  if(q.diff.th < 0) {
    stop("q.diff.th must be a non-negative value.")
  }
  if(de.score.th < 0) {
    stop("de.score.th must be a non-negative value.")
  }
  if(min.cells < 1) {
    stop("min.cells must be an integer with a value greater than 0.")
  }
  if(min.genes < 1) {
    stop("min.genes must be an integer with a value greater than 0.")
  }
  
  list(low.th = low.th, 
       padj.th = padj.th, 
       lfc.th = lfc.th, 
       q1.th = q1.th, 
       q2.th = q2.th, 
       q.diff.th = q.diff.th, 
       de.score.th = de.score.th, 
       min.cells = min.cells, 
       min.genes = min.genes)
  
  
}

#' Vectorized Chi-squared tests for differential gene detection
#' 
#' This function uses vectors of the number of samples in two sets that have detection of 
#' a set of genes and the total number of cells in each set to compute Chi-quared tests with 1 DOF for 
#' differential detection.
#' 
#' @param x an integer vector with the number of cells in group \emph{x} with detection of each gene.
#' @param x.total an integer value with the total number of cells in group \emph{x}.
#' @param y an integer vector with the number of cells in group \emph{y} with detection of each gene.
#' @param y.total an integer value with the total number of cells in group \emph{y}.
#' 
#' @return a data.frame with the following result for each gene:
#' \itemize{
#' \item{stats: The value of the chi-squared test statistic}
#' \item{pval: The p-value as reported by pchisq}
#' \item{logFC: The log2(fold change) in detection frequency between samples (x / y)}
#' \item{diff: The difference in proportions between the samples (x - y)}
#' }
#' 
#' @export
vec_chisq_test <- function(x, 
                           x.total, 
                           y, 
                           y.total) {
  
  total <- x.total + y.total
  present <- x + y
  absent <- total - x - y
  
  o <- cbind(x, 
             x.total - x, 
             y, 
             y.total - y)
  
  e <- cbind(present * x.total, 
             absent * x.total, 
             present * y.total, 
             absent * y.total)
  
  e <- e / as.vector(total)
  
  stat <- rowSums(pmax(0, abs(o - e) - 0.5) ^ 2 / e)
  
  results <- data.frame(stats = stat, 
                        pval = pchisq(stat, 1, lower.tail = FALSE), 
                        logFC = log2( (x * y.total) / (y * x.total)), 
                        diff = x / x.total - y / y.total)
  results
}

#' Perform pairwise DE tests using limma for a single pair of clusters
#' 
#' 
#' @param pair a numeric vector of length 2 specifying which clusters to compare
#' @param cl.present a data.frame of gene detection proportions (genes x clusters)
#' @param cl.means a data.frame of normalized mean gene expression values (genes x clusters)
#' @param design a limma design object
#' @param fit a limma fit object
#' @param genes the genes to use for pairwise comparisons
#' 
#' @return a data.frame with DE statistics:
#' \itemize{
#' \item{padj} P-values adjusted using the Holm (1979) method (\code{p.adjust()} default).
#' \item{pval} P-values reported by the \code{limma::eBayes()} function.
#' \item{lfc} Log fold change of mean expression values between the pair of clusters.
#' \item{meanA} Normalized mean expression value for the first cluster in the pair.
#' \item{meanB} Normalized mean expression value for the second cluster in the pair.
#' \item{q1} Proportion of cells expressing each gene for the first cluster in the pair.
#' \item{q2} Proportion of cells expressing each gene for the second cluster in the pair.
#' }
#' 
#' @export
#' 
de_pair_limma <- function(pair,
                             cl.present,
                             cl.means,
                             design,
                             fit,
                             genes) {
  
  x <- as.character(pair[1])
  y <- as.character(pair[2])
  
  ctr <- paste(paste0("cl", x), "-", paste0("cl", y))
  
  contrasts.matrix <- limma::makeContrasts(contrasts = ctr, 
                                           levels = design)
  
  fit2 <- limma::contrasts.fit(fit = fit, 
                               contrasts = contrasts.matrix)
  
  fit2 <- limma::eBayes(fit = fit2)
  
  pval <- fit2$p.value[, 1]
  padj <- p.adjust(pval)
  lfc <- coef(fit2)[, 1]
  # fc <- cl.means[, x] / cl.means[, y]
  # lfc <- log2(fc)
  # lfc[is.na(lfc)] <- 0
  
  results <- data.frame(padj = padj,
                        pval = pval,
                        lfc = lfc,
                        meanA = cl.means[,x], 
                        meanB = cl.means[,y],
                        q1 = cl.present[,x], 
                        q2 = cl.present[,y])
  
  row.names(results) <- genes
  
  return(results)
}


#' Perform pairwise differential detection tests using Chi-Squared for a single pair of clusters
#'
#' @param pair a numeric vector of length 2 specifying which clusters to compare
#' @param cl.present a data.frame of gene detection proportions (genes x clusters)
#' @param cl.means a data.frame of normalized mean gene expression values (genes x clusters)
#' @param cl.size a named numeric vector of cluster sizes
#' @param genes the genes to use for pairwise comparisons
#'
#' @return a data.frame with DE statistics:
#' \itemize{
#' \item{padj} P-values adjusted using the Holm (1979) method (\code{p.adjust()} default).
#' \item{pval} P-values reported by the \code{vec_chisq_test()} function.
#' \item{lfc} Log fold change of mean expression values between the pair of clusters.
#' \item{meanA} Mean expression value for the first cluster in the pair.
#' \item{meanB} Mean expression value for the second cluster in the pair.
#' \item{q1} Proportion of cells expressing each gene for the first cluster in the pair.
#' \item{q2} Proportion of cells expressing each gene for the second cluster in the pair.
#' }
#' 
#' @export
#'
de_pair_chisq <- function(pair,
                             cl.present,
                             cl.means,
                             cl.size,
                             genes) {
  
  x <- as.character(pair[1])
  y <- as.character(pair[2])
  
  chisq_results <- vec_chisq_test(cl.present[, x] * cl.size[[x]], 
                                  cl.size[x], 
                                  cl.present[, y] * cl.size[[y]], 
                                  cl.size[y])
  
  chisq_results$pval[is.na(chisq_results$pval)] <- 1
  
  pval <- chisq_results[, "pval"]
  padj <- p.adjust(pval)
  
  lfc <- cl.means[, x] - cl.means[, y]
  # lfc <- log2(fc)
  # lfc[is.na(lfc)] <- 0
  
  results <- data.frame(padj = padj,
                        pval = pval,
                        lfc = lfc,
                        meanA = cl.means[,x], 
                        meanB = cl.means[,y],
                        q1 = cl.present[,x], 
                        q2 = cl.present[,y])
  
  row.names(results) <- genes
  
  return(results)
  
}

#' Perform pairwise differential gene expression tests between main pairs of clusters in parallel
#' 
#' @param norm.dat a normalized data matrix for data.
#' @param cl a cluster factor object.
#' @param pairs A 2-column matrix of cluster pairs.
#' @param method Either "limma" or "chisq".
#' @param low.th The minimum expression value used to filter for expressed genes.
#' @param min.cells The minimum number of cells with expression above \code{low.th} to allow a gene to be used.
#' @param cl.present A matrix of proportions of cells in each cluster with gene detection. Can be generated with \code{get_cl_props()}. Default is NULL (will be generated).
#' @param use.voom Logical, whether or not to use \code{voom()} for \code{limma} calculations. Default is FALSE.
#' @param counts A matrix of raw count data for each cell. Required if \code{use.voom} is TRUE. Default is NULL.
#' @param mc.cores A number indicating how many processor cores to use for parallelization.
#' 
#' @return 
#' @export
de_selected_pairs <- function(norm.dat, 
                              cl,
                              pairs, 
                              method = "limma", 
                              low.th = 1, 
                              min.cells = 4, 
                              cl.present = NULL, 
                              use.voom = FALSE, 
                              counts = NULL,
                              mc.cores = 1) {
  
  method <- match.arg(method,
                      choices = c("limma", "chisq"))
  
  if(use.voom & is.null(counts)) {
    stop("The use.voom = TRUE parameter requires a raw count matrix via the counts parameter.")
  }
  
  # Sample filtering based on selected clusters
  select.cl <- unique(c(pairs[,1], pairs[,2]))
  
  cl <- cl[cl %in% select.cl]
  
  norm.dat <- norm.dat[, names(cl)]
  
  # Gene filtering based on low.th and min.cells thresholds
  if(length(low.th) == 1) {
    genes_above_low.th <- Matrix::rowSums(norm.dat >= low.th)
  } else {
    genes_above_low.th <- Matrix::rowSums(norm.dat >= low.th[row.names(norm.dat)])
  }
  
  genes_above_min.cells <- genes_above_low.th >= min.cells
  
  select.genes <- row.names(norm.dat)[genes_above_min.cells]
  
  norm.dat <- as.matrix(norm.dat[genes_above_min.cells, ])
  
  # Mean computation
  cl.means <- as.data.frame(get_cl_means(norm.dat, 
                                         cl))
  
  # Set thresholds per gene
  if(length(low.th) == 1) {
    low.th <- setNames(rep(low.th, nrow(norm.dat)), row.names(norm.dat))      
  }
  
  # Compute fraction of cells in each cluster with expression >= low.th
  if(is.null(cl.present)){
    cl.present <- as.data.frame(get_cl_means(norm.dat >= low.th[row.names(norm.dat)],
                                             cl))
  }
  
  cl.size <- table(cl)
  
  if(method=="limma"){
    cl <- setNames(as.factor(paste0("cl",cl)),names(cl))
    design <- model.matrix(~0 + cl)
    colnames(design) <- levels(as.factor(cl))
    
    if(use.voom & !is.null(counts)){
      v <- limma::voom(counts = as.matrix(counts[row.names(norm.dat), names(cl)]), 
                       design = design)
      
      fit <- limma::lmFit(object = v, 
                          design = design)		
    } else {
      fit <- limma::lmFit(object = norm.dat[, names(cl)], 
                          design = design)
    }
  }
  
  if(mc.cores == 1) {
    de_list <- list()
    
    for(i in 1:nrow(pairs)) {
      
      pair <- paste(pairs[i, 1], pairs[i, 2], sep = "_")
      if(method == "limma") {
        de_list[[pair]] <- de_pair_limma(pair = pairs[i,],
                                            cl.present = cl.present,
                                            cl.means = cl.means,
                                            design = design,
                                            fit = fit,
                                            genes = select.genes)
        
      } else if(method == "chisq") {
        de_list[[pair]] <- de_pair_chisq(pair = pairs[i,],
                                            cl.present = cl.present,
                                            cl.means = cl.means,
                                            cl.size = cl.size,
                                            genes = select.genes)
      }
      
    }
  } else {
    library(foreach)
    
    cluster <- parallel::makeCluster(mc.cores)
    doParallel::registerDoParallel(cluster)
    
    if(method == "limma") {
      de_list <- foreach::foreach(i = 1:nrow(pairs), 
                                  .combine='c') %dopar% 
        list(de_pair_limma(pair = pairs[i,],
                              cl.present = cl.present,
                              cl.means = cl.means,
                              design = design,
                              fit = fit,
                              genes = select.genes))
    } else if(method == "chisq") {
      de_list <- foreach::foreach(i = 1:nrow(pairs), 
                                  .combine='c') %dopar% 
        list(de_pair_chisq(pair = pairs[i,],
                              cl.present = cl.present,
                              cl.means = cl.means,
                              cl.size = cl.size,
                              genes = select.genes))
    }
    
    parallel::stopCluster(cluster)
    
    names(de_list) <- paste(pairs[,1],pairs[,2],sep="_")
  }
  
  return(de_list)
}


####Make sure dat and cl has the same dimension, and cells are in the same order

#' Perform all pairwise differential expression comparison between clusters
#' 
#' @param norm.dat
#' @param cl
#' @param ... Additional parameters passed to DE_genes_pairs()
#' 
#' @seealso \link{DE_genes_pairs}
#' 
#' @return a list containing DE results for every pair of clusters
#' 
#' @export
de_all_pairs <- function(norm.dat,
                         cl, 
                         ...) {
  
  if(sum(names(cl) %in% colnames(norm.dat)) != length(cl)) {
    stop("Missing data for some cells in cl.")
  }
  
  cn <- as.character(sort(unique(cl)))
  cl.n <- length(cn)
  pairs <- cbind(rep(cn, rep(cl.n,cl.n)), 
                 rep(cn, cl.n))
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  
  de_selected_pairs(norm.dat = norm.dat,
                    cl = cl,
                    pairs = pairs, 
                    ...)
  
}


#' Compute differential expression summary statistics based on a differential results data.frame and de_param().
#' 
#' @param df A data.frame of pairwise differential expression results (i.e. from \code{score_selected_pairs()}).
#' @param de.param A list of differential gene expression parameters from \code{de_param()}
#' @param cl.size1 Optional: The number of samples in the first/high cluster
#' @param cl.size2 Optional: The number of samples in the second/low cluster
#' 
#' @results A list of filtered differential expression results containing:
#' \itemize{
#' \item{score} The deScore value, equal to the sum of the -log10(p-values) of differentially expressed genes, with a cap of 20 per gene.
#' \item{up.score} The deScore value for up-regulated genes.
#' \item{down.score} The deScore value for down-regulated genes.
#' \item{num} The number of differentially expressed genes
#' \item{up.num} The number of up-regulated genes
#' \item{down.num} The number of down-regulated genes
#' \item{genes} Gene symbols for differentially expressed genes.
#' \item{up.genes} Gene symbols for up-regulated genes.
#' \item{down.genes} Gene symbols for down-regulated genes.
#' \item{de.df} The df used as input, filtered for differentially expressed genes.
#' }
#' 
de_stats_pair <- function(df,
                             de.param = de_param(), 
                             cl.size1 = NULL, 
                             cl.size2 = NULL) {
  
  df <- df[order(df$pval, -abs(df$lfc)), ]
  
  select <- with(df, which(padj < de.param$padj.th & abs(lfc) > de.param$lfc.th))
  select <- row.names(df)[select]
  
  if(is.null(select) | length(select) == 0){
    return(list())
  }
  
  up <- select[df[select, "lfc"] > 0]
  down <- select[df[select, "lfc"] < 0]
  df <- df[select, ]
  
  if(!is.null(de.param$q.diff.th) & is.null(df$q.diff)) {
    
    df$q.diff <- with(df, abs(q1 - q2) / pmax(q1, q2))
    df$q.diff[is.na(df$q.diff)] <- 0
    
  }
  
  if(!is.null(de.param$q1.th)) {
    up <- with(df[up, , drop = FALSE], up[q1 > de.param$q1.th])
    
    if(!is.null(cl.size1)){
      up <- with(df[up, , drop = FALSE], up[q1 * cl.size1 >= de.param$min.cells])
    }
    
    down <- with(df[down, , drop = FALSE], down[q2 > de.param$q1.th])
    
    if(!is.null(cl.size2)) {
      down <- with(df[down, , drop = FALSE], down[q2 * cl.size2 >= de.param$min.cells])
    }
  }
  
  if(!is.null(de.param$q2.th)) {
    up <- with(df[up, , drop = FALSE], up[q2 < de.param$q2.th])
    down <- with(df[down, , drop = FALSE], down[q1 < de.param$q2.th])
  }
  
  if(!is.null(de.param$q.diff.th)){
    up <- with(df[up, , drop = FALSE], up[abs(q.diff) > de.param$q.diff.th])
    down <- with(df[down, , drop = FALSE], down[abs(q.diff) > de.param$q.diff.th])
  }
  
  select <- c(up, down)
  
  if(length(select) == 0){
    return(list())
  } else {
    
    df$padj[df$padj < 1e-20] <- 1e-20
    up.score <- sum(-log10(df[up,"padj"]))
    down.score <- sum(-log10(df[down,"padj"]))
    
    if(length(up) == 0) { up.score <- 0 }
    if(length(down) == 0) { down.score <- 0 }
    
    list(score = up.score + down.score,
         up.score = up.score,
         down.score = down.score,
         num = length(select),
         up.num = length(up),
         down.num = length(down),
         genes = select,
         up.genes = up,
         down.genes = down, 
         de.df = df[df$padj < de.param$padj.th, ])
  }
  
}


de_stats_selected_pairs<- function(norm.dat, 
                           cl, 
                           pairs, 
                           de.df = NULL, 
                           de.param = de_param(), 
                           method = "limma", 
                           mc.cores = 1) {
  
  row.names(pairs) <- paste(pairs[, 1], pairs[, 2], sep = "_")
  
  select.cl <- unique(c(pairs[, 1],pairs[, 2]))
  cl <- cl[cl %in% select.cl]
  
  norm.dat <- as.matrix(norm.dat[, names(cl)])
  
  if(is.factor(cl)) { cl <- droplevels(cl) }
  
  cl.size <- table(cl)
  cl.n <- names(cl.size)
  
  cl.small <- cl.n[cl.size < de.param$min.cells]
  cl.big <- setdiff(cl.n,cl.small)
  
  select.pair <- pairs[, 1] %in% cl.big & pairs[, 2] %in% cl.big
  
  de.genes <- list()
  
  if(sum(select.pair) > 0) {
    
    cl <- cl[cl %in% c(pairs[select.pair, 1], pairs[select.pair, 2])]
    select.cells <- names(cl)
    low.th <- de.param$low.th
    
    if(length(low.th) == 1){
      low.th <- setNames(rep(low.th, nrow(norm.dat)),
                         row.names(norm.dat))      
    }
    
    if(is.null(de.df)){
      de.df <- de_selected_pairs(norm.dat, 
                              cl[select.cells], 
                              pairs[select.pair, , drop = FALSE], 
                              low.th = low.th,
                              min.cells = de.param$min.cells, 
                              method = method, 
                              mc.cores = mc.cores)
    }
    
    de.genes <- sapply(names(de.df), 
                       function(x) {
                         if(is.null(de.df[[x]])){
                           return(list())
                         }
                         
                         df = de.df[[x]]
                         if(!is.null(de.param$min.cells)) {
                           cl.size1 <- cl.size[as.character(pairs[x, 1])]
                           cl.size2 <- cl.size[as.character(pairs[x,2])]
                         } else {
                           cl.size1 <- NULL
                           cl.size2 <- NULL
                         }
                         score_pair_stats(df, 
                                 de.param = de.param, 
                                 cl.size1, 
                                 cl.size2)
                       },
                       simplify = FALSE)
  } else {
    de.df <- list()
  }
  
  for(i in which(!select.pair)) {
    pair <- paste(pairs[i, 1], pairs[i, 2], sep = "_")
    de.genes[[pair]] <- list()
    de.df[[pair]] <- list()
  }
  
  list(de.df = de.df, 
       de.genes = de.genes)
}


#' Title
#'
#' @param norm.dat log transformed normalized gene expression matrix
#' @param cl  Cluster membership. 
#' @param min.cells Minimal number of cells in a cluster
#' @param de.param  DE gene criteria set by de_param function. 
#' @param method    Default "limma".
#' @param de.genes  If DE genes already computed for some pairs of clusters, use them without recomputation. 
#'
#' @return a list with DE genes from all pairs of clusters. 
#' @export
#'
de_stats_all_pairs <- function(norm.dat, 
                     cl,  
                     de.param = de_param(), 
                     method = "limma", 
                     de.genes = NULL,
                     ...) {
  
  if(is.factor(cl)){
    cl <- droplevels(cl)
  }
  
  cn <- as.character(sort(unique(cl)))
  cl.n <- length(cn)	
  pairs <- cbind(rep(cn, rep(cl.n, cl.n)), rep(cn, cl.n))
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  row.names(pairs) <- paste(pairs[, 1], pairs[, 2], sep = "_")
  
  if(!is.null(de.genes)){
    pairs <- pairs[!row.names(pairs) %in% names(de.genes), , drop = FALSE]
  }
  
  de.result <- de_score_pairs(norm.dat, 
                              cl = cl, 
                              pairs = pairs, 
                              de.param = de.param, 
                              method = method,
                              ...)
  
  de.genes <- c(de.genes, de.result$de.genes)
  
  return(de.genes)
}


get_de_matrix <- function(de.genes, 
                          directed = FALSE, 
                          field = "num") {
  
  pairs <- get_pairs(names(de.genes))
  
  if(directed){
    f <- paste("up", field, sep = ".")
    up.de.num <- sapply(de.genes, function(x) { x[[f]] } )
    
    f <- paste("down", field,sep=".")
    down.de.num <- sapply(de.genes, function(x) { x[[f]] } )
    
    # head(sort(pmin(up.de.num, down.de.num)))
    
    names(down.de.num) <- with(pairs[names(down.de.num), ], paste0(P2, "_", P1))
    de.num <- c(up.de.num, down.de.num)
    
  } else {
    de.num <- sapply(de.genes, function(x) { x[[field]] })
  }
  
  de.matrix <- convert_pair_matrix(de.num, 
                                   directed = directed)
  
  return(de.matrix)    
}

plot_de_num <- function(de.genes, 
                        dend, 
                        cl.label = NULL, 
                        directed = FALSE, 
                        file = "log10.de.num.pdf", 
                        ...) {
  
  label <- as.hclust(dend)$label
  de.num.matrix <- get_de_matrix(de.genes, 
                                 directed = directed)
  de.num.matrix <- de.num.matrix[label, label]
  
  breaks <- c(-1, seq(0.2, 4, length.out = 100))
  
  if(!is.null(cl.label)) {
    colnames(de.num.matrix) <- cl.label[row.names(de.num.matrix)]
    row.names(de.num.matrix) <- cl.label[row.names(de.num.matrix)]
  }
  
  tmp.dat <- log10(de.num.matrix + 1)
  
  pdf(file, ...)
  heatmap.3(tmp.dat, 
            col = jet.colors(100), 
            breaks = breaks,
            trace = "none",
            Colv = dend, 
            Rowv = dend,
            dendrogram = "row",
            cexRow = 0.3,
            cexCol = 0.3)
  dev.off()
}


DE_genes_cat_by_cl <- function(norm.dat, 
                               cl, 
                               binary.cat, 
                               ...) {
  
  cl <- droplevels(as.factor(cl))
  cl.cat <- setNames(paste0(cl, binary.cat[names(cl)]), 
                     names(cl))
  
  tmp <- levels(cl)
  cl.cat.pairs <- data.frame(cat1 = paste0(tmp, binary.cat[1]), 
                             cat2 = paste0(tmp, binary.cat[2]), 
                             stringsAsFactors = FALSE)
  
  cl.cat.de.genes <- deScore.pairs(norm.dat[,names(cl.cat)], 
                                   cl = cl.cat, 
                                   pairs = cl.cat.pairs, 
                                   ...)
  
  cat1.de.num <- sapply(cl.cat.de.genes,
                        function(x) {
                          if(length(x) == 0) { return(0) }
                          length(x$up.genes)
                        })
  
  cat2.de.num <- sapply(cl.cat.de.genes,
                        function(x) {
                          if(length(x) == 0) { return(0) }
                          length(x$down.genes)
                        })
  
  cat1.de.genes <- sapply(cl.cat.de.genes,
                          function(x) {
                            if(is.null(x)) { return("") }
                            
                            return(paste(head(x$up.genes, 8),
                                         collapse = " "))
                          })
  
  cat2.de.genes <- sapply(cl.cat.de.genes,
                          function(x) {
                            if(is.null(x)) { return("") }
                            
                            return(paste(head(x$down.genes, 8),
                                         collapse = " "))
                          })
  
  cl.cat.de.df <- data.frame(cat1.de.num, 
                             cat2.de.num, 
                             cat1.de.genes, 
                             cat2.de.genes)
  
  return(list(cl.cat.de.df = cl.cat.de.df,
              cl.cat.de.genes = cl.cat.de.genes))
}



plot_de_lfc_num <- function(de.genes, 
                            top.n = 100, 
                            select.pair = NULL, 
                            cl.label = NULL) {
  
  tmp <- sapply(de.genes, length)
  de.genes <- de.genes[tmp!=0]
  
  de.score <- sapply(de.genes, 
                     function(x) {
                       x$score
                     })
  
  de.up.score <- sapply(de.genes, 
                        function(x) {
                          x$up.score
                        })
  
  de.down.score <- sapply(de.genes, 
                          function(x) {
                            x$down.score
                          })
  
  de.num <- sapply(de.genes, 
                   function(x) {
                     x$num
                   })
  
  de.lfc <- sapply(de.genes, 
                   function(x) {
                     top.genes <- head(x$genes[order(x$de.df[x$genes, "pval"])],
                                       top.n)
                     
                     mean(abs(x$de.df[top.genes, "lfc"]))
                   })
  
  de.q.diff <- sapply(de.genes, 
                      function(x) {
                        top.genes <- head(x$genes[order(x$de.df[x$genes, "pval"])],
                                          top.n)
                        q.diff <- with(x$de.df[top.genes, ], abs(q1 - q2) / pmax(q1, q2))
                        mean(q.diff)
                      })
  
  de.summary <- data.frame(de.num, 
                           de.lfc, 
                           de.q.diff, 
                           de.score, 
                           de.up.score, 
                           de.down.score)
  
  row.names(de.summary) <- names(de.genes)
  
  tmp <- do.call("rbind", strsplit(row.names(de.summary), "_"))
  de.summary$cl1 <- tmp[, 1]
  de.summary$cl2 <- tmp[, 2]
  
  g <- ggplot2::ggplot(de.summary, 
                       ggplot2::aes(de.num,
                                    de.lfc,
                                    color = de.q.diff)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_color_gradient2(midpoint = 0.85) + 
    ggplot2::scale_x_log10()
  
  if(!is.null(select.pair)){
    select.df <- de.summary[select.pair, ]
    
    if(!is.null(cl.label)){
      select.df$pair.label <- with(select.df, 
                                   paste(cl.label[as.character(cl1)], 
                                         cl.label[as.character(cl2)],
                                         sep = ":"))
    }
    
    g <- g + 
      geom_text(data = select.df, 
                aes(de.num - 0.02, 
                    de.lfc, 
                    label = pair.label),
                size = 2,
                color = "black") +
      geom_point(data = select.df, 
                 aes(de.num, 
                     de.lfc),
                 color = "red",
                 pch = 1)
  }
  
  g <- g + 
    xlab("Number of DE genes") + 
    ylab("Mean log2(FC) of top 100 DE.genes")
  
  return(list(g = g, 
              de.summary = de.summary))
}


