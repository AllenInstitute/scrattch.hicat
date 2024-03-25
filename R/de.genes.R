# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# de_param()
#
# vec_chisq_test()
#
#
#   vec_chisq_test() de.genes.R
#
# de_selected_pairs() Tests for multiple pairs of clusters. Was DE_genes_pairs().
#   get_cl_means() util.R
#   score_pair_limma() de.genes.R
#   score_pair_chisq() de.genes.R
#
#
# compute_pair_deScore() Compute deScores based on score_pair_X() results and de_param(). Was de_pairs().
# 
#

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
                          fit) {
  
  x <- as.character(pair[1])
  y <- as.character(pair[2])
  
  ctr <- paste(paste0("cl", x), "-", paste0("cl", y))
   
  contrasts.matrix <- limma::makeContrasts(contrasts = ctr, 
                                           levels = design)
  
  fit2 <- limma::contrasts.fit(fit = fit, 
                               contrasts = contrasts.matrix)
  
  # Using suppressWarnings here due to this message from genes with 0 detection in both pairs:
  # Zero sample variances detected, have been offset away from zero
  fit2 <- suppressWarnings(limma::eBayes(fit = fit2))
  
  pval <- fit2$p.value[, 1]
  padj <- p.adjust(pval)
  lfc <- coef(fit2)[, 1]
  # Note: Above depends on the data being log2 scaled already.
  # If we change this expectation, we may need a more generalized calculation.
  # fc <- cl.means[, x] / cl.means[, y]
  # lfc <- log2(fc)
  # lfc[is.na(lfc)] <- 0
  
  results <- data.frame(padj = padj,
                        pval = pval,
                        lfc = lfc,
                        meanA = cl.means[[x]],
                        meanB = cl.means[[y]],
                        q1 = cl.present[[x]],
                        q2 = cl.present[[y]])

  row.names(results) <- row.names(cl.means)
  
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
                          cl.size)
{
  
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
  # Note: Above depends on the data being log2 scaled already.
  # If we change this expectation, we may need a more generalized calculation.
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
  row.names(results) <- row.names(cl.means)
  
  return(results)
  
}


#' Perform pairwise differential detection tests using t.test for a single pair of clusters
#'
#' @param pair a numeric vector of length 2 specifying which clusters to compare
#' @param cl.means a data.frame of normalized mean gene expression values (genes x clusters)
#' @param cl.vars a data.frame of normalized variance gene expression values (genes x clusters)
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
de_pair_t.test <- function(pair,
                          cl.means,
                          cl.present,
                          cl.vars,
                          cl.size)
{  
  x <- as.character(pair[1])
  y <- as.character(pair[2])
  m1 = cl.means[[x]]
  m2 = cl.means[[y]]
  v1 = cl.vars[[x]]
  v2 = cl.vars[[y]]
  n1 = cl.size[[x]]
  n2 = cl.size[[y]]
  sd = sqrt( v1/n1 + v2/n2) 
  t.stats = (m1 - m2) / sd
  df = sd^4 / ((v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1))
  pval = pt(abs(t.stats), df, lower.tail=FALSE)  * 2
  padj <- p.adjust(pval)
  
  lfc <- m1 - m2
  # Note: Above depends on the data being log2 scaled already.
  # If we change this expectation, we may need a more generalized calculation.
  # fc <- cl.means[, x] / cl.means[, y]
  # lfc <- log2(fc)
  # lfc[is.na(lfc)] <- 0
  
  results <- data.frame(padj = padj,
                        pval = pval,
                        lfc = lfc,
                        meanA = m1,
                        meanB = m2,
                        q1 = cl.present[,x], 
                        q2 = cl.present[,y])
  row.names(results) <- row.names(cl.means)  
  return(results)
  
}


get_cl_sigma <- function(mat,cl, cl.means=NULL, cl.sqr.means = NULL)
  {
    cl.size = table(cl)
    cl.size = setNames(as.vector(cl.size),names(cl.size))
    if(is.null(cl.means)){
      cl.means = get_cl_means(mat,cl)
    }
    if(is.null(cl.sqr.means)){
      cl.sqr.means =  get_cl_sqr_means(mat,cl)
    }
    df = sum(cl.size) - length(cl.size)
    sigma = sqrt(colSums( t(cl.sqr.means - cl.means^2) * cl.size[colnames(cl.sqr.means)]) / df)
  }



simple_lmFit <- function(norm.dat, cl, cl.means = NULL, cl.sqr.means = NULL)
  {
    if(is.null(cl.means)){
      cl.means = get_cl_means(norm.dat, cl)
    }
    if(is.null(cl.sqr.means)){
      cl.sqr.means = get_cl_sqr_means(norm.dat, cl)
    }
    cl.size = table(cl)
    cl.means = cl.means[,names(cl.size)]
    cl.sqr.means = cl.sqr.means[,names(cl.size)]
    fit = list()
    fit$coefficients = cl.means
    fit$rank = ncol(cl.means)
    fit$df.residual = length(cl) - fit$rank       
    fit$sigma = get_cl_sigma(norm.dat, cl, cl.means = cl.means, cl.sqr.means= cl.sqr.means)
    fit$stdev.unscaled = 1/sqrt(cl.size)
    return(fit)
  }
 


simple_ebayes <- function(fit,proportion=0.01,stdev.coef.lim=c(0.1,4),trend=FALSE,robust=FALSE,winsor.tail.p=c(0.05,0.1))
#	Empirical Bayes statistics to select differentially expressed genes
#	Gordon Smyth
#	8 Sept 2002.  Last revised 1 May 2013.
#	Made a non-exported function 18 Feb 2018.
{
  require(limma)
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
  if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
  if(all(!is.finite(sigma))) stop("No finite residual standard deviations")
  if(trend) {
    covariate <- fit$Amean
    if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
  } else {
    covariate <- NULL
  }
  
                                        #	Moderated t-statistic
  out <- squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust, winsor.tail.p=winsor.tail.p)
  out$s2.prior <- out$var.prior
  out$s2.post <- out$var.post
  out$var.prior <- out$var.post <- NULL
  out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
  df.total <- df.residual + out$df.prior
  df.pooled <- sum(df.residual,na.rm=TRUE)
  df.total <- pmin(df.total,df.pooled)
  out$df.total <- df.total
  out$p.value <- 2*pt(-abs(out$t),df=df.total)
  return(out)
}


de_pair_fast_limma <- function(pair,
                               fit,
                               cl.means,
                               cl.present)
{
  x <- as.character(pair[1])
  y <- as.character(pair[2])
  fit2 = fit
  coef = fit$coefficients
  fit2$coefficients = coef[[x]] - coef[[y]]
  stdev.unscaled = fit$stdev.unscaled
  fit2$stdev.unscaled = sqrt(sum(stdev.unscaled[c(x,y)]^2))
    
  fit2 <- simple_ebayes(fit = fit2)
  m1 = cl.means[[x]]
  m2 = cl.means[[y]]  


  lfc <- m1 - m2
  pval <- fit2$p.value
  padj <- p.adjust(pval) 
  
   # Note: Above depends on the data being log2 scaled already.
   # If we change this expectation, we may need a more generalized calculation.
   # fc <- cl.means[, x] / cl.means[, y]
   # lfc <- log2(fc)
   # lfc[is.na(lfc)] <- 0

  results <- data.frame(padj = padj,
                         pval = pval,
                         lfc = lfc,
                         meanA = m1,
                         meanB = m2,
                         q1 = cl.present[[x]],
                         q2 = cl.present[[y]])

  row.names(results) <- row.names(cl.means)

  return(results)

 }


 #' Perform pairwise differential gene expression tests between main pairs of clusters in parallel
 #' 
 #' @param norm.dat a normalized data matrix for data.
 #' @param cl a cluster factor object.
 #' @param pairs A 2-column matrix of cluster pairs.
 #' @param method Either "limma" or "chisq".
 #' @param low.th The minimum expression value used to filter for expressed genes.
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
                               cl.size=NULL,
                               de.param = de_parm(),
                               method = "fast_limma", 
                               cl.means = NULL,
                               cl.present = NULL,
                               cl.sqr.means = NULL,
                               use.voom = FALSE, 
                               counts = NULL,
                               mc.cores = 1,
                               block.size = 10000,
                               out.dir = NULL,
                               summary.dir = NULL,
                               top.n=500,
                               overwrite=FALSE,                               
                               return.df = FALSE,
                               return.summary=FALSE) {
   
   library(arrow)
   method <- match.arg(method,
                       choices = c("fast_limma", "limma","chisq", "t.test"))
   require(parallel)                                      
   if(use.voom & is.null(counts)) {
     stop("The use.voom = TRUE parameter requires a raw count matrix via the counts parameter.")
   }

   # Sample filtering based on selected clusters
   if(is.null(cl.size)){
     cl.size <- table(cl)
     cl.size = setNames(as.integer(cl.size), names(cl.size))
   }
   pairs.fn=NULL
   if(length(pairs)==1){
     pairs.fn = pairs
     pairs = open_dataset(pairs.fn)
   }
   else{     
     pairs =as.data.frame(pairs)
     if(is.null(pairs$pair)){
       pairs$pair = row.names(pairs)
     }
     if(is.null(pairs$pair_id)){
       pairs$pair_id = 1:nrow(pairs)
     }
   }

   select.cl <- unique(c(pairs %>% pull(P1), pairs %>% pull(P2)))
   select.cl <- intersect(select.cl, names(cl.size)[cl.size >= de.param$min.cells])
   cl <- cl[cl %in% select.cl]
   if(is.factor(cl)){
     cl = droplevels(cl)
   }
   pairs = pairs %>% filter(P1 %in% select.cl & P2 %in% select.cl)
   if(is.null(pairs$pair_bin)){
     pairs$pair_bin = ceiling(pairs$pair_id/block.size)
   }
   cl.size = cl.size[select.cl]   
   nbin = max(pairs$pair_bin)
   
  # Gene filtering based on low.th and min.cells thresholds
  # This was removed recently by Zizhen, as this can be computationally expensive
  # and can cause some inconsistent results based on which pairs are selected.
  # if(length(low.th) == 1) {
  #   genes_above_low.th <- Matrix::rowSums(norm.dat >= low.th)
  # } else {
  #   genes_above_low.th <- Matrix::rowSums(norm.dat >= low.th[row.names(norm.dat)])
  # }
  # 
  # genes_above_min.cells <- genes_above_low.th >= min.cells
  # 
  # select.genes <- row.names(norm.dat)[genes_above_min.cells]
  # 
  # norm.dat <- as.matrix(norm.dat[genes_above_min.cells, ])
  
  # Mean computation
   if(is.null(cl.means)) {
     cl.means <- as.data.frame(get_cl_means(norm.dat, cl))
   } else {
     cl.means <- as.data.frame(cl.means)     
   }
   
   # Compute fraction of cells in each cluster with expression >= low.th
   if(is.null(cl.present)){
     cl.present <- as.data.frame(get_cl_present(norm.dat, cl, de.param$low.th))
   } else{
     cl.present <- as.data.frame(cl.present)   
   }
   
   if(is.null(cl.sqr.means)){
     cl.sqr.means <- as.data.frame(get_cl_sqr_means(norm.dat, cl))
   } else{
     cl.sqr.means <- as.data.frame(cl.sqr.means)
   }
   
   if(method == "limma"){
     require("limma")    
     norm.dat <- as.matrix(norm.dat[, names(cl)])
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
   else if (method == "fast_limma"){
     fit = simple_lmFit(norm.dat, cl=cl, cl.means= cl.means, cl.sqr.means= cl.sqr.means)
   }
   else if (method == "t.test"){
     cl.vars <- as.data.frame(get_cl_vars(norm.dat, cl, cl.means = cl.means))
   }
 
   require(doMC)
   require(foreach)
   mc.cores = min(mc.cores, nbin)
   registerDoMC(cores=mc.cores)

   de_combine <- function(result.1, result.2)
     {
      library(data.table)
      de.genes = c(result.1$de.genes, result.2$de.genes)
      if(!is.null(result.1$de.summary)){
        de.summary = rbindlist(result.1$de.summary, result.2$de.summary)
        return(list(de.genes=de.genes, de.summary=de.summary))
      }
      else{
        return(list(de.genes=de.genes))
      }         
    }

   if(!is.null(out.dir)){
     if(!dir.exists(out.dir)){
       dir.create(out.dir)
     }
   }

   if(!is.null(summary.dir)){
     if(!dir.exists(summary.dir)){
       dir.create(summary.dir)
     }
   }
   mcoptions <- list(preschedule = FALSE)   
   de_list = foreach::foreach(bin=1:nbin,.combine="de_combine",.errorhandling="pass",.options.multicore=mcoptions) %dopar% {
     library(dplyr)
     library(arrow)     
     library(data.table)     
     if(!is.null(out.dir)){
       tmp.dir = file.path(out.dir,bin)     
       if(!dir.exists(tmp.dir)){
         dir.create(tmp.dir)
       }
       fn = file.path(tmp.dir,"data.parquet")     
       if(!overwrite){
         if(file.exists(fn)){
           return(NULL)
         }
       }
     }
     tmp.pair = pairs %>% filter(pair_bin==bin) %>% collect()
     x = tmp.pair %>% pull(pair_id)
     de.genes=sapply(x, function(i){
       pair = tmp.pair %>% filter(pair_id==i) 
       pair = unlist(pair[,1:2])
       if(method == "limma") {
         require("limma")
         df= de_pair_limma(pair = pair,
           cl.present = cl.present,
           cl.means = cl.means,
           design = design,
           fit = fit)
       }
       else if(method == "fast_limma") {
         df= de_pair_fast_limma(pair = pair,
           cl.present = cl.present,
           cl.means = cl.means,
           fit = fit)
       }
       else if(method =="t.test"){
         df = de_pair_t.test(pair = pair,
           cl.present = cl.present,
           cl.means = cl.means,
           cl.vars = cl.vars,
           cl.size = cl.size)                    
       }
       else if (method == "chisq"){
         df = de_pair_chisq(pair = pair,
           cl.present = cl.present,
           cl.means = cl.means,
           cl.size = cl.size)
       }      
       if(!is.null(de.param$min.cells)) {
         cl.size1 <- cl.size[as.character(pair[1])]
         cl.size2 <- cl.size[as.character(pair[2])]
       } else {
         cl.size1 <- NULL
         cl.size2 <- NULL
       }
       stats= de_stats_pair(df, 
         de.param = de.param, 
         cl.size1, 
         cl.size2,
         return.df = return.df)
     },simplify=F)
     pair = tmp.pair %>% pull(pair)
     names(de.genes) = pair
     de.summary = NULL
     if(return.summary){
       cat(bin,"compute summary\n")
       de.summary = de_pair_summary(de.genes, out.dir= summary.dir,block.size = block.size,blockstart= bin - 1, return.df = is.null(summary.dir))
     }     
     if(!is.null(out.dir)){
       cat(bin, "export de\n")
       result=export_de_genes(de.genes, cl.means, out.dir=out.dir, block.size = block.size,blockstart=bin - 1, mc.cores=1, top.n = top.n)
       de.genes= NULL
     }
     list(de.genes=de.genes, de.summary=de.summary)     
   }
   return(de_list)
 }

# Add docs and implement within functions

#' Title
#'
#' @param cn 
#' @param direction 
#' @param include.self 
#'
#' @return
#' @export
#'
#' @examples
create_pairs <- function(cn1, cn2=cn1,direction="nondirectional", include.self = FALSE)
  {
    cn1=as.character(cn1)
    cn2=as.character(cn2)
    cl.n1 = length(cn1)
    cl.n2 = length(cn2)	
    pairs = cbind(rep(cn1, rep(cl.n2,cl.n1)), rep(cn2, cl.n1))
    if(!identical(cn1,cn2)){
      down.pairs = cbind(rep(cn2, rep(cl.n1,cl.n2)), rep(cn1, cl.n2))
      pairs = rbind(pairs, down.pairs)
    }
    if(direction=="nondirectional"){
      pairs = pairs[pairs[,1]<=pairs[,2],,drop=F]
    }
    if(!include.self){
      pairs = pairs[pairs[,1]!=pairs[,2],,drop=F]
    }
    colnames(pairs)=c("P1","P2")
    row.names(pairs) = paste0(pairs[,1],"_",pairs[,2])
    return(pairs)
  }


null_de <- function()
  {
    tmp=sapply(c("score","up.score","down.score","num","up.num","down.num"), function(x)0)
    tmp = c(tmp, list(up.genes=NULL, down.genes=NULL))
  }


get_de_truncate_score_sum <- function(gene.score, th=20)
  {
    tmp = gene.score
    tmp[tmp > 20] = 20
    return(sum(tmp))
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
                          cl.size2 = NULL,
                          select.genes = NULL,
                          return.df = FALSE) {  
  df <- df[order(df$pval, -abs(df$lfc)), ]
  
  select <- with(df, which(padj < de.param$padj.th & abs(lfc) > de.param$lfc.th))
  select <- row.names(df)[select]
  
  if(!is.null(select.genes)){
    select <- select[select %in% select.genes]
  }
  
  if(is.null(select) | length(select) == 0){
    return(null_de())
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
    return(null_de())
  } else {

    up.genes = setNames(-log10(df[up,"padj"]), up)
    down.genes = setNames(-log10(df[down,"padj"]), down)
    
    tmp = up.genes
    tmp[tmp > 20] = 20
    up.score <- sum(tmp)
    tmp = down.genes
    tmp[tmp > 20] = 20
    down.score <- sum(tmp)    
   
    result=list(
      up.genes=up.genes,
      down.genes=down.genes,
      up.score = up.score,
      down.score = down.score,
      score = up.score + down.score,
      up.num = length(up.genes),
      down.num = length(down.genes),
      num = length(up.genes) + length(down.genes)
      )
    

    if(return.df){
      result$de.df = df[select,]
    }
    return(result)
  } 
}


#' Compute differential expression summary statistics for all pairs of clusters based on de_param()
#'
#' @param norm.dat a normalized data matrix for data.
#' @param cl a cluster factor object.
#' @param de.param A list of differential gene expression parameters from \code{de_param()}
#' @param method If de.df is NULL, use "limma" or "chisq" to compute differentially expressed genes.
#' @param de.df Optional. Pre-computed results from \code{de_all_pairs()} or \code{de_selected_pairs}. Default = NULL.
#' @param ... Additional parameters passed to \code{de_selected_pairs()}
#'
#' @return a character vector of all differentially expressed genes. 
#' @export
#'
de_all_pairs <- function(norm.dat, 
                         cl,
                         de.param = de_param(), 
                         method = "fast_limma", 
                         mc.cores=1,
                         pairs.fn = "pairs.parquet",
                         ...) {

  cn <- as.character(sort(unique(cl)))
  pairs= create_pairs(cn)
  if(nrow(pairs)>1000000){
    pairs = as.data.frame(pairs)
    pairs$pair = row.names(pairs)
    pairs$pair_id = 1:nrow(pairs)
    library(arrow)
    write_parquet(pairs, sink=pairs.fn)
    pairs = pairs.fn
  }
  de.result=de_selected_pairs(norm.dat,
    cl = cl,
    pairs = pairs,
    de.param = de.param,
    method = method,
    mc.cores=mc.cores,
    ...)  
  return(de.result$de.genes)
}


#' Generate a matrix of pairwise DE results
#'
#' @param de.genes Output from \code{de_score} or \code{de_stats_all_pairs}.
#' @param directed Logical indicating whether to select results based on all DE genes (Default, FALSE) or up/down regulated genes (see Details).
#' @param field The result to retrieve from de.results. Either "score" or "num". Default is "num".
#'
#' @details When directed = TRUE and field = "num", the minimum value from up or down-regulated genes is returned for each pair. When field = "score", the 
#' minimum deScore is returned.
#'
#' @return a matrix with clusters as rows and columns, and pairwise DE results as values.
#' @export
#'
get_de_matrix <- function(de.genes, 
                          directed = FALSE, 
                          field = "num") {
  
  field <- match.arg(field,
                     choices = c("num","score"))
    
  pairs <- get_pairs(names(de.genes))
  
  if(directed){
    # If directed, take the minimum of up or down-regulated
    f <- paste("up", field, sep = ".")
    up.de.num <- unlist(sapply(de.genes, function(x) { x[[f]] } ))
    
    f <- paste("down", field,sep=".")
    down.de.num <- unlist(sapply(de.genes, function(x) { x[[f]] } ))

    pairs = get_pairs(names(down.de.num))
    names(down.de.num) = paste(pairs[,2],pairs[,1], sep="_")
    de.num = c(up.de.num, down.de.num)
  } else {
    de.num <- unlist(sapply(de.genes, function(x) { x[[field]] }))
  }

  de.matrix <- convert_pair_matrix(de.num, 
                                   directed = directed)
  
  return(de.matrix)
}








#' Title
#'
#' @param de.genes 
#' @param dend 
#' @param cl.label 
#' @param directed 
#' @param file 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
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





#' Title
#'
#' @param de.genes 
#' @param top.n 
#' @param select.pair 
#' @param cl.label 
#'
#' @return
#' @export
#'
#' @examples
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




# New plotting functions from Zizhen - to be checked


#' Title
#'
#' @param pair.num 
#' @param file 
#' @param directed 
#' @param dend 
#' @param col 
#' @param cl.label 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_pair_matrix <- function(pair.num, file, directed=FALSE, dend=NULL, col=jet.colors(100), cl.label=NULL,...)
  {
    pair.matrix <- convert_pair_matrix(pair.num, directed = directed)
    if(!is.null(cl.label)){
      colnames(pair.matrix) = row.names(pair.matrix) = cl.label[row.names(pair.matrix)]
    }
    breaks = c(min(pair.num)-0.1, quantile(pair.num, seq(0.05,0.95,length.out=100)), max(pair.num)+0.1)
    pdf(file, ...)
    heatmap.3(pair.matrix, col = col, 
              trace = "none", Colv = dend, Rowv = dend, dendrogram = "row", 
              cexRow = 0.3, cexCol = 0.3)
    dev.off()     
  }




export_de_genes<- function(de.genes, cl.means, out.dir="de_parquet", block.size = 10000,blockstart=0, mc.cores=5, top.n = 1000, overwrite=FALSE)
  {
    library(data.table)
    library(arrow)
    if(!dir.exists(out.dir)){
      dir.create(out.dir)      
    }
    bins = split(names(de.genes), ceiling(1:length(de.genes)/block.size))
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    tmp=foreach::foreach(i=1:length(bins),.combine="c")%dopar% {
    #foreach::foreach(i=missing)%dopar% {
      library(dplyr)    
      library(arrow)
      library(data.table)      
      tmp.dir = file.path(out.dir,i+blockstart)
      if(!dir.exists(tmp.dir)){
        dir.create(tmp.dir)      
      }
      fn = file.path(tmp.dir,"data.parquet")
      if(!overwrite){
        if(file.exists(fn)){
          return(NULL)
        }
      }
      pairs=bins[[i]]      
      tmp = lapply(pairs, function(p){
        if(is.null(de.genes[[p]])|de.genes[[p]]$num==0){
          return(NULL)
        }
        up = de.genes[[p]]$up.genes
        down = de.genes[[p]]$down.genes
        up = head(up, top.n)
        down = head(down, top.n)
        pair = strsplit(p, "_")[[1]]
        p1 = pair[1]
        p2 = pair[2]
        gene = c(names(up),names(down))
        logPval = c(up, down)
        rank = c(seq_len(length(up)),seq_len(length(down)))        
        lfc =  abs(cl.means[gene, p1] - cl.means[gene, p2])
        sign = rep(c("up","down"),c(length(up),length(down)))
        df = data.frame(gene=gene, logPval=logPval,sign=factor(sign,c("up","down")), rank=rank)
        df$pair = p
        df$lfc = lfc       
        df       
      })
      df = data.table::rbindlist(tmp)
      write_parquet(df, sink=fn)
      return(NULL)
    }
  }


de_pair_summary <- function(de.genes, pairs=names(de.genes), mc.cores=1, blockstart = 0, block.size=10000, out.dir="de_summary",return.df = FALSE)
  {
    library(data.table)
    library(arrow)
    if(!is.null(out.dir) & !dir.exists(out.dir)){
      dir.create(out.dir)
    }
    cols = c("num","up.num","down.num","score", "up.score","down.score")
    bin = ceiling(1:length(pairs)/block.size)
    nbin = max(bin)
    require(doMC)
    require(foreach)
    mc.cores = min(mc.cores, nbin)
    registerDoMC(cores=mc.cores)

    de.df = foreach::foreach(i=1:nbin,.combine="c")%dopar% {
      bins = pairs[bin==i]
      tmp = sapply(cols, function(col){
        df=unlist(sapply(bins, function(p){
          de.genes[[p]][col]
        }))
      },simplify=F)
      df = do.call("data.frame",tmp)
      df$pair = bins
      if(!is.null(out.dir)){
        tmp.dir = file.path(out.dir,i + blockstart)
        if(!dir.exists(tmp.dir)){
          dir.create(tmp.dir)
        }
        fn = file.path(tmp.dir,"data.parquet")
        write_parquet(df, sink=fn)
      }
      if(return.df){
        list(df)
      }
      else{
        NULL
      }
    }
    if(return.df){
      de.df = rbindlist(de.df)
      return(de.df)
    }
    NULL
  }



