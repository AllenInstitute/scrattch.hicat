# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# find_vg()
#   compute_vg() vg.R
#   gene_loess_fit() vg.R
#   plot_vg() vg.R
#
# compute_vg()
#   rescale_samples() vg.R
#   gene_means() vg.R
#   gene_vars() vg.R
#   gene_dispersion() vg.R
#   gene_z() vg.R
#   gene_loess_z() vg.R
#
# rescale_samples()
#
# gene_means()
#   rescale_samples() vg.R
#
# gene_vars()
#   rescale_samples() vg.R
#   gene_means() vg.R
#
# gene_dispersion()
#   rescale_samples() vg.R
#   gene_means() vg.R
#   gene_vars() vg.R
#
# gene_z()
#   gene_dispersion() vg.R
#
# gene_loess_fit()
#   rescale_samples() vg.R
#   gene_means() vg.R
#   gene_dispersion() vg.R
#
# gene_loess_fit_z()
#
# plot_vg()
#



#' Identify high variability genes using Brennecke's method
#' 
#' This function is a wrapper around compute_vg_stats(), which computes variability statistics, 
#' and plot_vg(), which generates multiple plots for those results.
#' 
#' @param dat A matrix or dgCMatrix with samples as columns and genes as rows
#' @param rescaled logical - whether the data matrix should be normalized by columns). Assume input data have already been normalized
#' @param plot_file (optional) A file to use to output diagnostic plots. Default is NULL.
#' @param return_type What kind of objects to return. Default is "data", which will return only 
#' a data.frame with results for each gene. Other options are "plots", which will return a list 
#' containing plot objects, and "both" which will return a list with the statistics (gene_var_stats) 
#' and plots.
#' @param verbose logical - whether or not to display status messages
#' 
#' @return Either a data.frame or a list (see return_type parameter).
#' 
#' @export
#' 
find_vg <- function(dat,
                    rescaled = FALSE,
                    plot_file = NULL, 
                    return_type = "data",
                    verbose = FALSE) {
  
  return_type <- match.arg(return_type,
                           choices = c("data","plots","both"))
  
  gene_var_stats <- compute_vg_stats(dat,
                                     rescaled=rescaled,
                                     verbose = verbose)
  
  if(return_type == "data") {
    return(gene_var_stats)
  }
  
  if(!is.null(plot_file) | return_type %in% c("plots","both")){
    if(verbose) { cat("Generating plots\n") }
    loess_fit <- gene_loess_fit(dat = NULL,
                                dispersions = gene_var_stats$dispersion,
                                means = gene_var_stats$g.means,
                                rescale = FALSE)
    
    gene_var_plots <- plot_vg(gene_var_stats,
                              plots = "all",
                              loess_fit = loess_fit)
    
    if(!is.null(plot_file)) {
      pdf(plot_file)
      for(i in 1:length(gene_var_plots)) {
        print(gene_var_plots[[i]])
      }
      dev.off()
    }
    
    if(return_type == "plots") {
      return(gene_var_plots)
    } else if(return_type == "both") {
      return(list(gene_var_stats = gene_var_stats, 
                  plots = gene_var_plots))
    }
    
  }
  
}


#' Compute gene expression variance statistics
#' 
#' @param dat A matrix or dgCMatrix with samples as columns and genes as rows
#' 
#' @return a data.frame with gene symbols, mean, variance, dispersion, z score, pvalues, and loess statistics
#' @param verbose logical - whether or not to display status messages
#'
#' @export
#' 
compute_vg_stats <- function(dat,
                             rescaled = FALSE,
                             verbose = FALSE) {
  
  if(verbose) { cat("Rescaling data\n") }
  if(rescaled){
    scaled_data <- rescale_samples(dat)
  }
  else{
    scaled_data = dat
  }
  if(verbose) { cat("Computing scaled means\n") }
  means <- gene_means(scaled_data)
  
  if(verbose) { cat("Computing gene variance\n") }
  vars <- gene_vars(scaled_data,
                    means = means)
  
  if(verbose) { cat("Computing dispersions\n") }
  dispersion <- gene_dispersion(scaled_data,
                                means = means,
                                vars = vars)
  
  if(verbose) { cat("Computing z scores\n") }
  z <- gene_z(scaled_data,
              dispersions = dispersion)
  
  ###loess regression
  if(verbose) { cat("Performing Loess fit\n") }
  loess_fit <- gene_loess_fit(scaled_data,
                              dispersions = dispersion,
                              means = means)
  
  if(verbose) { cat("Computing Loess z\n") }
  loess_z <- gene_loess_fit_z(loess_fit,
                              dispersions = dispersion)
  
  gene_var_stats <- data.frame(gene = rownames(dat),
                               g.means = means, 
                               g.vars = vars, 
                               dispersion = dispersion,
                               z = z,
                               pval = 1 - pnorm(z),
                               padj = p.adjust(1 - pnorm(z), method = "fdr"),
                               loess.z = loess_z,
                               loess.pval = 1 - pnorm(loess_z),
                               loess.padj = p.adjust(1 - pnorm(loess_z), method = "fdr"))
  
  return(gene_var_stats)
}

#' Rescale a gene x sample matrix
#' 
#' This function scales each sample/column based on the median of the total of gene expression values for each sample.
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' 
#' @return a rescaled matrix of expression values
#' 
#' @export
#' 
rescale_samples <- function(dat) {
  
  if(is.matrix(dat)) {
    sample_totals <- colSums(dat)
  } else {
    sample_totals <- Matrix::colSums(dat)
    names(sample_totals) <- NULL
  }
  
  scaling_factor <- sample_totals / median(sample_totals)
  
  scaled_data <- dat / scaling_factor[col(dat)]
  
  scaled_data
}

#' Compute means for each gene in a gene x sample matrix
#' 
#' This is a convenient wrapper around rowMeans that is compatible with
#' multiple matrix types (including dgCMatrix)
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' @param rescale a logical indicating whether or not to rescale using rescale_samples(). Default is FALSE.
#' 
#' @return a vector of means for each gene/row
#' 
#' @export
#' 
gene_means <- function(dat,
                       rescale = FALSE) {
  
  if(rescale) {
    dat <- rescale_samples(dat)
  }
  
  if(is.matrix(dat)) {
    rowMeans(dat)
  } else {
    Matrix::rowMeans(dat)
  }
  
}

#' Compute variance for each gene in a gene x sample matrix
#'
#' This function uses mean(values ^ 2) - mean(values) ^ 2
#' 
#' Results will differ from the var() function.
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' @param means (optional) a vector of means for each gene/row. If NULL (default), will be computed using gene_means().
#' @param rescale a logical indicating whether or not to rescale using rescale_samples(). Default is FALSE.
#'
#' @return a numeric vector of variance values for each gene
#' 
#' @export
#' 
gene_vars <- function(dat,
                      means = NULL,
                      rescale = FALSE) {
  
  if(rescale) {
    dat <- rescale_samples(dat)
  }
  get_row_vars(dat, means=means)
}

#' Compute dispersion for each gene in a gene x sample matrix
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' @param means (optional) a vector of means for each gene/row. If NULL (default), will be computed using gene_means().
#' @param vars (optional) a vector of variances for each gene/row. If NULL (default), will be computed using gene_vars().
#' @param rescale a logical indicating whether or not to rescale using rescale_samples(). Default is FALSE.
#'
#' @return a numeric vector of dispersion values for each gene
#' 
#' @export
#' 
gene_dispersion <- function(dat,
                            means = NULL,
                            vars = NULL,
                            rescale = FALSE) {
  
  if(rescale) {
    dat <- rescale_samples(dat)
  }
  
  if(is.null(means)) {
    means <- gene_means(dat)
  }
  
  if(is.null(vars)) {
    vars <- gene_vars(dat, means = means)
  }
  
  log10(vars / means)
}

#' Compute dispersion z-scores for each gene in a gene x sample matrix
#' 
#' This version of z scoring will differ from base R's scale()
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' @param dispersions (optional) a vector of dispersions for each gene/row. If NULL (default), will be computed using gene_dispersion().
#' @param rescale a logical indicating whether or not to rescale using rescale_samples(). Default is FALSE.
#'
#' @return a numeric vector of z-scores for each gene
#'
#' @export
#' 
gene_z <- function(dat,
                   dispersions = NULL,
                   rescale = FALSE) {
  
  ###fit normal with 25% to 75%
  
  if(is.null(dispersions)) {
    dispersions <- gene_dispersion(dat, 
                                   rescale = rescale)
  }
  
  IQR <- quantile(dispersions, 
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
  
  (dispersions  - m) / delta 
}

#' Compute dispersion-mean Loess-fit for each gene in a gene x sample matrix
#' 
#' @param dat a matrix or dg[CT]Matrix of expression values with genes as rows and samples as columns.
#' @param dispersions (optional) a vector of dispersions for each gene/row. If NULL (default), will be computed using gene_dispersion().
#' @param means (optional) a vector of means for each gene/row. If NULL (default), will be computed using gene_means().
#' @param rescale a logical indicating whether or not to rescale using rescale_samples(). Default is FALSE.
#'
#' @return a loess fit object returned by the loess() function.
#'
#' @export
#' 
gene_loess_fit <- function(dat,
                           dispersions = NULL,
                           means = NULL,
                           rescale = FALSE) {
  
  if(rescale) {
    dat <- rescale_samples(dat)
  }
  
  if(is.null(means)) {
    means <- gene_means(dat)
  }
  
  if(is.null(dispersions)) {
    dispersions <- gene_dispersion(dat,
                                   means = means)
  }
  
  select <- !is.na(dispersions) & dispersions > 0
  
  selected_dispersions <- dispersions[select]
  selected_means <- means[select]
  
  fit <- limma::loessFit(selected_dispersions, log10(selected_means))
  fit$x = as.matrix(log10(selected_means))
  return(fit)
  
}

#' Compute dispersion z scores for the residuals of a loess fit
#' 
#' @param loess_fit an object output by the loess() function, with class "loess".
#' @param dispersions gene dispersion values as generated by the gene_dispersion() function.
#' 
#' @return a numeric vector of z score values
#' 
#' @export
#' 
gene_loess_fit_z <- function(loess_fit,
                             dispersions) {
  
  select <- !is.na(dispersions) & dispersions > 0
  
  residual <- resid(loess_fit)
  base <- min(loess_fit$fitted)
  
  diff <- dispersions - base
  diff[select] <- residual
  
  z <- gene_z(dat = NULL,
              dispersions = diff)
  
  return(z)
}

#' Generate gene variance plots
#' 
#' @param gene_var_stats a data.frame of statistics generated by compute_vg_stats()
#' @param plots a character vector specifying one or more of the following options (Default is "all"):
#' \itemize{
#'   \item qq_z_plot
#'   \item qq_loess.z_plot
#'   \item z_density_plot
#'   \item loess.z_density_plot
#'   \item dispersion_fit_plot
#'   \item all
#' }
#' @param loess_fit a loess fit object as generated by gene_loess_fit(). This is required only for the dispersion_fit_plot.
#' 
#' @return a list object with multiple ggplot2 plot objects depending on the plots parameter
#' 
#' @export
#' 
plot_vg <- function(gene_var_stats,
                    plots = "all",
                    loess_fit = NULL) {
  
  plots <- match.arg(plots,
                     choices = c("all","qq_z_plot","qq_loess.z_plot",
                                 "z_density_plot", "loess.z_density_plot",
                                 "dispersion_fit_plot"),
                     several.ok = TRUE)
  
  if(!is.null(loess_fit)) {
    fit.df <- data.frame(x = loess_fit$x[,1], 
                         y = loess_fit$fitted)
    fit.df <- fit.df[order(fit.df$x), ]
    rownames(fit.df) <- NULL
  }
  
  out_list <- list()
  
  if("all" %in% plots | "qq_z_plot" %in% plots) {
    out_list$qq_z_plot <- ggplot2::ggplot(gene_var_stats,
                                          ggplot2::aes(sample = z)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::ggtitle("z-score Q-Q Plot")
  }
  row.names(gene_var_stats) = gene_var_stats$gene
  
  if("all" %in% plots | "qq_loess.z_plot" %in% plots) {
    out_list$qq_loess.z_plot <- ggplot2::ggplot(gene_var_stats,
                                                ggplot2::aes(sample = loess.z)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::ggtitle("Loess z-score Q-Q Plot")
  }
  
  if("all" %in% plots | "z_density_plot" %in% plots) {
    out_list$z_density_plot <- ggplot2::ggplot(gene_var_stats, 
                                               ggplot2::aes(z)) + 
      ggplot2::geom_density() +
      ggplot2::ggtitle("z Density Plot")
  }
  
  if("all" %in% plots | "loess.z_density_plot" %in% plots) {
    out_list$loess.z_density_plot <- ggplot2::ggplot(gene_var_stats, 
                                                     ggplot2::aes(loess.z)) + 
      ggplot2::geom_density() +
      ggplot2::ggtitle("Loess z Density Plot")
  }
  
  if("all" %in% plots | "dispersion_fit_plot" %in% plots) {
    if(is.null(loess_fit)) {
      warning("No loess_fit provided. Skipping dispersion_fit_plot.")
    } else {
      
      out_list$dispersion_fit_plot <- ggplot2::ggplot() + 
        ggplot2::geom_point(data = gene_var_stats, 
                            ggplot2::aes(x = log10(g.means), 
                                         y = dispersion)) + 
        ggplot2::geom_line(data = fit.df, 
                           ggplot2::aes(x = x,
                                        y = y),
                           color = "blue") +
        ggplot2::ggtitle("Dispersion Fit Plot")
    }
    
    
  }
  
  return(out_list)
}

