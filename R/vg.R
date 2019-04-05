# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# findVG()
#   compute_vg() vg.R
#   plot_vg() vg.R
#
# compute_vg()
#
# plot_vg()
#



#' Identify high variability genes using Brennecke's method
#' 
#' This function is a wrapper around compute_vg(), which computes variability statistics, 
#' and plot_vg(), which generates multiple plots for those results.
#' 
#' @param dat A matrix or dgCMatrix with samples as columns and genes as rows
#' @param plot_file (optional) A file to use to output diagnostic plots. Default is NULL.
#' @param return_type What kind of objects to return. Default is "data", which will return only 
#' a data.frame with results for each gene. Other options are "plots", which will return a list 
#' containing plot objects, and "both" which will return a list with the statistics (gene_var_stats) 
#' and plots.
#' 
#' @return Either a data.frame or a list (see return_type parameter).
#' 
#' @export
#' 
findVG <- function(dat, 
                   plot_file = NULL, 
                   return_type = "data") {
  
  gene_var_stats <- compute_vg(dat)
  
  if(return_type == "data") {
    return(gene_var_stats)
  }
  
  if(!is.null(plot_file) | return_type %in% c("plots","both")){
    
    gene_var_plots <- plot_vg(gene_var_stats,
                              plots = "all")
    
    if(!is.null(plot_file)) {
      pdf(plot_file)
      gene_var_plots
      dev.off()
    }
    
    if(return_type == "plots") {
      return(gene_var_plots)
    } else if(return_type == "both") {
      return(c(list(gene_var_stats), gene_var_plots))
    }
    
  }
  
}

#' Compute gene expression variance statistics
#' 
#' @param dat A matrix or dgCMatrix with samples as columns and genes as rows
#' 
#' @return a data.frame with gene symbols, mean, variance, dispersion, z score, pvalues, and loess statistics
#' 
#' @export
#' 
compute_vg <- function(dat) {
  
  scaled_data <- rescale_vg(dat)
  
  means <- gene_means(dat)
  
  vars <- gene_vars(dat,
                    means = gene_means)
  
  dispersion <- gene_dispersion(dat,
                                means = means,
                                vars = vars)
  
  z <- gene_z(dat,
              dispersion = dispersion)
  
  ###loess regression
  
  loess_z <- gene_loess_z(dat,
                          dispersion = dispersion,
                          means = means)
  
  gene_var_stats <- data.frame(gene = rownames(dat),
                               g.means = means, 
                               g.vars = vars, 
                               dispersion = dispersion,
                               z = z,
                               pval = 1 - pnorm(z),
                               padj = p.adjust(1 - pnorm(z), method = "fdr"),
                               loess.z = loess_z,
                               loess.pval = 1 - pnorm(loess_z),
                               loess.pad = p.adjust(1 - pnorm(loess_z), method = "fdr"))
  
  return(gene_var_stats)
}

rescale_vg <- function(dat) {
  
  if(!is.matrix(dat)) {
    sample_totals <- Matrix::colSums(dat)
  } else {
    sample_totals <- colSums(dat)
  }
  
  scaling_factor <- sample_totals / median(sample_totals)
  
  scaled_data <- dat / scaling_factor[col(dat)]
  
  scaled_data
}

gene_means <- function(dat) {
  
  if(is.matrix(dat)) {
    rowMeans(scaled_data)
  } else {
    Matrix::rowMeans(scaled_data)
  }
  
}

gene_vars <- function(dat,
                     means = NULL) {
  
  
  if(is.null(means)) {
    means <- gene_means(dat)
  }
  
  squared_dat <- dat^2
  
  if(is.matrix(dat)) {
    rowMeans(squared_dat) - means ^ 2
  } else {
    Matrix::rowMeans(squared_dat) - means ^ 2
  }
  
}

gene_dispersion <- function(dat,
                            means = NULL,
                            vars = NULL) {
  
  if(is.null(means)) {
    means <- gene_means(dat)
  }
  
  if(is.null(vars)) {
    vars <- gene_vars(dat, means = means)
  }
  
  log10(vars / means)
}

gene_z <- function(dat,
                   dispersion = NULL) {
  
  #####test samples####
  ###fit normal with 25% to 75%
  
  if(is.null(dispersion)) {
    dispersion <- gene_dispersion(dat)
  }
  
  IQR <- quantile(dispersion, 
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
  
  (dispersion  - m) / delta 
}

gene_loess_z <- function(dat,
                         dispersion = NULL,
                         means = NULL) {
  
  if(is.null(means)) {
    means <- gene_means(dat)
  }
  
  if(is.null(dispersion)) {
    dispersion <- gene_dispersion(dat,
                                  means = means)
  }
  
  select <- !is.na(dispersion) & dispersion > 0
  
  selected_dispersion <- dispersion[select]
  selected_means <- means[select]
  
  fit <- loess(selected_dispersion ~ log10(selected_means))
  
  residual <- resid(fit)
  base <- min(predict(fit))
  
  diff <- dispersion - base
  diff[select] <- residual
  
  IQR <- quantile(diff, 
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
  
  (diff - m) / delta
}

#' Generate gene variance plots
#' 
#' @param gene_var_stats a data.frame of statistics generated by compute_vg()
#' @param plots a character vector specifying one or more of the following options (Default is "all"):
#' \itemize{
#'   \item qq_z_plot
#'   \item qq_loess.z_plot
#'   \item z_density_plot
#'   \item loess.z_density_plot
#'   \item dispersion_fit_plot
#'   \item all
#' }
#' 
#' @return a list object with multiple ggplot2 plot objects depending on the plots parameter
#' 
#' @export
#' 
plot_vg <- function(gene_var_stats,
                    plots = "all") {
  
  plots <- match.arg(plots,
                     choices = c("all","qq_z_plot","qq_loess.z_plot",
                                 "z_density_plot", "loess.z_density_plot",
                                 "dispersion_fit_plot"),
                     several.ok = TRUE)
  
  out_list <- list()
  
  if("all" %in% plots | "qq_z_plot" %in% plots) {
    out_list$qq_z_plot <- ggplot2::ggplot(gene_var_stats,
                                          ggplot2::aes(sample = z)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line()
  }
  
  if("all" %in% plots | "qq_loess.z_plot" %in% plots) {
    out_list$qq_loess.z_plot <- ggplot2::ggplot(gene_var_stats,
                                                ggplot2::aes(sample = loess.z)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line()
  }
  
  if("all" %in% plots | "z_density_plot" %in% plots) {
    out_list$z_density_plot <- ggplot2::ggplot(gene_var_stats, 
                                               ggplot2::aes(z)) + 
      ggplot2::geom_density()
  }
  
  if("all" %in% plots | "loess.z_density_plot" %in% plots) {
    out_list$loess.z_density_plot <- ggplot2::ggplot(gene_var_stats, 
                                                     ggplot2::aes(loess.z)) + 
      ggplot2::geom_density()
  }
  
  if("all" %in% plots | "dispersion_fit_plot" %in% plots) {
    fit.df <- data.frame(x = fit$x[,1], y = fit$fitted)
    fit.df <- fit.df[order(fit.df$x), ]
    
    out_list$dispersion_fit_plot <- ggplot2::ggplot() + 
      ggplot2::geom_point(data = gene_var_stats, 
                          ggplot2::aes(x = log10(g.means), 
                                       y = dispersion)) + 
      ggplot2::geom_line(data = fit.df, 
                         aes(x = x,
                             y = y),
                         color = "blue")
  }
  
  return(out_list)
}
