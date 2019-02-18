#' Identify high variability genes using Brennecke's method
#' 
#' @param dat A matrix or dgCMatrix with samples as columns and genes as rows
#' @param plot_file (optional) A file to use to output diagnostic plots. Default is NULL.
#' @param return_type What kind of objects to return. Default is "data", which will return only 
#' a data.frame with results for each gene. Other options are "plots", which will return a list 
#' containing plot objects, and "both" which will return a list with the statistics (gene_var_stats) 
#' and plots.
#' 
#' @return Either a data.frame or a list (see return_type parameter).
findVG <- function(dat, 
                   plot_file = NULL, 
                   return_type = "data") {
  
  require(Matrix)
  require(ggplot2)
  
  if(!is.matrix(dat)) {
    sample_totals <- Matrix::colSums(dat)
  } else {
    sample_totals <- colSums(dat)
  }
  
  scaling_factor <- sample_totals / median(sample_totals)
  
  scaled_data <- dat / scaling_factor[col(dat)]
  
  if(is.matrix(dat)) {
    gene_means <- rowMeans(scaled_data)
    gene_vars <-  rowMeans(scaled_data^2) - gene_means^2
  } else {
    gene_means <- Matrix::rowMeans(scaled_data)
    gene_vars <- Matrix::rowMeans(scaled_data^2) - gene_means^2
  }
  
  dispersion <- log10(gene_vars/gene_means)
  
  #####test samples####
  ###fit normal with 25% to 75%
  
  IQR <- quantile(dispersion, 
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
  z <- (dispersion  - m) / delta        
  
  gene_var_stats <- data.frame(gene = rownames(dat),
                               g.means = gene_means, 
                               g.vars = gene_vars, 
                               dispersion = dispersion,
                               z = z,
                               pval = 1 - pnorm(z),
                               padj = p.adjust(1 - pnorm(z), method = "fdr"))
  
  ###loess regression
  
  select<- !is.na(gene_var_stats$dispersion) & gene_var_stats$dispersion > 0
  
  fit <- with(gene_var_stats, loess(dispersion ~ log10(g.means), subset = select))
  
  residual <- resid(fit)
  base <- min(predict(fit))
  diff <- gene_var_stats$dispersion - base
  diff[select] <- residual
  
  IQR <- quantile(diff, 
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
  
  gene_var_stats$loess.z <- (diff - m) / delta
  gene_var_stats$loess.pval = 1 - pnorm(gene_var_stats$loess.z)
  gene_var_stats$loess.padj = p.adjust(gene_var_stats$loess.pval, method="fdr")
  
  if(!is.null(plot_file) | return_type %in% c("plots","both")){
    
    qq_z_plot <- ggplot(gene_var_stats,
                        aes(sample = z)) +
      stat_qq() +
      stat_qq_line()
    
    qq_loess.z_plot <- ggplot(gene_var_stats,
                              aes(sample = loess.z)) +
      stat_qq() +
      stat_qq_line()
    
    z_density_plot <- ggplot(gene_var_stats, aes(z)) + 
      geom_density()
    
    loess.z_density_plot <- ggplot(gene_var_stats, aes(loess.z)) + 
      geom_density()
    
    fit.df <- data.frame(x = fit$x[,1], y = fit$fitted)
    fit.df <- fit.df[order(fit.df$x), ]
    
    dispersion_fit_plot <- ggplot() + 
      geom_point(data = gene_var_stats, 
                 aes(x = log10(g.means), 
                     y = dispersion)) + 
      geom_line(data = fit.df, 
                aes(x = x,
                    y = y),
                color = "blue")
    
  }
  
  if(!is.null(plot_file)) {
    pdf(plot_file)
    qq_z_plot
    qq_loess.z_plot
    z_density_plot
    loess.z_density_plot
    dispersion_fit_plot
    dev.off()
  }
  
  if(return_type == "data") {
    return(gene_var_stats)
  } else if(return_type == "plots") {
    return(list(qq_z_plot = qq_z_plot,
                qq_loess.z_plot = qq_loess.z_plot,
                z_density_plot = z_density_plot,
                loess.z_density_plot = loess.z_density_plot,
                dispersion_fit_plot = dispersion_fit_plot))
  } else if(return_type == "both") {
    return(list(gene_var_stats = gene_var_stats,
                plots = list(qq_z_plot = qq_z_plot,
                             qq_loess.z_plot = qq_loess.z_plot,
                             z_density_plot = z_density_plot,
                             loess.z_density_plot = loess.z_density_plot,
                             dispersion_fit_plot = dispersion_fit_plot)))
  }
  
}
