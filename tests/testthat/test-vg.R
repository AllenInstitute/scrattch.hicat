context("test-vg")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_data <- log2(tasic_2016_counts[, glial_cells$sample_name] + 1)

glial_data_sparse <- as(glial_data, "dgCMatrix")

check_genes <- c("Opalin", "Hspa8", "Mbp")

check_idx <- match(check_genes, rownames(glial_data))

## rescale_samples() tests
test_that(
  "rescale_samples() rescales data for both matrix and dgCMatrix",
  {
    rescaled_data <- rescale_samples(glial_data)
    
    rescaled_data_sparse <- rescale_samples(glial_data_sparse)
    
    test_sums <- apply(glial_data, 2, sum)
    test_med <- median(test_sums)
    test_scales <- test_sums / test_med
    test_rescaled <- glial_data / test_scales[col(glial_data)]
    
    test_rescaled_sparse <- as(test_rescaled, "dgCMatrix")
    
    expect_is(rescaled_data, "matrix")
    expect_equal(rescaled_data, test_rescaled)
    
    expect_is(rescaled_data_sparse, "dgCMatrix")
    expect_equal(rescaled_data_sparse@x, test_rescaled_sparse@x)
    
    expect_equal(rescaled_data, as(rescaled_data_sparse, "matrix"))
  }
)

## gene_means() tests
test_that(
  "gene_means() computes means for both matrix and dgCMatrix",
  {
    # no rescaling
    mean_values <- gene_means(glial_data,
                              rescale = FALSE)
    
    mean_values_sparse <- gene_means(glial_data_sparse,
                                     rescale = FALSE)
    
    test_means <- apply(glial_data[check_genes,], 1, mean)
    
    expect_is(mean_values, "numeric")
    expect_equal(mean_values[check_idx], test_means)
    
    expect_is(mean_values_sparse, "numeric")
    expect_equal(mean_values_sparse[check_idx], test_means)
    
    expect_equal(mean_values, mean_values_sparse)
    
    # with rescaling
    mean_values <- gene_means(glial_data,
                              rescale = TRUE)
    
    mean_values_sparse <- gene_means(glial_data_sparse,
                                     rescale = TRUE)
    
    rescaled_data <- rescale_samples(glial_data)
    
    test_means <- gene_means(rescaled_data[check_idx,],
                             rescale = FALSE)
    
    expect_is(mean_values, "numeric")
    expect_equal(mean_values[check_idx], test_means)
    
    expect_is(mean_values_sparse, "numeric")
    expect_equal(mean_values_sparse[check_idx], test_means)
    
    expect_equal(mean_values, mean_values_sparse)
  }
)

## gene_vars() tests
test_that(
  "gene_vars() computes variances for both matrix and dgCMatrix",
  {
    # no rescaling, no means
    var_values1 <- gene_vars(glial_data,
                             means = NULL,
                             rescale = FALSE)
    
    var_values_sparse1 <- gene_vars(glial_data_sparse,
                                    means = NULL,
                                    rescale = FALSE)
    
    test_means <- apply(glial_data[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data[check_genes,] ^ 2) - test_means ^ 2
    
    expect_is(var_values1, "numeric")
    expect_equal(var_values1[check_idx], test_vars)
    
    expect_is(var_values_sparse1, "numeric")
    expect_equal(var_values_sparse1[check_idx], test_vars)
    
    expect_equal(var_values1, var_values_sparse1)
    
    # no rescaling, provided means
    mean_values <- gene_means(glial_data)
    
    var_values2 <- gene_vars(glial_data,
                             means = mean_values,
                             rescale = FALSE)
    
    var_values_sparse2 <- gene_vars(glial_data_sparse,
                                    means = mean_values,
                                    rescale = FALSE)
    
    test_means <- apply(glial_data[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data[check_genes,] ^ 2) - test_means ^ 2
    
    expect_is(var_values2, "numeric")
    expect_equal(var_values2[check_idx], test_vars)
    
    expect_is(var_values_sparse2, "numeric")
    expect_equal(var_values_sparse2[check_idx], test_vars)
    
    expect_equal(var_values2, var_values_sparse2)
    
    # with rescaling, no means
    var_values3 <- gene_vars(glial_data,
                             means = NULL,
                             rescale = TRUE)
    
    var_values_sparse3 <- gene_vars(glial_data_sparse,
                                    means = NULL,
                                    rescale = TRUE)
    
    glial_data_rescaled <- rescale_samples(glial_data)
    
    test_means <- apply(glial_data_rescaled[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data_rescaled[check_genes,] ^ 2) - test_means ^ 2
    
    expect_is(var_values3, "numeric")
    expect_equal(var_values3[check_idx], test_vars)
    
    expect_is(var_values_sparse3, "numeric")
    expect_equal(var_values_sparse3[check_idx], test_vars)
    
    expect_equal(var_values3, var_values_sparse3)
    
    # with rescaling, provided means
    mean_values <- gene_means(glial_data,
                              rescale = TRUE)
    
    var_values4 <- gene_vars(glial_data,
                             means = mean_values,
                             rescale = TRUE)
    
    var_values_sparse4 <- gene_vars(glial_data_sparse,
                                    means = mean_values,
                                    rescale = TRUE)
    
    glial_data_rescaled <- rescale_samples(glial_data)
    
    test_means <- apply(glial_data_rescaled[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data_rescaled[check_genes,] ^ 2) - test_means ^ 2
    
    expect_is(var_values4, "numeric")
    expect_equal(var_values4[check_idx], test_vars)
    
    expect_is(var_values_sparse4, "numeric")
    expect_equal(var_values_sparse4[check_idx], test_vars)
    
    expect_equal(var_values4, var_values_sparse4)
    
    # cross conditions
    expect_equal(var_values1, var_values2)
    expect_equal(var_values3, var_values4)
  }
)

## gene_dispersion() tests
test_that(
  "gene_dispersion() computes dispersion values using matrix and dgCMatrix inputs",
  {
    # no means or vars; no rescaling
    disp_values1 <- gene_dispersion(glial_data,
                                    means = NULL,
                                    vars = NULL,
                                    rescale = FALSE)
    
    disp_values_sparse1 <- gene_dispersion(glial_data_sparse,
                                           means = NULL,
                                           vars = NULL,
                                           rescale = FALSE)
    
    test_means <- apply(glial_data[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data[check_genes,] ^ 2) - test_means ^ 2
    test_disp <- log10(test_vars / test_means)
    
    expect_is(disp_values1, "numeric")
    expect_equal(disp_values1[check_idx], test_disp)
    
    expect_is(disp_values_sparse1, "numeric")
    expect_equal(disp_values_sparse1[check_idx], test_disp)
    
    expect_equal(disp_values1, disp_values_sparse1)
    
    # no means or vars; with rescaling
    disp_values2 <- gene_dispersion(glial_data,
                                    means = NULL,
                                    vars = NULL,
                                    rescale = TRUE)
    
    disp_values_sparse2 <- gene_dispersion(glial_data_sparse,
                                           means = NULL,
                                           vars = NULL,
                                           rescale = TRUE)
    
    glial_data_rescaled <- rescale_samples(glial_data)
    test_means <- apply(glial_data_rescaled[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data_rescaled[check_genes,] ^ 2) - test_means ^ 2
    test_disp <- log10(test_vars / test_means)
    
    expect_is(disp_values2, "numeric")
    expect_equal(disp_values2[check_idx], test_disp)
    
    expect_is(disp_values_sparse2, "numeric")
    expect_equal(disp_values_sparse2[check_idx], test_disp)
    
    expect_equal(disp_values2, disp_values_sparse2)
    
    # supplied means and vars; no rescaling
    mean_values <- gene_means(glial_data)
    var_values <- gene_vars(glial_data,
                            means = mean_values)
    
    disp_values3 <- gene_dispersion(glial_data,
                                    means = mean_values,
                                    vars = var_values,
                                    rescale = FALSE)
    
    disp_values_sparse3 <- gene_dispersion(glial_data_sparse,
                                           means = mean_values,
                                           vars = var_values,
                                           rescale = FALSE)
    
    test_means <- apply(glial_data[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data[check_genes,] ^ 2) - test_means ^ 2
    test_disp <- log10(test_vars / test_means)
    
    expect_is(disp_values3, "numeric")
    expect_equal(disp_values3[check_idx], test_disp)
    
    expect_is(disp_values_sparse3, "numeric")
    expect_equal(disp_values_sparse3[check_idx], test_disp)
    
    expect_equal(disp_values3, disp_values_sparse3)
    
    # supplied means and vars; with rescaling
    mean_values <- gene_means(glial_data,
                              rescale = TRUE)
    var_values <- gene_vars(glial_data,
                            means = mean_values,
                            rescale = TRUE)
    
    disp_values4 <- gene_dispersion(glial_data,
                                    means = mean_values,
                                    vars = var_values,
                                    rescale = TRUE)
    
    disp_values_sparse4 <- gene_dispersion(glial_data_sparse,
                                           means = mean_values,
                                           vars = var_values,
                                           rescale = TRUE)
    
    test_means <- apply(glial_data_rescaled[check_genes,], 1, mean)
    test_vars <- rowMeans(glial_data_rescaled[check_genes,] ^ 2) - test_means ^ 2
    test_disp <- log10(test_vars / test_means)
    
    expect_is(disp_values4, "numeric")
    expect_equal(disp_values4[check_idx], test_disp)
    
    expect_is(disp_values_sparse4, "numeric")
    expect_equal(disp_values_sparse4[check_idx], test_disp)
    
    expect_equal(disp_values4, disp_values_sparse4)
    
    # cross-condition comparisons
    # only conditions with matching rescaling will match
    expect_equal(disp_values1, disp_values3)
    expect_equal(disp_values2, disp_values4)
    
  }
)

## gene_z() tests
test_that(
  "gene_z() returns z values. No tests for accuracy.",
  {
    # no rescaling, no dispersions
    z_values1 <- gene_z(glial_data,
                        dispersions = NULL,
                        rescale = FALSE)
    
    z_values_sparse1 <- gene_z(glial_data_sparse,
                               dispersions = NULL,
                               rescale = FALSE)
    
    expect_is(z_values1, "numeric")
    
    expect_is(z_values_sparse1, "numeric")
    
    expect_equal(z_values1, z_values_sparse1)
    
    # no rescaling, provided dispersions
    disp_values <- gene_dispersion(glial_data)
    
    z_values2 <- gene_z(glial_data,
                        dispersions = disp_values,
                        rescale = FALSE)
    
    z_values_sparse2 <- gene_z(glial_data_sparse,
                               dispersions = disp_values,
                               rescale = FALSE)
    
    expect_is(z_values2, "numeric")
    
    expect_is(z_values_sparse2, "numeric")
    
    expect_equal(z_values2, z_values_sparse2)
    
    # with rescaling, no dispersions
    z_values3 <- gene_z(glial_data,
                        dispersions = NULL,
                        rescale = TRUE)
    
    z_values_sparse3 <- gene_z(glial_data_sparse,
                               dispersions = NULL,
                               rescale = TRUE)
    
    
    
    expect_is(z_values3, "numeric")
    
    expect_is(z_values_sparse3, "numeric")
    
    expect_equal(z_values3, z_values_sparse3)
    
    # with rescaling, provided dispersions
    disp_values <- gene_dispersion(glial_data,
                                   rescale = TRUE)
    
    z_values4 <- gene_z(glial_data,
                        dispersions = disp_values,
                        rescale = TRUE)
    
    z_values_sparse4 <- gene_z(glial_data_sparse,
                               dispersions = disp_values,
                               rescale = TRUE)
    
    expect_is(z_values4, "numeric")
    
    expect_is(z_values_sparse4, "numeric")
    
    expect_equal(z_values4, z_values_sparse4)
    
    # cross condition comparisons
    expect_equal(z_values1, z_values2)
    expect_equal(z_values3, z_values4)
  }
)

## gene_loess_z() tests
test_that(
  "gene_loess_fit() computes loess fit for dispersions. No tests for accuracy.",
  {
    # no rescaling, no dispersions
    loess_fit_values1 <- gene_loess_fit(glial_data,
                                        dispersions = NULL,
                                        means = NULL,
                                        rescale = FALSE)
    
    loess_fit_values_sparse1 <- gene_loess_fit(glial_data_sparse,
                                               dispersions = NULL,
                                               means = NULL,
                                               rescale = FALSE)
    
    expect_is(loess_fit_values1, "loess")
    
    expect_is(loess_fit_values_sparse1, "loess")
    
    expect_equal(loess_fit_values1, loess_fit_values_sparse1)
    
    # no rescaling, provided dispersions
    disp_values <- gene_dispersion(glial_data)
    mean_values <- gene_means(glial_data)
    
    loess_fit_values2 <- gene_loess_fit(glial_data,
                                        dispersions = disp_values,
                                        means = mean_values,
                                        rescale = FALSE)
    
    loess_fit_values_sparse2 <- gene_loess_fit(glial_data_sparse,
                                               dispersions = disp_values,
                                               means = mean_values,
                                               rescale = FALSE)
    
    expect_is(loess_fit_values2, "loess")
    
    expect_is(loess_fit_values_sparse2, "loess")
    
    expect_equal(loess_fit_values2, loess_fit_values_sparse2)
    
    # with rescaling, no dispersions
    loess_fit_values3 <- gene_loess_fit(glial_data,
                                        dispersions = NULL,
                                        means = NULL,
                                        rescale = TRUE)
    
    loess_fit_values_sparse3 <- gene_loess_fit(glial_data_sparse,
                                               dispersions = NULL,
                                               means = NULL,
                                               rescale = TRUE)
    
    
    
    expect_is(loess_fit_values3, "loess")
    
    expect_is(loess_fit_values_sparse3, "loess")
    
    expect_equal(loess_fit_values3, loess_fit_values_sparse3)
    
    # with rescaling, provided dispersions
    mean_values <- gene_means(glial_data,
                              rescale = TRUE)
    disp_values <- gene_dispersion(glial_data,
                                   rescale = TRUE)
    
    
    loess_fit_values4 <- gene_loess_fit(glial_data,
                                        dispersions = disp_values,
                                        means = mean_values,
                                        rescale = TRUE)
    
    loess_fit_values_sparse4 <- gene_loess_fit(glial_data_sparse,
                                               dispersions = disp_values,
                                               means = mean_values,
                                               rescale = TRUE)
    
    expect_is(loess_fit_values4, "loess")
    
    expect_is(loess_fit_values_sparse4, "loess")
    
    expect_equal(loess_fit_values4, loess_fit_values_sparse4)
    
    # cross condition comparisons
    expect_equal(loess_fit_values1, loess_fit_values2)
    expect_equal(loess_fit_values3, loess_fit_values4)
  }
)

## gene_loess_fit_z() tests
test_that(
  "gene_loess_fit_z() computes z values based on a loess fit and dispersion values.",
  {
    disp_values <- gene_dispersion(glial_data,
                                   means = NULL,
                                   vars = NULL,
                                   rescale = TRUE)
    
    loess_fit <- gene_loess_fit(glial_data,
                                dispersions = disp_values,
                                means = NULL,
                                rescale = TRUE)
    
    loess_z_values <- gene_loess_fit_z(loess_fit,
                                       dispersions = disp_values)
    
    expect_is(loess_z_values, "numeric")
    
    expect_equal(length(loess_z_values), nrow(glial_data))
    
  }
)

## compute_vg() tests
test_that(
  "compute_vg_stats() calculates gene variance statistics.",
  {
    vg_results <- compute_vg_stats(glial_data)
    
    expect_is(vg_results, "data.frame")
    
    sparse_data <- as(glial_data, "dgCMatrix")
    
    sparse_vg_results <- compute_vg_stats(sparse_data)
    
    expect_is(sparse_vg_results, "data.frame")
    
    expect_equal(sparse_vg_results, vg_results)
    
  }
)

## plot_vg() tests
test_that(
  "plot_vg() generates plots of gene variance statistics.",
  {
    glial_stats <- compute_vg_stats(glial_data)
    glial_loess <- gene_loess_fit(glial_data)
    
    expect_warning(plot_vg(glial_stats,
                           plots = "all",
                           loess_fit = NULL))
    
    all_plots <- suppressWarnings(plot_vg(glial_stats,
                                          plots = "all",
                                          loess_fit = NULL))
    
    expect_is(all_plots, "list")
    expect_equal(length(all_plots), 4)
    expect_equal(names(all_plots), c("qq_z_plot","qq_loess.z_plot","z_density_plot","loess.z_density_plot"))
    expect_is(all_plots$qq_z_plot, "ggplot")
    expect_is(all_plots$qq_loess.z_plot, "ggplot")
    expect_is(all_plots$z_density_plot, "ggplot")
    expect_is(all_plots$loess.z_density_plot, "ggplot")
    
    all_plots <- plot_vg(glial_stats,
                         plots = "all",
                         loess_fit = glial_loess)
    
    expect_is(all_plots, "list")
    expect_equal(length(all_plots), 5)
    expect_equal(names(all_plots), c("qq_z_plot","qq_loess.z_plot","z_density_plot","loess.z_density_plot","dispersion_fit_plot"))
    expect_is(all_plots$qq_z_plot, "ggplot")
    expect_is(all_plots$qq_loess.z_plot, "ggplot")
    expect_is(all_plots$z_density_plot, "ggplot")
    expect_is(all_plots$loess.z_density_plot, "ggplot")
    expect_is(all_plots$dispersion_fit_plot, "ggplot")
    
    one_plot <- plot_vg(glial_stats,
                        plots = "qq_z_plot")
    
    expect_is(one_plot, "list")
    expect_equal(length(one_plot), 1)
    expect_equal(names(one_plot), "qq_z_plot")
    expect_is(one_plot$qq_z_plot, "ggplot")
    
  }
)


list_remove <- function(x,
                        rem) {
  list_names <- names(x)
  x[list_names != rem]
}

## Remove plotenv - necessary for comparing plots.
strip_plot_env <- function(x) {
  if("list" %in% class(x)) {
    for(i in 1:length(x)) {
      x[[i]] <- list_remove(x[[i]], "plot_env")
    }
  } else {
    x <- list_remove(x, "plot_env")
  }
  x
}

## find_vg() tests
test_that(
  "find_vg() combines compute_vg() and plot_vg() as necessary.",
  {
    glial_stats <- compute_vg_stats(glial_data)
    glial_loess <- gene_loess_fit(glial_data,
                                  rescale = TRUE)
    
    all_plots <- plot_vg(glial_stats,
                         plots = "all",
                         loess_fit = glial_loess)
    
    all_plots <- strip_plot_env(all_plots)
    
    glial_both <- list(gene_var_stats = glial_stats,
                       plots = all_plots)
    
    find_vg_results <- find_vg(dat = glial_data, 
                               plot_file = NULL, 
                               return_type = "data")
    
    expect_is(find_vg_results, "data.frame")
    
    expect_equal(find_vg_results, glial_stats)
    
    find_vg_results2 <- find_vg(dat = glial_data,
                                plot_file = NULL,
                                return_type = "plots")
    
    find_vg_results2 <- strip_plot_env(find_vg_results2)
    
    expect_is(find_vg_results2, "list")
    
    expect_equal(find_vg_results2, all_plots)
    
    find_vg_results3 <- find_vg(dat = glial_data,
                                plot_file = NULL,
                                return_type = "both")
    
    find_vg_results3$plots <- strip_plot_env(find_vg_results3$plots)
    
    expect_is(find_vg_results3, "list")
    
    expect_equal(length(find_vg_results3), 2)
    
    expect_equal(find_vg_results3, glial_both)
  }
)

## Note: I can't get testthat to play nice with file writing
## -Lucas
# target_file <- "temp.pdf"
# 
# test_that(
#   "find_vg() can output multiple plots to a PDF",
#   {
#     
#     find_vg_results <- suppressWarnings(find_vg(dat = glial_data,
#                                                 plot_file = target_file,
#                                                 return_type = "both"))
#     
#     out_file_exists <- file.exists(target_file)
#     out_file_has_size <- file.size(target_file) > 0
#     
#     expect_true(out_file_exists)
#     expect_true(out_file_has_size)
#     
#   }
# )
# 
# file.remove(target_file)


## Old findVG for testing consistency between new and old version
original_findVG <- function(dat, 
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
  row.names(gene_var_stats) = gene_var_stats$gene
  
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

## Comparisons between find_vg() and original_findVG()
test_that(
  "find_vg() results match the original findVG() function",
  {
    new_results <- find_vg(glial_data,
                           plot_file = NULL,
                           return_type = "data")
    
    old_results <- original_findVG(glial_data,
                                   plot_file = NULL,
                                   return_type = "data")
    
    expect_equal(new_results$gene_var_stats, old_results$gene_var_stats)
  }
)
