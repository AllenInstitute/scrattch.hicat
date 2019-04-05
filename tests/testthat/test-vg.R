context("test-vg")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_data <- log2(tasic_2016_counts[, glial_cells$sample_name] + 1)

## compute_vg() tests

test_that(
  "compute_vg() calculates gene variance statistics.",
  {
    vg_results <- compute_vg(glial_data)
    
    expect_is(vg_results, "data.frame")
    
    sparse_data <- as(glial_data, "dgCMatrix")
    
    sparse_vg_results <- compute_vg(sparse_data)
    
    expect_is(sparse_vg_results, "data.frame")
    
    expect_equal(sparse_vg_results, vg_results)
    
    ## Spot Checks
    check_genes <- c("Opalin", "Hspa8", "Mbp")
    glial_sample_totals <- colSums(glial_data)
    glial_scaling <- glial_sample_totals / median(glial_sample_totals)
    glial_data_scaled <- glial_data / glial_scaling[col(glial_data)]
    
    # g.means
    check_means <- rowMeans(glial_data_scaled[check_genes,])
    names(check_means) <- NULL
    expect_equal(vg_results$g.means[match(check_genes, vg_results$gene)],
                 check_means)
    
    # g.vars
    check_vars <- apply(glial_data_scaled[check_genes,], 1, 
                        function(x) {
                          # correct for differences in var calculation
                          # var() uses the N-1 denominator version; we use N.
                          var(x) * (length(x) - 1) / length(x) 
                        })
    names(check_vars) <- NULL
    expect_equal(vg_results$g.vars[match(check_genes, vg_results$gene)],
                 check_vars)
    
    # dispersion
    check_dispersion <- log10(check_vars / check_means)
    
    expect_equal(vg_results$dispersion[match(check_genes, vg_results$gene)],
                 check_dispersion)
    
    # z
    d <- vg_results$dispersion
    IQR <- quantile(d, probs = c(0.25, 0.75), na.rm = T)
    m <- mean(IQR)
    delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))
    check_z <- (d - m) / delta 
    
    expect_equal(vg_results$z[match(check_genes, vg_results$gene)],
                 check_z[match(check_genes, vg_results$gene)])
  }
)

## plot_vg() tests
test_that(
  "plot_vg() generates plots of gene variance statistics.",
  {
    
  }
)

## findVG() tests
test_that(
  "findVG() combines compute_vg() and plot_vg() as necessary.",
  {
    
  }
)
