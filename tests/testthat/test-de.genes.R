context("test-de.genes")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_data <- log2(tasic_2016_counts[, glial_cells$sample_name] + 1)

glial_data_sparse <- as(glial_data, "dgCMatrix")

glial_cl <- as.factor(glial_cells$primary_type_id)
names(glial_cl) <- glial_cells$sample_name

glial_var <- apply(glial_data, 1, var)
glial_var <- glial_var[order(glial_var, decreasing = TRUE)]
glial_hv_genes <- names(glial_var[1:1000])

## de_param() tests
test_that(
  "de_param() returns a list of parameters.",
  {
    blank_param <- de_param()
    
    expect_is(blank_param, "list")
    expect_equal(length(blank_param), 9)
    
    test_param <- de_param(low.th = 1,
                           padj.th = 0.01,
                           lfc.th = 1,
                           q1.th = 0.5,
                           q2.th = NULL,
                           q.diff.th = 0.7,
                           de.score.th = 150,
                           min.cells = 4,
                           min.genes = 5)
    
    expect_is(test_param, "list")
    expect_equal(test_param, blank_param)
    
  }
)

test_that(
  "de_param() provides errors to guide users to usable parameters.",
  {
    expect_error(de_param(low.th = -1))
    
    expect_error(de_param(padj.th = -1))
    
    expect_error(de_param(padj.th = 2))
    expect_error(de_param(lfc.th = -1))
    
    expect_error(de_param(q1.th = -1))
    expect_error(de_param(q1.th = 2))
    
    expect_error(de_param(q2.th = -1))
    expect_error(de_param(q2.th = 2))
    
    expect_error(de_param(q.diff.th = -1))
    
    expect_error(de_param(de.score.th = -1))
    
    expect_error(de_param(min.cells = -1))
    expect_error(de_param(min.cells = 0))
    
    expect_error(de_param(min.genes = -1))
    expect_error(de_param(min.genes = 0))
  }
)

## vec_chisq_test() tests
test_that(
  "vec_chisq_test() correctly computes Chi-squared tests.",
  {
    n_x <- 17
    n_y <- 14
    n_x_success <- 5
    n_x_fail <- n_x - n_x_success
    n_y_success <- 7
    n_y_fail <- n_y - n_y_success
    
    test_matrix <- matrix(c(n_x_fail, n_x_success,
                            n_y_fail, n_y_success),
                          ncol = 2)
    
    chisq.test_result1 <- chisq.test(test_matrix)
    
    vec_chisq_test_result <- vec_chisq_test(x = n_x_success,
                                            x.total = n_x,
                                            y = n_y_success,
                                            y.total = n_y)
    
    expect_is(vec_chisq_test_result, "data.frame")
    
    expect_equal(vec_chisq_test_result$pval[1],
                 chisq.test_result1$p.value)
    
    expect_equal(vec_chisq_test_result$stats[1],
                 unname(chisq.test_result1$statistic))
    
    
    n_x <- 423
    n_y <- 1567
    n_x_success <- 132
    n_x_fail <- n_x - n_x_success
    n_y_success <- 353
    n_y_fail <- n_y - n_y_success
    
    test_matrix <- matrix(c(n_x_fail, n_x_success,
                            n_y_fail, n_y_success),
                          ncol = 2)
    
    chisq.test_result2 <- chisq.test(test_matrix)
    
    vec_chisq_test_result <- vec_chisq_test(x = n_x_success,
                                            x.total = n_x,
                                            y = n_y_success,
                                            y.total = n_y)
    
    expect_is(vec_chisq_test_result, "data.frame")
    
    expect_equal(vec_chisq_test_result$pval[1],
                 chisq.test_result2$p.value)
    
    expect_equal(vec_chisq_test_result$stats[1],
                 unname(chisq.test_result2$statistic))
    
    
  }
)

## de_pair_limma() tests
test_that(
  "de_pair_limma() uses limma to run DEGene tests.",
  {
    limma_data <- glial_data[glial_hv_genes,]
    
    limma_cl <- setNames(as.factor(paste0("cl",glial_cl)),names(glial_cl))
    design <- model.matrix(~0 + limma_cl)
    colnames(design) <- levels(as.factor(limma_cl))
    
    fit <- limma::lmFit(object = limma_data[, names(limma_cl)], 
                        design = design)
    
    glial_means <- get_cl_means(limma_data,
                                glial_cl)
    
    glial_props <- get_cl_prop(limma_data,
                               glial_cl)
    
    spl_result <- de_pair_limma(pair = c("43","46"),
                                   cl.present = glial_props,
                                   cl.means = glial_means,
                                   design = design,
                                   fit = fit,
                                   genes = glial_hv_genes)
    
    expect_is(spl_result, "data.frame")
    expect_equal(rownames(spl_result), glial_hv_genes)
    expect_equal(sum(spl_result$padj <= 1), nrow(spl_result))
    expect_equal(sum(spl_result$pval <= 1), nrow(spl_result))
    expect_equal(sum(spl_result$lfc > 0), sum(spl_result$meanA > spl_result$meanB))
    
    # Exact result here to test for changes affecting results
    expect_equal(sum(spl_result$padj < 0.01), 290)
    
    # Spot check for means
    check_genes <- sample(glial_hv_genes, 10)
    expect_equal(spl_result[check_genes, "meanA"], unname(glial_means[check_genes, "43"]))
    expect_equal(spl_result[check_genes, "meanB"], unname(glial_means[check_genes, "46"]))
    
    # Spot check for props
    expect_equal(spl_result[check_genes, "q1"], unname(glial_props[check_genes, "43"]))
    expect_equal(spl_result[check_genes, "q2"], unname(glial_props[check_genes, "46"]))
    
    # Spot check for log fold change
    test_lfc <- glial_means[check_genes, "43"] - glial_means[check_genes, "46"]

    expect_equal(spl_result[check_genes, "lfc"], unname(test_lfc))
  }
)

## de_pair_chisq() tests
test_that(
  "de_pair_chisq() performs Chi-squared tests for a pair of clusters.",
  {
    
    hv_data <- glial_data[glial_hv_genes,]
    glial_means <- get_cl_means(hv_data,
                                glial_cl)
    
    glial_props <- get_cl_prop(hv_data,
                               glial_cl)
    
    glial_cl_size <- table(glial_cl)
    
    chi_result <- de_pair_chisq(pair = c("43","46"),
                                   cl.present = glial_props,
                                   cl.means = glial_means,
                                   cl.size = glial_cl_size,
                                   genes = glial_hv_genes)
    
    expect_is(chi_result, "data.frame")
    expect_equal(rownames(chi_result), glial_hv_genes)
    expect_equal(sum(chi_result$padj <= 1), nrow(chi_result))
    expect_equal(sum(chi_result$pval <= 1), nrow(chi_result))
    expect_equal(sum(chi_result$lfc > 0), sum(chi_result$meanA > chi_result$meanB))
    
    # Exact result here to test for changes affecting results
    expect_equal(sum(chi_result$padj < 0.01), 156)
    
    # Spot check for means
    check_genes <- sample(glial_hv_genes, 10)
    expect_equal(chi_result[check_genes, "meanA"], unname(glial_means[check_genes, "43"]))
    expect_equal(chi_result[check_genes, "meanB"], unname(glial_means[check_genes, "46"]))
    
    # Spot check for props
    expect_equal(chi_result[check_genes, "q1"], unname(glial_props[check_genes, "43"]))
    expect_equal(chi_result[check_genes, "q2"], unname(glial_props[check_genes, "46"]))
    
    # Spot check for log fold change
    test_lfc <- glial_means[check_genes, "43"] - glial_means[check_genes, "46"]

    expect_equal(chi_result[check_genes, "lfc"], unname(test_lfc))
  }
)

## de_selected_pairs() tests
test_that(
  "de_selected_pairs() performs pairwise DE gene tests using chisq.",
  {
    test_pairs <- matrix(c(43, 43, 44, 46), ncol = 2)
    
    de_results <- de_selected_pairs(norm.dat = glial_data, 
                                    cl = glial_cl,
                                    pairs = test_pairs, 
                                    method = "chisq", 
                                    low.th = 1, 
                                    min.cells = 4, 
                                    cl.present = NULL, 
                                    use.voom = FALSE, 
                                    counts = NULL,
                                    mc.cores = 1)
    
    expect_is(de_results, "list")
    expect_equal(length(de_results), nrow(test_pairs))
    
  }
)

test_that(
  "de_selected_pairs() performs pairwise DE gene tests using limma without voom.",
  {
    test_pairs <- matrix(c(43, 43, 44, 46), ncol = 2)
    
    de_results <- de_selected_pairs(norm.dat = glial_data, 
                                    cl = glial_cl,
                                    pairs = test_pairs, 
                                    method = "limma", 
                                    low.th = 1, 
                                    min.cells = 4, 
                                    cl.present = NULL, 
                                    use.voom = FALSE, 
                                    counts = NULL,
                                    mc.cores = 1)
    
    expect_is(de_results, "list")
    expect_equal(length(de_results), nrow(test_pairs))
    
  }
)

test_that(
  "de_selected_pairs() performs pairwise DE gene tests using limma with voom.",
  {
    test_pairs <- matrix(c(43, 43, 44, 46), ncol = 2)
    
    count_data <- tasic_2016_counts[, glial_cells$sample_name]
    
    de_results <- de_selected_pairs(norm.dat = glial_data, 
                                    cl = glial_cl,
                                    pairs = test_pairs, 
                                    method = "limma", 
                                    low.th = 1, 
                                    min.cells = 4, 
                                    cl.present = NULL, 
                                    use.voom = TRUE, 
                                    counts = count_data,
                                    mc.cores = 1)
    
    expect_is(de_results, "list")
    expect_equal(length(de_results), nrow(test_pairs))
    
  }
)
## de_all_pairs() tests
test_that(
  "de_all_pairs() performs all pairwise tests for a set of clusters.",
  {
    de_results <- de_all_pairs(norm.dat = glial_data,
                               cl = glial_cl)
    
    n_clusters <- length(levels(glial_cl))
    n_combos <- sum(1:(n_clusters-1))
    
    expect_is(de_results, "list")
    expect_equal(length(de_results), n_combos)
    expect_is(de_results[[1]], "data.frame")
    
    # can we pass parameters?
    de_results2 <- de_all_pairs(norm.dat = glial_data,
                                cl = glial_cl,
                                method = "chisq")
    
    expect_is(de_results, "list")
    expect_equal(length(de_results), n_combos)
    expect_is(de_results[[1]], "data.frame")
  }
)

## de_stats_pair() tests
test_that(
  "de_stats_pair() needs tests.",
  {
    test_pairs <- matrix(c(43, 43, 44, 46), ncol = 2)
    
    de_results <- de_selected_pairs(norm.dat = glial_data, 
                                    cl = glial_cl,
                                    pairs = test_pairs, 
                                    method = "chisq", 
                                    low.th = 1, 
                                    min.cells = 4, 
                                    cl.present = NULL, 
                                    use.voom = FALSE, 
                                    counts = NULL,
                                    mc.cores = 1)
    
    de_result <- de_results[[1]]
    
    results <- de_stats_pair(df = de_result)
    
    expect_is(results, "list")
    expect_length(results, 10)
    
    expect_is(results$score, "numeric")
    expect_is(results$up.score, "numeric")
    expect_is(results$down.score, "numeric")
    expect_equal(results$down.score + results$up.score, results$score)
    
    expect_is(results$num,"integer")
    expect_is(results$up.num, "integer")
    expect_is(results$down.num, "integer")
    expect_equal(results$up.num + results$down.num, results$num)
    
    expect_is(results$genes, "character")
    expect_is(results$up.genes, "character")
    expect_is(results$down.genes, "character")
    expect_equal(sum(results$up.genes %in% results$genes), length(results$up.genes))
    expect_equal(sum(results$down.genes %in% results$genes), length(results$down.genes))
    expect_equal(length(results$genes), results$num)
    expect_equal(length(results$up.genes), results$up.num)
    expect_equal(length(results$down.genes), results$down.num)
  }
)

## de_stats_selected_pairs() tests
test_that(
  "de_score() needs tests.",
  {
    
  }
)
## de_stats_all_pairs() tests
test_that(
  "de_score_pairs() needs tests.",
  {
    
  }
)
## get_de_matrix() tests
test_that(
  "get_de_matrix() needs tests.",
  {
    
  }
)
## plot_de_num() tests
test_that(
  "plot_de_num() needs tests.",
  {
    
  }
)
## DE_genes_cat_by_cl() tests
test_that(
  "DE_genes_cat_by_cl() needs tests.",
  {
    
  }
)
## plot_de_lfc_num() tests
test_that(
  "plot_de_lfc_num() needs tests.",
  {
    
  }
)
