context("test-cluster")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)
library(Matrix)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_counts <- tasic_2016_counts[, glial_cells$sample_name]

glial_counts_sparse <- as(glial_counts, "dgCMatrix")

glial_data <- logCPM(glial_counts)

glial_data_sparse <- as(glial_data, "dgCMatrix")

## jaccard() tests
test_that(
  "jaccard() computes jaccard distances between cells using a binary sparse matrix.",
  {
    glial_subset <- Matrix::t(glial_data_sparse[1:5000,])
    glial_subset@x <- rep(1, length(glial_subset@x))
    
    results <- jaccard(glial_subset)
    
    expect_is(results, "dgTMatrix")
    expect_equal(dim(results), c(nrow(glial_cells),nrow(glial_cells)))
    
    # spot checks
    check_cells <- c(34, 37)
    
    cell1 <- glial_subset[check_cells[1],]
    cell2 <- glial_subset[check_cells[2],]
    
    ol <- sum(cell1 > 0 & cell2 > 0)
    denom <- sum(cell1) + sum(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
    
    check_cells <- c(93, 52)
    
    cell1 <- glial_subset[check_cells[1],]
    cell2 <- glial_subset[check_cells[2],]
    
    ol <- sum(cell1 > 0 & cell2 > 0)
    denom <- sum(cell1) + sum(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
  }
)

test_that(
  "jaccard() computes jaccard distances between cells using a binary matrix.",
  {
    glial_subset <- t(glial_data[1:5000,])
    glial_subset[glial_subset != 0] <- 1
    
    results <- jaccard(glial_subset)
    
    expect_is(results, "dgTMatrix")
    expect_equal(dim(results), c(nrow(glial_cells),nrow(glial_cells)))
    
    # spot checks
    check_cells <- c(34, 37)
    
    cell1 <- glial_subset[check_cells[1],]
    cell2 <- glial_subset[check_cells[2],]
    
    ol <- sum(cell1 > 0 & cell2 > 0)
    denom <- sum(cell1) + sum(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
    
    check_cells <- c(93, 52)
    
    cell1 <- glial_subset[check_cells[1],]
    cell2 <- glial_subset[check_cells[2],]
    
    ol <- sum(cell1 > 0 & cell2 > 0)
    denom <- sum(cell1) + sum(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
  }
)

## knn_jaccard() tests
test_that(
  "knn_jaccard() computes jaccard distances between cells using the output of RANN::nn2().",
  {
    glial_subset <- t(glial_data[1:5000,])
    glial_knn <- RANN::nn2(glial_subset, k = 10)[[1]]
    
    results <- knn_jaccard(glial_knn)
    
    expect_is(results, "dgTMatrix")
    expect_equal(dim(results), c(nrow(glial_cells),nrow(glial_cells)))
    
    # spot checks
    check_cells <- c(34, 37)
    
    cell1 <- glial_knn[check_cells[1],]
    cell2 <- glial_knn[check_cells[2],]
    
    ol <- length(intersect(cell1, cell2))
    denom <- length(cell1) + length(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
    
    check_cells <- c(93, 52)
    
    cell1 <- glial_knn[check_cells[1],]
    cell2 <- glial_knn[check_cells[2],]
    
    ol <- length(intersect(cell1, cell2))
    denom <- length(cell1) + length(cell2) - ol
    js <- ol / denom
    
    expect_equal(js, results[check_cells[1],check_cells[2]])
  }
)

## pass_louvain() tests
test_that(
  "pass_jaccard() needs tests.",
  {
    
  }
)

## jaccard_louvain.RANN() tests
test_that(
  "jaccard_louvain.RANN() needs tests.",
  {
    
  }
)

## phenograph_jaccard_coeff() tests
test_that(
  "phenograph_jaccard_coeff() needs tests.",
  {
    
  }
)

## phenograph() tests
test_that(
  "phenograph() performs PhenoGraph clustering identical to Rphenograph",
  {
    pca <- rd_PCA(glial_data_sparse)
    pcs <- pca$rd.dat
    
    ref_results <- Rphenograph::Rphenograph(pcs,
                                        k = 10)
    results <- phenograph(pcs,
                          k = 10,
                          verbose = TRUE)
    
    # results[[1]] and ref_results[[1]] will
    # differ due to a random id assigned by igraph.
    
    expect_is(results, "list")
    expect_equal(length(results), 2)
    
    expect_identical(results[[2]],
                     ref_results[[2]])
  }
)

## jaccard_louvain() tests
test_that(
  "jaccard_louvain() needs tests.",
  {
    
  }
)

## onestep_clust() tests
test_that(
  "onestep_clust() louvain clustering.",
  {
    de.param <- de_param(low.th = 1,
                         padj.th = 0.01,
                         lfc.th = 1,
                         q1.th = 0.5,
                         q2.th = NULL,
                         q.diff.th = 0.7,
                         de.score.th = 150,
                         min.cells = 4,
                         min.genes = 5)
    
    # PCA for dimensionality reduction with louvain clustering
    onestep_result1 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts_sparse,
                                     vg.padj.th = 0.5,
                                     method = "louvain",
                                     dim.method = "pca",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result1, "list")
    expect_equal(length(onestep_result1), 2)
    
    expect_is(onestep_result1$cl, "factor")
    expect_equal(length(onestep_result1$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result1$cl), glial_cells$sample_name)
    
    # WGCNA for dimensionality reduction with louvain clustering
    onestep_result2 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts,
                                     vg.padj.th = 0.5,
                                     method = "louvain",
                                     dim.method = "WGCNA",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result2, "list")
    expect_equal(length(onestep_result2), 2)
    
    expect_is(onestep_result2$cl, "factor")
    expect_equal(length(onestep_result2$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result2$cl), glial_cells$sample_name)
    
  }
)

test_that(
  "onestep_clust() Rphenograph clustering.",
  {
    de.param <- de_param(low.th = 1,
                         padj.th = 0.01,
                         lfc.th = 1,
                         q1.th = 0.5,
                         q2.th = NULL,
                         q.diff.th = 0.7,
                         de.score.th = 150,
                         min.cells = 4,
                         min.genes = 5)
    
    # PCA for dimensionality reduction with Rphenograph clustering
    onestep_result1 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts_sparse,
                                     vg.padj.th = 0.5,
                                     method = "Rphenograph",
                                     dim.method = "pca",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result1, "list")
    expect_equal(length(onestep_result1), 2)
    
    expect_is(onestep_result1$cl, "factor")
    expect_equal(length(onestep_result1$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result1$cl), glial_cells$sample_name)
    
    # WGCNA for dimensionality reduction with Rphenograph clustering
    onestep_result2 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts,
                                     vg.padj.th = 0.5,
                                     method = "Rphenograph",
                                     dim.method = "WGCNA",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result2, "list")
    expect_equal(length(onestep_result2), 2)
    
    expect_is(onestep_result2$cl, "factor")
    expect_equal(length(onestep_result2$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result2$cl), glial_cells$sample_name)
    
  }
)

test_that(
  "onestep_clust() hclust clustering.",
  {
    de.param <- de_param(low.th = 1,
                         padj.th = 0.01,
                         lfc.th = 1,
                         q1.th = 0.5,
                         q2.th = NULL,
                         q.diff.th = 0.7,
                         de.score.th = 150,
                         min.cells = 4,
                         min.genes = 5)
    
    # PCA for dimensionality reduction with louvain clustering
    onestep_result1 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts_sparse,
                                     vg.padj.th = 0.5,
                                     method = "ward.D",
                                     dim.method = "pca",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result1, "list")
    expect_equal(length(onestep_result1), 2)
    
    expect_is(onestep_result1$cl, "factor")
    expect_equal(length(onestep_result1$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result1$cl), glial_cells$sample_name)
    
    # WGCNA for dimensionality reduction with louvain clustering
    onestep_result2 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts,
                                     vg.padj.th = 0.5,
                                     method = "ward.D",
                                     dim.method = "WGCNA",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result2, "list")
    expect_equal(length(onestep_result2), 2)
    
    expect_is(onestep_result2$cl, "factor")
    expect_equal(length(onestep_result2$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result2$cl), glial_cells$sample_name)
    
  }
)

test_that(
  "onestep_clust() kmeans clustering.",
  {
    de.param <- de_param(low.th = 1,
                         padj.th = 0.01,
                         lfc.th = 1,
                         q1.th = 0.5,
                         q2.th = NULL,
                         q.diff.th = 0.7,
                         de.score.th = 150,
                         min.cells = 4,
                         min.genes = 5)
    
    # PCA for dimensionality reduction with louvain clustering
    onestep_result1 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts_sparse,
                                     vg.padj.th = 0.5,
                                     method = "kmeans",
                                     dim.method = "pca",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result1, "list")
    expect_equal(length(onestep_result1), 2)
    
    expect_is(onestep_result1$cl, "factor")
    expect_equal(length(onestep_result1$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result1$cl), glial_cells$sample_name)
    
    # WGCNA for dimensionality reduction with louvain clustering
    onestep_result2 <- onestep_clust(glial_data_sparse, 
                                     select.cells = glial_cells$sample_name, 
                                     counts = glial_counts,
                                     vg.padj.th = 0.5,
                                     method = "kmeans",
                                     dim.method = "WGCNA",
                                     max.dim = 20,
                                     rm.eigen = NULL,
                                     rm.th = 0.7,
                                     de.param = de.param,
                                     merge.type = "undirectional",
                                     max.genes = 3000,
                                     sample.size = 4000,
                                     max.cl.size = 300,
                                     k.nn = 15,
                                     prefix = NULL,
                                     verbose = FALSE,
                                     regress.x = NULL)
    
    expect_is(onestep_result2, "list")
    expect_equal(length(onestep_result2), 2)
    
    expect_is(onestep_result2$cl, "factor")
    expect_equal(length(onestep_result2$cl), length(glial_cells$sample_name))
    expect_identical(names(onestep_result2$cl), glial_cells$sample_name)
    
  }
)

## iter_clust() tests
test_that(
  "iter_clust() needs tests.",
  {
    
  }
)

## reorder_cl() tests
test_that(
  "reorder_cl() needs tests.",
  {
    
  }
)

## combine_finer_split() tests
test_that(
  "combine_finer_split() needs tests.",
  {
    
  }
)
