context("test-cluster")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_data <- log2(tasic_2016_counts[, glial_cells$sample_name] + 1)

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
  "knn_jaccard() needs tests.",
  {
    
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

## jaccard_louvain() tests
test_that(
  "jaccard_louvain() needs tests.",
  {
    
  }
)

## onestep_clust() tests
test_that(
  "onestep_clust() needs tests.",
  {
    
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
