context("test-utils")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

set.seed(42)
glial_test_cell_sample <- sample(1:nrow(glial_cells), floor(nrow(glial_cells) / 4))

glial_train_cells <- glial_cells[-glial_test_cell_sample, ]
glial_test_cells <- glial_cells[glial_test_cell_sample, ]

glial_train_data <- log2(tasic_2016_counts[, glial_train_cells$sample_name] + 1)
glial_test_data <- log2(tasic_2016_counts[, glial_test_cells$sample_name] + 1)

glial_var <- apply(glial_train_data, 1, var)
glial_var <- glial_var[order(glial_var, decreasing = TRUE)]
glial_hv_genes <- names(glial_var[1:1000])

glial_train_cl <- as.factor(glial_train_cells$primary_type_id)
names(glial_train_cl) <- glial_train_cells$sample_name
glial_test_cl <- as.factor(glial_test_cells$primary_type_id)
names(glial_test_cl) <- glial_test_cells$sample_name

train_cl.df <- unique(glial_train_cells[, grepl("primary_type", names(glial_train_cells))])
rownames(train_cl.df) <- train_cl.df$primary_type_id
names(train_cl.df) <- c("cluster_id","cluster_label","cluster_color")


test_that(
  "get_pair_matrix_coor() resolves vector positions using integers or names", 
  {
    row_idx <- 101:105
    col_idx <- 55:51
    
    int_coor <- get_pair_matrix_coor(m = tasic_2016_counts, 
                                     rows = row_idx, 
                                     cols = col_idx)
    
    expect_is(int_coor, "integer")
    
    expect_equal(length(int_coor), length(row_idx))
    
    row_names <- rownames(tasic_2016_counts)[row_idx]
    col_names <- colnames(tasic_2016_counts)[col_idx]
    
    chr_coor <- get_pair_matrix_coor(m = tasic_2016_counts,
                                     rows = row_idx,
                                     cols = col_idx)
    
    expect_is(chr_coor, "integer")
    
    expect_equal(length(chr_coor), length(row_names))
    
    expect_identical(int_coor, chr_coor)
    
  }
)

test_that(
  "get_pair_matrix() gets a vector of results for row and col positions.", 
  {
    
    row_idx <- 101:105
    col_idx <- 55:51
    
    int_vector <- get_pair_matrix(m = tasic_2016_counts, 
                                  rows = row_idx, 
                                  cols = col_idx)
    
    expect_is(int_vector, "numeric")
    
    expect_equal(length(int_vector), length(row_idx))
    
    expect_identical(int_vector, 
                     c(tasic_2016_counts[101,55],
                       tasic_2016_counts[102,54],
                       tasic_2016_counts[103,53],
                       tasic_2016_counts[104,52],
                       tasic_2016_counts[105,51]))
    
    row_names <- rownames(tasic_2016_counts)[row_idx]
    col_names <- colnames(tasic_2016_counts)[col_idx]
    
    chr_vector <- get_pair_matrix(m = tasic_2016_counts,
                                  rows = row_names,
                                  cols = col_names)
    
    expect_is(chr_vector, "numeric")
    
    expect_equal(length(chr_vector), length(row_idx))
    
    expect_identical(chr_vector, 
                     c(tasic_2016_counts[rownames(tasic_2016_counts)[101],colnames(tasic_2016_counts)[55]],
                       tasic_2016_counts[rownames(tasic_2016_counts)[102],colnames(tasic_2016_counts)[54]],
                       tasic_2016_counts[rownames(tasic_2016_counts)[103],colnames(tasic_2016_counts)[53]],
                       tasic_2016_counts[rownames(tasic_2016_counts)[104],colnames(tasic_2016_counts)[52]],
                       tasic_2016_counts[rownames(tasic_2016_counts)[105],colnames(tasic_2016_counts)[51]]))
    
    
    expect_identical(int_vector, chr_vector)
  }
)

test_that(
  "set_pair_matrix() updates a matrix using a vector of values and coordinates",
  {
    
    original_mat <- tasic_2016_counts[1:10, 1:10]
    
    row_idx <- 1:5
    col_idx <- 5:1
    new_values <- 11:15
    
    updated_mat_idx <- set_pair_matrix(m = original_mat, 
                                       rows = row_idx, 
                                       cols = col_idx, 
                                       vals = new_values)
    
    expect_is(updated_mat_idx, "matrix")
    
    expect_equal(sum(updated_mat_idx != original_mat), 5)
    expect_equal(new_values, c(updated_mat_idx[1,5],
                               updated_mat_idx[2,4],
                               updated_mat_idx[3,3],
                               updated_mat_idx[4,2],
                               updated_mat_idx[5,1]))
    
    row_names <- rownames(tasic_2016_counts)[row_idx]
    col_names <- colnames(tasic_2016_counts)[col_idx]
    
    
    updated_mat_chr <- set_pair_matrix(m = original_mat, 
                                       rows = row_names, 
                                       cols = col_names, 
                                       vals = new_values)
    
    expect_is(updated_mat_chr, "matrix")
    
    expect_equal(sum(updated_mat_chr != original_mat), 5)
    expect_equal(new_values, c(updated_mat_chr[1,5],
                               updated_mat_chr[2,4],
                               updated_mat_chr[3,3],
                               updated_mat_chr[4,2],
                               updated_mat_chr[5,1]))
    
    expect_identical(updated_mat_idx, updated_mat_chr)
    
    
  }
)

test_that(
  "get_pairs() splits underscore_separated names",
  {
    p1 <- c("a","b","x")
    p2 <- c("z","d","e")
    p1_p2 <- paste(p1, p2, sep = "_")
    
    pairs_df <- get_pairs(pairs.str = p1_p2)
    
    expect_is(pairs_df, "data.frame")
    expect_equal(nrow(pairs_df), length(p1_p2))
    expect_equal(rownames(pairs_df), p1_p2)
    expect_equal(pairs_df$P1, p1)
    expect_equal(pairs_df$P2, p2)
  }
)

test_that(
  "convert_pair_matrix() converts a vector of values to a matrix based on paired names of the values.",
  {
    
    p1 <- c("a","b","x")
    p2 <- c("z","d","e")
    p1_p2 <- paste(p1, p2, sep = "_")
    pair_values <- c(10:12)
    named_pair_values <- pair_values
    names(named_pair_values) <- p1_p2
    
    pair_matrix <- convert_pair_matrix(pair.num = named_pair_values, 
                                       l = NULL,
                                       directed = FALSE)
    
    expect_is(pair_matrix, "matrix")
    expect_true(is.numeric(pair_matrix))
    
    expect_equal(nrow(pair_matrix), length(unique(c(p1, p2))))
    expect_equal(ncol(pair_matrix), length(unique(c(p1, p2))))
    
    expect_equal(pair_values, c(pair_matrix[p1[1], p2[1]],
                                pair_matrix[p1[2], p2[2]],
                                pair_matrix[p1[3], p2[3]]))
    
    expect_equal(pair_values, c(pair_matrix[p2[1], p1[1]],
                                pair_matrix[p2[2], p1[2]],
                                pair_matrix[p2[3], p1[3]]))
    
    expect_identical(pair_matrix, t(pair_matrix))
    
    pair_matrix_directed <- convert_pair_matrix(pair.num = named_pair_values,
                                                l = NULL,
                                                directed = TRUE)
    
    expect_is(pair_matrix_directed, "matrix")
    expect_true(is.numeric(pair_matrix_directed))
    
    expect_equal(nrow(pair_matrix_directed), length(unique(c(p1, p2))))
    expect_equal(ncol(pair_matrix_directed), length(unique(c(p1, p2))))
    
    expect_equal(pair_values, c(pair_matrix_directed[p1[1], p2[1]],
                                pair_matrix_directed[p1[2], p2[2]],
                                pair_matrix_directed[p1[3], p2[3]]))
    
    expect_equal(rep(0, 3), c(pair_matrix_directed[p2[1], p1[1]],
                              pair_matrix_directed[p2[2], p1[2]],
                              pair_matrix_directed[p2[3], p1[3]]))
    
    expect_false(identical(pair_matrix, t(pair_matrix_directed)))
    expect_false(identical(pair_matrix_directed, t(pair_matrix_directed)))
    
  }
)

test_that(
  "convert_pair_matrix_str() converts a character vector of values to a matrix based on paired names of the values.",
  {
    
    p1 <- c("a","b","x")
    p2 <- c("z","d","e")
    p1_p2 <- paste(p1, p2, sep = "_")
    pair_values <- c("yes","no","maybe")
    named_pair_values <- pair_values
    names(named_pair_values) <- p1_p2
    
    pair_matrix <- convert_pair_matrix_str(pair.str = named_pair_values, 
                                           l = NULL,
                                           directed = FALSE)
    
    expect_is(pair_matrix, "matrix")
    expect_true(is.character(pair_matrix))
    
    expect_equal(nrow(pair_matrix), length(unique(c(p1, p2))))
    expect_equal(ncol(pair_matrix), length(unique(c(p1, p2))))
    
    expect_equal(pair_values, c(pair_matrix[p1[1], p2[1]],
                                pair_matrix[p1[2], p2[2]],
                                pair_matrix[p1[3], p2[3]]))
    
    expect_equal(pair_values, c(pair_matrix[p2[1], p1[1]],
                                pair_matrix[p2[2], p1[2]],
                                pair_matrix[p2[3], p1[3]]))
    
    expect_identical(pair_matrix, t(pair_matrix))
    
    pair_matrix_directed <- convert_pair_matrix_str(pair.str = named_pair_values,
                                                    l = NULL,
                                                    directed = TRUE)
    
    expect_is(pair_matrix_directed, "matrix")
    expect_true(is.character(pair_matrix_directed))
    
    expect_equal(nrow(pair_matrix_directed), length(unique(c(p1, p2))))
    expect_equal(ncol(pair_matrix_directed), length(unique(c(p1, p2))))
    
    expect_equal(pair_values, c(pair_matrix_directed[p1[1], p2[1]],
                                pair_matrix_directed[p1[2], p2[2]],
                                pair_matrix_directed[p1[3], p2[3]]))
    
    expect_equal(rep("", 3), c(pair_matrix_directed[p2[1], p1[1]],
                               pair_matrix_directed[p2[2], p1[2]],
                               pair_matrix_directed[p2[3], p1[3]]))
    
    expect_false(identical(pair_matrix, t(pair_matrix_directed)))
    expect_false(identical(pair_matrix_directed, t(pair_matrix_directed)))
    
  }
)

test_that(
  "get_cl_mat() generates a one-hot matrix for a set of clusters.",
  {
    cl_mat <- get_cl_mat(cl = glial_train_cl)
    
    expect_is(cl_mat, "dgCMatrix")
    
    expect_equal(ncol(cl_mat), length(levels(glial_train_cl)))
    expect_equal(nrow(cl_mat), length(glial_train_cl))
    expect_equal(sum(cl_mat), length(glial_train_cl))
    
    row_values <- Matrix::rowSums(cl_mat)
    names(row_values) <- NULL
    
    expect_equal(row_values, rep(1, length(glial_train_cl)))
    
    expect_equal(colnames(cl_mat), levels(glial_train_cl))
    expect_equal(rownames(cl_mat), names(glial_train_cl))
  }
)

test_that(
  "get_cl_sums() computes per-cluster sums from a matrix and cluster designations",
  {
    
    genes <- c("Opalin","Hspa8","Mbp")
    
    test_mat <- glial_test_data[genes,]
    
    cl_sums <- get_cl_sums(mat = test_mat,
                           cl = glial_test_cl)
    
    expect_is(cl_sums, "matrix")
    expect_equal(ncol(cl_sums), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_sums), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_sums["Opalin", "44"],
                 sum(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_sums["Hspa8", "46"],
                 sum(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
    # Using a sparse matrix
    test_mat_sparse <- as(test_mat, "dgCMatrix")
    
    cl_sums_sparse <- get_cl_sums(mat = test_mat_sparse,
                                  cl = glial_test_cl)
    
    expect_is(cl_sums_sparse, "matrix")
    expect_equal(ncol(cl_sums_sparse), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_sums_sparse), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_sums_sparse["Opalin", "44"],
                 sum(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_sums_sparse["Hspa8", "46"],
                 sum(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
    
  }
)

test_that(
  "get_cl_means() computes per-cluster means from a matrix and cluster designations",
  {
    
    genes <- c("Opalin","Hspa8","Mbp")
    
    test_mat <- glial_test_data[genes,]
    
    cl_means <- get_cl_means(mat = test_mat,
                             cl = glial_test_cl)
    
    expect_is(cl_means, "matrix")
    expect_equal(ncol(cl_means), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_means), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_means["Opalin", "44"],
                 mean(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_means["Hspa8", "46"],
                 mean(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
    # Using a sparse matrix
    test_mat_sparse <- as(test_mat, "dgCMatrix")
    
    cl_means_sparse <- get_cl_means(mat = test_mat_sparse,
                                    cl = glial_test_cl)
    
    expect_is(cl_means_sparse, "matrix")
    expect_equal(ncol(cl_means_sparse), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_means_sparse), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_means_sparse["Opalin", "44"],
                 mean(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_means_sparse["Hspa8", "46"],
                 mean(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
  }
)

test_that(
  "get_cl_medians() computes per-cluster medians from a matrix and cluster designations",
  {
    
    genes <- c("Opalin","Hspa8","Mbp")
    
    test_mat <- glial_test_data[genes,]
    
    cl_medians <- get_cl_medians(mat = test_mat,
                                 cl = glial_test_cl)
    
    expect_is(cl_medians, "matrix")
    expect_equal(ncol(cl_medians), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_medians), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_medians["Opalin", "44"],
                 median(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_medians["Hspa8", "46"],
                 median(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
    test_mat_sparse <- as(test_mat, "dgCMatrix")
    
    
    cl_medians_sparse <- get_cl_medians(mat = test_mat_sparse,
                                        cl = glial_test_cl)
    
    expect_is(cl_medians_sparse, "matrix")
    expect_equal(ncol(cl_medians_sparse), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_medians_sparse), nrow(test_mat))
    
    # Spot Checks
    expect_equal(cl_medians_sparse["Opalin", "44"],
                 median(test_mat["Opalin", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))]))
    
    expect_equal(cl_medians_sparse["Hspa8", "46"],
                 median(test_mat["Hspa8", which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))]))
    
  }
)


test_that(
  "get_cl_prop() computes per-cluster proportions from a matrix and cluster designations",
  {
    
    genes <- c("Opalin","Hspa8","Mbp")
    
    test_mat <- glial_test_data[genes,]
    
    test_thresh <- 1
    
    cl_prop <- get_cl_prop(mat = test_mat,
                           cl = glial_test_cl,
                           threshold = test_thresh)
    
    expect_is(cl_prop, "matrix")
    expect_equal(ncol(cl_prop), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_prop), nrow(test_mat))
    
    # Spot Checks
    cells_44 <- which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))
    expect_equal(cl_prop["Opalin", "44"],
                 sum(test_mat["Opalin", cells_44] > test_thresh) / length(cells_44))
    
    cells_46 <- which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))
    expect_equal(cl_prop["Hspa8", "46"],
                 sum(test_mat["Hspa8", cells_46] > test_thresh) / length(cells_46))
    
    # Using a sparse matrix
    test_mat_sparse <- as(test_mat, "dgCMatrix")
    
    cl_prop_sparse <- get_cl_prop(mat = test_mat_sparse,
                                  cl = glial_test_cl,
                                  threshold = test_thresh)
    
    expect_is(cl_prop_sparse, "matrix")
    expect_equal(ncol(cl_prop_sparse), length(levels(glial_test_cl)))
    expect_equal(nrow(cl_prop_sparse), nrow(test_mat))
    
    # Spot Checks
    cells_44 <- which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "44"]))
    expect_equal(cl_prop_sparse["Opalin", "44"],
                 sum(test_mat["Opalin", cells_44] > test_thresh) / length(cells_44))
    
    cells_46 <- which(colnames(test_mat) %in% names(glial_test_cl[glial_test_cl == "46"]))
    expect_equal(cl_prop_sparse["Hspa8", "46"],
                 sum(test_mat["Hspa8", cells_46] > test_thresh) / length(cells_46))
  }
)

test_that(
  "sparse_cor() computes column-wise correlations for a sparse matrix",
  {
    
    sparse_mat <- as(glial_test_data, "dgCMatrix")
    
    cor_results <- sparse_cor(sparse_mat)
    
    expect_is(cor_results, "matrix")
    expect_equal(nrow(cor_results), ncol(cor_results))
    expect_equal(ncol(cor_results), ncol(sparse_mat))
    expect_equal(rownames(cor_results), colnames(sparse_mat))
    expect_equal(colnames(cor_results), colnames(sparse_mat))
    
    expect_equal(cor_results, t(cor_results))
    
    # Spot checks
    
    vals_5 <- sparse_mat[,5]
    vals_13 <- sparse_mat[,13]
    vals_21 <- sparse_mat[,21]
    
    cor_5_13 <- as.vector(cor(vals_5, vals_13))
    expect_equal(cor_results[5,13], cor_5_13)
    cor_5_21 <- as.vector(cor(vals_5, vals_21))
    expect_equal(cor_results[5,21], cor_5_21)
    
  }
)

test_that(
  "calc_tau() computes tau scores for a set of genes",
  {
    tau_values <- calc_tau(m = glial_test_data,
                           byRow = TRUE)
    
    expect_is(tau_values, "numeric")
    expect_equal(length(tau_values), nrow(glial_test_data))
    
    sparse_mat <- as(glial_test_data, "dgCMatrix")
    
    tau_values_sparse <- calc_tau(m = sparse_mat,
                                  byRow = TRUE)
    
    expect_is(tau_values_sparse, "numeric")
    expect_equal(length(tau_values), nrow(sparse_mat))
    
    expect_equal(tau_values, tau_values_sparse)
  }
)

test_that(
  "sample_cells() equally downsamples cells from multiple clusters",
  {
    
    sample_size <- 5
    
    # NULL seed 1
    sampled_cells_random_1 <- sample_cells(cl = glial_train_cl, 
                                           sample.size = sample_size, 
                                           weights = NULL,
                                           seed = NULL)
    
    expect_is(sampled_cells_random_1, "character")
    
    sampled_cl_random_1 <- glial_train_cl[sampled_cells_random_1]
    
    expect_equal(length(sampled_cl_random_1), length(levels(glial_train_cl)) * sample_size)
    expect_true(sum(table(sampled_cl_random_1) == sample_size) == length(levels(glial_train_cl)))
    
    # NULL seed 2
    sampled_cells_random_2 <- sample_cells(cl = glial_train_cl, 
                                           sample.size = sample_size, 
                                           weights = NULL,
                                           seed = NULL)
    
    expect_is(sampled_cells_random_2, "character")
    
    sampled_cl_random_2 <- glial_train_cl[sampled_cells_random_2]
    
    expect_equal(length(sampled_cl_random_2), length(levels(glial_train_cl)) * sample_size)
    expect_true(sum(table(sampled_cl_random_2) == sample_size) == length(levels(glial_train_cl)))
    
    expect_true(length(setdiff(sampled_cells_random_1, sampled_cells_random_2)) > 1)
    expect_true(length(intersect(sampled_cells_random_1, sampled_cells_random_2)) < length(sampled_cells_random_1))
    
    # defined seed 1
    seed <- 42
    
    sampled_cells_set_1 <- sample_cells(cl = glial_train_cl, 
                                        sample.size = sample_size, 
                                        weights = NULL,
                                        seed = seed)
    
    expect_is(sampled_cells_set_1, "character")
    
    sampled_cl_set_1 <- glial_train_cl[sampled_cells_set_1]
    
    expect_equal(length(sampled_cl_set_1), length(levels(glial_train_cl)) * sample_size)
    expect_true(sum(table(sampled_cl_set_1) == sample_size) == length(levels(glial_train_cl)))
    
    # defined seed 2
    seed <- 42
    
    sampled_cells_set_2 <- sample_cells(cl = glial_train_cl, 
                                        sample.size = sample_size, 
                                        weights = NULL,
                                        seed = seed)
    
    expect_identical(sampled_cells_set_1, sampled_cells_set_2)
    
    # sampling more cells than are available
    # should just return all cells, not replace
    
    sample_size <- 50
    seed <- 42
    
    oversampled_cells <- sample_cells(cl = glial_test_cl,
                                      sample.size = sample_size,
                                      weights = NULL,
                                      seed = seed)
    
    expect_is(oversampled_cells, "character")
    
    expect_equal(length(oversampled_cells), length(glial_test_cl))
  }
)

test_that(
  "cpm() returns equivalent results for different type of matrices", 
  {
    library(Matrix)
    
    counts_mat <- tasic_2016_counts[1:100, 1:100]
    counts_dgt <- as(counts_mat, "dgTMatrix")
    counts_dgc <- as(counts_mat, "dgCMatrix")
    
    results_mat <- cpm(counts_mat)
    results_dgc <- cpm(counts_dgc)
    results_dgt <- cpm(counts_dgt)
    
    expect_is(results_mat, "matrix")
    expect_is(results_dgc, "dgCMatrix")
    expect_is(results_dgt, "dgTMatrix")
    
    values_mat <- as.vector(results_mat)
    values_mat <- values_mat[values_mat != 0]
    values_dgc <- results_dgc@x
    names(values_dgc) <- NULL
    values_dgt <- results_dgt@x
    
    expect_equal(values_mat, values_dgc)
    expect_equal(values_mat, values_dgt)
    expect_equal(values_dgc, values_dgt)
    
  }
)

test_that(
  "logCPM() returns equivalent results for different type of matrices", 
  {
    library(Matrix)
    
    counts_mat <- tasic_2016_counts[1:100, 1:100]
    counts_dgt <- as(counts_mat, "dgTMatrix")
    counts_dgc <- as(counts_mat, "dgCMatrix")
    
    results_mat <- logCPM(counts_mat)
    results_dgc <- logCPM(counts_dgc)
    results_dgt <- logCPM(counts_dgt)
    
    expect_is(results_mat, "matrix")
    expect_is(results_dgc, "dgCMatrix")
    expect_is(results_dgt, "dgTMatrix")
    
    values_mat <- as.vector(results_mat)
    values_mat <- values_mat[values_mat != 0]
    values_dgc <- results_dgc@x
    names(values_dgc) <- NULL
    values_dgt <- results_dgt@x
    
    expect_equal(values_mat, values_dgc)
    expect_equal(values_mat, values_dgt)
    expect_equal(values_dgc, values_dgt)
    
  }
)

test_that(
  "pair_cor() fails for mismatched matrices.",
  {
    mat_1 <- matrix(0, nrow = 20, ncol = 10)
    mat_2 <- matrix(0, nrow = 10, ncol = 20)
    
    expect_error(pair_cor(mat1 = mat_1,
                          mat2 = mat_2,
                          margin = 1))
    
    expect_error(pair_cor(mat1 = mat_1,
                          mat2 = mat_2,
                          margin = 2))
    
    expect_error(pair_cor(mat1 = mat_1,
                          mat2 = mat_2,
                          margin = 5))
  }
)


test_that(
  "pair_cor() runs paired correlations between matched matrices.",
  {
    mat_1 <- matrix(runif(200), nrow = 20, ncol = 10)
    mat_2 <- matrix(runif(200), nrow = 20, ncol = 10)
    mat_1_duplicate <- mat_1
    
    cor_results_1_duplicate <- pair_cor(mat1 = mat_1,
                                        mat2 = mat_1_duplicate,
                                        margin = 1)
    
    expect_is(cor_results_1_duplicate, "numeric")
    expect_equal(sum(cor_results_1_duplicate > 0.999), length(cor_results_1_duplicate))
    
    cor_results_1_2_rows <- pair_cor(mat1 = mat_1,
                                     mat2 = mat_2,
                                     margin = 1)
    
    expect_is(cor_results_1_2_rows, "numeric")
    expect_equal(length(cor_results_1_2_rows), nrow(mat_1))
    
    
    cor_results_1_2_cols <- pair_cor(mat1 = mat_1,
                                     mat2 = mat_2,
                                     margin = 2)
    
    expect_is(cor_results_1_2_cols, "numeric")
    expect_equal(length(cor_results_1_2_cols), ncol(mat_1))
    
    cor_results_t_1_2_cols <- pair_cor(mat1 = t(mat_1),
                                       mat2 = t(mat_2),
                                       margin = 2)
    
    expect_is(cor_results_t_1_2_cols, "numeric")
    expect_equal(length(cor_results_t_1_2_cols), nrow(mat_1))
    
    expect_equal(cor_results_t_1_2_cols, cor_results_1_2_rows)
    
  }
)
