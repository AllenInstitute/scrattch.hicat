context("test-annotate")
library(scrattch.hicat)
library(devtools)
# Load glial test data
devtools::install_github("AllenInstitute/tasic2016data")
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

## map_by_cor() tests
test_that(
  "map_by_cor() performs mapping to a reference using medians",
  {
    glial_mapping <- map_by_cor(train.dat = glial_train_data, 
                                train.cl = glial_train_cl, 
                                test.dat = glial_test_data,
                                method = "median")
    
    mapped_cluster_ids <- as.numeric(as.character(glial_mapping$pred.df$pred.cl))
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$pred.df, "data.frame")
    expect_is(glial_mapping$cor.matrix, "matrix")
    
    # expect most of them to match, but allow a few differences
    expect_gt(sum(mapped_cluster_ids == glial_test_cells$primary_type_id) / length(mapped_cluster_ids),
              0.95)
    
  }
)

test_that(
  "map_by_cor() performs mapping to a reference using means",
  {
    glial_mapping <- map_by_cor(train.dat = glial_train_data, 
                                train.cl = glial_train_cl, 
                                test.dat = glial_test_data,
                                method = "mean")
    
    mapped_cluster_ids <- as.numeric(as.character(glial_mapping$pred.df$pred.cl))
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$pred.df, "data.frame")
    expect_is(glial_mapping$cor.matrix, "matrix")
    
    # expect most of them to match, but allow a few differences
    expect_gt(sum(mapped_cluster_ids == glial_test_cells$primary_type_id) / length(mapped_cluster_ids),
              0.95)
    
  }
)

## map_cl_summary() tests
test_that(
  "map_cl_summary() performs mapping comparison using medians",
  {
    glial_test_cell_broad_cl <- factor(glial_test_cells$broad_type)
    names(glial_test_cell_broad_cl) <- glial_test_cells$sample_name
    
    glial_mapping <- map_cl_summary(ref.dat = glial_train_data,
                                    ref.cl = glial_train_cl,
                                    map.dat = glial_test_data,
                                    map.cl = glial_test_cell_broad_cl,
                                    method = "median")
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$map.df, "data.frame")
    expect_is(glial_mapping$cl.map.df, "data.frame")
    
  }
)

test_that(
  "map_cl_summary() performs mapping comparison using means",
  {
    glial_test_cell_broad_cl <- factor(glial_test_cells$broad_type)
    names(glial_test_cell_broad_cl) <- glial_test_cells$sample_name
    
    glial_mapping <- map_cl_summary(ref.dat = glial_train_data,
                                    ref.cl = glial_train_cl,
                                    map.dat = glial_test_data,
                                    map.cl = glial_test_cell_broad_cl,
                                    method = "mean")
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$map.df, "data.frame")
    expect_is(glial_mapping$cl.map.df, "data.frame")
    
  }
)

## compare_annotate() tests
test_that(
  "compare_annotate() fails if cl and ref.cl have to overlap",
  {
    glial_test_cell_broad_cl <- factor(glial_test_cells$broad_type)
    names(glial_test_cell_broad_cl) <- glial_test_cells$sample_name
    
    expect_error(compare_annotate(cl = glial_test_cell_broad_cl, 
                                         ref.cl = glial_train_cl, 
                                         ref.cl.df = train_cl.df, 
                                         reorder = TRUE,
                                         rename = TRUE))
  }
)

test_that(
  "compare_annotate() reports comparisons between cluster sets",
  {
    glial_train_cell_broad_cl <- factor(glial_train_cells$broad_type)
    names(glial_train_cell_broad_cl) <- glial_train_cells$sample_name
    
    glial_comparison <- compare_annotate(cl = glial_train_cell_broad_cl, 
                                         ref.cl = glial_train_cl, 
                                         ref.cl.df = train_cl.df, 
                                         reorder = TRUE,
                                         rename = TRUE)
    
    expect_is(glial_comparison, "list")
    expect_equal(length(glial_comparison), 6)
    expect_is(glial_comparison$cl, "factor")
    expect_is(glial_comparison$cl.df, "data.frame")
    expect_is(glial_comparison$g, "ggplot")
    expect_is(glial_comparison$tb.df, "data.frame")
    expect_is(glial_comparison$cl.id.map, "data.frame")
    
  }
)

## predict_annotate_cor() tests
test_that(
  "predict_annotate_cor() performs comparisons using medians and reports results",
  {
    glial_train_cell_broad_cl <- factor(glial_train_cells$broad_type)
    names(glial_train_cell_broad_cl) <- glial_train_cells$sample_name
    
    glial_mapping <- predict_annotate_cor(cl = glial_train_cell_broad_cl, 
                                          norm.dat = glial_train_data, 
                                          ref.markers = glial_hv_genes, 
                                          ref.cl = glial_train_cl, 
                                          ref.cl.df = train_cl.df, 
                                          ref.norm.dat = glial_train_data, 
                                          method = "median", 
                                          reorder = TRUE)
    
    expect_is(glial_mapping, "list")
    expect_equal(length(glial_mapping), 3)
    
    expect_is(glial_mapping$pred.df, "data.frame")
    expect_is(glial_mapping$cor.matrix, "matrix")
    expect_is(glial_mapping$annotate, "list")
    
    expect_equal(length(glial_mapping$annotate), 6)
  }
)

test_that(
  "predict_annotate_cor() performs comparisons using means and reports results",
  {
    glial_train_cell_broad_cl <- factor(glial_train_cells$broad_type)
    names(glial_train_cell_broad_cl) <- glial_train_cells$sample_name
    
    glial_mapping <- predict_annotate_cor(cl = glial_train_cell_broad_cl, 
                                          norm.dat = glial_train_data, 
                                          ref.markers = glial_hv_genes, 
                                          ref.cl = glial_train_cl, 
                                          ref.cl.df = train_cl.df, 
                                          ref.norm.dat = glial_train_data, 
                                          method = "mean", 
                                          reorder = TRUE)
    
    expect_is(glial_mapping, "list")
    expect_equal(length(glial_mapping), 3)
    
    expect_is(glial_mapping$pred.df, "data.frame")
    expect_is(glial_mapping$cor.matrix, "matrix")
    expect_is(glial_mapping$annotate, "list")
    
    expect_equal(length(glial_mapping$annotate), 6)
  }
)

## map_sampling() tests
test_that(
  "map_sampling() performs bootstrapped mapping using medians.",
  {
    glial_mapping <- map_sampling(train.dat = glial_train_data, 
                                  train.cl = glial_train_cl, 
                                  test.dat = glial_test_data, 
                                  markers = glial_hv_genes, 
                                  markers.perc = 0.8, 
                                  iter = 10, 
                                  method = "median",
                                  verbose = TRUE)
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$map.df, "data.frame")
    expect_is(glial_mapping$map.freq, "matrix")
  }
)

test_that(
  "map_sampling() performs bootstrapped mapping using means",
  {
    glial_mapping <- map_sampling(train.dat = glial_train_data, 
                                  train.cl = glial_train_cl, 
                                  test.dat = glial_test_data, 
                                  markers = glial_hv_genes, 
                                  markers.perc = 0.8, 
                                  iter = 10, 
                                  method = "mean",
                                  verbose = TRUE)
    
    expect_is(glial_mapping, "list")
    
    expect_equal(length(glial_mapping), 2)
    
    expect_is(glial_mapping$map.df, "data.frame")
    expect_is(glial_mapping$map.freq, "matrix")
  }
)

## map_cv() tests
test_that(
  "map_cv() performs cross-validation.",
  {
    glial_mapping <- map_cv(norm.dat = glial_train_data, 
                            cl = glial_train_cl, 
                            markers = glial_hv_genes, 
                            n.bin = 5,
                            g.perc = 1, 
                            method = "median",
                            verbose = TRUE)
    
    expect_is(glial_mapping, "character")
    expect_equal(length(glial_mapping), length(glial_train_cl))
  }
)

## adjust_color() tests
test_that(
  "adjust_color() removes duplicate hex colors.",
  {
    test_colors <- c("#00FF00","#00FF00","#FF0000","#FF0000","#00FF00","#000000")
    
    fixed_colors <- adjust_color(test_colors)
    
    expect_equal(length(test_colors), length(fixed_colors))
    expect_equal(length(unique(fixed_colors)), length(test_colors))
    expect_true(length(unique(test_colors)) < length(unique(fixed_colors)))
    expect_equal(sum(unique(test_colors) %in% fixed_colors), length(unique(test_colors)))
  }
)

test_that(
  "adjust_color() removes duplicate R colors.",
  {
    test_colors <- c("blue","blue","red","red","blue","black")
    
    fixed_colors <- adjust_color(test_colors)
    
    test_rgb <- col2rgb(test_colors)
    test_hex <- rgb(test_rgb[1,], test_rgb[2,], test_rgb[3,], maxColorValue = 255)
    
    expect_equal(length(test_colors), length(fixed_colors))
    expect_equal(length(unique(fixed_colors)), length(test_colors))
    expect_true(length(unique(test_colors)) < length(unique(fixed_colors)))
    expect_equal(sum(unique(test_hex) %in% fixed_colors), length(unique(test_colors)))
  }
)

test_that(
  "adjust_color() removes duplicate mixed R and hex colors.",
  {
    test_colors <- c("blue","#0000FF","red","#FF0000","blue","black")
    
    fixed_colors <- adjust_color(test_colors)
    
    test_rgb <- col2rgb(test_colors)
    test_hex <- rgb(test_rgb[1,], test_rgb[2,], test_rgb[3,], maxColorValue = 255)
    
    expect_equal(length(test_colors), length(fixed_colors))
    expect_equal(length(unique(fixed_colors)), length(test_colors))
    expect_true(length(unique(test_colors)) < length(unique(fixed_colors)))
    expect_equal(sum(unique(test_hex) %in% fixed_colors), length(unique(test_hex)))
  }
)

## get_cl_df() tests
test_that(
  "get_cl_df() initializes a cl.df data.frame.",
  {
    glial_cl.df <- get_cl_df(glial_train_cl)
    
    expect_is(glial_cl.df, "data.frame")
    expect_equal(nrow(glial_cl.df), length(levels(glial_train_cl)))
  }
)

## match_cl() tests
test_that(
  "match_cl() performs correlations between annotations.",
  {
    glial_test_cell_broad_cl <- factor(glial_test_cells$broad_type)
    names(glial_test_cell_broad_cl) <- glial_test_cells$sample_name
    
    glial_mapping <- match_cl(cl = glial_test_cell_broad_cl, 
                              dat = glial_test_data, 
                              ref.cl = glial_train_cl, 
                              ref.cl.df = train_cl.df, 
                              ref.dat = glial_train_data, 
                              rename = TRUE)
    
    expect_is(glial_mapping, "list")
    expect_equal(length(glial_mapping), 3)
    
    expect_is(glial_mapping$cl, "factor")
    expect_is(glial_mapping$cl.df, "data.frame")
    expect_is(glial_mapping$cor, "matrix")
    
  }
)



