context("test-annotate")
library(scrattch.hicat)

# Load glial test data
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                  "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glia_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

set.seed(42)
glial_test_cell_sample <- sample(1:nrow(glial_cells), floor(nrow(glial_cells) / 4))

glial_train_cells <- glial_cells[-glial_test_cell_sample, ]
glial_test_cells <- glial_cells[glial_test_cell_sample, ]

glial_train_data <- tasic_2016_counts[, glial_train_cells$sample_name]
glial_test_data <- tasic_2016_counts[, glial_test_cells$sample_name]

glial_train_cl <- as.factor(glial_train_cells$primary_type_id)
names(glial_train_cl) <- glial_train_cells$sample_name

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

## predict_annotate_cor() tests

## map_sampling() tests

## map_cv() tests

## compare_annotate() tests

## adjust_color() tests

## get_cl_df() tests

## match_cl() tests

## find_low_quality_cl() tests

## plot_low_qc() tests
