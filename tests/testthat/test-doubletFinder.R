context("test-doubletFinder")
library(scrattch.hicat)
library(tasic2016data)

glial_classes <- c("Astrocyte", "Endothelial Cell", "Microglia", 
                   "Oligodendrocyte", "Oligodendrocyte Precursor Cell")
glial_cells <- tasic_2016_anno[tasic_2016_anno$broad_type %in% glial_classes, ]
glial_cells <- glial_cells[glial_cells$secondary_type_id == 0, ]

glial_counts <- tasic_2016_counts[, glial_cells$sample_name]

glial_counts_sparse <- as(glial_counts, "dgCMatrix")

glial_data <- log2(glial_counts_sparse + 1)

glial_vg <- find_vg(glial_data,
                    plot_file = NULL,
                    return_type = "data")

glial_vg <- glial_vg[order(glial_vg$g.vars,
                           decreasing = TRUE),]

glial_genes <- as.character(glial_vg$gene[1:200])

## doubletFinder() tests

test_that(
  "doubletFinder generates artificial doublets and returns doublet scores",
  {
    df_results <- doubletFinder(glial_counts_sparse,
                                select.genes = glial_genes,
                                proportion.artificial = 0.2,
                                k = 10,
                                plot = FALSE)
    
    expect_type(df_results, "double")
    expect_equal(names(df_results), colnames(glial_counts_sparse))
    expect_equal(length(df_results), ncol(glial_counts_sparse))
    
  }
)