context("test-de.genes")
library(scrattch.hicat)

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
## DE_genes_pairs() tests
test_that(
  "DE_genes_pairs() needs tests.",
  {
    
  }
)
## DE_genes_pw() tests
test_that(
  "DE_genes_pw() needs tests.",
  {
    
  }
)
## de_pair() tests
test_that(
  "de_pair() needs tests.",
  {
    
  }
)
## de_score() tests
test_that(
  "de_score() needs tests.",
  {
    
  }
)
## de_score_pairs() tests
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
