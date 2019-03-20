context("Testing Human Functions")
# 1. First ever test
test_that("test check_qc function", {
  result_inside_test_1 <- check_qc(macaque.anno.feather$Genes.Detected.CPM, qc.iqr.mult = 3)
  result_inside_test_2 <- check_qc(macaque.anno.feather$percent_reads_aligned_total, qc.iqr.mult = 3)
  expect_known_output(result_inside_test_1, system.file("testdata", filename = "check_qc_result_1.rda", package = "scrattch.hicat")) 
  expect_known_output(result_inside_test_2, system.file("testdata", filename = "check_qc_result_2.rda", package = "scrattch.hicat"))
  })


# 2. check_qc sanity test
test_that("test check_qc sanity", {
  # length of output
  check_qc_result_1 <- check_qc(macaque.anno.feather$Genes.Detected.CPM, qc.iqr.mult = 3)
  check_qc_result_2 <- check_qc(macaque.anno.feather$percent_reads_aligned_total, qc.iqr.mult = 3)
  check_qc_result_length <- dim(macaque.anno.feather)[1]
  expect_length(check_qc_result_1, check_qc_result_length)
  expect_length(check_qc_result_2, check_qc_result_length)
  # type of output
  expect_that(check_qc_result_1, is_a("numeric") )
  expect_that(check_qc_result_2, is_a("numeric") )
  
  })




