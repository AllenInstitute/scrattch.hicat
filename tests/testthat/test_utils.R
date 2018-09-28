context("test-utils")
library(scrattch.hicat)


test_that("checking get pair matrix", {
  
  coor <- get_pair_matrix(tasic16.dat, row.names(tasic16.dat)[1], colnames(tasic16.dat)[2])
  
  expect_that(coor, equals(tasic16.dat[1, 2]) )
})


test_that("checking get pair matrix coor 1", {
  
  coor <- get_pair_matrix_coor(tasic16.dat, 1, 2)
  
  expect_that(coor, is_a("numeric") )
  expect_that(length(coor), equals(1) )
})


test_that("checking get pair matrix coor 2", {
  
  coor <- get_pair_matrix_coor(tasic16.dat, row.names(tasic16.dat)[1], colnames(tasic16.dat)[1])
  
  expect_that(coor, is_a("numeric") )
  expect_that(length(coor), equals(1) )
})


test_that("cpm is working for different type of matrices", {
  
  counts <- tasic_2016_counts[1:5, 1:5]
  counts_dgt <- as(counts, "dgTMatrix")
  counts_dgc <- as(counts, "dgCMatrix")
  
  results1 = cpm(counts)
  results2 = cpm(counts_dgc)
  results3 = cpm(counts_dgt)
  
  expect_that(all(results1[ which(!results1 == 0)] == results2@x), equals(TRUE))
  expect_that(all(results3@x == results2@x), equals(TRUE))
  
})

