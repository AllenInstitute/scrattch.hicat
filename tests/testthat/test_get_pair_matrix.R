test_that("compare the function output with the real data", {
  
  coor <- get_pair_matrix(tasic16.dat, row.names(tasic16.dat)[1], colnames(tasic16.dat)[2])
  
  expect_that(coor, equals(tasic16.dat[1, 2]) )
})
