test_that("check when input is numeric", {
  
  coor <- get_pair_matrix_coor(tasic16.dat, 1, 2)

  expect_that(coor, is_a("numeric") )
  expect_that(length(coor), equals(1) )
})

test_that("check when input is character", {
  
  coor <- get_pair_matrix_coor(tasic16.dat, row.names(tasic16.dat)[1], colnames(tasic16.dat)[1])
  
  expect_that(coor, is_a("numeric") )
  expect_that(length(coor), equals(1) )
})
