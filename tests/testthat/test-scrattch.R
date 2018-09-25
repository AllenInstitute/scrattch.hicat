context("test-scrattch")
library(scrattch.hicat)

test_WGCNA_louvain_consistent <- function()
{
  require(mclust)
  de.param = de_param(q1.th=0.5, de.score.th=40)
  result = onestep_clust(tasic16.dat,dim.method="WGCNA", de.param = de.param)
  adj.rand.index=adjustedRandIndex(result$cl, tasic16.cl[names(result$cl)])
  print(adj.rand.index)
  adj.rand.index
}

test_PCA_louvain_consistent <- function()
{
  require(mclust)
  de.param = de_param(q1.th=0.5, de.score.th=40)
  result = onestep_clust(tasic16.dat,dim.method="PCA", de.param = de.param)
  adj.rand.index=adjustedRandIndex(result$cl, tasic16.cl[names(result$cl)])
  print(adj.rand.index)
  adj.rand.index
}

test_WGCNA_ward_consistent <- function()
{
  require(mclust)
  de.param = de_param(q1.th=0.5, de.score.th=40)
  result = onestep_clust(tasic16.dat, dim.method="WGCNA", method="ward.D", de.param = de.param)
  adj.rand.index=adjustedRandIndex(result$cl, tasic16.cl[names(result$cl)])
  print(adj.rand.index)
  adj.rand.index
}

test_PCA_ward_consistent <- function()
{
  require(mclust)
  de.param = de_param(q1.th=0.5, de.score.th=40)
  result = onestep_clust(tasic16.dat,dim.method="PCA", method="ward.D", de.param = de.param)
  adj.rand.index=adjustedRandIndex(result$cl, tasic16.cl[names(result$cl)])
  print(adj.rand.index)
  adj.rand.index
}

test_markers <- function()
{
  display.result= display_cl(tasic16.cl, tasic16.dat, de.param = de.param)
  return(length(display.result$markers))
}


test_that("one step clustering consistent", {
  expect_gt(test_WGCNA_louvain_consistent(), 0.3)
  expect_gt(test_PCA_louvain_consistent(), 0.3)
  expect_gt(test_WGCNA_ward_consistent(), 0.3)
  expect_gt(test_PCA_ward_consistent(), 0.3)
  expect_gt(test_markers(), 100)
})




