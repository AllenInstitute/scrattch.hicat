context("Testing Human Functions")
#options(stringsAsFactors = FALSE) 
#exons    <- read.csv("human_MTG_2018-06-14_exon-matrix.csv",row.names = 1)
#introns  <- read.csv("human_MTG_2018-06-14_intron-matrix.csv",row.names = 1)
#geneInfo <- read.csv("human_MTG_2018-06-14_genes-rows.csv",row.names = 1)
#sampInfo <- read.csv("human_MTG_2018-06-14_samples-columns.csv",row.names = 1)

#ibrary("Matrix") 
#kpSamp   <- sampInfo$class == "Non-neuronal"
#exons2   <- as(as.matrix(exons[,kpSamp]),"dgCMatrix")
#introns2 <- as(as.matrix(introns[,kpSamp]),"dgCMatrix")
#anno     <- sampInfo[kpSamp,]
#rownames(exons2) <- rownames(introns2) <- rownames(geneInfo)

#norm.dat   <- Matrix(cpm(exons2 + introns2), sparse = TRUE)
#norm.dat@x <- log2(norm.dat@x + 1)

#save(anno, file =  "glia_anno.Rdata")
#save(norm.dat, file = "glia_data.Rdata")

#-------------------------------------------------------------------------------------------------------------------
load(system.file("testdata", "glia_anno.Rdata", package = "scrattch.hicat"))

# 1. First ever test on check_qc function
test_that("test check_qc function", {
  
  # loading input and output
  result_1 <- readRDS(system.file("testdata", "check_qc_result_1.RData", package = "scrattch.hicat"))
  result_2 <- readRDS(system.file("testdata", "check_qc_result_2.RData", package = "scrattch.hicat"))
  
  # running function
  result_inside_test_1 <- check_qc(anno$genes_detected_cpm_criterion, qc.iqr.mult = 3)
  result_inside_test_2 <- check_qc(anno$percent_aligned_reads_total, qc.iqr.mult = 3)
  
  # testing
  expect_equal(result_inside_test_1, result_1) 
  expect_equal(result_inside_test_2, result_2) 
})


# 2. check_qc sanity test
test_that("test check_qc sanity", {
 
  # running function
  result_inside_test_1 <- check_qc(anno$genes_detected_cpm_criterion, qc.iqr.mult = 3)
  result_inside_test_2 <- check_qc(anno$percent_aligned_reads_total, qc.iqr.mult = 3)
  check_qc_result_length <- dim(anno)[1]

  # length of output  
  expect_length(result_inside_test_1, check_qc_result_length)
  expect_length(result_inside_test_2, check_qc_result_length)
  
  # type of output
  expect_that(result_inside_test_1, is_a("numeric") )
  expect_that(result_inside_test_2, is_a("numeric") )
  
})

#-------------------------------------------------------------------------------------------------------------------

# 3. check_neun function test

test_that("test check_neun() function needs review", {
  
  # # loading input and output
  # result <- readRDS(system.file("testdata", "check_neun_result.RData", package = "scrattch.hicat"))
  # 
  # # running function
  # cluster_with_names <- setNames(anno$cluster, anno$sample_id)
  # result_inside_test <- check_neun(anno, cluster_with_names, 0.5, "facs_sort_criteria")
  # 
  # # testing
  # expect_equal(result_inside_test, result) 
})

# 4. check_neun sanity test

test_that("test check_neun sanity", {
  
  # running function
  cluster_with_names <- setNames(anno$cluster, anno$sample_id)
  result_inside_test <-  check_neun(anno, cluster_with_names, 0.5, "facs_sort_criteria")
  
  check_neun_result_length <- length(unique(anno$cluster))
  
  # length of output test
  expect_length(result_inside_test, check_neun_result_length)
  
  # type of output test
  expect_is(result_inside_test, "array" )
  
})

#-------------------------------------------------------------------------------------------------------------------


