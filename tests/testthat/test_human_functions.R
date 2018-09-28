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

# 1. First ever test
test_that("test check_qc function", {
  
  # loading input
  load(system.file("testdata", "glia_anno.Rdata", package = "scrattch.hicat"))
  load(system.file("testdata", "check_qc_result_1.RData", package = "scrattch.hicat"))
  load(system.file("testdata", "check_qc_result_2.RData", package = "scrattch.hicat")) 
  
  # running fnction
  result_inside_test_1 <- check_qc(anno$genes_detected_cpm_criterion, qc.iqr.mult = 3)
  result_inside_test_2 <- check_qc(anno$percent_aligned_reads_total, qc.iqr.mult = 3)
  
  # testing
  expect_known_output(result_inside_test_1, check_qc_result_1) 
  expect_known_output(result_inside_test_2, check_qc_result_2) 
})



