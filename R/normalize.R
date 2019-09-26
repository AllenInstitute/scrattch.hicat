# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# lm_matrix()
#
# lm_normalize()
#   lm_matrix() normalize.R
# 

#' Perform a linear model fit on a matrix of data
#'
#' @param dat a matrix
#' @param x a vector
#'
#' @return a list with fit and residuals
#' @export
#'
lm_matrix <- function(dat, x) {
  
  coef <- colSums(t(dat) * x) / sum(x^2)
  
  fit <- matrix(coef, ncol = 1) %*% matrix(x, nrow = 1)
  
  resid <- dat - fit
  
  se_r <- rowVars(dat)
  se_f <- rowVars(resid)
  
  R_2 <- 1 - se_f / se_r
  R_2[is.na(R_2)] <- 0
  R_2 <- setNames(R_2, row.names(dat))
  
  n <- ncol(dat)
  
  F_stats <- (se_r - se_f) / (se_f / (n - 2))
  F_stats <- setNames(F_stats, row.names(dat))
  
  p.value <- df(F_stats, 1, n - 2)
  p.value[is.na(p.value)] <- 1
  padj <- p.adjust(p.value)
  
  fit <- data.frame(coef, R_2, F_stats, p.value, padj)
  
  return(list(fit = fit, 
              resid = resid))
}

# To do: Add batch substracted median

#' Normalize using a linear model
#'
#' @param dat a data matrix
#' @param x a vector
#' @param R_2.th R^2 threshold. Default is 0.2.
#' @param padj.th p-value threshold. Default is 0.01.
#' @param min.genes Minimum genes. Default is 5.
#'
#' @return
#' @export
#' 
lm_normalize <- function(dat, 
                         x, 
                         R_2.th = 0.2, 
                         padj.th = 0.01,
                         min.genes = 5) {
    
    m <- rowMedians(dat)
    m <- setNames(m, row.names(dat))
    
    q.max <- rowMaxs(dat)
    q.max <- setNames(q.max, row.names(dat))
    
    dat <- dat - m
    
    x <- x - median(x)
    
    norm.result <- lm_matrix(dat, x)
    fit <- norm.result$fit
    resid <- norm.result$resid
    
    select.genes <- row.names(fit)[fit$padj < padj.th & fit$R_2 > R_2.th]
    
    if(length(select.genes) >= min.genes) {
      dat[select.genes,] <- resid[select.genes,]
      dat <- dat + m
      
      # Avoid over-correction. Keep all values in the original range. 
      dat[dat < 0] <- 0
      dat[select.genes,] <- apply(dat[select.genes,], 
                                  2, 
                                  function(x){
                                    select <- x > q.max[select.genes]
                                    x[select] <- q.max[select]
                                    return(x)
                                  })
    } else {
      dat <- dat + m
    }
    
    return(list(dat, fit))
  
  }

