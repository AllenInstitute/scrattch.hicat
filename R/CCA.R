#' Canon cor
#'
#' @param mat1 
#' @param mat2 
#' @param standardize 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
CanonCor <- function(mat1, mat2, standardize = TRUE, k = 20) {
  set.seed(seed = 42)
  if (standardize) {
    mat1 <- Standardize(mat = mat1, display_progress = FALSE)
    mat2 <- Standardize(mat = mat2, display_progress = FALSE)
  }
  mat3 <- crossprod(mat1,mat2)
  cca.svd <- irlba(A = mat3, nv = k)
  return(list(u = cca.svd$u, v = cca.svd$v, d = cca.svd$d))
}



##Extract from Seurat Package
#' CCA
#'
#' @param mat1 
#' @param mat2 
#' @param k 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
CCA <- function(mat1, mat2, k=20, verbose=FALSE)
  {
    library(irlba)
    cca.results <- CanonCor(mat1, mat2, standardize = FALSE, k=k)
    cca.data <- rbind(cca.results$u, cca.results$v)
    colnames(x = cca.data) <- paste0("CC", 1:num.cc)
    rownames(cca.data) <- c(colnames(data.use1), colnames(data.use2))
    cca.data <- apply(cca.data, MARGIN = 2, function(x){
      if(sign(x[1]) == -1) {
        x <- x * -1
      }
      return(x)
    })

    
    object = list(u, v)
    ###alignment
    num.groups=length(objects)
    
    for (cc.use in k) {
      for (g in 2:num.groups){
        if (verbose) {
          cat(paste0("Aligning dimension ", cc.use, "\n"), file = stderr())
        }
        genes.rank <- data.frame(rank(x = abs(x = cc.loadings[[1]][, cc.use])),
                                 rank(x = abs(x = cc.loadings[[g]][, cc.use])),
                                 cc.loadings[[1]][, cc.use],
                                 cc.loadings[[g]][, cc.use]
                                 )
        genes.rank$min <- apply(X = genes.rank[,1:2], MARGIN = 1, FUN = min)
        genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
        genes.top <- rownames(x = genes.rank)[1:min(num.possible.genes, nrow(genes.rank))]
        bicors <- list()
        for (i in c(1, g)) {
          cc.vals <- cc.embeds[[i]][, cc.use]
          if(verbose) {
            bicors[[i]] <- pbsapply(
                                    X = genes.top,
                                    FUN = function(x) {
                                      return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
                                    }
                                    )
          } else {
            bicors[[i]] <- sapply(
                                  X = genes.top,
                                  FUN = function(x) {
                                    return(BiweightMidcor(x = cc.vals, y = scaled.data[[i]][x, ]))
                                  }
                                  )
          }
        }
        genes.rank <- data.frame(
                                 rank(x = abs(x = bicors[[1]])),
                                 rank(x = abs(x = bicors[[g]])),
                                 bicors[[1]],
                                 bicors[[g]]
                                 )
        genes.rank$min <- apply(X = abs(x = genes.rank[, 1:2]), MARGIN = 1, FUN = min)
                                        # genes must be correlated in same direction in both datasets
        genes.rank <- genes.rank[sign(genes.rank[,3]) == sign(genes.rank[,4]), ]
        genes.rank <- genes.rank[order(genes.rank$min, decreasing = TRUE), ]
        genes.use <- rownames(x = genes.rank)[1:min(num.genes, nrow(genes.rank))]
        if(length(genes.use) == 0) {
          stop("Can't align group ", g, " for dimension ", cc.use)
        }
        metagenes <- list()
        multvar.data <- list()
        for (i in c(1, g)) {
          scaled.use <- sweep(
                              x = scaled.data[[i]][genes.use, ],
                              MARGIN = 1,
                              STATS = sign(x = genes.rank[genes.use, which(c(1, g) == i) + 2]),
                              FUN = "*"
                              )
          scaled.use <- scaled.use[, names(x = sort(x = cc.embeds[[i]][, cc.use]))]
          metagenes[[i]] <- (
                             cc.loadings[[i]][genes.use, cc.use] %*% scaled.data[[i]][genes.use, ]
                             )[1, colnames(x = scaled.use)]
        }
        mean.difference <- mean(x = ReferenceRange(x = metagenes[[g]])) -
          mean(x = ReferenceRange(x = metagenes[[1]]))
        align.1 <- ReferenceRange(x = metagenes[[g]])
        align.2 <- ReferenceRange(x = metagenes[[1]])
        a1q <- quantile(x = align.1, probs =  seq(from = 0, to = 1, by = 0.001))
        a2q <- quantile(x = align.2, probs =  seq(from = 0, to = 1, by = 0.001))
        iqr <- (a1q - a2q)[100:900]
        iqr.x <- which.min(x = abs(x = iqr))
        iqrmin <- iqr[iqr.x]
        if (show.plots) {
          print(iqrmin)
        }
        align.2 <- align.2 + iqrmin
        alignment <- dtw(x = align.1,
                         y = align.2,
                         keep.internals = TRUE)
        alignment.map <- data.frame(alignment$index1, alignment$index2)
        alignment.map$cc_data1 <- sort(cc.embeds[[g]][, cc.use])[alignment$index1]
        alignment.map$cc_data2 <- sort(cc.embeds[[1]][, cc.use])[alignment$index2]
        alignment.map.orig <- alignment.map
        alignment.map$dups <- duplicated(x = alignment.map$alignment.index1) |
        duplicated(x = alignment.map$alignment.index1, fromLast = TRUE)
        alignment.map %>% group_by(alignment.index1) %>% mutate(cc_data1_mapped = ifelse(dups, mean(cc_data2), cc_data2)) -> alignment.map
        alignment.map <- alignment.map[! duplicated(x = alignment.map$alignment.index1), ]
        cc.embeds.all[names(x = sort(x = cc.embeds[[g]][, cc.use])), cc.use] <- alignment.map$cc_data1_mapped
        if (show.plots) {
          par(mfrow = c(3, 2))
          plot(x = ReferenceRange(x = metagenes[[1]]), main = cc.use)
          plot(x = ReferenceRange(x = metagenes[[g]]))
          plot(
               x = ReferenceRange(x = metagenes[[1]])[(alignment.map.orig$alignment.index2)],
               pch = 16
               )
          points(
                 x = ReferenceRange(metagenes[[g]])[(alignment.map.orig$alignment.index1)],
                 col = "red",
                 pch = 16,
                 cex = 0.4
                 )
          plot(x = density(x = alignment.map$cc_data1_mapped))
          lines(x = density(x = sort(x = cc.embeds[[1]][, cc.use])), col = "red")
          plot(x = alignment.map.orig$cc_data1)
          points(x = alignment.map.orig$cc_data2, col = "red")
        }
      }
    }
    
    
  }
