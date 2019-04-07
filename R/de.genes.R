# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# de_param()
#
# vec_chisq_test()
#
# DE_genes_pairs()
#   get_cl_means() util.R
#   vec_chisq_test() de.genes.R
#
# DE_genes_pw()
#   DE_genes_pairs() de.genes.R
#
# de_pair()
#   
# de_score()
#   de_score_pairs() de.genes.R
#
# de_score_pairs()
#   DE_genes_pairs de.genes.R
#
# get_de_matrix()
#   get_pairs() de.genes.R
#   convert_pair_matrix() util.R
#
# plot_de_num()
#   get_de_matrix() de.genes.R
#   heatmap.3() heatmap.R
#
#  DE_genes_cat_by_cl
#    deScore.pairs() ??
#
# plot_de_lfc_num()

#' Set differential expression (DE) threshold for genes and clusters. 
#'
#' @param padj.th Threshold on adjusted Pvalue for limma. Default vaue is 0.01. 
#' @param lfc.th  Threshold on log2 fold change for limma. Default value is 1. 
#' @param low.th  Threshold on normalized gene expression to determiend whether a gene is detected in a cell. Default value 1. Users can specifiy different thresholds for different genes if necessary by using a vector with genes as names.
#' @param q1.th For every pair of clusters (one as foreground, and the other as background), we define q1, and q2 as the proportion of cells detecting the gene in the foregound and background cluster respectively. Up reguated genes should have q1 greater than q1.th
#' @param q2.th  Up regualted genes should have q2 smaller than q2.th. q2.th is NULL by default.  
#' @param q.diff.th Threshold on scaled proportion difference, defined as abs(q1 - q2)/max(q1, q2). DE genes should have differences greater than q.diff.th. Default 0.7 
#' @param de.score.th To determine whether two clusters are seperable based on DE genes, we define de.score as the sum of -log10(adjusted Pvalub) for all DE genes. Each gene contributes at most 20 towards the sum. All clusters should have pairwise de.score greater than de.score.th. 
#'
#' @return returns a list of parameters for reuse
#' @export
#'
de_param <- function(low.th = 1, 
                     padj.th = 0.01, 
                     lfc.th = 1, 
                     q1.th = 0.5, 
                     q2.th = NULL,
                     q.diff.th = 0.7, 
                     de.score.th = 100, 
                     min.cells = 4, 
                     min.genes = 5) {
  list(low.th = low.th, 
       padj.th = padj.th, 
       lfc.th = lfc.th, 
       q1.th = q1.th, 
       q2.th = q2.th, 
       q.diff.th = q.diff.th, 
       de.score.th = de.score.th, 
       min.cells = min.cells, 
       min.genes = min.genes)
}


vec_chisq_test <- function(x, 
                           x.total, 
                           y, 
                           y.total) {
  
  total <- x.total + y.total
  present <- x + y
  absent <- total - x - y
  
  o <- cbind(x, 
             x.total - x, 
             y, 
             y.total - y)
  
  e <- cbind(present * x.total, 
             absent * x.total, 
             present * y.total, 
             absent * y.total)
  
  e <- e / as.vector(total)
  
  stat <- rowSums(pmax(0, abs(o - e) - 0.5) ^ 2 / e)
  
  results <- cbind(stats = stat, 
                   pval = pchisq(stat, 1, lower.tail = FALSE), 
                   logFC = log2( (x * y.total) / (y * x.total)), 
                   diff = x / x.total - y / y.total)    
  results
}



DE_genes_pairs <- function(norm.dat, 
                           cl,
                           pairs, 
                           method = "limma", 
                           low.th = 1, 
                           min.cells = 4, 
                           cl.present = NULL, 
                           use.voom = FALSE, 
                           counts = NULL,
                           mc.cores = 1) {
  
  select.cl <- unique(c(pairs[,1], pairs[,2]))
  
  cl <- cl[cl %in% select.cl]
  
  select.genes <- row.names(norm.dat)[Matrix::rowSums(norm.dat[,names(cl)] >= low.th[row.names(norm.dat)]) >= min.cells]
  
  norm.dat <- as.matrix(norm.dat[select.genes,names(cl)])
  
  cl.means <- as.data.frame(get_cl_means(norm.dat, 
                                         cl))
  
  if(length(low.th) == 1) {
    low.th <- setNames(rep(low.th, nrow(norm.dat)), row.names(norm.dat))      
  }
  
  if(is.null(cl.present)){
    cl.present <- as.data.frame(get_cl_means(norm.dat >= low.th[row.names(norm.dat)],
                                             cl))
  }
  
  cl.size <- table(cl)
  
  fit <- NULL
  
  if(method=="limma"){
    cl <- setNames(as.factor(paste0("cl",cl)),names(cl))
    design <- model.matrix(~0 + cl)
    colnames(design) <- levels(as.factor(cl))
    
    if(use.voom & !is.null(counts)){
      v <- limma::voom(counts = as.matrix(counts[row.names(norm.dat), names(cl)]), 
                       design = design)
      
      fit <- lmFit(object = v, 
                   design = design)		
    } else {
      fit <- lmFit(object = norm.dat[, names(cl)], 
                   design = design)
    }
  }
  
  score_pair <- function(pair) {
    x <- as.character(pair[1])
    y <- as.character(pair[2])
    
    if(method == "limma") {
      ctr <<- paste(paste0("cl", x), "-", paste0("cl", y))
      
      contrasts.matrix <- limma::makeContrasts(contrasts = ctr, 
                                               levels = design)
      
      fit2 <- limma::contrasts.fit(fit = fit, 
                                   contrasts = contrasts.matrix)
      fit2 <- limma::eBayes(fit = fit2)
      
      pval <- fit2$p.value[, 1]
      padj <- p.adjust(pval)
      lfc <- coef(fit2)[, 1]
      
    } else if(method == "chisq") {
      
      chisq_results <- vec_chisq_test(cl.present[, x] * cl.size[[x]], 
                                      cl.size[[x]], 
                                      cl.present[, y] * cl.size[[y]], 
                                      cl.size[[y]])
      
      pval <- chisq_results[, "pval"]
      padj <- p.adjust(pval)
      lfc <- cl.means[, x] - cl.means[, y]
      
    }
    
    results <- data.frame(padj = padj,
                          pval = pval,
                          lfc = lfc,
                          meanA = cl.means[[x]], 
                          meanB = cl.means[[y]],
                          q1 = cl.present[[x]], 
                          q2 = cl.present[[y]])
    
    row.names(results) = row.names(norm.dat)
    
    return(results)
  }
  
  if(mc.cores == 1) {
    de.df <- list()
    
    for(i in 1:nrow(pairs)) {
      
      pair <- paste(pairs[i, 1], pairs[i, 2], sep = "_")    
      de.df[[pair]] <- score_pair(pairs[i,])
      
    }
  } else {
    library(foreach)
    
    cluster <- parallel::makeCluster(mc.cores)
    doParallel::registerDoParallel(cluster)
    
    de.df <- foreach::foreach(i = 1:nrow(pairs), 
                              .combine='c') %dopar% 
      list(score_pair(pairs[i,]))
    
    parallel::stopCluster(cluster)
    
    names(de.df) <- paste(pairs[,1],pairs[,2],sep="_")
  }
  
  return(de.df)
}

####Make sure dat and cl has the same dimension, and cells are in the same order
DE_genes_pw <- function(norm.dat,
                        cl, 
                        ...) {
  
  cn <- as.character(sort(unique(cl)))
  cl.n <- length(cn)	
  pairs <- cbind(rep(cn, rep(cl.n,cl.n)), 
                 rep(cn, cl.n))
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  
  de.df <- DE_genes_pairs(norm.dat = norm.dat,
                          cl = cl,
                          pairs = pairs, 
                          ...)
  return(de.df)
}


de_pair <- function(df,
                    de.param = de_param(), 
                    cl.size1 = NULL, 
                    cl.size2 = NULL) {
  
  df <- df[order(df$pval, -abs(df$lfc)), ]
  
  select <- with(df, which(padj < de.param$padj.th & abs(lfc) > de.param$lfc.th))
  select <- row.names(df)[select]
  
  if(is.null(select) | length(select) == 0){
    return(list())
  }
  
  up <- select[df[select, "lfc"] > 0]
  down <- select[df[select, "lfc"] < 0]
  df <- df[select, ]
  
  if(!is.null(de.param$q.diff.th) & is.null(df$q.diff)) {
    
    df$q.diff <- with(df, abs(q1 - q2) / pmax(q1, q2))
    df$q.diff[is.na(df$q.diff)] <- 0
    
  }
  
  if(!is.null(de.param$q1.th)) {
    up <- with(df[up, , drop = FALSE], up[q1 > de.param$q1.th])
    
    if(!is.null(cl.size1)){
      up <- with(df[up, , drop = FALSE], up[q1 * cl.size1 >= de.param$min.cells])
    }
    
    down <- with(df[down, , drop = FALSE], down[q2 > de.param$q1.th])
    
    if(!is.null(cl.size2)) {
      down <- with(df[down, , drop = FALSE], down[q2 * cl.size2 >= de.param$min.cells])
    }
  }
  
  if(!is.null(de.param$q2.th)) {
    up <- with(df[up, , drop = FALSE], up[q2 < de.param$q2.th])
    down <- with(df[down, , drop = FALSE], down[q1 < de.param$q2.th])
  }
  
  if(!is.null(de.param$q.diff.th)){
    up <- with(df[up, , drop = FALSE], up[abs(q.diff) > de.param$q.diff.th])
    down <- with(df[down, , drop = FALSE], down[abs(q.diff) > de.param$q.diff.th])
  }
  
  select <- c(up, down)
  
  if(length(select) == 0){
    return(list())
  } else {
    
    df$padj[df$padj < 1e-20] <- 1e-20
    up.score <- sum(-log10(df[up,"padj"]))
    down.score <- sum(-log10(df[down,"padj"]))
    
    if(length(up) == 0) { up.score <- 0 }
    if(length(down) == 0) { down.score <- 0 }
    
    list(score = up.score + down.score,
         up.score = up.score,
         down.score = down.score,
         num = length(select),
         up.num = length(up),
         down.num = length(down),
         genes = select,
         up.genes = up,
         down.genes = down, 
         de.df = df[df$padj < de.param$padj.th, ])
  }
  
}



#' Title
#'
#' @param norm.dat log transformed normalized gene expression matrix
#' @param cl  Cluster membership. 
#' @param min.cells Minimal number of cells in a cluster
#' @param de.param  DE gene criteria set by de_param function. 
#' @param method    Default "limma".
#' @param de.genes  If DE genes already computed for some pairs of clusters, use them without recomputation. 
#'
#' @return a list with DE genes from all pairs of clusters. 
#' @export
#'
de_score <- function(norm.dat, 
                     cl,  
                     de.param = de_param(), 
                     method = "limma", 
                     de.genes = NULL,
                     ...) {
  
  if(is.factor(cl)){
    cl <- droplevels(cl)
  }
  
  cn <- as.character(sort(unique(cl)))
  cl.n <- length(cn)	
  pairs <- cbind(rep(cn, rep(cl.n, cl.n)), rep(cn, cl.n))
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  row.names(pairs) <- paste(pairs[, 1], pairs[, 2], sep = "_")
  
  if(!is.null(de.genes)){
    pairs <- pairs[!row.names(pairs) %in% names(de.genes), , drop = FALSE]
  }
  
  de.result <- de_score_pairs(norm.dat, 
                              cl = cl, 
                              pairs = pairs, 
                              de.param = de.param, 
                              method = method,
                              ...)
  
  de.genes <- c(de.genes, de.result$de.genes)
  
  return(de.genes)
}

de_score_pairs <- function(norm.dat, 
                           cl, 
                           pairs, 
                           de.df = NULL, 
                           de.param = de_param(), 
                           method = "limma", 
                           mc.cores = 1) {
  
  row.names(pairs) <- paste(pairs[, 1], pairs[, 2], sep = "_")
  
  select.cl <- unique(c(pairs[, 1],pairs[, 2]))
  cl <- cl[cl %in% select.cl]
  
  norm.dat <- as.matrix(norm.dat[, names(cl)])
  
  if(is.factor(cl)) { cl <- droplevels(cl) }
  
  cl.size <- table(cl)
  cl.n <- names(cl.size)
  
  cl.small <- cl.n[cl.size < de.param$min.cells]
  cl.big <- setdiff(cl.n,cl.small)
  
  select.pair <- pairs[, 1] %in% cl.big & pairs[, 2] %in% cl.big
  
  de.genes <- list()
  
  if(sum(select.pair) > 0) {
    
    cl <- cl[cl %in% c(pairs[select.pair, 1], pairs[select.pair, 2])]
    select.cells <- names(cl)
    low.th <- de.param$low.th
    
    if(length(low.th) == 1){
      low.th <- setNames(rep(low.th, nrow(norm.dat)),
                         row.names(norm.dat))      
    }
    
    if(is.null(de.df)){
      de.df <- DE_genes_pairs(norm.dat, 
                              cl[select.cells], 
                              pairs[select.pair, , drop = FALSE], 
                              low.th = low.th,
                              min.cells = de.param$min.cells, 
                              method = method, 
                              mc.cores = mc.cores)
    }
    
    de.genes <- sapply(names(de.df), 
                       function(x) {
                         if(is.null(de.df[[x]])){
                           return(list())
                         }
                         
                         df = de.df[[x]]
                         if(!is.null(de.param$min.cells)) {
                           cl.size1 <- cl.size[as.character(pairs[x, 1])]
                           cl.size2 <- cl.size[as.character(pairs[x,2])]
                         } else {
                           cl.size1 <- NULL
                           cl.size2 <- NULL
                         }
                         de_pair(df, 
                                 de.param = de.param, 
                                 cl.size1, 
                                 cl.size2)
                       },
                       simplify = FALSE)
  } else {
    de.df <- list()
  }
  
  for(i in which(!select.pair)) {
    pair <- paste(pairs[i, 1], pairs[i, 2], sep = "_")
    de.genes[[pair]] <- list()
    de.df[[pair]] <- list()
  }
  
  list(de.df = de.df, 
       de.genes = de.genes)
}


get_de_matrix <- function(de.genes, 
                          directed = FALSE, 
                          field = "num") {
  
  pairs <- get_pairs(names(de.genes))
  
  if(directed){
    f <- paste("up", field, sep = ".")
    up.de.num <- sapply(de.genes, function(x) { x[[f]] } )
    
    f <- paste("down", field,sep=".")
    down.de.num <- sapply(de.genes, function(x) { x[[f]] } )
    
    # head(sort(pmin(up.de.num, down.de.num)))
    
    names(down.de.num) <- with(pairs[names(down.de.num), ], paste0(P2, "_", P1))
    de.num <- c(up.de.num, down.de.num)
    
  } else {
    de.num <- sapply(de.genes, function(x) { x[[field]] })
  }
  
  de.matrix <- convert_pair_matrix(de.num, 
                                   directed = directed)
  
  return(de.matrix)    
}

plot_de_num <- function(de.genes, 
                        dend, 
                        cl.label = NULL, 
                        directed = FALSE, 
                        file = "log10.de.num.pdf", 
                        ...) {
  
  label <- as.hclust(dend)$label
  de.num.matrix <- get_de_matrix(de.genes, 
                                 directed = directed)
  de.num.matrix <- de.num.matrix[label, label]
  
  breaks <- c(-1, seq(0.2, 4, length.out = 100))
  
  if(!is.null(cl.label)) {
    colnames(de.num.matrix) <- cl.label[row.names(de.num.matrix)]
    row.names(de.num.matrix) <- cl.label[row.names(de.num.matrix)]
  }
  
  tmp.dat <- log10(de.num.matrix + 1)
  
  pdf(file, ...)
  heatmap.3(tmp.dat, 
            col = jet.colors(100), 
            breaks = breaks,
            trace = "none",
            Colv = dend, 
            Rowv = dend,
            dendrogram = "row",
            cexRow = 0.3,
            cexCol = 0.3)
  dev.off()
}


DE_genes_cat_by_cl <- function(norm.dat, 
                               cl, 
                               binary.cat, 
                               ...) {
  
  cl <- droplevels(as.factor(cl))
  cl.cat <- setNames(paste0(cl, binary.cat[names(cl)]), 
                     names(cl))
  
  tmp <- levels(cl)
  cl.cat.pairs <- data.frame(cat1 = paste0(tmp, binary.cat[1]), 
                             cat2 = paste0(tmp, binary.cat[2]), 
                             stringsAsFactors = FALSE)
  
  cl.cat.de.genes <- deScore.pairs(norm.dat[,names(cl.cat)], 
                                   cl = cl.cat, 
                                   pairs = cl.cat.pairs, 
                                   ...)
  
  cat1.de.num <- sapply(cl.cat.de.genes,
                        function(x) {
                          if(length(x) == 0) { return(0) }
                          length(x$up.genes)
                        })
  
  cat2.de.num <- sapply(cl.cat.de.genes,
                        function(x) {
                          if(length(x) == 0) { return(0) }
                          length(x$down.genes)
                        })
  
  cat1.de.genes <- sapply(cl.cat.de.genes,
                          function(x) {
                            if(is.null(x)) { return("") }
                            
                            return(paste(head(x$up.genes, 8),
                                         collapse = " "))
                          })
  
  cat2.de.genes <- sapply(cl.cat.de.genes,
                          function(x) {
                            if(is.null(x)) { return("") }
                            
                            return(paste(head(x$down.genes, 8),
                                         collapse = " "))
                          })
  
  cl.cat.de.df <- data.frame(cat1.de.num, 
                             cat2.de.num, 
                             cat1.de.genes, 
                             cat2.de.genes)
  
  return(list(cl.cat.de.df = cl.cat.de.df,
              cl.cat.de.genes = cl.cat.de.genes))
}



plot_de_lfc_num <- function(de.genes, 
                            top.n = 100, 
                            select.pair = NULL, 
                            cl.label = NULL) {
  
  tmp <- sapply(de.genes, length)
  de.genes <- de.genes[tmp!=0]
  
  de.score <- sapply(de.genes, 
                     function(x) {
                       x$score
                     })
  
  de.up.score <- sapply(de.genes, 
                        function(x) {
    x$up.score
  })
  
  de.down.score <- sapply(de.genes, 
                          function(x) {
    x$down.score
  })
  
  de.num <- sapply(de.genes, 
                   function(x) {
    x$num
  })
  
  de.lfc <- sapply(de.genes, 
                   function(x) {
                     top.genes <- head(x$genes[order(x$de.df[x$genes, "pval"])],
                                       top.n)
                     
                     mean(abs(x$de.df[top.genes, "lfc"]))
                   })
  
  de.q.diff <- sapply(de.genes, 
                      function(x) {
                        top.genes <- head(x$genes[order(x$de.df[x$genes, "pval"])],
                                          top.n)
                        q.diff <- with(x$de.df[top.genes, ], abs(q1 - q2) / pmax(q1, q2))
                        mean(q.diff)
                      })
  
  de.summary <- data.frame(de.num, 
                           de.lfc, 
                           de.q.diff, 
                           de.score, 
                           de.up.score, 
                           de.down.score)
  
  row.names(de.summary) <- names(de.genes)
  
  tmp <- do.call("rbind", strsplit(row.names(de.summary), "_"))
  de.summary$cl1 <- tmp[, 1]
  de.summary$cl2 <- tmp[, 2]
  
  g <- ggplot2::ggplot(de.summary, 
                       ggplot2::aes(de.num,
                                    de.lfc,
                                    color = de.q.diff)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_color_gradient2(midpoint = 0.85) + 
    ggplot2::scale_x_log10()
  
  if(!is.null(select.pair)){
    select.df <- de.summary[select.pair, ]
    
    if(!is.null(cl.label)){
      select.df$pair.label <- with(select.df, 
                                   paste(cl.label[as.character(cl1)], 
                                         cl.label[as.character(cl2)],
                                         sep = ":"))
    }
    
    g <- g + 
      geom_text(data = select.df, 
                aes(de.num - 0.02, 
                    de.lfc, 
                    label = pair.label),
                size = 2,
                color = "black") +
      geom_point(data = select.df, 
                 aes(de.num, 
                     de.lfc),
                 color = "red",
                 pch = 1)
  }
  
  g <- g + 
    xlab("Number of DE genes") + 
    ylab("Mean log2(FC) of top 100 DE.genes")
  
  return(list(g = g, 
              de.summary = de.summary))
}


