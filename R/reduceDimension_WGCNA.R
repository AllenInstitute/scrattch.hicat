# Function call map
# function_1()
#   called_by_function_1() called_function_file.R
#
# score_gene_mod()
#   get_eigen() reduceDimension_WGCNA.R
#   jaccard_louvain() cluster.R
#   de_stats_all_pairs() de.genes.R
#   
# filter_gene_mod()
#   get_eigen() reduceDimension_WGCNA.R
#   score_gene_mod() reduceDimension_WGCNA.R
#
# get_eigen()
#   heatmap.3() heatmap.R
#   
# rd_WGCNA()
#   get_eigen() reduceDimension_WGCNA.R
#   filter_gene_mod() reduceDimension_WGCNA.R
#


#' Compute scores for a set of gene modules
#' 
#' @param norm.dat
#' @param select.cells
#' @param gene.mod
#' @param eigen
#' @param method
#' @param de.param
#' @param max.cl.size
#' 
#' @return 
#' @export
#' 
score_gene_mod <-  function(norm.dat, 
                            select.cells, 
                            gene.mod, 
                            eigen = NULL,
                            method = "average", 
                            de.param = de_param(), 
                            max.cl.size = NULL) {
  
  if(is.null(select.cells)) {
    select.cells <- colnames(norm.dat)
  }
  
  if(length(gene.mod) == 0){
    return(NULL)
  }
  
  if(is.null(eigen)){
    eigen <- get_eigen(gene.mod, 
                       norm.dat[,names(select.cells)])
  }
  colnames(eigen) <- names(gene.mod)
  
  gene.mod.val <- sapply(names(gene.mod), 
                         function(gene.mod.name) {
                           gene.mod.genes <- gene.mod[[gene.mod.name]]
                           
                           if(is.null(max.cl.size)) {
                             # If no max.cl.size, use all cells
                             gene.mod.dat <- norm.dat[gene.mod.genes, select.cells]
                           } else {
                             # Otherwise, use cells at the ends of the eigengene
                             v <- eigen[select.cells, gene.mod.name]
                             ord <- order(v)
                             
                             top.cells <- select.cells[head(ord, max.cl.size)]
                             bot.cells <- select.cells[tail(ord, max.cl.size)]
                             
                             select.cells <- unique(c(top.cells, bot.cells))
                             
                             gene.mod.dat <- norm.dat[gene.mod.genes, select.cells]
                           }
                           
                           gene.mod.dat <- as.matrix(gene.mod.dat)
                           
                           # Re-center each gene based on means
                           gene.mod.dat <- gene.mod.dat - rowMeans(gene.mod.dat)
                           
                           if(method == "average") {
                             gene.mod.dist <- dist(t(gene.mod.dat))
                             gene.mod.hclust <- hclust(gene.mod.dist,
                                                       method = "average")
                             gene.mod.cl <- cutree(gene.mod.hclust,
                                                   k = 2)
                           } else if(method == "ward.D") {
                             gene.mod.dist <- dist(t(gene.mod.dat))
                             gene.mod.hclust <- hclust(gene.mod.dist,
                                                       method = "ward.D")
                             gene.mod.cl <- cutree(gene.mod.hclust,
                                                   k = 2)
                           } else if(method == "kmeans") {
                             gene.mod.kmeans <- kmeans(t(gene.mod.dat),
                                                       centers = 2)
                             gene.mod.cl <- gene.mod.kmeans$cluster
                           } else if(method == "louvain") {
                             gene.mod.jl <- jaccard_louvain(t(gene.mod.dat), 
                                                            k = 15)
                             if(is.null(gene.mod.jl)) {
                               return(list(c(0,0),
                                           NULL))
                             }
                             
                             gene.mod.cl <- gene.mod.jl$cl
                           } else {
                             stop(paste("Unknown method", method))
                           }
                           
                           de.gene.stats <- de_stats_all_pairs(as.matrix(norm.dat[,select.cells]), 
                                                               cl = gene.mod.cl, 
                                                               de.param = de.param)
                           de.genes <- de.gene.stats$de.genes
                           
                           de.genes <- de.genes[sapply(de.genes, length) > 1]
                           
                           if(length(de.genes) > 0) {
                             sc <- max(sapply(de.genes, 
                                              function(x) {
                                                x$score
                                              }))
                             gene.num <- max(sapply(de.genes, 
                                                    function(x) {
                                                      x$num
                                                    }))
                           } else {
                             sc <- 0
                             gene.num <- 0
                           }
                           
                           return(list(c(sc=sc, 
                                         gene.num = gene.num),
                                       de.genes = de.genes))
                         }, simplify = FALSE)
  
  return(gene.mod.val)
}


filter_gene_mod <- function(norm.dat, 
                            select.cells, 
                            gene.mod, 
                            min.module.size = 10, 
                            min.deScore = 40, 
                            de.param = de_param(), 
                            max.cl.size = NULL,
                            rm.eigen = NULL, 
                            rm.th = 0.7, 
                            max.size = 200, 
                            prefix = "cl", 
                            max.mod = NULL) {
  
  eigen <- get_eigen(gene.mod = gene.mod, 
                     norm.dat = norm.dat,
                     select.cells = select.cells)[[1]]
  
  if(!is.null(rm.eigen)) {
    rm.cor <- cor(eigen, 
                  rm.eigen[select.cells,])
    
    rm.cor[is.na(rm.cor)] <- 0
    
    rm.score <- setNames(rowMaxs(abs(rm.cor)), 
                         colnames(eigen))
    print(rm.score)
    select <- rm.score < rm.th
    
    if(sum(!select)) {
      print("Remove module")
      print(rm.score[!select, drop = FALSE])
    }
    
    eigen <- eigen[, select, drop = FALSE]
    gene.mod <- gene.mod[select]
  }
  
  if(length(select.cells) > 4000){
    method <- "louvain"
  } else {
    method <- c("ward.D","kmeans")
  }
  
  nmod <- min(20, length(gene.mod))
  if(nmod == 0) {
    return(NULL)
  }
  
  gene.mod <- head(gene.mod, nmod)
  eigen <- eigen[, 1:nmod, drop = FALSE]
  
  #print("Score modules")
  mod.score <- setNames(rep(0, length(gene.mod)), 
                        names(gene.mod))
  not.selected <- 1:length(gene.mod)
  
  for(m in method) {
    tmp <- score_gene_mod(norm.dat, 
                          select.cells, 
                          gene.mod = gene.mod[not.selected], 
                          eigen = eigen[select.cells, not.selected, drop = FALSE], 
                          method = m, 
                          de.param = de.param,
                          max.cl.size = max.cl.size)
    
    x <- do.call("cbind", 
                 sapply(tmp, 
                        function(x) {
                          x[[1]]
                        },
                        simplify = FALSE))
    
    tmp <- x["sc",] > mod.score[not.selected]
    mod.score[not.selected[tmp]] <- x["sc",tmp]
    tmp <- x["sc",] > min.deScore 
    not.selected <- not.selected[!tmp]
  }
  
  ord <- order(mod.score,
               decreasing = TRUE)
  
  gene.mod <- gene.mod[ord]
  mod.score <- mod.score[ord]
  eigen <- eigen[, ord, drop = FALSE]
  select.mod <- head(which(mod.score > min.deScore), 
                     max.mod)
  
  if(length(select.mod) == 0) {
    select.mod <- head(which(mod.score > min.deScore / 2), 2)
  }
  
  if(length(select.mod) > 0) {
    gene.mod <- gene.mod[select.mod]
    mod.score <- mod.score[select.mod]
    eigen <- eigen[, select.mod, drop = FALSE]
    
    return(list(gene.mod = gene.mod,
                eigen = eigen, 
                gene.mod.val = mod.score))
  } else {
    return(NULL)
  }
}



#' Compute module eigengenes 
#' 
#' @param gene.mod A list of gene modules. 
#' @param norm.dat log transformed normalized data matrix. 
#' @param select.cells Cells used to compute module eigen genes. 
#' @param prefix Default NULL. If not NULL, a heatmap of the gene module eigen genes will be produced with "prefix" as the prefix for the pdf file. 
#' @param method Default "ward.D". Used by hclust method to create the cell dendrogram for the heatmap display.
#' @param hc Precomputed cell dendrogram for heatmap display. Default NULL.  
#' @param ... Other plotting parameters passed to the heatmap function. 
#' 
#' @return A list with two elements: module eigen genes, and if prefix is not NULL, dendrogram for selected cells. 
#' 
get_eigen <- function(gene.mod, 
                      norm.dat, 
                      select.cells = colnames(norm.dat), 
                      prefix = NULL,
                      method = "ward.D",
                      hc = NULL,
                      ...) {
  
  #gene.vector = setNames(rep(names(gene.mod), sapply(gene.mod, length)), unlist(gene.mod))
  #eigen = moduleEigengenes(t(norm.dat[names(gene.vector),select.cells]), gene.vector)[[1]]
  tmp.dat <- as.matrix(norm.dat[unlist(gene.mod), select.cells, drop = FALSE])
  
  eigen <- sapply(gene.mod, 
                  function(x) {
                    tmp.num <- sum(rowSums(tmp.dat[x, select.cells] > 0) > 0)
                    
                    if(tmp.num < 3){
                      return(rep(0, length(select.cells)))
                    }
                    
                    pr.result <- prcomp(t(tmp.dat[x,select.cells]),
                                        tol = 0.8)
                    
                    pc1 <- pr.result$x[,1]
                    rot <- pr.result$rotation[,1, drop = FALSE]
                    
                    if(sum(rot > 0) < length(x) / 2) {
                      pc1 <- -pc1
                    }
                    
                    pc1
                  })
  
  colnames(eigen) <- paste0("ME", names(gene.mod))
  rownames(eigen) <- select.cells
  
  kME <- cor(t(as.matrix(norm.dat[unlist(gene.mod), select.cells])), 
             eigen)
  
  hub <- sapply(names(gene.mod), 
                function(i) {
                  x <- gene.mod[[i]]
                  x[which.max(kME[x, paste0("ME",i)])]
                })
  
  hub <- setNames(hub, paste0("ME", names(gene.mod)))
  
  colnames(eigen) <- paste(colnames(eigen), hub[colnames(eigen)])
  rownames(eigen) <- select.cells
  
  if(!is.null(prefix) & ncol(eigen) > 1){
    if(is.null(hc)) {
      hc <- hclust(dist(eigen), 
                   method = method)
    }
    
    colv <- as.dendrogram(hc)
    
    pdf(paste0(prefix, ".eigen.pdf"),
        height = 6,
        width = 7)
    heatmap.3(t(eigen), 
              Colv = colv, 
              col = blue.red(100), 
              trace = "none", 
              dendrogram = "column",
              ...)
    dev.off()
  }
  
  return(list(eigen = eigen,
              hc = hc))
}

#' Compute reduced dimensions using WGCNA
#' 
#' @param norm.dat
#' @param select.genes
#' @param select.cells
#' @param sampled.cells
#' @param min.module.size
#' @param cutHeight
#' @param type Character, type parameter passed to \code{WGCNA::adjacency()}. Must be one of "unsigned", "signed", "signed hybrid", or "distance". Default is "unsigned".
#' @param soft.power Numeric, power parameter passed to \code{WGCNA::adjacency()}
#' @param rm.gene.mod
#' @param rm.eigen
#' @param ... Additional parameters passed to \code{filter_gene_mod()}
#' 
#' @return a list object with reduced dimensions (rd.dat) and gene modules (gm)
#' @export
#' 
rd_WGCNA <- function(norm.dat, 
                     select.genes, 
                     select.cells, 
                     sampled.cells = select.cells,
                     min.module.size = 10, 
                     cutHeight = 0.99,
                     type = "unsigned",
                     soft.power = 4,
                     rm.gene.mod = NULL,
                     rm.eigen = NULL,
                     ...) {
  
  type <- match.arg(type,
                    choices = c("unsigned", "signed", "signed hybrid", "distance"))
  
  dat <- as.matrix(norm.dat[select.genes, sampled.cells])
  
  adj <- WGCNA::adjacency(t(dat), 
                          power = soft.power,
                          type = type)
  adj[is.na(adj)] <- 0
  
  TOM <- WGCNA::TOMsimilarity(adj,
                              TOMType = type,
                              verbose = 0)
  
  dissTOM <- as.matrix(1 - TOM)
  rownames(dissTOM) <- rownames(dat)
  colnames(dissTOM) <- rownames(dat)
  
  rm(dat)
  gc()
  
  geneTree <- hclust(as.dist(dissTOM), 
                     method = "average")
  
  dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, 
                                               distM = dissTOM, 
                                               cutHeight = cutHeight,
                                               deepSplit = 2, 
                                               pamRespectsDendro = FALSE,
                                               minClusterSize = min.module.size)
  
  gene.mod <- split(row.names(dissTOM), dynamicMods)
  gene.mod <- gene.mod[setdiff(names(gene.mod), "0")]
  
  if (!is.null(rm.gene.mod)) {                                              
    rm.eigen <- get_eigen(rm.gene.mod, 
                          norm.dat, 
                          select.cells)[[1]]
  } else {
    rm.eigen <- NULL
  }
  
  if(is.null(gene.mod) | length(gene.mod) == 0) {
    return(NULL)
  }
  
  gm <- filter_gene_mod(norm.dat, 
                        select.cells, 
                        gene.mod, 
                        min.module.size = min.module.size, 
                        rm.eigen = rm.eigen,
                        ...)
  if(is.null(gm)) {
    return(NULL)
  }
  
  rd.dat <- gm$eigen
  
  return(list(rd.dat = rd.dat, 
              gm = gm))
}
