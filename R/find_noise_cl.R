#' Get de score
#'
#' @param de.df 
#' @param top.genes 
#' @param upper 
#'
#' @return
#' @export
#'
#' @examples
get_de_score <- function(de.df, top.genes, upper=20)
  {
    gene.score = -log10(de.df[top.genes,"padj"])
    gene.score[gene.score > upper] = upper
    gene.score = sum(gene.score)
    return(gene.score)
  }

#' Get de genes sym
#'
#' @param de.genes 
#'
#' @return
#' @export
#'
#' @examples
get_de_genes_sym <- function(de.genes)
  {
    de.genes.sym = de.genes
    pairs= names(de.genes)
    pairs = as.data.frame(do.call("rbind",strsplit(pairs,"_")))
    row.names(pairs) = names(de.genes)
    pairs$tmp = paste0(pairs[,2],"_",pairs[,1])
    de.genes.sym[pairs$tmp] = de.genes[row.names(pairs)]
    for(i in 1:nrow(pairs)){
      de.genes.sym[[pairs[i,"tmp"]]]$up.genes = de.genes[[row.names(pairs)[i]]]$down.genes
      de.genes.sym[[pairs[i,"tmp"]]]$down.genes = de.genes[[row.names(pairs)[i]]]$up.genes 
    }
    return(de.genes.sym)
  }


#' Get de pair
#'
#' @param de.genes 
#' @param cl1 
#' @param cl2 
#'
#' @return
#' @export
#'
#' @examples
get_de_pair<- function(de.genes, cl1, cl2)
  {
    pair = paste0(cl1, "_", cl2)
    if(pair %in% names(de.genes)){
      return(de.genes[[pair]])
    }
    else{
      pair = paste0(cl2, "_", cl1)
      de = de.genes[[pair]]
      tmp  = de$up.genes
      de$up.genes = de$down.genes
      de$down.genes = tmp
      return(de)
    }
  }


#' Find doublet
#'
#' @param cl.df 
#' @param cl.sim 
#' @param cl.good 
#' @param de.genes 
#'
#' @return
#' @export
#'
#' @examples
find_doublet <- function(cl.df, cl.sim, cl.good, de.genes=NULL)
  {    
    diag(cl.sim)=0
    if(is.null(de.genes)){
      stop("Need to specify de.genes")
    }
    nn = setNames(cl.good[apply(cl.sim[, cl.good],1, function(x){
      which.max(x)
    })],row.names(cl.sim))
    nn.df= data.frame(cl_label = cl.df[names(nn), "cluster_label"],cl.gene.counts = cl.df[names(nn), "gene.counts"],
      nn.cl=nn, nn.cl_label = cl.df[as.character(nn),"cluster_label"],nn.cl.gene.counts = cl.df[as.character(nn), "gene.counts"],stringsAsFactors=FALSE)
    nn.df$pair = paste0(names(nn), "_",nn)
    nn.df$pair.sim = get_pair_matrix(cl.sim, row.names(nn.df), nn.df$nn.cl)
    nn.df$up.genes=0
    nn.df$up.genes.score = 0
    nn.df$max.olap.cl = NA
    nn.df$max.olap.genes1 = 0
    nn.df$max.olap.score1 = 0
    nn.df$max.olap.ratio1 = 0

    nn.df$max.olap.genes2 = 0
    nn.df$max.olap.score2 = 0
    nn.df$max.olap.ratio2 = 0

     
    for(i in 1:nrow(nn.df)){
      cl1= row.names(nn.df)[i]
      cl2=as.character(nn.df$nn.cl[i])
      de = get_de_pair(de.genes, cl1, cl2)
      up.genes = head(de$up.genes, 50)
      nn.df$up.genes[i] = length(up.genes)
      nn.df$up.genes.score[i] = get_de_score(de$de.df, up.genes)
      
      up.genes.olap =sapply( setdiff(cl.good, c(cl1,cl2)), function(k){
        p = paste0(k,"_",cl2)
        tmp.de = get_de_pair(de.genes, k, cl2)
        olap.genes= intersect(tmp.de$up.genes, up.genes)
        olap.num = length(olap.genes)
        olap.score= get_de_score(de$de.df, olap.genes)
        c(olap.num, olap.score)
      })
      cl3 = names(which.max(up.genes.olap[1,]))
      nn.df$max.olap.cl[i] = cl3
      nn.df$max.olap.genes1[i] = up.genes.olap[1,cl3]
      nn.df$max.olap.score1[i] = up.genes.olap[2,cl3]
      nn.df$max.olap.ratio1[i] = up.genes.olap[2,cl3]/nn.df$up.genes.score[i]
      
      de1 = get_de_pair(de.genes, cl1, cl3)
      up.genes = head(de1$up.genes, 50)
      up.gene.score= get_de_score(de1$de.df, up.genes)
      de2 = get_de_pair(de.genes, cl2, cl3)
      olap.genes = intersect(de2$up.genes, up.genes)  
      nn.df$max.olap.genes2[i] = length(olap.genes)
      nn.df$max.olap.score2[i] = get_de_score(de1$de.df, olap.genes)
      nn.df$max.olap.ratio2[i] = nn.df$max.olap.score2[i]/ up.gene.score
    }
    nn.df$max.olap.cl_label = cl.df[as.character(nn.df$max.olap.cl),"cluster_label"]
    nn.df$olap.cl.sim = get_pair_matrix(cl.sim, nn.df$max.olap.cl, nn.df$nn.cl)
    return(nn.df)
  }


#' Plot doublet
#'
#' @param norm.dat 
#' @param cl 
#' @param nn.df 
#' @param de.genes 
#' @param all.col 
#'
#' @return
#' @export
#'
#' @examples
plot_doublet <- function(norm.dat, cl, nn.df, de.genes, all.col)
  {
    for(i in 1:nrow(nn.df)){
        tmp.cl = droplevels(cl[cl %in% as.character(c(row.names(nn.df)[i],nn.df[i,"nn.cl"],nn.df[i,"max.olap.cl"]))])
        tmp=display_cl(tmp.cl, norm.dat, prefix=paste(levels(tmp.cl), collapse="_"), col=all.col, max.cl.size=100, de.genes=de.genes)
      }
  }


#' Plot cl low
#'
#' @param norm.dat 
#' @param cl 
#' @param low.df 
#' @param de.genes 
#' @param all.col 
#'
#' @return
#' @export
#'
#' @examples
plot_cl_low <- function(norm.dat, cl, low.df, de.genes, all.col)
  {
    for(i in 1:nrow(low.df)){
        tmp.cl = droplevels(cl[cl %in% as.character(unlist(low.df[i,1:2]))])
        tmp=display_cl(tmp.cl, norm.dat, prefix=paste(levels(tmp.cl), collapse="_"), col=all.col, max.cl.size=100, de.genes=de.genes)
      }
  }


#' Use DEGenes to identify clusters with few differentially expressed genes
#'
#' To use this, you'll need to perform pairwise gene expression tests with de_score(), and 
#' a set of clusters with high confidence (cl.good).
#' 
#' Low-quality clusters will be from all clusters in cl.df that are not in cl.good.
#'
#' @param cl.df A data.frame with cluster annotations that includes a "cluster_label" column 
#' @param cl.good a cluster factor object for cells from high-quality clusters.
#' @parade.genesm de.score.mat an optional matrix of pairwise degene scores generated by get_de_matrix().
#' @param de.genes differential gene expression results like those generated by de_score().
#'
#' @return a data.frame with the closest "good" cluster for each non-good/low cluster.
#' 
#' @export
#'
find_low_quality_cl <- function(cl.df, 
                                cl.good, 
                                de.score.mat = NULL, 
                                de.genes) {
  
  if(is.null(de.score.mat)){
    de.score.mat <- get_de_matrix(de.genes, 
                                  directed = TRUE, 
                                  field = "num")
  }
  
  diag(de.score.mat) <- max(de.score.mat)
  
  tmp.mat <- de.score.mat[row.names(cl.df), cl.good]
  ###Only search the good clusters that have bigger size than the tested cluster
  not.select = sapply(cl.good, function(x) cl.df[x, "size"] <= cl.df$size)
  tmp.mat[not.select] = max(tmp.mat)+1
  low.pair <- data.frame(cl.low = row.names(tmp.mat), 
                         cl.good = colnames(tmp.mat)[apply(tmp.mat, 1, which.min)],
                        stringsAsFactors = FALSE)
  
  low.pair$low.gene.counts <- cl.df[low.pair[,1], "gene.counts"]
  low.pair$high.gene.counts <- cl.df[low.pair[,2], "gene.counts"]
  
  low.pair$low.size <- cl.df[low.pair[,1], "size"]
  low.pair$high.size <- cl.df[low.pair[,2], "size"]
  
  low.pair$up.genes <- get_pair_matrix(de.score.mat, 
                                       low.pair$cl.low, 
                                       low.pair$cl.good)
  
  low.pair$cl.low.label <- cl.df[as.character(low.pair$cl.low), "cluster_label"]
  low.pair$cl.good.label <- cl.df[as.character(low.pair$cl.good), "cluster_label"]
  
  row.names(low.pair) <- low.pair$cl.low
  
  return(low.pair)
}

#' Plot low-quality clusters
#'
#' @param norm.data normalized data matrix for clustered cells.
#' @param cl A cluster factor object for all cells.
#' @param low.df a data.frame with results from find_low_quality_cl().
#' @param nn.df 
#' @param de.genes differential gene expression results like those generated by de_score().
#' @param all.col 
#'
#' @return
#'
plot_low_qc <- function(norm.dat, 
                        cl, 
                        low.df, 
                        nn.df, 
                        de.genes, 
                        all.col) {
  
  for(i in 1:nrow(low.df)) {
    x <- low.df[i, "cl.low"]
    y <- low.df[i, "cl.good"]
    i <- nn.df[x, "nn.cl"]
    j <- nn.df[y, "nn.cl"]
    
    tmp.cl <- droplevels(cl[cl %in% c(x, y, i, j)])
    tmp.cl <- tmp.cl[names(tmp.cl) %in% colnames(norm.dat)]
    
    tmp <- display_cl(tmp.cl, 
                      norm.dat, 
                      prefix = paste(levels(tmp.cl), 
                                     collapse = "_"), 
                      col = all.col, 
                      max.cl.size = 100, 
                      de.genes = de.genes)
  }  
}
