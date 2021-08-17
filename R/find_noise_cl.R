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
      up.genes = head(names(de$up.genes), 50)
      nn.df$up.genes[i] = length(up.genes)
      nn.df$up.genes.score[i] = sum(pmin(-log10(de$up.genes[up.genes]),20))
      
      up.genes.olap =sapply( setdiff(cl.good, c(cl1,cl2)), function(k){
        p = paste0(k,"_",cl2)
        tmp.de = get_de_pair(de.genes, k, cl2)
        olap.genes= intersect(names(tmp.de$up.genes), up.genes)
        olap.num = length(olap.genes)
        olap.score= sum(pmin(-log10(de$up.genes[olap.genes]), 20))
        c(olap.num, olap.score)
      })
      cl3 = names(which.max(up.genes.olap[1,]))
      nn.df$max.olap.cl[i] = cl3
      nn.df$max.olap.genes1[i] = up.genes.olap[1,cl3]
      nn.df$max.olap.score1[i] = up.genes.olap[2,cl3]
      nn.df$max.olap.ratio1[i] = up.genes.olap[2,cl3]/nn.df$up.genes.score[i]
      
      de1 = get_de_pair(de.genes, cl1, cl3)
      up.genes = head(names(de1$up.genes), 50)
      up.gene.score= sum(pmin(-log10(de1$up.genes[up.genes]),20))
      de2 = get_de_pair(de.genes, cl2, cl3)
      olap.genes = intersect(names(de2$up.genes), up.genes)  
      nn.df$max.olap.genes2[i] = length(olap.genes)
      nn.df$max.olap.score2[i] = sum(pmin(-log10(de1$up.genes[olap.genes]),20))
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
#' @param doublet.df 
#' @param de.genes 
#' @param all.col 
#'
#' @return
#' @export
#'
#' @examples
plot_doublet <- function(norm.dat, cl, doublet.df, de.genes, all.col)
  {
    for(i in 1:nrow(doublet.df)){                                  
      x = as.character(doublet.df[i, "cl"])
      y = as.character(doublet.df[i, "cl1"])
      z = as.character(doublet.df[i, "cl2"])
      tmp.cl = cl[cl %in% c(x, y, z)]
      tmp.cl = setNames(factor(as.character(tmp.cl), c(x,y,z)), names(tmp.cl))
      tmp=display_cl(tmp.cl, norm.dat, prefix=paste0("doublet.",paste(levels(tmp.cl), collapse="_")), col=all.col, max.cl.size=100, de.genes=de.genes)
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






find_doublet_all <- function(de.genes, mc.cores=5, min.genes=100)
  {    
    require(parallel)
    if(is.null(de.genes)){
      stop("Need to specify de.genes")
    }
    pairs.df = get_pairs(names(de.genes))
        
    cl = union(pairs.df[,1],pairs.df[,2])
    result.list=  parallel::pvec(sample(row.names(pairs.df)), function(pairs){
      result.list= sapply(pairs, function(p){
        de = de.genes[[p]]
        if(length(de)==0){
          return(NULL)
        }
        if(de$up.num < min.genes | de$down.num < min.genes){
          return(NULL)
        }
        cl1 = pairs.df[p,1]
        cl2 = pairs.df[p,2]
        
        up.genes.score = head(de$up.genes, 50)
        down.genes.score = head(de$down.genes,50)
        up.genes.score[up.genes.score > 20] = 20
        down.genes.score[down.genes.score > 20] = 20
        up.genes = names(up.genes.score)
        down.genes = names(down.genes.score)
        up.genes.score=sum(up.genes.score)
        down.genes.score = sum(down.genes.score)
        
        results = sapply(setdiff(as.character(cl),c(cl1,cl2)), function(cl3){
          tmp1.de = get_de_pair(de.genes, cl1, cl3)
          tmp2.de = get_de_pair(de.genes, cl3, cl2)
          if(tmp1.de$num < 20 | tmp2.de$num < 20 ){
            return(NULL)
          }
          olap.up.genes1 = intersect(names(tmp2.de$up.genes), up.genes)
          olap.up.num1 = length(olap.up.genes1)
          olap.up.score1 = get_de_truncate_score_sum(de$up.genes[olap.up.genes1])
          
          olap.up.ratio1 = olap.up.score1 / up.genes.score
          
          olap.down.genes1 = intersect(names(tmp1.de$down.genes), down.genes)
          olap.down.num1 = length(olap.down.genes1)
          
          olap.down.score1 = get_de_truncate_score_sum(de$down.genes[olap.down.genes1])
          olap.down.ratio1 = olap.down.score1 / down.genes.score
          
          up.genes2 = head(names(tmp1.de$up.genes), 50)
          up.genes.score2 = get_de_truncate_score_sum(tmp1.de$up.genes[up.genes2])          
          olap.up.genes2 = intersect(names(up.genes2),de$up.genes)
          olap.up.num2 = length(olap.up.genes2)
          olap.up.score2 = get_de_truncate_score_sum(tmp1.de$up.genes[olap.up.genes2])
          olap.up.ratio2 = olap.up.score2 /up.genes.score2
          
          
          down.genes2 = head(names(tmp2.de$down.genes), 50)
          down.genes.score2 = get_de_truncate_score_sum(tmp2.de$down.genes[down.genes2])
          olap.down.genes2 = intersect(down.genes2,names(de$down.genes))
          olap.down.num2 = length(olap.down.genes2)
          olap.down.score2 = get_de_truncate_score_sum(tmp2.de$down.genes[olap.down.genes2])
          olap.down.ratio2 = olap.down.score2 /down.genes.score2
          
          result = list(
            cl1=cl1,
            cl2=cl2,
            up.num = length(up.genes),
            down.num = length(down.genes),
            olap.num=c(olap.up.num1, olap.down.num1, olap.up.num2, olap.down.num2),
            olap.ratio = c(olap.up.ratio1, olap.down.ratio1, olap.up.ratio2, olap.down.ratio2),
            olap.score = c(olap.up.score1, olap.down.score1, olap.up.score2, olap.down.score2)
            )
          result$score = sum(result$olap.score) / sum(c(up.genes.score, down.genes.score, up.genes.score2, down.genes.score2))
          return(result)
        },simplify=F)
        if(is.null(results) | length(results)==0){
          return(NULL)
        }
        test.score=unlist(sapply(results, function(x)x$score))
        tmp = names(which.max(test.score))
        result = results[[tmp]]
        result$cl = tmp
        return(result)                                        
      },simplify=F)
    },mc.cores=mc.cores)

    result.list = result.list[!sapply(result.list, is.null)]
    if(is.null(result.list)| length(result.list)==0){
      return(NULL)
    }
    cl = sapply(result.list, function(x)x$cl)
    cl1 = sapply(result.list, function(x)x$cl1)
    cl2 = sapply(result.list, function(x)x$cl2)
    
    up.num = sapply(result.list, function(x)x$up.num)
    down.num = sapply(result.list, function(x)x$down.num)
    olap.num.df = t(sapply(result.list, function(x)x$olap.num))
    colnames(olap.num.df) = paste0("olap.num.",c("up.1","down.1","up.2","down.2"))
              
    olap.ratio.df = t(sapply(result.list, function(x)x$olap.ratio))
    colnames(olap.ratio.df) = paste0("olap.ratio.",c("up.1","down.1","up.2","down.2"))
    
    score = sapply(result.list, function(x)x$score)
    result.df = data.frame(cl=cl, cl1=cl1, cl2=cl2, score=score, up.num=up.num, down.num=down.num,  as.data.frame(olap.num.df), as.data.frame(olap.ratio.df))
    return(result.df)
  }



find_low_quality_all <- function(de.genes=NULL, de.score.mat=NULL,low.th = 2)
  {
    library(dplyr)
    if(is.null(de.score.mat)){
      de.score.mat <- get_de_matrix(de.genes, 
                                    directed = TRUE, 
                                    field = "num")
    }
    de.score.mat.df = as.data.frame(as.table(de.score.mat))
    de.score.mat.df = de.score.mat.df %>% filter(Freq < low.th & Var1!=Var2)
    colnames(de.score.mat.df) = c("cl.low","cl", "up.genes")
    return(de.score.mat.df)
  }

