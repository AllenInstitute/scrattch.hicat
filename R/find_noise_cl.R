get_de_score <- function(de.df, top.genes, upper=20)
  {
    gene.score = -log10(de.df[top.genes,"padj"])
    gene.score[gene.score > upper] = upper
    gene.score = sum(gene.score)
    return(gene.score)
  }

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


plot_doublet <- function(norm.dat, cl, nn.df, de.genes, all.col)
  {
    for(i in 1:nrow(nn.df)){
        tmp.cl = droplevels(cl[cl %in% as.character(c(row.names(nn.df)[i],nn.df[i,"nn.cl"],nn.df[i,"max.olap.cl"]))])
        tmp=display_cl(tmp.cl, norm.dat, prefix=paste(levels(tmp.cl), collapse="_"), col=all.col, max.cl.size=100, de.genes=de.genes)
      }
  }


plot_cl_low <- function(norm.dat, cl, low.df, de.genes, all.col)
  {
    for(i in 1:nrow(low.df)){
        tmp.cl = droplevels(cl[cl %in% as.character(unlist(low.df[i,1:2]))])
        tmp=display_cl(tmp.cl, norm.dat, prefix=paste(levels(tmp.cl), collapse="_"), col=all.col, max.cl.size=100, de.genes=de.genes)
      }
  }

