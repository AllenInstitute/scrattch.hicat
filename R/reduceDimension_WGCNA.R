score_gene_mod <-  function(norm.dat, select.cells, gene.mod, eigen=NULL,method="average",max.cl.size=NULL,de.param=de_param()){
    if(length(gene.mod)==0){
      return(NULL)
    }
    if(is.null(eigen)){
      eigen = get_eigen(gene.mod, norm.dat[,names(select.cells)])
    }
    colnames(eigen) = names(gene.mod)
    gene.mod.val=sapply(names(gene.mod), function(y){
      x= gene.mod[[y]]
      if(is.null(max.cl.size)){
        tmp.dat = norm.dat[x,select.cells]
      }
      else{
        v = eigen[select.cells,y]
        ord = order(v)
        tmp.cells = unique(select.cells[c(head(ord, max.cl.size),tail(ord, max.cl.size))])
        tmp.dat = norm.dat[x,tmp.cells]
      }
      tmp.dat = as.matrix(tmp.dat)
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      if(method=="average"){
        tmp.cl = cutree(hclust(dist(t(tmp.dat)),method="average"),2)
      }
      else if(method=="ward"){
        tmp.cl = cutree(hclust(dist(t(tmp.dat)),method="ward"),2)
      }
      else if(method=="kmeans"){
        tmp.cl = kmeans(t(tmp.dat), 2)$cluster
      }
      else if(method=="louvain"){
        tmp = jaccard_louvain(t(tmp.dat), 15)
        if(is.null(tmp)){
          return(list(c(0,0),NULL))
        }
        tmp.cl = tmp$cl
      }
      else if(method=="mclust"){
        require(mclust)
        tmp = Mclust(t(tmp.dat), G=2)
        if(!is.na(tmp[[1]])){
          tmp.cl = tmp$classification
        }
        else{
          return(list(c(0,0),NULL))
        }
      }
      
      tmp = de_score(as.matrix(norm.dat[,names(tmp.cl)]), cl=tmp.cl, select.cells=names(tmp.cl), ...)
      de.genes = tmp[[2]]
      de.genes = de.genes[sapply(de.genes, length)>1]
      if(length(de.genes) > 0){
        sc =max(sapply(de.genes, function(x)x$score))
        gene.num = max(sapply(de.genes, function(x)x$num))
      }
      else{
        sc=0
        gene.num=0
      }
      return(list(c(sc=sc, gene.num = gene.num),de.genes=de.genes))
    },simplify=F)
    return(gene.mod.val)
  }



filter_gene_mod <- function(norm.dat, select.cells, gene.mod, minModuleSize=10,min.cells=10, min.padj=40, de.param = de_param(), max.cl.size=NULL,rm.eigen=NULL, rm.th = 0.6, maxSize=200, prefix="cl", max.mod=NULL)
  {
    eigen = get_eigen(gene.mod, norm.dat,select.cells)[[1]]
    if(!is.null(rm.eigen)){
      rm.cor=cor(eigen, rm.eigen[select.cells,])
      rm.cor[is.na(rm.cor)]=0
      rm.score = setNames(rowMaxs(abs(rm.cor)), colnames(eigen))
      print(rm.score)
      select1 = rm.score < rm.th
      ###Check of overlapping genes in rm.gene.mod
      select2 = sapply(gene.mod, function(x){
        all(sapply(rm.gene.mod, function(y){
          m=length(intersect(x,y))
          max(m/length(x),m/length(y)) < 0.5
        }))
      })
      #select = select1 & select2
      select = select1 
      if(sum(!select)){
        print("Remove module")
        print(rm.score[!select,drop=F])
      }
      eigen= eigen[,select,drop=F]
      gene.mod = gene.mod[select]
    }
    if(length(select.cells) > 4000){
      method="louvain"
    }
    else{
      method=c("ward","kmeans")
    }
    nmod = min(20, length(gene.mod))
    if(nmod==0){
      return(NULL)
    }
    gene.mod = head(gene.mod, nmod)
    eigen = eigen[,1:nmod,drop=F]
    print("Score modules")
    mod.score = setNames(rep(0, length(gene.mod)), names(gene.mod))
    not.selected=1:length(gene.mod)
    for(m in method){
      tmp=score_gene_mod(norm.dat, select.cells, gene.mod=gene.mod[not.selected],eigen = eigen[select.cells,not.selected,drop=F], min.cells=min.cells, method=method, de.param=de.param,max.cl.size=max.cl.size)
      x = do.call("cbind", sapply(tmp, function(x)x[[1]],simplify=F))
      tmp= x["sc",] > min.padj 
      mod.score[not.selected[tmp]]= x["sc",tmp]
      not.selected = not.selected[!tmp]
      print(mod.score[mod.score>0,drop=F])
    }
    select.mod = mod.score > 0
    if(sum(select.mod)>0){
      gene.mod = gene.mod[select.mod]
      eigen=eigen[,select.mod,drop=F]
      mod.score=mod.score[select.mod]
      ord = order(mod.score,decreasing=T)
      if(!is.null(max.mod)){
        ord = head(ord, max.mod)
      }
      gene.mod = gene.mod[ord]
      mod.score = mod.score[ord,drop=F]
      eigen = eigen[,ord,drop=F]
     
      for(x in 1:length(gene.mod)){
        g = gene.mod[[x]]
        if(length(g) > maxSize){
          kME=cor(as.matrix(t(norm.dat[g,select.cells])), eigen[select.cells,x])
          g = head(g[order(abs(kME),decreasing=T)], maxSize)
          gene.mod[[x]] <- g
        }
      }
      return(list(gene.mod=gene.mod,eigen=eigen, gene.mod.val=mod.score))
    }
    return(NULL)
  }


    
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param gene.mod 
##' @param norm.dat 
##' @param select.cells 
##' @param prefix 
##' @param method 
##' @param hc 
##' @param ... 
##' @return 
##' @author 
get_eigen <- function(gene.mod, norm.dat, select.cells, prefix=NULL,method="ward",hc=NULL,...)
  {
    #gene.vector = setNames(rep(names(gene.mod), sapply(gene.mod, length)), unlist(gene.mod))
    #eigen = moduleEigengenes(t(norm.dat[names(gene.vector),select.cells]), gene.vector)[[1]]
    tmp.dat= as.matrix(norm.dat[unlist(gene.mod),select.cells, drop=F])
    eigen = sapply(gene.mod, function(x){
      tmp.num=sum(rowSums(tmp.dat[x, select.cells] > 0) > 0)
      if(tmp.num < 3){
        return(rep(0, length(select.cells)))
      }
      pr.result = prcomp(t(tmp.dat[x,select.cells]),tol=0.8)
      pc1=pr.result$x[,1]
      rot  =  pr.result$rotatio[,1,drop=F]
      if(sum(rot>0) < length(x)/2){
        pc1 = -pc1
      }
      pc1
    })
    colnames(eigen) = paste0("ME", names(gene.mod))
    row.names(eigen)= select.cells
    kME=cor(t(as.matrix(norm.dat[unlist(gene.mod),select.cells])), eigen)
    hub=sapply(names(gene.mod), function(i){
      x = gene.mod[[i]]
      x[which.max(kME[x, paste0("ME",i)])]
    })
    hub = setNames(hub, paste0("ME",names(gene.mod)))
    colnames(eigen)=paste(colnames(eigen),hub[colnames(eigen)])
    row.names(eigen)=select.cells
    if(!is.null(prefix) & ncol(eigen)>1){
      if(is.null(hc)){
        hc = hclust(dist(eigen), method=method)
      }
      colv = as.dendrogram(hc)
      pdf(paste0(prefix, ".eigen.pdf"),height=6,width=7)
      heatmap.3(t(eigen), Colv= colv, col=blue.red(100), trace="none", dendrogram="column",...)
      dev.off()
    }
    return(list(eigen,hc))
  }



####Extraced from WGCNA to reduce package dependency
adjacency <- function (datExpr, selectCols = NULL, type = "unsigned", power = if (type == "distance") 1 else 6, corFnc = "cor", corOptions = "use = 'p'", distFnc = "dist", distOptions = "method = 'euclidean'") 
{
  intType = charmatch(type, .adjacencyTypes)
  if (is.na(intType)) 
    stop(paste("Unrecognized 'type'. Recognized values are", 
               paste(.adjacencyTypes, collapse = ", ")))
  if (intType < 4) {
    if (is.null(selectCols)) {
      corExpr = parse(text = paste(corFnc, "(datExpr ", 
                        prepComma(corOptions), ")"))
      cor_mat = eval(corExpr)
    }
    else {
      corExpr = parse(text = paste(corFnc, "(datExpr, datExpr[, selectCols] ", 
                        prepComma(corOptions), ")"))
      cor_mat = eval(corExpr)
    }
  }
  else {
    if (!is.null(selectCols)) 
      stop("The argument 'selectCols' cannot be used for distance adjacency.")
        corExpr = parse(text = paste(distFnc, "(t(datExpr) ", 
                          prepComma(distOptions), ")"))
    d = eval(corExpr)
    if (any(d < 0)) 
      warning("Function WGCNA::adjacency: Distance function returned (some) negative values.")
    cor_mat = 1 - as.matrix((d/max(d, na.rm = TRUE))^2)
  }
  if (intType == 1) {
    cor_mat = abs(cor_mat)
  }
  else if (intType == 2) {
    cor_mat = (1 + cor_mat)/2
  }
  else if (intType == 3) {
    cor_mat[cor_mat < 0] = 0
  }
  cor_mat^power
}

MyTOM <- function(adj)
{
  require(Matrix)
  adj[adj < sparse.cutoff] = 0
  adj = Matrix(adj, sparse=T)
  print("matrix multiplication")
  olap= crossprod(adj)
  print("normalize")
  adj.counts = rowSums(adj)-1
  nsize = nrow(adj)
  tmp1 = matrix(rep(adj.counts, ncol(olap)), nrow=length(adj.counts))
  N = pmin(tmp1, t(tmp1))
  TOM = (olap - adj)/(N + 1 - adj)
  diag(TOM)=1
  TOM
}


rd_WGCNA <- function(norm.dat, select.genes, select.cells, sampled.cells=select.cells,minModuleSize=10, cutHeight=0.99,type="unsigned",softPower=4,rm.gene.mod=NULL,rm.eigen=NULL,...)
  {
    dat =norm.dat[select.genes,sampled.cells]
    adj=adjacency(t(dat), power = softPower,type=type)
    adj[is.na(adj)]=0
    TOM = TOMsimilarity(adj,TOMType=type)
    dissTOM = as.matrix(1-TOM)
    row.names(dissTOM)= colnames(dissTOM) = row.names(dat)
    rm(dat)
    gc()
    geneTree = flashClust(as.dist(dissTOM), method = "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight=cutHeight,
      deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
    gene.mod = split(row.names(dissTOM), dynamicMods)
    gene.mod = gene.mod[setdiff(names(gene.mod),"0")]
    if (!is.null(rm.gene.mod)) {                                              
      rm.eigen = get_eigen(rm.gene.mod, norm.dat, select.cells)[[1]]
    }
    else{
      rm.eigen=NULL
    }
    if(is.null(gene.mod)|length(gene.mod)==0){return(NULL)}
    gm= filter_gene_mod(norm.dat, select.cells, gene.mod, minModuleSize=minModuleSize, rm.eigen=rm.eigen,...)
    if(is.null(gm)){
      return(NULL)
    }
    rd.dat = gm$eigen    
    return(rd.dat)
  }
