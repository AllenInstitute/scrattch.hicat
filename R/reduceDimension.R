rdPCA <- function(norm.dat, select.genes, select.cells,sampled.cells=select.cells, max.pca=10, w=NULL)
{
  if(is.null(w)){
    pca = prcomp(t(as.matrix(norm.dat[select.genes,sampled.cells])),tol=0.01)
    pca.importance = summary(pca)$importance
    v = pca.importance[2,]
  }
  else{
    dat <- scale(norm.dat[select.genes, select.cells], center = TRUE, scale = TRUE)
    w = w[select.genes, select.cells]
    
    for (x in 1:nsamp) {
      for (y in 1:nsamp) {
        wt1 <- w[, x] * w[, y]
        cov1 <- cov.wt(e[, c(x, y)], wt = wt1, center = FALSE)$cov
        e.cov[x, y] <- cov1[1, 2]
      }
    }
        
    nsamp <- ncol(e)
    wt <- crossprod(w)
    cov = cov.wt(dat, wt, center=FALSE)$cov
    eig1 <- eigen(cov, symmetric = TRUE)
    eig.val <- eig1$values
    eig.val[is.na(eig.val) | eig.val < 0] <- 0
    eig.vec <- eig1$vectors
    dimnames(eig.vec) <- list(colnames(e), paste0("PC", 1:ncol(eig.vec)))
    pca1 <- list()
    pca1$sdev <- sqrt(eig.val)
    pca1$rotation <- eig.vec
    pca1$x <- e %*% eig.vec
  }

  select= which((v - mean(v))/sd(v)>2) 
  tmp = head(select,max.pca)
  if(length(sampled.cells)< length(select.cells)){
    rot  =  pca$rotatio[,tmp]
    rd.dat = as.matrix(t(norm.dat[row.names(rot),select.cells])  %*% rot)
  }
  else{
    rd.dat=pca$x[,tmp,drop=F]
  }
  return(rd.dat)
}
  
rdWGCNA <- function(norm.dat, select.genes, select.cells, sampled.cells=select.cells,minModuleSize=10, cutHeight=0.99,type="unsigned",softPower=4,rm.gene.mod=NULL,rm.eigen=NULL,...)
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
      rm.eigen = getEigen(rm.gene.mod, norm.dat, select.cells)[[1]]
    }
    else{
      rm.eigen=NULL
    }
    if(is.null(gene.mod)|length(gene.mod)==0){return(NULL)}
    gm= filterGeneMod(norm.dat, select.cells, gene.mod, minModuleSize=minModuleSize, rm.eigen=rm.eigen,...)
    if(is.null(gm)){
      return(NULL)
    }
    rd.dat = gm$eigen    
    return(rd.dat)
  }


scoreGeneMod <-  function(norm.dat, select.cells, gene.mod, eigen=NULL,method="average",max.cl.size=NULL,...){
    if(length(gene.mod)==0){
      return(NULL)
    }
    if(is.null(eigen)){
      eigen = getEigen(gene.mod, norm.dat[,names(select.cells)])
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
      
      tmp = deScore(as.matrix(norm.dat[,names(tmp.cl)]), cl=tmp.cl, select.cells=names(tmp.cl), ...)
      gc()
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



filterGeneMod <- function(norm.dat, select.cells, gene.mod, minModuleSize=10,min.cells=10, padj.th= 0.01, lfc.th = 1, min.padj=40, q1.th=NULL, q2.th=NULL, q.diff.th=NULL, max.cl.size=NULL,rm.eigen=NULL, rm.th = 0.6, maxSize=200, prefix="cl", max.mod=NULL,...)
  {
    eigen = getEigen(gene.mod, norm.dat,select.cells)[[1]]
    if(!is.null(rm.eigen)){
      rm.cor=cor(eigen, rm.eigen[select.cells,])
      rm.cor[is.na(rm.cor)]=0
      rm.score = setNames(rowMaxs(abs(rm.cor)), colnames(eigen))
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
      tmp=scoreGeneMod(norm.dat, select.cells, gene.mod=gene.mod[not.selected],eigen = eigen[select.cells,not.selected,drop=F], min.cells=min.cells, padj.th=padj.th, lfc.th=lfc.th,method=m, q1.th=q1.th, q2.th=q2.th, q.diff.th=q.diff.th,max.cl.size=max.cl.size)
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



plotGeneMod <- function(norm.dat, select.cells,gene.mod,prefix="cl",hc=NULL,method="ward",...)
  { 
    pdf(paste(prefix,"gene.mod.pdf",sep="."), height=9, width=9)
    for(x in 1:length(gene.mod)){
      g = gene.mod[[x]]
      tmp.dat = norm.dat[g,select.cells]
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      if(is.null(hc)){
        hc = hclust(dist(t(tmp.dat[g,])), method=method)
      }
      gene.hc = hclust(dist(tmp.dat[g,]), method=method)
      cexCol = min(70/ncol(tmp.dat),1)
      cexRow = min(60/length(g),1)
      heatmap.3(tmp.dat[g,],  Colv=as.dendrogram(hc), Rowv=as.dendrogram(gene.hc),col=blue.red(100), trace="none",  dendrogram="column", cexCol=cexCol,cexRow=cexRow,...)
    }
    dev.off()
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
getEigen <- function(gene.mod, norm.dat, select.cells, prefix=NULL,method="ward",hc=NULL,...)
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

