library(WGCNA)
library(flashClust)
library(limma)
library(matrixStats)
library(Matrix)
library(IRanges)
require(igraph)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))

jaccard <- function(m) {
  require(Matrix)
  ## common values:
  A =  m %*% t(m)
  ## indexes for non-zero common values
  im = Matrix::which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = Matrix::rowSums(m)  
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
    )  
  return( J )
}

###rows are cells, columns are feathers
jaccard_louvain <- function(dat, k=10)
  {
    require(FNN)
    require(igraph)
    knn.result= get.knn(dat, k)
    knn.matrix = knn.result[[1]]
    p = as.vector(t(knn.matrix))
    edge = cbind(rep(1:nrow(knn.matrix),rep(k, nrow(knn.matrix))),p)
    edge.unique = cbind(rowMins(edge), rowMaxs(edge))
    edge.unique=unique(edge.unique)
    knn.gr = graph(t(edge))
    knn.matrix = get.adjacency(knn.gr)
    
    jaccard.adj  = jaccard(knn.matrix)
    jaccard.gr = graph.adjacency(jaccard.adj, mode="undirected",weighted=TRUE)
    louvain.result= cluster_louvain(jaccard.gr)
    mod.sc = modularity(louvain.result)
    if(pass_louvain(mod.sc,jaccard.adj)){
      cl  = setNames(louvain.result$membership,row.names(dat))
      return(list(cl=cl, result=louvain.result))
    }
    else{
      return(NULL)
    }
  }


cluster_one_round <- function(norm.dat, counts=NULL, select.cells, prefix, rm.gene.mod=NULL,  rm.eigen=NULL, rm.th=0.6, col=NULL, vg.padj.th=0.5, min.cells=10,min.padj=100,de.padj.th=0.05, lfc.th=1, max.dim=20,maxGenes=3000,sampleSize=4000,type="undirectional",display=FALSE,method="louvain",dim.method="pca",low.th=1,min.genes=5,max.cl.size=300,use.voom=FALSE,verbose=FALSE,...)
  {
    if(length(select.cells)>sampleSize){
      tmp.cells = sample(select.cells, pmin(length(select.cells),sampleSize))
    }
    else{
      tmp.cells = select.cells
    }
    ###Find high variance genes
    if(is.matrix(norm.dat)){
      select.genes = row.names(norm.dat)[which(rowSums(norm.dat[,select.cells] > low.th) >= min.cells)]
    }
    else{
      select.genes = row.names(norm.dat)[which(Matrix::rowSums(norm.dat[,select.cells] > low.th) >= min.cells)]
    }
    ###Find high variance genes.
    if(is.null(counts)){
      counts = 2^(norm.dat[select.genes, tmp.cells])-1
    }
    plot.fig=NULL
    if(verbose){
      plot.fig=paste0(prefix,".vg.pdf")
    }
    vg = findVG(counts[select.genes,tmp.cells],plot.fig=plot.fig)
    print("Dimension reduction")
    print(dim.method)
    if(dim.method=="WGCNA"){
      ###Ignore vg.padj.th for WGCNA
      select.genes = row.names(vg)[which(vg$loess.padj < 1)]
      select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
      rd.dat = rd_WGCNA(norm.dat, select.genes=select.genes,select.cells=select.cells,sampled.cells=tmp.cells,min.cells=min.cells,min.padj=min.padj,max.mod=max.dim,max.cl.size=max.cl.size,...)
    }
    else{
       ###If most genes are differentially expressed, then use absolute dispersion value
      select.genes = row.names(vg)[which(vg$loess.padj < vg.padj.th| vg$dispersion >3)]
      select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "z"])],maxGenes)
      if(verbose){
        cat("Num high variance genes:",length(select.genes),"\n")
      }
      if(length(select.genes)< min.genes){
        return(NULL)
      }
      rd.dat = rd_PCA(norm.dat,select.genes, select.cells, sampled.cells=tmp.cells, max.pca = max.dim)
    }
    if(is.null(rd.dat)||ncol(rd.dat)==0){
      return(NULL)
    }
    if(!is.null(rm.eigen)){
      rm.cor=cor(rd.dat, rm.eigen[row.names(rd.dat),])
      rm.score = rowMaxs(abs(rm.cor))
      print(rm.score)
      select = rm.score < rm.th
      if(sum(!select)>0){
        print("Remove dimension:")
        print(rm.score[!select])
      }
      if(sum(select)==0){
        return(NULL)
      }
      rd.dat = rd.dat[,select,drop=F]
    }
    print(method)
    max.cl = ncol(rd.dat)*2 + 1
    if(method=="apclust"){
      require(apclust)
      sim = negDistMat(rd.dat, r=2)
      ap.result=apcluster(sim)
      tmp=unlist(sapply(1:length(ap.result), function(i){
        setNames(rep(i, length(ap.result[[i]])),names(ap.result[[i]]))
      },simplify=F))
      tmp.means =sapply(1:length(ap.result), function(i){
        colMeans(rd.dat[ap.result[[i]],,drop=F])
      })
      tmp.hc = hclust(dist(t(tmp.means)), method="average")
      tmp.cl= cutree(tmp.hc, max.cl)
      cl = setNames(ap.cl[tmp], names(tmp))
    }
    else if(method=="louvain"){
      tmp = jaccard_louvain(rd.dat, 15)
      if(is.null(tmp)){
        return(NULL)
      }
      cl = tmp$cl
      if(length(unique(cl))>max.cl){
        tmp.means =do.call("cbind",tapply(names(cl),cl, function(x){
          colMeans(rd.dat[x,,drop=F])
        },simplify=F))
        tmp.hc = hclust(dist(t(tmp.means)), method="average")
        tmp.cl= cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
        cl = setNames(tmp.cl[as.character(cl)], names(cl))
      }
    }
    else if(method=="hclust"){
      hc = hclust(dist(rd.dat),method="ward")
      print("Cluster cells")
      cl = cutree(hc, max.cl)
    }
    else if(method=="kmeans"){
      cl = kmeans(rd.dat, max.cl)$cluster
    }
    print(table(cl))
    merge.result=merge_cl(norm.dat, cl=cl, rd.dat=rd.dat, type=type, min.padj=min.padj,min.genes=min.genes, min.cells=min.cells,max.cl.size=max.cl.size,counts=counts, use.voom=use.voom,...)
    gc()
    if(is.null(merge.result))return(NULL)
    #save(merge.result, file=paste0(prefix, ".merge.rda"))
    sc = merge.result$sc
    #print(sc)
    cl = merge.result$cl
    if(length(unique(cl))>1){
      de.genes = merge.result$de.genes
      if(ncol(rd.dat)>1){
        rd.de = DE_genes_pw(t(rd.dat),cl)
        select = rowMins(sapply(rd.de, function(x)x$padj)) < de.padj.th
        rd.dat = rd.dat[,select,drop=F]
      }
      ##ordering clusters 
      cl.dat = do.call("cbind",as.list(tapply(names(cl),cl, function(x)colMeans(rd.dat[x,,drop=F]) )))
      cl.hc = hclust(dist(t(cl.dat)),method="average")
      cl = setNames(factor(as.character(cl), levels= colnames(cl.dat)[cl.hc$order]), names(cl))
      markers = select_markers(norm.dat, cl, de.df = NULL, de.genes=de.genes,n.markers=10)[[1]]
      
      if(verbose){
        cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
        tmp.col =rbind(col[,select.cells],cl.col)
        displayResult(cl, norm.dat,prefix, col=NULL, max.cl.size=NULL,markers=NULL,low.th=1,de.df=NULL,de.genes=NULL, main="",method="limma",...)
          
        tmp.cells = unlist(tapply(names(cl),cl, function(x){
          if(length(x)>max.cl.size){
            x= sample(x, max.cl.size)
          }
          x
        },simplify=FALSE))
        plot_cl_heatmap(norm.dat, cl[tmp.cells],markers=markers, ColSideColors=tmp.col[,tmp.cells, drop=F], prefix=prefix, by.cl=TRUE)
      }
      levels(cl) = 1:length(levels(cl))
      result=list(cl=cl, markers=markers)
      return(result)
    }
    return(NULL)
  }


recursive_cluster<-function(norm.dat, select.cells,prefix, split.size = 10, result=NULL,method="auto",...)
  {
    print(prefix)
    if(method=="auto"){
      if(length(select.cells)>3000){
        select.method="louvain"
      }
      else{
        select.method="hclust"
      }
    }
    else{
      select.method=method
    }
    if(length(select.cells)<3000){
      if(!is.matrix(norm.dat)){
        norm.dat = as.matrix(norm.dat[,select.cells])
      }
    }
    if(is.null(result)){        
      result=cluster_one_round(norm.dat, select.cells=select.cells, prefix=prefix,method=select.method,...)
      gc()
    }
    if(!is.null(result)){
      #save(result, file=paste0(prefix,".rda"))
      cl = result$cl[select.cells]
      gene.mod = result$gene.mod
      markers=result$markers
      cl = setNames(as.integer(cl),names(cl))
      new.cl =cl
      cl.size = table(cl)
      print(cl.size)
      to.split = names(cl.size)[cl.size >=split.size]
      if(length(to.split)>0){
        n.cl = 1
        for(x in sort(unique(cl))){
          tmp.cells = names(cl)[cl==x]
          if(!x %in% to.split){
            new.cl[tmp.cells]=n.cl
          }
          else{
            tmp.prefix = paste(prefix, x, sep=".")
            tmp.result=recursive_cluster(split.size=split.size, norm.dat=norm.dat, select.cells=tmp.cells, prefix=tmp.prefix,method= method,...)
            gc()
            if(is.null(tmp.result)){
              new.cl[tmp.cells]=n.cl
            }
            else{
              tmp.cl = tmp.result$cl
              if(length(unique(tmp.cl)>1)){
                cat("Expand",tmp.prefix, "\n")
                print(table(tmp.cl))
                new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
                markers=union(markers, tmp.result$markers)
              }
            }
          }
          n.cl = max(new.cl)+1
        }
        cl = new.cl
      }
      result=list(cl=cl, markers=markers)
      #save(result, file=paste0(prefix,".rda"))
      return(result)
    }
    return(NULL)
  }



tsneGroup<- function(prefix, cl, select.cl,cl.col,cl.label,tsne.df=NULL,markers=NULL,cex=1,...)
  {
    if(is.null(tsne.df)){
      select.cells = names(cl)[cl %in% select.cl]
      tmp.cl= droplevels(cl[select.cells])
      if(is.null(markers)){
        de.df = DE.genes.pw(norm.dat[, select.cells],paste0("cl",tmp.cl))
        markers =selectMarkers(norm.dat, tmp.cl, de.df=de.df, ...)$markers
      }
      require(Rtsne)
      tsne.result=Rtsne(t(norm.dat[markers, select.cells]))$Y
      save(tsne.result, file=paste0(prefix, ".tsne.result.rda"))
      tsne.df = data.frame(Lim1=tsne.result[,1], Lim2=tsne.result[,2], cl = droplevels(cl[select.cells]))
      tsne.df$cl_label = factor(cl.label[as.character(tsne.df$cl)], levels=cl.label[levels(tsne.df$cl)])
    }
    library(ggplot2)
    cl.center=do.call("rbind",tapply(1:nrow(tsne.df), tsne.df$cl, function(x){
      c(x=median(tsne.df[x,1]), y= median(tsne.df[x,2]))
    }))
    row.names(cl.center)= cl.label[row.names(cl.center)]
    shape = setNames(1:length(levels(tsne.df$cl)) %% 20 + 1,levels(tsne.df$cl))
    p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=cl_label,shape=cl_label),size=cex)
    
    
    p = p+ scale_color_manual(values=as.vector(cl.col[levels(tsne.df$cl)]))+ scale_shape_manual(values=as.vector(shape[levels(tsne.df$cl)]))
    p = p +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(tsne.df$cl)])),ncol=1)
    p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    for(i in 1:nrow(cl.center)){
      p = p +  annotate("text", label=row.names(cl.center)[i], x=cl.center[i,1], y=cl.center[i,2],size=2,color="black")
      }
    pdf(paste0(prefix, ".tsne.pdf"))
    print(p)
    dev.off()
    return(list(markers, tsne.df,p))
  }


selectMarkersPairGroup <- function(norm.dat,cl, g1,g2,de.genes,top.n=50,max.num=1000,n.markers=20,up.gene.score=NULL, down.gene.score=NULL)
  {
    pairs = do.call("rbind",strsplit(names(de.genes), "_"))
    pairs = gsub("cl", "",pairs)
    row.names(pairs)= names(de.genes)
    up.pairs = row.names(pairs)[pairs[,1] %in% g1 & pairs[,2] %in% g2]
    down.pairs = row.names(pairs)[pairs[,1] %in% g2 & pairs[,2] %in% g1]
    select.pairs = c(up.pairs, down.pairs)
    if(is.null(up.gene.score)){
      tmp=get_gene_score(de.genes,top.n=top.n, max.num=max.num,bin.th=4)
      up.gene.score=tmp$up.gene.score
      down.gene.score=tmp$down.gene.score
      row.names(down.gene.score)= all.genes
    }
    all.genes = row.names(up.gene.score)
    tmp.up.gene.score = cbind(up.gene.score[,up.pairs,drop=F], down.gene.score[,down.pairs,drop=F])
    tmp.down.gene.score = cbind(down.gene.score[,up.pairs,drop=F], up.gene.score[,down.pairs,drop=F])
    
    up.genes = row.names(tmp.up.gene.score)[head(order(rowSums(tmp.up.gene.score)), n.markers)]
    down.genes = row.names(tmp.down.gene.score)[head(order(rowSums(tmp.down.gene.score)), n.markers)]
        
    up.num = colSums(tmp.up.gene.score[up.genes,,drop=F] < max.num)
    down.num = colSums(tmp.down.gene.score[down.genes,,drop=F] < max.num)
    total.num = up.num + down.num
    add.genes = setNames(rep(n.markers, ncol(tmp.up.gene.score)), colnames(tmp.up.gene.score)) - total.num
    add.genes = add.genes[add.genes > 0]
    
    up.genes = up.genes[rowMins(tmp.up.gene.score[up.genes,,drop=F]) < max.num]
    down.genes = down.genes[rowMins(tmp.down.gene.score[down.genes,,drop=F]) < max.num]
    genes = union(up.genes, down.genes)
    if(length(add.genes)>0){
      tmp=selectMarkersPair(norm.dat, add.genes= add.genes,de.genes= de.genes, gene.score=pmin(tmp.up.gene.score, tmp.down.gene.score), rm.genes=c(up.genes, down.genes),top.n=top.n)
      genes=union(genes, unlist(tmp))
    }
    return(genes)
  }


run_tSNE <- function(norm.dat, select.genes, cl, cl.df, tsne.result = NULL,...)
  {
    library(Rtsne)
    if(is.null(tsne.result)){
      tsne.result = Rtsne(t(norm.dat[select.genes,names(cl)]),...)$Y
      row.names(tsne.result)=names(cl)
    }
    tsne.df = as.data.frame(tsne.result[names(cl),])
    tsne.df$cl = cl
    tsne.df$cl_label = factor(cl.df[as.character(tsne.df$cl),"cluster_label"], levels=as.character(cl.df$cluster_label))
    tsne.df$cl_label= droplevels(tsne.df$cl_label)
    colnames(tsne.df)[1:2]=c("Lim1","Lim2")
   
    cl.center=do.call("rbind",tapply(1:nrow(tsne.df), tsne.df$cl, function(x){
      center  = c(median(tsne.df[x,1]), median(tsne.df[x,2]))
      #i = min()
      #center= x[which.max(cell.cl.co.ratio[x, as.character(i)])]
      #c(x=median(tsne.df[center,1]), y= median(tsne.df[center,2]))
    }))
    row.names(cl.center)= cl.df[row.names(cl.center), "cluster_label"]
    cex=0.15
    cl.col = setNames(as.character(cl.df$cluster_color),cl.df$cluster_label)
    shape = setNames(1:length(levels(tsne.df$cl_label)) %% 20 + 1,levels(tsne.df$cl_label))
    p=ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color=cl_label,shape=cl_label),size=cex)
    p = p+ scale_color_manual(values=as.vector(cl.col[levels(tsne.df$cl_label)]))+ scale_shape_manual(values=as.vector(shape[levels(tsne.df$cl_label)]))
    for(i in 1:nrow(cl.center)){
      p = p +  annotate("text", label=row.names(cl.center)[i], x=cl.center[i,1], y=cl.center[i,2],size=2,color="black")
    }
    p = p +  guides(colour = guide_legend(override.aes = list(shape = shape[levels(tsne.df$cl_label)])),ncol=5)
    p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.position="bottom")
    
    return(list(tsne.df=tsne.df, p=p))    
  }


###Merge multiple pairs are the same type if their merging are independent. 
merge_cl<- function(norm.dat, cl, rd.dat, min.padj=100, type="undirectional", min.genes=5, min.cells=4,max.cl.size=300,de.method="limma",de.param = de_param(), 
                    eigen=NULL,rm.gene.mod=NULL,rm.eigen=NULL,rm.th=0.7)
  {
    cl = setNames(as.integer(as.character(cl)), names(cl))
    de.genes=list()
    de.df=list()
    select.cells = names(cl)
    if(is.matrix(norm.dat)){
      cell.gene.counts= colSums(norm.dat[,select.cells]>0)
    }
    else{
      cell.gene.counts= Matrix::colSums(norm.dat[,select.cells]>0)
    }
    cell.weights = cell.gene.counts - min(cell.gene.counts)+1
    ###Remove genes correlate with rm.eigen
    if(!is.null(rm.eigen)){
      tmp.cells=sample(select.cells, min(length(select.cells), 4000))
      rm.kME=cor(t(as.matrix(norm.dat[,tmp.cells])), rm.eigen[tmp.cells,])
      rm.kME[is.na(rm.kME)] =0
      rm.score=setNames(rowMaxs(abs(rm.kME)),row.names(rm.kME))
      select = rm.score < rm.th 
      if (!is.null(eigen)){
        kME=cor(t(tmp.dat), eigen)
        kME[is.na(kME)] =0          
        select = select | rowMaxs(abs(kME)) > rm.score
      }
      norm.dat = norm.dat[select,]
    }
    pairs=NULL
    while(length(unique(cl)) > 1){
      cl.size = table(cl)
      ##Down sample cells for efficiency   
      if(!is.null(max.cl.size)){
        tmp.cells = unlist(tapply(names(cl),cl, function(x){
          if(length(x)>max.cl.size){
            #x= sample(x, max.cl.size)
            x=head(x[order(cell.gene.counts[x],decreasing=TRUE)], min(length(select.cells), max.cl.size))
            #x= sample(x, max.cl.size,prob=cell.weights[x])
          }
          x
        },simplify=FALSE))
      }
      else{
        tmp.cells= names(cl)
      }
      tmp.dat = as.matrix(norm.dat[,tmp.cells])
      ###Merge small clusters with the closest neighbors first.
      cl.rd = t(get.cl.means(t(rd.dat),cl))
      #cl.rd <- do.call("rbind", tapply(names(cl),cl, function(x){
      #    colSums(rd.dat[x,,drop=F])
      #  },simplify=F))
      while(TRUE){
        cl.size = table(cl)
        ##Compute cluster similary on reduced dimension
        if(ncol(cl.rd)>2){
          cl.sim = cor(t(cl.rd))
        }
        else{
          cl.diff=as.matrix(dist(cl.rd/as.vector(cl.size[row.names(cl.rd)])))
          cl.sim = 1 - cl.diff/max(cl.diff)
        }
        cl.small = names(cl.size)[cl.size < min.cells]
        ###Merge small clusters with its closest neighbors.
        if(length(cl.small)>0){
          tmp=as.data.frame(as.table(cl.sim[cl.small,,drop=F]))
          tmp[,1]=as.integer(as.character(tmp[,1]))
          tmp[,2]=as.integer(as.character(tmp[,2]))
          tmp = tmp[tmp[,1]!=tmp[,2],,drop=F]
          closest.pair = which.max(tmp$Freq)
          x = tmp[closest.pair,1]
          y=  tmp[closest.pair,2]
          cat("Merge: ", x,y, "sim:", tmp[closest.pair,3],"\n")
          cl[cl==x]= y
          cl.rd[as.character(y),]= cl.rd[as.character(y),] + cl.rd[as.character(x),]
          cl.rd = cl.rd[row.names(cl.rd)!=x,,drop=F]
        }		
        else{
          break
        }	
      }
      if(length(cl.size)==1){
        return(NULL)
      }
      
      ####if the number of clusters is fewer than 10, compute all pairs. If no clusters are significant, quit
      if(length(cl.size)<10 & length(de.genes)==0){
        tmp.result=de_score(tmp.dat, cl=cl[tmp.cells], min.cells=min.cells, method=de.method, de.param= de.param)
        de.genes=tmp.result$de.genes
        de.df = tmp.result$de.df
        sc = sapply(de.genes, function(x){
          if(length(x)>0){x$score}
          else{0}
        })
        merge.pairs=do.call("rbind",strsplit(names(de.genes), "_"))
        merge.pairs = gsub("cl","", merge.pairs)
        row.names(merge.pairs)=names(de.genes)
        pairs = merge.pairs
      }
      else{
        ####Choose top 2 nearest neighbor based on correlation matrix.
        if(nrow(cl.rd)>2){
          k = 3
          knn.matrix=t(sapply(1:nrow(cl.sim), function(i){colnames(cl.sim)[-i][head(order(cl.sim[i,-i],decreasing=T), k)]}))
          row.names(knn.matrix)=row.names(cl.sim)
          merge.pairs = do.call("rbind",apply(knn.matrix, 2,function(x)data.frame(c1=row.names(knn.matrix),c2=x)))
          merge.pairs[,1] = as.integer(as.character(merge.pairs[,1]))
          merge.pairs[,2] = as.integer(as.character(merge.pairs[,2]))
          tmp1 =pmin(merge.pairs[,1],merge.pairs[,2])
          tmp2 =pmax(merge.pairs[,1],merge.pairs[,2])
          merge.pairs = cbind(tmp1,tmp2)        
        }
        else{
          merge.pairs = matrix(as.integer(row.names(cl.rd)), nrow=1)
        }
        
        row.names(merge.pairs) = paste(merge.pairs[,1],merge.pairs[,2],sep="_")
        merge.pairs= merge.pairs[!duplicated(row.names(merge.pairs)),,drop=F]
        ###Determine the de score for these pairs
        if(nrow(merge.pairs)==0){
          break
        }
        #####Check pairs already known but not yet merged yet.
        new.pairs = setdiff(row.names(merge.pairs),names(de.genes))
        pairs = rbind(pairs, merge.pairs[new.pairs,,drop=F])
        tmp.result=de_score_pairs(tmp.dat, counts = counts, cl=cl[tmp.cells], pairs=merge.pairs[new.pairs,,drop=F], min.cells=min.cells, method=method, use.voom=use.voom, ...)
        tmp.de.df = tmp.result[[1]]
        tmp.de.genes = tmp.result[[2]]
        de.genes[names(tmp.de.genes)] = tmp.de.genes
        de.df[names(tmp.de.df)] = tmp.de.df
        sc = sapply(de.genes[row.names(merge.pairs)], function(x){
          if(length(x)>0){x$score}
          else{0}
        })
      }
      sc = sort(sc)
      print(head(sc,10))      
      to.merge = sapply(names(sc), function(p){
        x = de.genes[[p]]
        if(length(x)==0){
          to.merge = TRUE
        }
        else if(type=="undirectional"){
          to.merge=x$score < min.padj | x$num < min.genes    
        }
        else{
          to.merge=x$up.score < min.padj | x$down.score < min.padj | x$up.num < min.genes | x$down.num < min.genes 
        }
        to.merge
      })
      if(sum(to.merge)==0){
        break
      }
      sc = sc[to.merge]
      to.merge= merge.pairs[names(sc),,drop=FALSE]
      merged =c()
      ###The first pair in to.merge always merge. For the remaining pairs, if both clusters have already enough cells,
      ###or independent of previus merging, then they can be merged directly merged as well, without re-assessing DE genes. 
      for(i in 1:nrow(to.merge)){
        p = to.merge[i,]
        if(i == 1 | sc[i] < min.padj/2  & length(intersect(p, merged))==0){
          cat("Merge ",p[1], p[2], sc[i],"\n")                   
          cl[cl==p[2]] = p[1]
          rm.pairs = row.names(pairs)[pairs[,1]%in% p | pairs[,2]%in% p]
          de.genes = de.genes[setdiff(names(de.genes),rm.pairs)]
        }
        merged = c(merged, p)
      }
      pairs = pairs[names(de.genes),,drop=F]
    }
    return(list(cl=cl, de.genes=de.genes,sc=sc))
  }

pass_louvain <- function(mod.sc, adj.mat)
  {
    p = mean(Matrix::colSums(adj.mat>0)-1)/ nrow(adj.mat)
    print(p)
    n = ncol(adj.mat)
    rand.mod1 <- 0.97 * sqrt((1 - p)/(p*n))
    rand.mod2 <- (1 - 2 / sqrt(n)) * (2 / (p * n))^(2/3)
    rand.mod.max <- max(rand.mod1, rand.mod2, na.rm=TRUE)
    cat("Modularity:",mod.sc, "threshold:",rand.mod.max, "\n")
    return(mod.sc > rand.mod.max)
  }

#####internal heterogenity#######
pca_hetero <- function(norm.dat,select.cells, perm.num=10,vg.padj=0.01,rm.gene.mod= NULL, rm.eigen=NULL, rm.th=0.6)
  {

    vg= findVG(2^norm.dat[,select.cells]-1)
    select.genes = row.names(vg)[vg$loess.padj < vg.padj]
    if(length(select.genes)==0){
      return(0)
    }
    if(length(select.genes)==1){
      return(1)
    }
    pca = prcomp(t(as.matrix(norm.dat[select.genes, select.cells])),tol=0.01)
    if(!is.null(rm.gene.mod)| !is.null(rm.eigen)){
      if(is.null(rm.eigen)){
        rm.eigen = get_eigen(rm.gene.mod, norm.dat, select.cells)[[1]]
      }
      rm.cor=cor(pca$x, rm.eigen)
      rm.score =rowMaxs(abs(rm.cor))
      select.pca = which(rm.score < rm.th)[1]
    }
    pca.var = summary(pca)$importance[2,select.pca]
    ###shuffle expression of every gene.
    perm.pca.var = sapply(1:10, function(x){
      perm.dat = norm.dat[select.genes, select.cells]
      for(i in 1:nrow(perm.dat)){
        perm.dat[i, ] = sample(perm.dat[i, ])
      }
      perm.pca = prcomp(t(as.matrix(perm.dat)),tol=0.5)
      summary(perm.pca)$importance[2,1]
    })
    vg = vg[select.genes,]
    vg = vg[order(vg$loess.padj),]
    list(pca.ratio=pca.var/mean(perm.pca.var), pca.var = pca.var, vg.num=length(select.genes), vg=vg)
  }


WGCNA_hetero <- function(norm.dat,select.cells, perm.num=10,vg.padj=0.5,rm.gene.mod=NULL)
  {
    vg= findVG(2^norm.dat[,select.cells]-1)
    select.genes = row.names(vg)[vg$loess.padj < vg.padj]
    dat = norm.dat[select.genes,select.cells]
    adj=adjacency(t(dat), power = 4,type="unsigned")
    adj[is.na(adj)]=0
    TOM = TOMsimilarity(adj,TOMType="unsigned")
    dissTOM = as.matrix(1-TOM)
    row.names(dissTOM)= colnames(dissTOM) = row.names(dat)
    geneTree = flashClust(as.dist(dissTOM), method =  "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight=0.99, deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = 5)
    gene.mod = split(row.names(dissTOM), dynamicMods)
    gene.mod = gene.mod[setdiff(names(gene.mod),"0")]
    if(!is.null(rm.gene.mod)| !is.null(rm.eigen)){
      if(is.null(rm.eigen)){
        rm.eigen = get_eigen(rm.gene.mod, norm.dat, select.cells)[[1]]
      }
    }
    else{
      rm.eigen=NULL
    }
    gm= filter_gene_mod(norm.dat, select.cells, gene.mod, minModuleSize=5, rm.eigen=rm.eigen,min.padj = 40, padj.th=0.05,min.cells=3, ...)
    if(length(select.genes)==0){
      return(0)
    }
    if(length(select.genes)==1){
      return(1)
    }
    rdWGCNA(norm.dat, select.genes, select.cells, minModuleSize=5,...)

    pca.var = summary(pca)$importance[2,1]
 
    ###shuffle expression of every gene.
    perm.pca.var = sapply(1:10, function(x){
      perm.dat = norm.dat[select.genes, select.cells]
      for(i in 1:nrow(perm.dat)){
        perm.dat[i, ] = sample(perm.dat[i, ])
      }
      perm.pca = prcomp(t(as.matrix(perm.dat)),tol=0.5)
      summary(perm.pca)$importance[2,1]
    })
    vg = vg[select.genes,]
    vg = vg[order(vg$loess.padj),]
    list(pca.ratio=pca.var/mean(perm.pca.var), pca.var = pca.var, vg.num=length(select.genes), vg=vg)
  }

