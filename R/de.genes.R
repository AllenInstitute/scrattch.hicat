library(limma)

vec_chisq_test <- function(x, x.total, y, y.total)
  {
    total <- x.total + y.total
    present = x + y
    absent = total - x - y
    o <- cbind(x, x.total - x, y, y.total - y)
    e <- cbind(present * x.total, absent * x.total, present * y.total, absent * y.total)
    e <- e/as.vector(total)
    stat <- rowSums(pmax(0, abs(o - e)-0.5)^2/e)    
    tmp <- cbind(stats=stat, pval=pchisq(stat, 1, lower.tail = FALSE), logFC = log2( (x * y.total)/(y * x.total)), diff = x/x.total - y/y.total)    
    tmp
  }


###If use voom mode, use counts, otherwise use norm.dat
DE_genes_pairs <- function(norm.dat, cl,pairs,counts=NULL, low.th=1,cl.present=NULL, method="limma", use.voom=FALSE){
  require(limma)
  select.cl = unique(c(pairs[,1],pairs[,2]))
  cl = cl[cl%in% select.cl]
  cl.means = as.data.frame(get_cl_means(norm.dat, cl))
  if(length(low.th)==1){
    low.th =setNames(rep(low.th, nrow(norm.dat)),row.names(norm.dat))      
  }
  if(is.null(cl.present)){
    cl.present = as.data.frame(get.cl.means(norm.dat > low.th[row.names(norm.dat)],cl))
  }
  cl.size = table(cl)
  de.df=list()
  fit = NULL
  if(method=="limma"){
    cl = setNames(as.factor(paste0("cl",cl)),names(cl))
    design=model.matrix(~0+ cl)
    colnames(design)=levels(as.factor(cl))
    if(use.voom & !is.null(counts)){
      v=voom(as.matrix(counts[row.names(norm.dat),names(cl)]), design)
      fit = lmFit(v, design)		
    }
    else{
      fit = lmFit(norm.dat[,names(cl)] , design=design)
    }
  }
  for(i in 1:nrow(pairs)){
    x = as.character(pairs[i,1])
    y = as.character(pairs[i,2])
    pair=paste(x,y,sep="_")
    if(i%%100==1){
      cat("Process pair", i, pair,"\n")
    }
    if(method=="limma"){
      ctr <<- paste(paste0("cl",x), "- ", paste0("cl",y))
      contrasts.matrix <- makeContrasts(ctr,  levels=design)
      fit2 = contrasts.fit(fit, contrasts.matrix)
      fit2 = eBayes(fit2)
      pval = fit2$p.value[,1]
      padj = p.adjust(pval)
      lfc = coef(fit2)[,1]
    }
    else if(method=="chisq"){
      lfc = cl.means[,x] - cl.means[,y]
      df=vec_chisq_test(cl.present[,x]*cl.size[[x]], cl.size[[x]], cl.present[,y]*cl.size[[y]], cl.size[[y]])
      pval = df[,"pval"]
      padj = p.adjust(pval)
    }
    
    de.df[[pair]]=data.frame(padj=padj,pval=pval,lfc=lfc,meanA=cl.means[[x]], meanB=cl.means[[y]],q1=cl.present[[x]], q2=cl.present[[y]])
    row.names(de.df[[pair]])= row.names(norm.dat)
  }
  return(de.df)
}

####Make sure dat and cl has the same dimension, and cells are in the same order
DE_genes_pw <- function(norm.dat,cl, ...)
{
	cn = as.character(sort(unique(cl)))
	cl.n = length(cn)	
  	pairs = cbind(rep(cn, rep(cl.n,cl.n)), rep(cn, cl.n))
  	pairs = pairs[pairs[,1]<pairs[,2],,drop=F]
  	de.df=DE_genes_pairs(norm.dat=norm.dat,cl=cl,pairs=pairs, ...)
  	return(de.df)
}

DE_genes_one_vs_other <- function(dat,cl, x.cl,low.th=setNames(rep(1,nrow(dat)),row.names(dat)),cl.present=NULL){
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  cl.means = tapply(1:ncol(dat), cl, function(x){
    rowMeans(dat[,x,drop=F])
  })
  cl.means = cl.means[tmp.cl]
  if(length(low.th)==1){
    low.th =setNames(rep(low.th, nrow(dat)),row.names(dat))      
  }
  if(is.null(cl.present)){
    cl.present = do.call("cbind",tapply(1:ncol(dat), cl, function(x){
      rowSums(dat[,x]>low.th[row.names(dat)])/length(x)  
    }))
  }
  de.df=list()
  for(y.cl in setdiff(tmp.cl, x.cl)){
    x=min(x.cl, y.cl)
    y=max(x.cl, y.cl)
    ctr <<- paste(x, "- ", y)
    contrasts.matrix <- makeContrasts(ctr,  levels=design)
    fit2 = contrasts.fit(fit, contrasts.matrix)
    fit2 = eBayes(fit2)
    padj = apply(fit2$p.value, 2, p.adjust)
    lfc = coef(fit2)
    pair=paste(x,y,sep="_")
    de.df[[pair]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1],meanA=cl.means[[x]], meanB=cl.means[[y]],q1=cl.present[,x], q2=cl.present[,y])
    row.names(de.df[[pair]])= row.names(dat)
    row.names(de.df[[pair]])= row.names(dat)
  }
  return(de.df)
}

####DESeq not used. too slow. 
DESeq_genes<- function(dat,cl,...)
  {
    require("DESeq2")
    require("parallel")
    df= data.frame(cl=cl)
    dds = DESeqDataSetFromMatrix(dat, colData=df, design=~ cl)
    dds=DESeq(dds,...)
    df = results(dds)
    colnames(df)[2] = "lfc"
    return(df)
  }

DESeq_genes_pw <- function(dat,cl, dds.file="dds.rda",mc.cores=4){
  require("DESeq2")
  require("parallel")
  if(is.null(dds.file) || !file.exists(dds.file)){
    df= data.frame(cl=cl)
    dds = DESeqDataSetFromMatrix(dat, colData=df, design=~ cl)
    dds=DESeq(dds)
    if(!is.null(dds.file)){
      save(dds, file=dds.file)
    }
  }
  tmp.cl = sort(unique(cl))
  pairs = data.frame(X=as.character(rep(tmp.cl,length(tmp.cl))), Y=as.character(rep(tmp.cl, rep(length(tmp.cl), length(tmp.cl)))), stringsAsFactors=F)
  pairs = pairs[pairs[,1]<pairs[,2],]
  de.df=sapply(1:nrow(pairs), function(i){
    x=pairs[i,1]
    y=pairs[i,2]
    res=results(dds, contrast = c("cl", x,y))
    colnames(res)[2]="lfc"
    res
  },simplify=F)
  names(de.df) = paste(pairs[,1],pairs[,2],sep="_")
  return(de.df)
}

de_param <- function(padj.th=0.01, lfc.th=1, q1.th=NULL, q2.th=NULL,q.diff.th=NULL)
  {
    list(padj.th=padj.th, lfc.th=lfc.th, q1.th=q1.th, q2.th=q2.th, q.diff.th = q.diff.th)
  }

de_pair <- function(df,  de.param = de_param())
  {
    df = df[order(df$pval,-abs(df$lfc)),]
    select=with(df, which(padj < padj.th & abs(lfc)>de.param$lfc.th))
    select=row.names(df)[select]
    if(is.null(select) | length(select)==0){
      return(list())
    }
    up = select[df[select, "lfc"]>0]
    down = select[df[select, "lfc"]<0]
    df = df[select,]
    if(!is.null(de.param$q.diff.th) & is.null(df$q.diff)){
      df$q.diff = with(df, abs(q1-q2)/pmax(q1,q2))
      df$q.diff[is.na(df$q.diff)]=0
    }
    if(!is.null(de.param$q1.th)){
      up = with(df[up,,drop=F], up[q1 > de.param$q1.th])
      down = with(df[down,,drop=F], down[q2 > de.param$q1.th])
    }
    if(!is.null(de.param$q2.th)){
      up = with(df[up,,drop=F], up[q2 < de.param$q2.th])
      down = with(df[down,,drop=F], down[q1 < de.param$q2.th])
    }
    if(!is.null(de.param$q.diff.th)){
      up = with(df[up,,drop=F], up[abs(q.diff) > de.param$q.diff.th])
      down = with(df[down,,drop=F], down[abs(q.diff) > de.param$q.diff.th])
    }
    select= c(up, down)    
    if(length(select)==0){
      return(list())
    }     
    else{
      df$padj[df$padj < 10^-20]=10^-20
      up.score=sum(-log10(df[up,"padj"]))
      down.score=sum(-log10(df[down,"padj"]))
      if(length(up)==0){up.score=0}
      if(length(down)==0){down.score=0}
      tmp=list(score=sum(-log10(df[select,"padj"])),
        up.score=up.score,
        down.score=down.score,
        num=length(select),
        up.num = length(up),
        down.num=length(down),
        genes= select,
        up.genes=up,
        down.genes=down, de.df=df[df$padj < padj.th,])
    }
}

de_score <- function(norm.dat, cl, min.cells=4,method="limma", low.th=1, de.param= de_param())
{
   if(is.factor(cl)){cl = droplevels(cl)}
   cn = as.character(sort(unique(cl)))
   cl.n = length(cn)	
   pairs = cbind(rep(cn, rep(cl.n,cl.n)), rep(cn, cl.n))
   pairs = pairs[pairs[,1]<pairs[,2],,drop=F]
   de_score_pairs(norm.dat, cl=cl, pairs=pairs,min.cells=min.cells, method=method, de.param= de.param)
}

de_score_pairs <- function(norm.dat, cl, pairs,  de.df=NULL, min.cells=4, method="limma", low.th=1, de.param=de_param())
{
  print(method)
  select.cl = unique(c(pairs[,1],pairs[,2]))
  cl = cl[cl %in% select.cl]
  if(is.factor(cl)){cl = droplevels(cl)}
  cl.size = table(cl)
  cl.n = names(cl.size)
  cl.small = cl.n[cl.size < min.cells]
  cl.big =  setdiff(cl.n,cl.small)
  select.pair = pairs[,1] %in% cl.big & pairs[,2] %in% cl.big
  de.genes=list()
  if(sum(select.pair)>0){
    cl = cl[cl %in% c(pairs[select.pair,1],pairs[select.pair,2])]
    select.cells = names(cl)
    select.genes= row.names(norm.dat)[rowSums(norm.dat[,select.cells] >low.th[row.names(norm.dat)]) > min.cells]
    if(is.null(de.df)){
      de.df = DE_genes_pairs(norm.dat[select.genes,select.cells], cl[select.cells], pairs[select.pair,,drop=F],low.th=low.th,method=method)
    }
    de.genes = sapply(names(de.df), function(x){
      if(is.null(de.df[[x]])){
        return(list())
      }
      df = de.df[[x]]
      de_pair(df,  de.param = de_param)
    },simplify=F)
  }
  else{
    de.df=list()
  }
  for(i in which(!select.pair)){
    pair=paste(pairs[i,1],pairs[i,2],sep="_")
    de.genes[[pair]]=list()
    de.df[[pair]]=list()
  }
  list(de.df=de.df, de.genes=de.genes)
}

de_matrix <- function(de.genes)
  {
    de.num = sapply(de.genes, function(x)length(x$genes))
    up.de.num = sapply(de.genes, function(x)length(x$up.genes))
    down.de.num = sapply(de.genes, function(x)length(x$down.genes))
    de.num.mat  = convert_pair_matrix(de.num)
    up.de.num.mat  = convert_pair_matrix(up.de.num)
    down.de.num.mat  = convert_pair_matrix(down.de.num)
    return(list(de.num.mat, up.de.num.mat, down.de.num.mat))
  }

