map_by_cor <- function(train.dat, train.cl, test.dat,method="median")
  {
    if(method=="median"){
      require(matrixStats)
      cl.dat= do.call("cbind",tapply(names(train.cl), train.cl, function(x){
        rowMedians(as.matrix(train.dat[,x,drop=F]))
      }))
    }
    else{
      cl.dat = get_cl_means(train.dat, train.cl)
    }
    row.names(cl.dat)=row.names(train.dat)
    test.cl.cor = cor(as.matrix(test.dat), cl.dat)
    test.cl.cor[is.na(test.cl.cor)]=0
    pred.cl= setNames(colnames(test.cl.cor)[apply(test.cl.cor, 1, which.max)], row.names(test.cl.cor))
    pred.score = apply(test.cl.cor, 1, max)
    if(is.factor(train.cl)){
      pred.cl = setNames(factor(pred.cl, levels=levels(train.cl)), names(pred.cl))
    }
    return(list(pred.df= data.frame(pred.cl=pred.cl,pred.score=pred.score), cor.matrix=test.cl.cor))
  }


predict_annotate_cor <- function(cl, ref.markers, ref.cl, ref.cl.df,norm.dat)
  {
    common.cells= intersect(names(ref.cl),colnames(norm.dat))
    tmp = map_by_cor(norm.dat[ref.markers,common.cells], ref.cl[common.cells], norm.dat[ref.markers,names(cl)])
    pred.cl= setNames(factor(as.character(tmp$pred.df$pred.cl),levels=row.names(ref.cl.df)), row.names(tmp$pred.df))
    compare_annotate(cl, pred.cl, ref.cl.df)
  }

###cluster annotation ref.cl.df must include "cluster_label" column
compare_annotate<- function(cl, ref.cl, ref.cl.df, reorder=TRUE)
{
  common.cells=intersect(names(cl),names(ref.cl))
  ###compare predicted cluster member with the new clustering result 
  tb = table(cl[common.cells], ref.cl[common.cells])
  ###Reorder clusters
  tmp = apply(tb, 1, which.max)
  cl = setNames(factor(as.character(cl), levels=row.names(tb)[order(tmp)]), names(cl))
  cl.id.map=NULL
  if(reorder){
    cl.id.map <- data.frame(new=1:length(levels(cl)),old=levels(cl))
    levels(cl)=1:length(levels(cl))
  }
  ###Assign the best matching old cluster to each new cluster. 
  tb=table(cl=cl[common.cells],ref.cl=ref.cl[common.cells])
  tmp = colnames(tb)[apply(tb, 1, which.max)]
  
  cl.df = data.frame(ref.cl=tmp)
  cl.df = cbind(cl.df, ref.cl.df[tmp,])
  
  cl_label= rep("", nrow(cl.df))
  tmp = split(1:nrow(cl.df), cl.df$cluster_label)
  for(label in names(tmp)){
    x = tmp[[label]]
    if(length(x)>1){
      cl_label[x] <- paste(label, 1:length(x),sep="_")
    }
    else{
      cl_label[x]<-label
    }
  }
  cl.df$cluster_label = cl_label
  row.names(cl.df) = levels(cl)
  cl.size = table(cl)
  cl.df$size = cl.size[row.names(cl.df)]
  ###plot the mapping
  tb.df = as.data.frame(tb)
  tb.df = tb.df[tb.df$Freq > 0,]
  library(ggplot2)
  select.cells=names(cl)
  tb.df$jaccard=0
  for(i in 1:nrow(tb.df)){
    tb.df[i,"jaccard"] = tb.df[i,"Freq"]/length(union(names(cl)[cl==as.character(tb.df[i,1])],names(ref.cl)[ref.cl==as.character(tb.df[i,2])]))
  }
  tb.df$ref.cl.label = factor(ref.cl.df[as.character(tb.df$ref.cl),"cluster_label"], levels=ref.cl.df$cluster_label)
  g= ggplot(tb.df, aes(x=cl, y=ref.cl.label)) + geom_point(aes(size=sqrt(Freq),color=jaccard))
  g = g+ theme(axis.text.x=element_text(angle=90,size=7),axis.text.y=element_text(size=6)) + scale_color_gradient(low="yellow",high="darkblue")
  g= g+scale_size(range=c(0,3))
  return(list(cl=cl, cl.df=cl.df,g = g,tb.df=tb.df,cl.id.map=cl.id.map))
}



