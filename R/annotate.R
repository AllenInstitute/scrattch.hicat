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


predict_annotate_cor <- function(cl, old.markers, old.cl, old.cl.df,norm.dat)
  {
    tmp = mapByCor(norm.dat[old.markers,names(old.cl)], old.cl, norm.dat[old.markers,names(cl)])
    pred.score= tmp$pred.score
    pred.cl= tmp$pred.cl
    ###compare predicted cluster member with the new clustering result 
    tb = table(cl, pred.cl)
    ###Reorder clusters
    tmp = apply(tb, 1, which.max)
    cl = setNames(factor(as.character(cl), levels=row.names(tb)[order(tmp)]), names(cl))
    levels(cl)=1:length(levels(cl))
    ###Assign the best matching old cluster to each new cluster. 
    tb=table(cl,pred.cl=pred.cl[names(cl)])
    tmp = colnames(tb)[apply(tb, 1, which.max)]
    
    cl.df = data.frame(pred.cl=tmp)
    match.id = match(cl.df$pred.cl, row.names(old.cl.df))
    cl.df = cbind(cl.df, old.cl.df[match.id,])
    
    ###plot the mapping
    tb.df = as.data.frame(tb)
    tb.df = tb.df[tb.df$Freq > 0,]
    library(ggplot2)
    tb.df$pred.prob=0
    select.cells=names(cl)
    for(i in 1:nrow(tb.df)){
      tmp.cells=select.cells[cl== as.character(tb.df[i, 1])]
      tb.df[i,"pred.prob"] = mean(pred.score[tmp.cells, as.character(tb.df[i,2])])
    }
    
    g= ggplot(tb.df, aes(x=cl, y=pred.cl)) + geom_point(aes(size=sqrt(Freq),color=pred.prob))
    g = g+ theme(axis.text.x=element_text(angle=90,size=7),axis.text.y=element_text(size=6)) + scale_color_gradient(low="white",high="darkblue")
    g= g+scale_size(range=c(0,3))
    return(list(cl=cl, cl.df=cl.df,g = g,tb.df=tb.df))
  }


compare_annotate<- function(cl, old.cl, old.cl.df)
  {
    common.cells=intersect(names(cl),names(old.cl))
    ###compare predicted cluster member with the new clustering result 
    tb = table(cl[common.cells], old.cl[common.cells])
    ###Reorder clusters
    tmp = apply(tb, 1, which.max)
    cl = setNames(factor(as.character(cl), levels=row.names(tb)[order(tmp)]), names(cl))
    levels(cl)=1:length(levels(cl))
    ###Assign the best matching old cluster to each new cluster. 
    tb=table(cl=cl[common.cells],old.cl=old.cl[common.cells])
    tmp = colnames(tb)[apply(tb, 1, which.max)]
  
    cl.df = data.frame(old.cl=tmp)
    cl.df$top = old.cl.df[tmp,"top"]
    cl.df$first = old.cl.df[tmp,"first"]
    cl.df$cluster_color = as.character(old.cl.df[tmp,"cluster_color"])
    cl.df$cluster_label = as.character(old.cl.df[tmp,"cluster_label"])

    ###plot the mapping
    tb.df = as.data.frame(tb)
    tb.df = tb.df[tb.df$Freq > 0,]
    library(ggplot2)
    select.cells=names(cl)
    tb.df$jaccard=0
    for(i in 1:nrow(tb.df)){
      tb.df[i,"jaccard"] = tb.df[i,"Freq"]/length(union(names(cl)[cl==as.character(tb.df[i,1])],names(old.cl)[old.cl==as.character(tb.df[i,2])]))
    }
    tb.df$old.cl.label = factor(old.cl.df[as.character(tb.df$old.cl),"cluster_label"], levels=old.cl.df$cluster_label)
    g= ggplot(tb.df, aes(x=cl, y=old.cl.label)) + geom_point(aes(size=sqrt(Freq),color=jaccard))
    g = g+ theme(axis.text.x=element_text(angle=90,size=7),axis.text.y=element_text(size=6)) + scale_color_gradient(low="white",high="darkblue")
    g= g+scale_size(range=c(0,3))
    return(list(cl=cl, cl.df=cl.df,g = g,tb.df=tb.df))
  }


