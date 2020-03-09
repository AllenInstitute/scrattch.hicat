get_core_intermediate <- function(norm.dat, cl, select.markers, n.bin=5, n.iter=100, mc.cores=10, method="median")
  {
    require(parallel)
    cl.cv <- mclapply(1:n.iter, function(i){
      tmp=map_cv(norm.dat, cl, select.markers, n.bin=n.bin,method=method)
    }, mc.cores=mc.cores)
    
    cell.cl.cor.map = do.call("rbind",sapply(cl.cv, function(x){
      df = data.frame(cell=names(x),cl=x)
    },simplify=F))
    cell.cl.cor.map = table(cell.cl.cor.map[,1],cell.cl.cor.map[,2])
    cell.cl.cor.map = cell.cl.cor.map / rowSums(cell.cl.cor.map)

    cell.cl.map.df = data.frame(org.cl= as.character(cl[row.names(cell.cl.cor.map)]),best.score=rowMaxs(cell.cl.cor.map), best.cl = colnames(cell.cl.cor.map)[apply(cell.cl.cor.map, 1, which.max)], stringsAsFactors=FALSE)
    row.names(cell.cl.map.df) = row.names(cell.cl.cor.map)
    tmp=get_pair_matrix_coor(cell.cl.cor.map, row.names(cell.cl.map.df), as.character(cell.cl.map.df$best.cl))
    tmp1 = cell.cl.cor.map
    tmp1[tmp]= 0
    cell.cl.map.df$second.score = rowMaxs(tmp1)
    cell.cl.map.df$second.cl =colnames(tmp1)[apply(tmp1,1, which.max)]
    cell.cl.map.df$second.cl[cell.cl.map.df$second.score ==0] = NA
    
    cell.cl.map.df$transition.cl = NA
    tmp = with(cell.cl.map.df, org.cl!=best.cl | best.score < 0.9)
    cell.cl.map.df[tmp,"transition.cl"] = as.character(cell.cl.map.df[tmp,"best.cl"])
    tmp = with(cell.cl.map.df, which(org.cl==transition.cl))
    cell.cl.map.df$transition.cl[tmp] = as.character(cell.cl.map.df[tmp,"second.cl"])
    
    cl.med <- do.call("cbind",tapply(names(cl), cl, function(x){
      rowMedians(as.matrix(norm.dat[select.markers,x]))
    }))
    row.names(cl.med) = select.markers
    
    cell.cl.cor=cor(as.matrix(norm.dat[select.markers, row.names(cell.cl.map.df)]), cl.med[select.markers,])
    cell.cl.map.df$cor = with(cell.cl.map.df, get_pair_matrix(cell.cl.cor, row.names(cell.cl.map.df),as.character(org.cl)))
    cell.cl.map.df$core = is.na(cell.cl.map.df$transition.cl)

    ###compute the transition edges between clusters
    cl.transition.df = with(cell.cl.map.df, as.data.frame(table(org.cl, transition.cl)))
    colnames(cl.transition.df)[1:2] = c("cl1","cl2")
    cl.transition.df = cl.transition.df[cl.transition.df$Freq > 0,]
    cl.transition.df$org.cl = as.character(cl.transition.df$org.cl)
    cl.transition.df$transition.cl = as.character(cl.transition.df$transition.cl)

    ###combine edges from both directions
    cl.transition.df$cl.min = pmin(cl.transition.df$org.cl, cl.transition.df$transition.cl)
    cl.transition.df$cl.max = pmax(cl.transition.df$org.cl, cl.transition.df$transition.cl)
    cl.transition.df$cl.pair = paste(cl.transition.df$cl.min, cl.transition.df$cl.max)
    trainsition.df.bi= do.call("rbind",tapply(1:nrow(cl.transition.df),cl.transition.df$cl.pair, function(x){
      tmp = cl.transition.df[x,][1,]
      tmp$Freq = sum(cl.transition.df[x,"Freq"])
      tmp[,c(4,5,3)]
    }))
    cl.size = table(cl)
    cl.transition.df.bi$cl.min.size = cl.size[cl.transition.df.bi$cl.min]
    cl.transition.df.bi$cl.max.size = cl.size[cl.transition.df.bi$cl.max]
    cl.transition.df.bi$ratio = with(cl.transition.df.bi,Freq/pmin(cl.min.size,cl.max.size))

    cl.transition.df = cl.transition.df.bi
    
    return(list(cell.cl.map.df=cell.cl.map.df, cl.transition.df = cl.transition.df))
  }
           



 








