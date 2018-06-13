library(iterclust)
library(matrixStats)
 
test_cv_cor <- function(norm.dat, cl, markers, n.bin=5,g.perc=1){
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  pred.cl = setNames(rep(NA, length(cl)), names(cl))
  for(i in 1:n.bin){
    print(i)
    train.cells = names(cl)[bins!=i]
    test.cells =names(cl)[bins==i]
    select.markers=sample(markers, round(length(markers)*g.perc))
    map.result <- map_by_cor(norm.dat[select.markers,], cl[train.cells], norm.dat[select.markers, test.cells])$pred.df
    pred.cl[test.cells] = as.character(map.result[test.cells, "pred.cl"])
  }
  return(pred.cl)
}


get_core_transition <- function(norm.dat, cl, select.markers, n.bin=5, n.iter=100, mc.cores=10)
  {
    cl.cv <- mclapply(1:n.iter, function(i){
      tmp=test_cv_cor(norm.dat, cl, select.markers, n.bin=n.bin)
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
    return(cell.cl.map.df)
  }

