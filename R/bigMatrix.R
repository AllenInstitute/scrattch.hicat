#library(bit64)

bigMatrix <- setClass("bigMatrix", 
                      slots = c(x = "integer64", 
                                i = "integer64", 
                                p = "numeric", 
                                dim = "numeric", 
                                row_id = "ANY", 
                                col_id = "ANY"),
                      package = "scrattch.hicat")


combine_dat <- function(dat.list)
{
  library(bit64)
  total.size = sum(sapply(dat.list, function(dat)length(dat@x)))
  x = integer64(total.size)
  i = integer64(total.size)
  p = 0
  k=0
  for(dat in dat.list){
    n = length(dat@x)
    x[(k+1):(k+n)] = dat@x
    i[(k+1):(k+n)] = dat@i
    p = c(p, dat@p[-1] + k)
    k = k + n 
  }
  
  row_id = row.names(dat.list[[1]])
  col_id = unlist(lapply(dat.list, colnames))
  ncol = sum(sapply(dat.list, ncol))
  dim = c(nrow(dat.list[[1]]),ncol)
  comb.dat = bigMatrix(x=x, i=i, p=p, row_id=row_id, col_id=col_id, dim=dim)
  return(comb.dat)
}

get_cols <- function(big.dat, cols)
{
  p = big.dat@p
  if(is.character(cols)){
    cols = match(cols, big.dat@col_id)
  }
  select = lapply(cols, function(col){
    if(p[col]+1 <= p[col+1]){
      (p[col]+1):(p[col+1])
    }
    else{
      NULL
    }
  })
  select.index = do.call("c", select)
  l = sapply(select, length)
  p = c(0,cumsum(l))
  i= (big.dat@i)[select.index]
  x = (big.dat@x)[select.index]
  mat=sparseMatrix(i=as.integer(i+1), x=as.double(x), p=p, dims=c(big.dat@dim[1],length(l)))
  colnames(mat) = big.dat@col_id[cols]
  row.names(mat) = big.dat@row_id
  return(mat)
}
