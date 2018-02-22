get_pair_matrix <- function(m, rows, cols)
  {
    v = as.vector(m)
    if(!is.numeric(rows)){
      rows= match(rows, row.names(m))
    }
    if(!is.numeric(cols)){
      cols=match(cols, colnames(m))
    }
    coor = (cols - 1)* nrow(m) + rows
    return(v[coor])
  }



set_pair_matrix <- function(m, rows, cols,vals)
  {
    v = as.vector(m)
    if(!is.numeric(rows)){
      rows= match(rows, row.names(m))
    }
    if(!is.numeric(cols)){
      cols=match(cols, colnames(m))
    }
    coor = (cols - 1)* nrow(m) + rows
    m[coor] = vals
    return(m)
  }

convert_pair_matrix <- function(pair.num, l=NULL,directed=FALSE)
  {
    pairs = do.call("rbind",strsplit(names(pair.num),"_"))
    if(is.null(l)){
      l = sort(unique(as.vector(pairs)))
    }
    n.cl = length(l)
    pair.num.mat = matrix(0, nrow = n.cl, ncol=n.cl)
    row.names(pair.num.mat) = colnames(pair.num.mat) = l
    for(i in 1:nrow(pairs)){
      pair.num.mat[pairs[i,1], pairs[i,2]] = pair.num[i]
      if(!directed){
        pair.num.mat[pairs[i,2], pairs[i,1]] = pair.num[i]
      }
    }
    return(pair.num.mat)
  }


get_cl_sums <- function(mat, cl)
  {
    tmp.df= data.frame(cell=names(cl), cl=cl)
    tb=xtabs(~cl+cell, data=tmp.df)
    tb = Matrix(tb, sparse=TRUE)
    tmp=tcrossprod(mat[,colnames(tb)], tb)
    cl.sums = as.matrix(tmp)
    return(cl.sums)
  }

get_cl_means <- function(mat, cl)
  {
    cl.sums = get_cl_sums(mat, cl)
    cl.size = table(cl)
    cl.means = as.matrix(t(t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
    return(cl.means)
  }

sparse_cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans
  
  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
             crossprod(t(Ex))*(n-length(ii))
             )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}

calc_tau <- function(m, byRow=TRUE)
{
  if(!byRow){
    m = t(m)
  }
  m = m/rowMaxs(m)
  tau = rowSums(1 - m)/(ncol(m) - 1)
  tau[is.na(tau)]=0
  return(tau)
}


sample_cells <- function(cl,max.cl.size, weights=NULL)
{
  
  sampled.cells = unlist(tapply(names(cl),cl, function(x){
    if(length(x) > max.cl.size){
      if(!is.null(weights)){
        x= sample(x, max.cl.size, prob= weights[x])
      }
      else{
        x= sample(x, max.cl.size)
      }
    }
    x
  },simplify=FALSE))
}
  

sample_cells_by_genecounts <- function(cl, norm.dat, max.cl.size=200)
{
  select.cells=names(cl)
  if(is.matrix(norm.dat)){
    cell.gene.counts= colSums(norm.dat[,select.cells]>0)
  }
  else{
    cell.gene.counts= Matrix::colSums(norm.dat[,select.cells]>0)
  }
  cell.weights = cell.gene.counts - min(cell.gene.counts)+1
  sample_cells(cl, weights=cell.weights, max.cl.size=max.cl.size)  
}
