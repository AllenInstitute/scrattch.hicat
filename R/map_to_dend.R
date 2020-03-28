#' Title
#'
#' @param cl.g 
#' @param cl.med 
#' @param markers 
#' @param dat 
#' @param map.dat 
#' @param select.cells 
#' @param p 
#' @param low.th 
#'
#' @return
#' @export
#'
#' @examples
resolve_cl <-
  function(cl.g,
           cl.med,
           markers,
           dat,
           map.dat,
           select.cells,
           p = 0.7,
           low.th = 0.2)
  {
    ##
    genes = names(markers)[markers > 0]
    tmp.cl = unlist(cl.g)
    
    ###For each branch point, find the highest expression cluster.
    tmp.med = sapply(cl.g, function(g)
      rowMaxs(cl.med[genes, g, drop = F]))
    row.names(tmp.med) = genes
    ###Make sure the genes are discriminative between all the branches.
    genes = genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]
    
    ###Sample the markers based on the weigts.
    ##TO DO: randomforest sometimes give importance value of 0. adjust for that.
    genes = sample(genes, round(length(genes) * p), prob = markers[genes])
    
    ###Compute the correlation with the median cluster profile.
    ###add drop=F
    cl.cor = cor(as.matrix(map.dat[genes, select.cells, drop = F]), cl.med[genes, tmp.cl, drop =
                                                                             F])
    cl.cor[is.na(cl.cor)] = 0
    ###Compute the best match in each branch.
    tmp.score = do.call("cbind", sapply(cl.g, function(x)
      rowMaxs(cl.cor[, x, drop = F]), simplify = F))
    row.names(tmp.score) = row.names(cl.cor)
    ####Determine the best match.
    best.score = setNames(rowMaxs(tmp.score), row.names(tmp.score))
    ###determine the difference from the best match.
    diff.score = best.score - tmp.score
    
    ####Give up on cells can't be discriminated,choose one branch randomly.
    unresolved.cl = row.names(tmp.score)[rowSums(diff.score < low.th) ==
                                           ncol(diff.score)]
    mapped.cl = setNames(sample(colnames(tmp.score), length(unresolved.cl), replace =
                                  T), unresolved.cl)
    
    ###Cells mapped to one or more branches.
    mapped.cells = setdiff(row.names(cl.cor), unresolved.cl)
    ###For binary branch, done already
    if (length(cl.g) == 2) {
      mapped.cl = c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells, , drop =
                                                                                F], 1, which.min)], mapped.cells))
      return(mapped.cl)
    }
    ##The remaining options for mapped cells
    tmp.cl = sapply(mapped.cells, function(x)
      colnames(diff.score)[which(diff.score[x,] < low.th)], simplify = F)
    ###cells with multiple options
    resolve.cells = names(tmp.cl)[sapply(tmp.cl, length) > 1]
    ###cells with only one option. Not further job.
    mapped.cells = setdiff(mapped.cells, resolve.cells)
    if (length(mapped.cells) > 0) {
      mapped.cl = c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
    }
    ###Resolve further options.
    if (length(resolve.cells) > 0) {
      tmp.cat = sapply(tmp.cl[resolve.cells], function(x)
        paste(x, collapse = " "))
      for (cat in unique(tmp.cat)) {
        tmp.cl = unlist(strsplit(cat, " "))
        select.cells = names(tmp.cat)[tmp.cat == cat]
        mapped.cl = c(
          mapped.cl,
          resolve_cl(
            cl.g[tmp.cl],
            cl.med,
            markers,
            dat,
            map.dat,
            select.cells,
            p = p,
            low.th = low.th
          )
        )
      }
    }
    return(mapped.cl)
  }


#' Title
#'
#' @param dend 
#' @param cl 
#' @param cl.med 
#' @param dat 
#' @param map.dat 
#' @param select.cells 
#' @param p 
#' @param low.th 
#' @param default.markers 
#'
#' @return
#' @export
#'
#' @examples
map_dend <-
  function(dend,
           cl,
           cl.med,
           dat,
           map.dat,
           select.cells,
           p = 0.8,
           low.th = 0.2,
           default.markers = NULL)
  {
    final.cl = c(setNames(rep(
      attr(dend, "label"), length(select.cells)
    ), select.cells))
    if (length(dend) <= 1) {
      return(final.cl)
    }
    markers = attr(dend, "markers")
    markers = markers[names(markers) %in% row.names(map.dat)]
    cl.g = sapply(dend, labels, simplify = F)
    names(cl.g) = 1:length(cl.g)
    select.cl = cl[cl %in% unlist(cl.g)]
    ###Sampling the cells from the reference cluster
    cells = unlist(tapply(names(select.cl), select.cl, function(x)
      sample(x, round(length(
        x
      ) * p))))
    genes = names(markers)
    genes = union(genes, default.markers)
    mapped.cl = resolve_cl(cl.g,
                           cl.med,
                           markers,
                           dat,
                           map.dat,
                           select.cells,
                           p = p,
                           low.th = low.th)
    if (length(mapped.cl) > 0) {
      for (i in unique(mapped.cl)) {
        select.cells = names(mapped.cl)[mapped.cl == i]
        if (length(select.cells) > 0) {
          final.cl = c(
            final.cl,
            map_dend(
              dend[[as.integer(i)]],
              cl,
              cl.med,
              dat,
              map.dat,
              select.cells,
              p = p,
              low.th = low.th
            )
          )
        }
      }
      return(cl = final.cl)
    }
    
  }

#' Title
#'
#' @param markers.cl.list 
#' @param map.dat 
#'
#' @return
#' @export
#'
#' @examples
get.markers.num <- function(markers.cl.list, map.dat)
{
  all.markers = unique(unlist(markers.cl.list))
  markers.num = sapply(markers.cl.list, function(x) {
    colSums(map.dat[x, ] > 0.5)
  })
}

#' Title
#'
#' @param dend 
#' @param memb 
#' @param map.dat 
#' @param exp.th 
#' @param conf.th 
#' @param min.genes.ratio 
#' @param min.genes 
#'
#' @return
#' @export
#'
#' @examples
summarize_cl <-
  function(dend,
           memb,
           map.dat,
           exp.th = 1,
           conf.th = 0.7,
           min.genes.ratio = 0.3,
           min.genes = 3)
  { 
    require(dendextend)
    node.height = setNames(get_nodes_attr(dend, "height"),
                           get_nodes_attr(dend, "label"))
    dend.list = dend_list(dend)
    tmp = sapply(dend.list, length) > 1
    tmp.dend.list = dend.list[tmp]
    markers.cl.list = lapply(tmp.dend.list, function(x) {
      tmp = attr(x, "markers.byCl")
      if (!is.null(tmp))
        # So it doesn't crash on the leaves
        names(tmp) = sapply(x, function(y)
          attr(y, "label"))
      tmp
    })
    names(markers.cl.list) = NULL
    markers.cl.list = do.call("c", markers.cl.list)
    
    all.markers = unique(unlist(markers.cl.list))
    memb.th = lapply(row.names(memb), function(cell) {
      ###Check all the node with confidence > conf.th
      x = memb[cell,]
      mapped.node = colnames(memb)[which(x > conf.th)]
      
      ###Check for detected markers at the given cell
      det.genes = all.markers[map.dat[all.markers, cell] >= exp.th]
      
      ##compute detected markers at every branch point.
      gene.olap = sapply(mapped.node, function(i)
        intersect(markers.cl.list[[i]], det.genes), simplify = F)
      gene.olap.num = sapply(gene.olap, length)
      #TO DO: weight markers instead of using absolute counts/ratio.
      
      ####set the root, so that root always succeed.
      gene.olap.num[attr(dend, "label")] = min.genes
      gene.olap[attr(dend, "label")] = ""
      gene.olap.ratio = gene.olap.num / sapply(markers.cl.list[names(gene.olap.num)], length)
      gene.olap.ratio[is.na(gene.olap.ratio)] = 1
      
      ###mapped nodes not met the minimal gene number/ratio constraints
      fail.node = mapped.node[gene.olap.ratio < min.genes.ratio |
                                gene.olap.num < min.genes]
      if (length(fail.node) > 0) {
        ###choose the mapped nodes above any failed nodes
        mapped.node = mapped.node[node.height[mapped.node] > max(node.height[fail.node])]
      }
      ###Choose the deepest nodes that pass all the criteria.
      mapped.node = mapped.node[order(node.height[mapped.node])]
      best.node = mapped.node[1]
      ###Get the markers on every mapped nodes.
      gene.anno = sapply(mapped.node, function(x) {
        paste0(x, ":", paste0(gene.olap[[x]], collapse = " "))
      })
      c(
        cl = best.node,
        score = x[best.node],
        marker.num = gene.olap.num[best.node],
        markers = paste(gene.anno, collapse = ",")
      )
    })
    memb.th = do.call("rbind", memb.th)
    row.names(memb.th) = row.names(memb)
    memb.df = data.frame(
      cl = memb.th[, 1],
      score = as.numeric(memb.th[, 2]),
      marker.num = as.integer(memb.th[, 3]),
      stringsAsFactors = F
    )
    memb.df$resolution.index = 1 - (node.height[memb.df$cl] / attr(dend, "height"))
    memb.df$cl = factor(memb.df$cl, names(node.height))
    #tmp = t(t(memb) * (1-node.height[colnames(memb)]/max(node.height)))
    #memb.df$h.score = rowMaxs(tmp)
    memb.df$h.score = memb.df$resolution.index * memb.df$score
    memb.df$markers = memb.th[, 4]
    ord = order(-memb.df$resolution.index, memb.df$cl,-memb.df$h.score)
    memb.df = memb.df[ord,]
    
    # adding resolution.index.percentile
    qq <-
      sort(quantile(memb.df$resolution.index, probs = seq(0, 1, by = 0.01)))
    memb.df$resolution.index.percentile <-
      findInterval(memb.df$resolution.index, qq) - 1
    
    
    return(memb.df)
  }

#' Title
#'
#' @param dend 
#' @param cl 
#' @param cl.med 
#' @param dat 
#' @param map.dat 
#' @param map.cells 
#' @param mc.cores 
#' @param bs.num 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
map_dend_membership <-
  function(dend,
           cl,
           cl.med,
           dat,
           map.dat,
           map.cells,
           mc.cores = 1,
           bs.num = 100,
           ...)
  {
    if (mc.cores == 1) {
      mem = sapply(1:bs.num, function(i) {
        print(i)
        ###determine which branch to t
        
        tmp = map_dend(dend, cl, cl.med, dat, map.dat, map.cells, ...)
      }, simplify = F)
      memb = unlist(mem)
    }
    else{
      require(foreach)
      require(doParallel)
      #fcluster <-makeCluster(mc.cores)
      fcluster <- makeForkCluster(mc.cores)
      registerDoParallel(fcluster)
      #on.exit(stopCluster(fcluster))
      mem = foreach(i = 1:bs.num, .combine = 'cbind') %dopar% map_dend(dend, cl, cl.med, dat, map.dat, map.cells, ...)
      
      stopCluster(fcluster)
      memb = as.character(mem)
      names(memb) = rownames(mem)
    }
    
    memb = data.frame(cell = names(memb), cl = memb)
    memb = table(memb$cell, memb$cl)
    memb = memb / bs.num
    tmp = get_nodes_attr(dend, "label")
    tmp = tmp[tmp %in% colnames(memb)]
    memb = memb[, tmp]
    return(memb)
  }

#' Title
#'
#' @param dend 
#' @param cl.df 
#' @param cl 
#' @param norm.dat 
#' @param query.dat 
#' @param bp.collapse.th 
#' @param mc.cores 
#' @param bs.num 
#' @param p 
#' @param low.th 
#' @param conf.th 
#' @param min.genes 
#' @param min.genes.ratio 
#'
#' @return
#' @export
#'
#' @examples
mapping <-
  function(dend,
           cl.df,
           cl,
           norm.dat,
           query.dat,
           bp.collapse.th = NULL,mc.cores=1, bs.num=100, p=0.7, low.th=0.15, conf.th=0.7, min.genes=1, min.genes.ratio=0.3
           )
  { 
    rownames(cl.df)=cl.df$cluster_id
    cltmp=cl.df[as.character(cl),"cluster_label"]
    names(cltmp)=names(cl)
    cl=factor(cltmp)
    
        
    query.dat.norm = log2(as.matrix(query.dat+1))
    idx=match(rownames(norm.dat), rownames(query.dat.norm))
    query.dat.norm=query.dat.norm[idx,]
    
    ### Some more initializations
    query.dat.cells = colnames(query.dat.norm)
    cl.med = get_cl_means(norm.dat, cl)
    rownames(cl.med)=rownames(norm.dat)
    
    # The key algorithm : Run 100 iterations, in each one sample 70% of the cells.
    # low.th is the minimum differnce in Pearson correlation required to decide on which branch to map to. If the difference 
    # is lower than this threshold, a random branch is chosen.
    # The resulted memb object is a matrix where every row is a cell, and columns are the nodes of the dendrogram. The values are the bootstrap support for that node.
    memb = map_dend_membership(dend, cl,cl.med, norm.dat, query.dat.norm, query.dat.cells, mc.cores=mc.cores, bs.num=bs.num, p=p, low.th=low.th)
    # analyze the results to generate the final output, i.e. for every cell the deepest node with th=0.8 confidence.
    # Also apply constraints on the minimum number of genes (or genes ratio)
    mapping.df = summarize_cl(dend, memb, query.dat.norm, conf.th=conf.th, min.genes=min.genes, min.genes.ratio=min.genes.ratio)
    
    # save results
    save(mapping.df, file=file.path(paste0("mapping.df.rda")))
    save(memb, file=file.path(paste0("mapping.memb.rda")))
    write.csv(mapping.df, file=file.path(paste0("mapping.df.csv")))
    write.csv(memb, file=file.path( paste0("mapping.memb.csv")))
    return (list(memb, mapping.df))
  }