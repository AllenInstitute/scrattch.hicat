get_gene_score <- function(ds, to.add, genes, de=NULL, max.num=1000,all.pairs=NULL)
  {
    if(is.null(de)){
      up.pair = to.add %>% filter(sign=="up") %>% pull(pair)
      down.pair = to.add %>% filter(sign=="down") %>% pull(pair)        
      if(!is.null(all.pairs)){
        select.pair.bin = all.pairs %>% filter(pair %in% c(up.pair, down.pair)) %>% pull(pair_bin) %>% unique
        ds = ds %>% filter(pair_bin %in% select.pair.bin)        
      }
      de = ds %>% filter(((sign=="up" & pair %in% up.pair) | (sign=="down" & pair %in% down.pair)))      
    }
    de = de %>% filter(gene %in% genes) %>% group_by(gene) %>% collect()
    gene.score = de %>% summarize(score = sum(max.num - rank))  %>% arrange(-score)
    return(list(de=de, gene.score=gene.score))
  }

select_markers <- function(ds, pairs = NULL, top.n=20)
  {
    if(is.null(pairs)){
      select.markers = ds %>% filter(rank < top.n) %>% pull(gene) %>% collect() %>% unique
    }
    else{
      select.markers = ds %>% filter(pair %in% pairs & rank < top.n) %>% pull(gene) %>% collect() %>% unique
    }
  }



select_markers_pair_direction <- function(ds, add.num, genes.allowed, de=NULL, max.num=1000, cl.means=NULL, all.pairs=NULL)
  {
    select.genes=c()
    all.de=NULL
    while(nrow(add.num)>0 & length(genes.allowed)>0){
      tmp= get_gene_score(ds, to.add=add.num, genes=genes.allowed,de=de, max.num=max.num, all.pairs=all.pairs)      
      de = tmp$de
      if(nrow(de)==0){
        break
      }
      if(is.null(all.de)){
        all.de=de
      }
      gene.score = tmp$gene.score
      if(!is.null(cl.means)){
        ###Penalize genes expressed in too many clbusters
        penalty = rowSums(cl.means[gene.score$gene, ,drop=F]> 3)
        gene.score$score = gene.score$score - penalty[gene.score$gene]
        gene.score = gene.score %>% arrange(-score)
      }
      g = gene.score$gene[1]
      print(g)
      #print(head(gene.score,1))
      new.checked = de  %>% filter(gene==g) %>% group_by(pair,sign) %>% summarize(checked=n())
      print(sum(new.checked$checked))
      add.num$checked=NULL
      add.num = add.num %>% left_join(new.checked, by=c("sign","pair"))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      to.remove.up = add.num %>% filter(num <= 0 & sign=="up") %>% pull(pair)
      to.remove.down = add.num %>% filter(num <= 0 & sign=="down") %>% pull(pair)
      add.num = add.num %>% filter(num > 0)                  
      genes.allowed = setdiff(genes.allowed,g)
      select.genes=c(select.genes,g)
      de = de %>% filter(gene!=g & ((sign=="up" & !pair %in% to.remove.up)| (sign=="down"& !pair %in% to.remove.down)))      
    }
    return(list(markers=select.genes, de= all.de))
  }

#' Title
#'
#' @param cl 
#' @param g1 
#' @param g2 
#' @param ds
#' @param top.n 
#' @param max.num 
#' @param n.markers 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
# Always directory
select_markers_pair_group_top<- function(cl, g1,g2,ds, genes, select.sign=c("up","down"),max.num=1000,all.pairs=NULL, n.markers=20)
{
  require(matrixStats)
  require(data.table)
  require(arrow)
  require(dplyr)
  pairs.df = as.data.frame(create_pairs(g1,g2))
  pairs.df$pair =row.names(pairs.df)
  to.add = pairs.df
  to.add$sign = "up"
  to.add = rbindlist(list(to.add, to.add %>% mutate(sign="down")))
  flip.sign = c("up"="down","down"="up")
  to.add = to.add %>% mutate(group.sign = ifelse(P1 %in% g1,sign, flip.sign[sign]))
  to.add = to.add %>% filter(group.sign %in% select.sign)
  if(!is.null(all.pairs)){
    select.pair.bin = all.pairs %>% filter(pair %in% to.add$pair) %>% pull(pair_bin) %>% unique
    ds = ds %>% filter(pair_bin %in% select.pair.bin)
  }
  de = ds %>% filter(pair %in% to.add$pair & gene %in% genes & rank < max.num) %>% collect() %>% left_join(to.add) %>% filter(!is.na(group.sign))  
  gene.score = de %>% group_by(group.sign, gene) %>% summarize(score = sum(max.num - rank))  %>% arrange(-score)
  up.genes = head(gene.score %>% filter(group.sign=="up") %>% pull(gene), n.markers)
  down.genes = head(gene.score %>% filter(group.sign=="down") %>% pull(gene), n.markers)
  return(list(up.genes=up.genes, down.genes=down.genes,de=de,to.add=to.add))
}

select_markers_pair_group<- function(cl, g1,g2,ds, genes, select.sign=c("up", "down"), max.num=1000,n.markers=20,de=NULL,all.pairs=NULL)
  {
    tmp = select_markers_pair_group_top(cl, g1,g2,ds=ds,genes=genes, select.sign=select.sign, max.num=max.num, n.markers=n.markers,all.pairs=all.pairs)
    de=tmp$de
    default.markers=c(tmp$up.genes,tmp$down.genes)
    add.num = tmp$to.add
    add.num$num = n.markers
    result = select_N_markers(ds, pairs=NULL, genes=genes, add.num=add.num, de= de, default.markers=default.markers,all.pairs=all.pairs)
    result$markers= c(default.markers, result$markers)
    return(result)
  }

select_N_markers <- function(ds, pairs, genes, pair.num=1, add.num=NULL, de=NULL, default.markers=NULL,all.pairs=NULL)
  {
    if(is.null(add.num)){
      add.up=data.frame(pair=pairs, sign="up", num=pair.num)
      add.down= data.frame(pair=pairs, sign="down", num=pair.num)
      add.num = rbind(add.up, add.down)
    }    
    if(!is.null(default.markers)){      
      if(is.null(de)){
        if(!is.null(all.pairs)){
          select.pair.bin = all.pairs %>% filter(pair %in% to.add$pair) %>% pull(pair_bin) %>% unique
          ds = ds %>% filter(pair_bin %in% select.pair.bin)
        }
        de.checked = ds %>% filter(pair %in% add.num$pair &gene %in% default.markers)  %>% collect()  
      }
      else{
        de.checked = de %>% filter(pair %in% add.num$pair &gene %in% default.markers) 
      }    
      de.checked.num = de.checked %>% group_by(pair, sign) %>% summarize(checked=n())
      de.checked.num$sign = factor(de.checked.num$sign)
      add.num = add.num %>% left_join(de.checked.num)
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      to.add = add.num %>% filter(num > 0)
      genes = setdiff(genes, default.markers)      
    }
    else{
      to.add=add.num
    }
    marker.select.result <- select_markers_pair_direction(ds, to.add, genes.allowed=genes)
  }



#' Title
#'
#' @param de.genes 
#' @param cl 
#' @param n.markers 
#' @param default.markers 
#' @param rm.genes 
#' @param up.gene.score 
#' @param down.gene.score 
#'
#' @return
#' @export
#'
#' @examples
select_pos_markers <- function(ds, cl, select.cl, genes, n.markers=1,  mc.cores=1, max.num=1000, all.pairs=NULL)
  {
    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types

    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
    #cl.markers <- list()
    #for(x in select.cl){
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      marker.result <- select_markers_pair_group(cl, g1,g2, ds=ds, genes=genes, select.sign="up", max.num=max.num, all.pairs=all.pairs, n.marker=n.markers)
      list(marker.result$markers)
    }
    return(cl.markers)
  }


select_top_pos_markers <- function(ds, cl, select.cl, genes, n.markers=3, mc.cores=10, max.num=1000,all.pairs=NULL)
  {
    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types

    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      print(x)
      g1=x
      g2 = setdiff(cl, x)
      markers= select_markers_pair_group_top(cl, g1,g2,ds, genes=genes, select.sign="up",max.num=max.num,n.markers=n.markers,all.pairs=NULL)$up.genes
      list(markers)
    }
    return(cl.markers)
  }



###select_markers_groups 
#' Title
#'
#' @param de.genes 
#' @param cl.group Assignment of clusters to groups cluster as names, and group id as values.
#'
#' @return
#' @export
#'
#' @examples


select_markers_groups <- function(ds, cl.group, select.groups=names(cl.group), n.markers=3,mc.cores=1,...)
  {

    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types
    all.cl = unlist(cl.group)
    #group.markers <- foreach(x=select.groups, .combine="c") %dopar% {
    group.markers=list()
    for(x in select.groups){
      print(x)
      g1 = cl.group[[x]]
      g2 = setdiff(all.cl, g1)              
      markers=select_markers_pair_group_top(cl=all.cl, g1,g2,ds=ds, select.sign="up",n.markers=n.markers, ...)$up.genes
      group.markers[[x]] = markers
    }
    return(group.markers)
  }

  
