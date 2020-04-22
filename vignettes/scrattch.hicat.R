## ----tasic2016data------------------------------------------------------------
if(!"tasic2016data" %in% rownames(installed.packages())) {
  devtools::install_github("AllenInstitute/tasic2016data")
}
library(tasic2016data)

## ----Load libraries, messages = F, include=FALSE------------------------------
library(dendextend)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)

## ----Set up anno--------------------------------------------------------------
# Load sample annotations (anno)
anno <- tasic_2016_anno

# Make a data.frame of unique cluster id, type, color, and broad type
ref.cl.df <- as.data.frame(unique(anno[,c("primary_type_id", "primary_type_label", "primary_type_color", "broad_type")]))

#standardize cluster annoation with cluster_id, cluster_label and cluster_color. These are the required fields to visualize clusters properly.
colnames(ref.cl.df)[1:3] <- c("cluster_id", "cluster_label", "cluster_color")

# Sort by cluster_id
ref.cl.df <- ref.cl.df[order(ref.cl.df$cluster_id),]
row.names(ref.cl.df) <- ref.cl.df$cluster_id

ref.cl <- setNames(factor(anno$primary_type_id), anno$sample_name)

## ----Normalize data-----------------------------------------------------------
norm.dat <- log2(cpm(tasic_2016_counts)+1)

## ----Convert to sparse--------------------------------------------------------
norm.dat <- Matrix(cpm(tasic_2016_counts), sparse = TRUE)

norm.dat@x <- log2(norm.dat@x+1)

## ----Filter samples-----------------------------------------------------------
select.cells <- with(anno, sample_name[primary_type_label!="unclassified" & grepl("Igtp|Ndnf|Vip|Sncg|Smad3",primary_type_label)])

## ----Set Params---------------------------------------------------------------
de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.5, 
                     q.diff.th   = 0.7, 
                     de.score.th = 40)

## ----Run iter_clust, message=FALSE, warning=FALSE, echo=TRUE------------------
onestep.result <- onestep_clust(norm.dat, 
                               select.cells = select.cells, 
                               dim.method = "WGCNA", 
                               de.param = de_param(de.score.th=500))

## ---- fig.height=7, fig.width=7-----------------------------------------------
display.result = display_cl(onestep.result$cl, norm.dat, plot=TRUE, de.param=de.param)

## ---- message=FALSE, warning=FALSE, results="hide"----------------------------
WGCNA.clust.result <- iter_clust(norm.dat, 
                               select.cells = select.cells, 
                               dim.method = "WGCNA", 
                               de.param = de.param, 
                               result=onestep.result)

## ---- eval=FALSE--------------------------------------------------------------
#  WGCNA.clust.result <- iter_clust(norm.dat,
#                                 select.cells = select.cells,
#                                 dim.method = "WGCNA",
#                                 de.param = de.param)

## ---- message=FALSE, warning=FALSE, results="hide"----------------------------
gene.counts <- colSums(norm.dat > 0)
rm.eigen <- matrix(log2(gene.counts), ncol = 1)
row.names(rm.eigen) <- names(gene.counts)
colnames(rm.eigen) <- "log2GeneCounts"

## ---- eval=FALSE--------------------------------------------------------------
#  WGCNA.clust.result <- iter_clust(norm.dat,
#                                 select.cells = select.cells,
#                                 dim.method = "WGCNA",
#                                 de.param = de.param,
#                                 rm.eigen = rm.eigen)

## ----Merge clusters, message=FALSE, warning=FALSE, result="hide"--------------
WGCNA.merge.result <- merge_cl(norm.dat, 
                         cl = WGCNA.clust.result$cl, 
                         rd.dat = t(norm.dat[WGCNA.clust.result$markers, select.cells]),
                         de.param = de.param)

## ----Compare to benchmark, echo=FALSE, fig.height = 5, fig.width = 6.5--------
compare.result <- compare_annotate(WGCNA.merge.result$cl, ref.cl, ref.cl.df)
compare.result$g
cl <- compare.result$cl
cl.df <- compare.result$cl.df

## ---- fig.height=7, fig.width=7-----------------------------------------------
display.result = display_cl(cl, norm.dat, plot=TRUE, de.param=de.param, min.sep=4, n.markers=20)
de.genes= display.result$de.genes

## ----Drop cluster levels------------------------------------------------------
cl.clean <- droplevels(cl)

## ----Build dendrogram, warning=FALSE, message=FALSE, results="hide", fig.keep="all", fig.height = 4.5, fig.width = 7----
select.markers = select_markers(norm.dat, cl.clean, de.genes=de.genes,n.markers=50)$markers
cl.med <- get_cl_medians(norm.dat[select.markers,], cl)
##The prefered order for the leaf nodes.
l.rank <- setNames(1:nrow(cl.df), row.names(cl.df))
##Color of the leaf nodes.
l.color <- setNames(as.character(cl.df$cluster_color), row.names(cl.df))
dend.result <- build_dend(cl.med[,levels(cl.clean)],
                          l.rank, 
                          l.color,
                          nboot = 100)
dend <- dend.result$dend
###attach cluster labels to the leafs of the tree 
dend.labeled = dend
labels(dend.labeled) <- cl.df[labels(dend), "cluster_label"]
plot(dend.labeled)

## ----Reorder dendrogram-------------------------------------------------------
cl.clean <- setNames(factor(as.character(cl.clean), levels = labels(dend)), names(cl.clean))
cl.df.clean <- cl.df[levels(cl.clean),]

