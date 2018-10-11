##########################################################################################################
# Functions for formatting the region and layer calls.
# It may or may not make sense to include these in hicat...


#' Numeric layer for annotation object
#'
#' This function expects a specific format for input that may not make sense for most users.
#'
#' @param Samp.dat annotation data frame or tibble in standard format (for human data)
#' @param column name of column that contains layer information in format "REGION_L#"
#'
#' @return a numeric vector corresponding to cortical layer (or NA for subcortical regions)
#' @export
#'
makeLayerLabel <- function(Samp.dat,
                           column = "roi_label") {
  val <- as.character(as.matrix(Samp.dat[, column]))
  val <- gsub("2_3", "3", val)
  val <- as.character(sapply(val, function(x) return(strsplit(x, "_")[[1]][2])))
  val <- substr(val, 2, nchar(val))
  val[val == "ZZ_Missing"] <- NA
  val <- gsub("1-6", "0", val)
  val <- as.numeric(substr(val, 1, 1))
  val[val == 0] <- NA
  val
}

#' Split region from roi for annotation object
#'
#' This function expects a specific format for input that may not make sense for most users.
#'
#' @param Samp.dat annotation data frame or tibble in standard format (for human data)
#' @param column name of column that contains layer information in format "REGION_L#"
#'
#' @return a chracter vector corresponding to brain region
#' @export
#'
makeRegionLabel <- function(Samp.dat, column = "roi_label") {
  val <- as.character(as.matrix(Samp.dat[, column]))
  val[nchar(val) <= 4] <- paste0("LGN_", val[nchar(val) <= 4]) # To account for incorrect formatting of roi
  val <- gsub("-", "_", val)
  val <- as.character(sapply(val, function(x) return(strsplit(x, "_")[[1]][1])))
  val
}


##########################################################################################################
# Function to make sure colors are unique

#' Adjust color vector
#'
#' Note: This function is an internal function
#'
#' @param col a vector of colors (in a format compatable with col2rgb)
#' @param r integer value to increment red (-255 to 255)
#' @param g integer value to increment green (-255 to 255)
#' @param b integer value to increment blue (-255 to 255)
#'
#' @return a vector with all colors adjusted by the r, g, b integer
#'   values and in #xxxxxx format.
#'
incrementHex <- function(col, r = 0, g = 0, b = 0) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  apply(
    sapply(col, col2rgb) / 255, 2,
    function(x)
      rgb(
        pmin(1, pmax(0, x[1] + r / 255)),
        pmin(1, pmax(0, x[2] + g / 255)),
        pmin(1, pmax(0, x[3] + b / 255))
      )
  )
}

#' Make a vector of unique colors
#'
#' This function converts any vector of colors to a new vector where each
#'   element is unique but the elements overall are as close as possible to
#'   the original elements.
#'
#' @param colorVector a vector of colors (in a format compatable with col2rgb)
#' @param seed random seem for reproducibility
#'
#' @return a color vector very close to the initial colorVector, but
#'   where all color values are unique
#' @export
#'
makeColorsUnique <- function(colorVector, seed = 1) {
  set.seed(seed)
  while (max(table(colorVector)) > 1) {
    cl <- names(table(colorVector))[table(colorVector) > 1]
    for (i in which(is.element(colorVector, cl)))
      colorVector[i] <- as.character(incrementHex(
        colorVector[i], sample(-2:2, 1), sample(-2:2, 1), sample(-2:2, 1)
      ))
  }
  colorVector
}


##########################################################################################################
# Function for calculting top DEX genes.

#' Find specific marker gene for each cluster
#'
#' This function identifies the top marker gene (or "none") for each cluster
#'   based on the cluster medians and proportions of cells expressing a gene
#'   and based on various parameters.
#'
#' Note: this could potentially be replaced by an existing hicat function for
#'   identifying marker genes.  This function is fast and chooses 'reasonable" genes.
#'
#' @param propExpr matrix of proportions of cells expressing a gene in each cluster
#'   (genes=rows, clusters=columns)
#' @param medianExpr matrix of median expression per cluster (genes=rows, clusters=columns)
#' @param propDiff Must have difference in proportion higher than this value in "on" cluster
#'   compared with each other cluster
#' @param propMin Must have higher proportion in "on" cluster
#' @param medianFC Must have median fold change greater than this value in "on" group vs.
#'   each other cluster
#' @param excludeGenes Genes exlcuded from marker consideration (NULL by default)
#' @param sortByMedian Should genes passing all filters be prioritized by median fold
#'   change (TRUE, default) or by difference in proportion between clusters (FALSE)
#'
#' @return a vector of the top marker gene per cluster (or "none" for clusters with none)
#' @export
#'
getTopMarkersByPropNew <- function(
                                   propExpr,
                                   medianExpr,
                                   propDiff = 0,
                                   propMin  = 0.5,
                                   medianFC = 1,
                                   excludeGenes = NULL,
                                   sortByMedian = TRUE) {
  specGenes <- rep("none", dim(propExpr)[2])
  names(specGenes) <- colnames(propExpr)
  propSort  <- t(apply(propExpr, 1, function(x) return(-sort(-x))))
  propWhich <- t(apply(propExpr, 1, function(x, y) return(y[order(-x)]), colnames(propExpr)))
  medianDif <- apply(data.frame(propWhich[, 1], medianExpr), 1, function(x, y) {
    wIn  <- y == as.character(x[1])
    mIn  <- as.numeric(x[2:length(x)])[wIn]
    mOut <- max(as.numeric(x[2:length(x)])[!wIn])
    return(mIn - mOut)
  }, colnames(propExpr))
  keepProp  <- (propSort[, 1] >= propMin) & ((propSort[, 1] - propSort[, 2]) > propDiff) &
    (medianDif >= medianFC) & (!is.element(rownames(propExpr), excludeGenes))
  propSort  <- propSort[keepProp, ]
  propWhich <- propWhich[keepProp, ]
  ord <- order(
    -medianDif[keepProp] * ifelse(sortByMedian, 1, 0),
    propSort[, 2] - propSort[, 1]
  )
  propSort  <- propSort[ord, ]
  propWhich <- propWhich[ord, ]
  while (sum(keepProp) > 1) {
    keepProp <- !is.element(propWhich[, 1], names(specGenes)[specGenes != "none"])
    if (sum(keepProp) <= 1) break
    tmp <- propWhich[keepProp, ]
    specGenes[tmp[1, 1]] <- rownames(tmp)[1]
  }
  specGenes
}


##########################################################################################################
# Function for renaming and reordering the HCT clusters

#' Rename clusters using genes and metadata
#'
#' This function uses information from both the data (gene expression) and meta-data
#'   (annotation) objects to build a useful automated cluster name.
#'
#' When all options are selected, the outputed format is as follows: [cell class]_[layer
#'   range]_[broad marker gene]_[specific marker gene]_[brain region with most cells (and
#'   scaled fraction of cells)]_[best matched type from previous taxonomy]_[number of cells
#'   in cluster].  The output is a data frame with information about each cluster, including
#'   the new cluster names.  \code{updateSampDat} needs to be run after \code{renameAndOrderClusters}
#'   to apply the new cluster names to each sample.  If a dendrogram has already been
#'   created, the dendrogram labels will also need to be changed separately.
#'
#' @param sampleInfo Sample information with rows as samples and columns for annotations.
#'   All samples in sampleInfo are used for renaming (so subset prior to running this function
#'   if desired). Columns must include "cluster_id", "cluster_label", and "cluster_color".
#' @param classNameColumn Column name where class information is stored (e.g., inh/exc/glia),
#'   or NULL if you'd like it to be defined based on \code{classGenes}
#' @param classGenes Set of genes for defining classes (which is ignored in this context
#'   if \code{classNameColumn!=NULL}). Also used if \code{broadClass} gene is not expressed.
#' @param classLevels A vector of the levels for classes of the same length (and in the same
#'   order) as classGenes.  Either include all relevant levels or set to NA for none if using
#'   \code{classNameColumn}
#' @param layerNameColumn Column name where the (numeric) layer info is stored (NA if none)
#' @param regionNameColumn Column name where the (character) region info is stored (NA if none)
#' @param matchNameColumn Column name where the (character) comparison info stored (e.g.,
#'   closest mapping cell type for each cell pre-calculated against a previous taxonomy; NA if none)
#' @param newColorNameColumn Column name where the new cluster colors are found (e.g., color
#'   column corresponding to \code{matchNameColumn}). NA keeps the current colors.
#' @param otherColumns Other columns to transfer to the output variable.  Note that the value from
#'   a random sample in the cluster is returned, so this usually should be left as default (NULL).
#' @param propLayer Proportion of cells (relative to max) must be higher than this for a cluster
#'   to be considered as expressed in a particular layer (default is 0.3).
#' @param dend Dendrogram object, only used for ordering of clusters (NULL as default)
#' @param orderbyColumns column names indicating the outputted cluster order (not used unless 
#'   dend=NULL).  Must be some combination of "layer", "region", and "topMatch" in any order (or 
#'   NULL).  Default is first by "layer" than "region" then "topMatch".
#' @param includeClusterCounts Should the number of cells in each cluster be included in name?
#' @param includeBroadGenes Should broad genes be included in the name (if so, \code{broadGenes}
#'   must be provided)?
#' @param broadGenes List of broad genes, where the top median CPM in cluster is included in name
#' @param includeSpecificGenes Should specific genes be included in the name?  If TRUE, the next
#'   seven parameters are used to call \code{getTopMarkersByPropNew}.
#' @param propExpr matrix of proportions of cells expressing a gene in each cluster
#'   (genes=rows, clusters=columns)
#' @param medianExpr matrix of median expression per cluster (genes=rows, clusters=columns)
#' @param propDiff Must have difference in proportion higher than this value in "on" cluster
#'   compared with each other cluster
#' @param propMin Must have higher proportion in "on" cluster
#' @param medianFC Must have median fold change greater than this value in "on" group vs.
#'   each other cluster
#' @param excludeGenes Genes exlcuded from marker consideration (NULL by default)
#' @param sortByMedian Should genes passing all filters be prioritized by median fold
#'   change (TRUE, default) or by difference in proportion between clusters (FALSE)
#' @param sep  Separation character for renaming (default is "_")
#'
#' @return A data frame of cluster information, which includes the new and old names, the
#'   requested variables from sampleInfo, and all the specific components of the new name.
#'   This is the required input for \code{updateSampDat} in the appropriate format.
#' @export
#'
renameAndOrderClusters <- function(
                                   sampleInfo,

                                   # Variables for broad class call (two options)
                                   classNameColumn = "cluster_type_label",
                                   classGenes  = c("GAD1", "SLC17A7", "SLC1A3"),
                                   classLevels = c("inh", "exc", "glia"),

                                   # Variables for other metadata (e.g., layer, region, matching cell types with previous data sets)
                                   layerNameColumn  = "layer_label",
                                   regionNameColumn = "Region_label",
                                   matchNameColumn  = "cellmap_label",
                                   newColorNameColumn = "cellmap_color",
                                   otherColumns = NULL,
                                   propLayer = 0.3,

                                   # Other naming and ordering options
                                   dend = NULL,
                                   orderbyColumns = c("layer", "region", "topMatch"),
                                   includeClusterCounts = FALSE,

                                   # Variables for including broad genes
                                   includeBroadGenes = FALSE,
                                   broadGenes = NULL,

                                   # Variables for including specific genes
                                   includeSpecificGenes = FALSE,
                                   propExpr = NULL,
                                   medianExpr = NULL,
                                   propDiff = 0,
                                   propMin  = 0.5,
                                   medianFC = 1,
                                   excludeGenes = NULL,
                                   sortByMedian = TRUE,

                                   sep = "_") {
  sampleInfo  <- as.data.frame(sampleInfo)

  ## Define clusterInfo variable to store cluster-level info
  kpColumns   <- unique(c("cluster_id", "cluster_label", "cluster_color", classNameColumn, matchNameColumn, regionNameColumn, otherColumns))
  kpColumns   <- intersect(kpColumns, colnames(sampleInfo))
  clusterInfo <- t(sampleInfo[, kpColumns])
  colnames(clusterInfo) <- sampleInfo[,"cluster_label"]
  clusterInfo <- clusterInfo[, unique(colnames(clusterInfo))]
  clusterInfo <- t(clusterInfo)
  clusterInfo <- as.data.frame(clusterInfo)
  clusterInfo$old_cluster_label <- clusterInfo$cluster_label
  clusterInfo <- clusterInfo[order(clusterInfo[, "cluster_id"]), ]
  rownames(clusterInfo) <- 1:dim(clusterInfo)[1]


  ## Define class information
  if (!is.null(classNameColumn)) {
    if (!is.na(classLevels[1])) {
      clusterInfo[, classNameColumn] <-
        factor(clusterInfo[, classNameColumn], levels = classLevels)
    }
    classLabel <- clusterInfo[, classNameColumn]
    if (is.factor(classLabel)) classLabel <- droplevels(classLabel)
    clusterInfo$class <- classLabel
    if (is.null(classGenes)) classGenes <- NA
    if (is.na(classGenes[1])) {
      classLab <- rep("none", length(classLabel)) # Don't include any gene here. # NEW
    } else {
      classLab <- classGenes[apply(propExpr[classGenes, ], 2, which.max)] # NEW
    } 
  } else {
    classLab <- classGenes[apply(propExpr[classGenes, ], 2, which.max)]
    names(classLevels) <- classGenes
    clusterInfo$class  <- factor(classLevels[classLab], levels = classLevels)
    classNameColumn    <- "class"
  }

  ## Identify broad genes, if desired
  if (includeBroadGenes) {
    broadLab  <- broadGenes[apply(propExpr[broadGenes, ], 2, which.max)]
    broadProp <- apply(propExpr[broadGenes, ], 2, max)
    broadLab[broadProp < propMin] <- classLab[broadProp < propMin]
    names(broadLab) <- colnames(propExpr)
    clusterInfo$broadGene <- broadLab[clusterInfo$old_cluster_label]  # FIX
  }

  ###### (BEGIN: THIS PART COULD BE REPLACED BY OTHER HICAT MARKER GENE SELECTION STRATEGY) ######
  ## Identify specific genes, if desired
  if (includeSpecificGenes) {
    kpGn <- rep(TRUE, dim(propExpr)[1]) # betaScore>=minBeta
    specGenes <- getTopMarkersByPropNew(
      propExpr = propExpr[kpGn, ], 
      medianExpr = medianExpr[kpGn, ], 
      propDiff = propDiff, 
      propMin  = propMin,
      medianFC = medianFC, 
      excludeGenes = excludeGenes
    )
    specGenes0 <- getTopMarkersByPropNew(
      propExpr = propExpr[kpGn, ], 
      medianExpr = medianExpr[kpGn, ], 
      propDiff = 0, 
      propMin  = propMin,
      medianFC = 0, 
      excludeGenes = excludeGenes
    )
    for (s in colnames(propExpr)[(specGenes == "none")]) {
      if ((specGenes0[s] != "none") & (specGenes[s] == "none")) {
        specGenes[s] <- specGenes0[s]
      }
    }
    clusterInfo$specificGene <- specGenes[clusterInfo$old_cluster_label] # FIX
  }
  ###### (END: THIS PART COULD BE REPLACED BY OTHER HICAT MARKER GENE SELECTION STRATEGY) ######

  ## Add additional cluster information for renaming
  cl3 <- sampleInfo[, "cluster_label"] # FIX
  names(cl3) <- sampleInfo$sample_id # FIX
  cl3 <- factor(cl3, levels = clusterInfo$cluster_label) # FIX

  ## Match cluster color to desired cluster color column
  if (is.null(newColorNameColumn)) newColorNameColumn <- NA
  if (is.na(newColorNameColumn)) newColorNameColumn <- "cluster_color"
  colorVec <- as.character(tapply(names(cl3), cl3, function(x) {
    col <- as.factor(sampleInfo[, newColorNameColumn])
    names(col) <- sampleInfo$sample_id
    return(names(sort(-table(col[x])))[1])
  }))
  clusterInfo$cluster_color <- makeColorsUnique(colorVec)

  ## Include information about brain region in name, if desired
  if (is.null(regionNameColumn)) regionNameColumn <- NA
  if (!is.na(regionNameColumn)) {
    regionVec <- as.character(tapply(names(cl3), cl3, function(x) {
      rg <- as.factor(sampleInfo[, regionNameColumn])
      names(rg) <- sampleInfo$sample_id
      rg <- rg[x]
      rg <- table(rg) / table(sampleInfo[, regionNameColumn])
      rg <- -sort(-round(100 * rg / sum(rg)))[1]
      return(paste(names(rg), rg, sep = "~"))
    }))
    clusterInfo$region <- regionVec
  }

  ## Include information about layer, if desired
  if (!is.na(layerNameColumn)) {
    clLayer <- sampleInfo[, layerNameColumn]
    names(clLayer) <- names(cl3)
    layerVec <- (tapply(names(cl3), cl3, function(x) {
      lyy <- factor(clLayer)[x]
      if (mean(is.na(lyy)) >= 0.5) return(c(0, 0, 0, 0, 0, 0)) # Address non-cortical areas
      lyy <- lyy[!is.na(lyy)] # Address non-cortical areas
      layTab <- cbind(as.numeric(names(table(lyy))), table(lyy), table(clLayer))
      return(((layTab[, 2] / layTab[, 3]) / max(layTab[, 2] / layTab[, 3]))) # replace max with sum?
    }))
    rn <- names(layerVec)
    layerVec <- matrix(unlist(layerVec), ncol = 6, byrow = TRUE)
    rownames(layerVec) <- rn
    colnames(layerVec) <- 1:6
    layLab <- apply(layerVec, 1, function(x, y) {
      z <- as.numeric(colnames(layerVec)[x >= y])
      if (length(z) == 0) return("x") # Replace with whatever we want to call layers outside cortex
      if (length(z) == 1) return(z)
      return(paste(range(z), collapse = "-"))
    }, propLayer)
    layLab <- paste0("L", layLab)
    clusterInfo$layer <- layLab
  }

  ## Determine closest matching other column (e.g., MTG cluster) if desired
  if (is.null(matchNameColumn)) matchNameColumn <- NA
  if (!is.na(matchNameColumn)) {
    matchVec <- as.character(tapply(names(cl3), cl3, function(x) {
      y <- is.element(sampleInfo$sample_id, x)
      nm <- -sort(-table(sampleInfo[y, matchNameColumn]))
      return(names(nm)[1])
    }))
    clusterInfo$topMatch <- matchVec
  }

  ## Count the number of cells per cluster, if desired
  if (includeClusterCounts) {
    clusterInfo$cellCount <- table(factor(sampleInfo$cluster_id, levels = as.numeric(clusterInfo$cluster_id)))
  }

  ## Rename the clusters based on the above info
  for (i in 1:dim(clusterInfo)[1]) {
    lab <- NULL
    if (!is.null(clusterInfo$class)) lab <- paste(lab, clusterInfo$class[i], sep = sep)
    if (!is.null(clusterInfo$layer)) lab <- paste(lab, clusterInfo$layer[i], sep = sep)
    if (!is.null(clusterInfo$broadGene)) lab <- paste(lab, clusterInfo$broadGene[i], sep = sep)
    if (!is.null(clusterInfo$specificGene)) lab <- paste(lab, clusterInfo$specificGene[i], sep = sep)
    if (!is.null(clusterInfo$region)) lab <- paste(lab, clusterInfo$region[i], sep = sep)
    if (!is.null(clusterInfo$topMatch)) lab <- paste(lab, clusterInfo$topMatch[i], sep = sep)
    if (!is.null(clusterInfo$cellCount)) lab <- paste0(lab, sep, "N=", clusterInfo$cellCount[i])
    lab <- substr(lab, 2, nchar(lab))
    clusterInfo[i, "cluster_label"] <- lab
  }

  # Postprocessing things for cluster names, which probably should be removed
  clNames <- clusterInfo[, "cluster_label"]
  clNames <- gsub("Exc L1", "Exc L2", clNames) # Ensure that no excitatory clusters are included in layer 1
  clNames <- gsub("exc L1", "exc L2", clNames) # Ensure that no excitatory clusters are included in layer 1
  clNames <- gsub("L2-2", "L2", clNames) # Ensure that no excitatory clusters are included in layer 1
  clusterInfo[, "cluster_label"] <- clNames

  ## Determine a new optimal order based on inputted parameters (default broad class, then layer, then region)
  ordNew <- 1:dim(clusterInfo)[1]
  ordCols<- intersect(intersect(orderbyColumns,colnames(clusterInfo)),c("topMatch","layer","region"))
  ordVal <- "ordNew = order(clusterInfo[,classNameColumn],"
  if (length(ordCols)>=1) for (i in 1:length(ordCols)){
    ordVal <- paste0(ordVal,"clusterInfo[,\"",ordCols[i],"\"],")
  }
  ordVal <- paste0(ordVal,"clusterInfo[,\"cluster_id\"])")

  eval(parse(text = ordVal))
  clusterInfo <- clusterInfo[ordNew, ]

  ## Now reorder based on dendrogram order (if provided)
  if (!is.null(dend)) {
    lab <- c(labels(dend), setdiff(labels(dend), clusterInfo$old_cluster_label))
    clusterInfo <- clusterInfo[match(lab, clusterInfo$old_cluster_label), ] # DOUBLE CHECK THAT MATCH IS NOT BACKWARDS
    print("NOTE: You will need to regenerate the dendrogram with these new cluster labels now.")
  }

  ## Make new cluster ids that correspond to the current order of the data set
  rownames(clusterInfo) <- clusterInfo$lrank <- clusterInfo$cluster_id <- 1:dim(clusterInfo)[1]

  ## Return clusterInfo
  clusterInfo
}


#' Update cluster names in annotation data frame
#'
#' @param Samp.dat Sample information with rows as samples and columns for annotations.
#' @param clusterInfo Output cluster information data frame from \code{renameAndOrderClusters}
#'
#' @return Samp.dat variable but with updated cluster names
#' @export
#'
updateSampDat <- function(Samp.dat, 
                          clusterInfo) {
  lab   <- as.character(Samp.dat$cluster_label)
  for (i in 1:dim(clusterInfo)[1]) {
    kp <- lab == clusterInfo$old_cluster_label[i]
    Samp.dat$cluster_label[kp] <- clusterInfo$cluster_label[i]
    Samp.dat$cluster_id[kp]    <- clusterInfo$cluster_id[i]
    Samp.dat$cluster_color[kp] <- clusterInfo$cluster_color[i]
  }
  Samp.dat
}
