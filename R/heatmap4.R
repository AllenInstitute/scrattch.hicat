heatmap.4 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                       distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                        "row", "column", "none"), symm = FALSE, scale = c("none", 
                           "row", "column"), na.rm = TRUE, revC = identical(Colv, 
                           "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
                         scale != "none", col = "heat.colors", colsep, rowsep, 
                       sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
                       notecol = "cyan", na.color = par("bg"), trace = c("column", 
                         "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                       vline = median(breaks), linecol = tracecol, margins = c(5, 
                         5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                       cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
                       key = TRUE, keysize = 1, density.info=c("histogram", 
                        "density", "none"), denscol = tracecol, symkey = min(x < 
                          0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
                       xlab = NULL, ylab = NULL, lmat = NULL, lhei = 0.01, lwid = NULL, 
                       ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  
  if (exists("hcc")) {retval$hcc <- hcc}
  if (exists("hcr")) {retval$hcr <- hcr}
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
      1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <-  rbind(4:3,2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors)) #|| ncol(ColSideColors) != nc) 
        stop("'ColSideColors' must be a character ") #vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      if (is.vector(ColSideColors)) nnn=1 
      else nnn=nrow(ColSideColors)
      lhei <- c(lhei[1], nnn*0.1, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors)) #|| length(RowSideColors) != nr) 
        stop("'RowSideColors' must be a character ")
      if (is.vector(RowSideColors)) nnn=1 
      else nnn=ncol(RowSideColors)
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], nnn*0.1, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    
    if (is.vector(RowSideColors)) {
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } 
    if (is.matrix(RowSideColors)) {
      jk.row = RowSideColors
      jk.xy = matrix(which(jk.row != "0"), dim(jk.row))
      colnames(jk.xy) <- colnames(jk.row)
      #image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
      #grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
      image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
      #            axis(3, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
      axis(1, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
      #   grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
    }
    
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    if (is.vector(ColSideColors)) {
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } 
    if (is.matrix(ColSideColors)) {
      jk.col = ColSideColors
      jk.xy = matrix(which(jk.col != "0"), dim(jk.col))
      rownames(jk.xy) <- rownames(jk.col)
      jk.xy=t(jk.xy)
      image(jk.xy, col = jk.col[, colInd], axes = FALSE)
      axis(2, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, tick=FALSE)
    }
    
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
    cellnote <- t(cellnote)
  }
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  # if (!invalid(na.color) & any(is.na(x))) {
  #     mmat <- ifelse(is.na(x), 1, NA)
  #    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
  #         col = na.color, add = TRUE)
  # }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  #    axis(3, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
  #        cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 3, line = margins[1] - 1.25)
  
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1], 
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 3, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x))
      tmpbreaks[length(tmpbreaks)] <- max(abs(x))
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 2)
    else { 
      #mtext(side = 1, "Value", line = 2)
      mtext(side = 1, "", line = 2)
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("") ;#Color Key and Density Plot", cex=0.25)
      #title("Color Key and Density Plot", cex=0.25)
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
      #mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      #title("Color Key and Histogram", cex=0.25)
      title("", cex=0.25)
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
      #mtext(side = 2, "Count", line = 2)
    }
    else { title("", cex=0.25)
      #title("Color Key", cex=0.25)
    }
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



plot_cl_heatmap4 <- function(norm.dat, cl, markers, prefix=NULL,hc=NULL, gene.hc=NULL,centered=FALSE,labels=names(cl),sorted=FALSE,by.cl=TRUE,ColSideColors=NULL,maxValue=5,min.sep=4,main="", height=13, width=9, legend.dat=NULL)
{
  select.cells=names(cl)
  tmp.dat = as.matrix(norm.dat[markers,select.cells,drop=F])
  if(!is.null(ColSideColors)){
    ColSideColors=ColSideColors[, select.cells,drop=F]
  }
  if(centered){
    tmp.dat = tmp.dat - rowMeans(tmp.dat)
    breaks=c(min(min(tmp.dat)-0.1,-maxValue),  seq(-maxValue,maxValue, length.out=99), max(max(tmp.dat)+1))
  }
  else{
    tmp.dat = tmp.dat/pmax(rowMaxs(tmp.dat), 1)
    breaks=c(0, seq(0.05, 1, length.out=100))
  }
  colnames(tmp.dat)=labels
  cexCol = min(70/ncol(tmp.dat),1)
  cexRow = min(60/nrow(tmp.dat),1)
  if(is.null(gene.hc)){
    gene.hc = hclust(dist(tmp.dat), method="ward.D")
  }
  if(is.null(hc) & !sorted & length(select.cells)< 2000){
    hc = hclust(dist(t(tmp.dat)), method="ward.D")
  }
  col = blue.red(150)[51:150]
  if(!is.null(prefix)){
    pdf(paste(prefix,"pdf",sep="."), height=height, width=width)
  }
  if(by.cl){
    if(sorted){
      ord = 1:length(cl)
    }
    else{
      if(!is.null(hc)){
        ord = order(cl, order(hc$order))
      }
      else{
        ord = order(cl)
      }
    }
    sep = cl[ord]
    sep=which(sep[-1]!=sep[-length(sep)])
    sep = c(sep[1], sep[which(sep[-1] - sep[-length(sep)] >=min.sep)+1])
    heatmap.4(tmp.dat[,ord],Rowv=as.dendrogram(gene.hc), Colv=NULL, col=col, trace="none", dendrogram="none", cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors[,ord],breaks=breaks,colsep=sep, sepcolor="black",main=main)
    cells.order=colnames(tmp.dat)[ord]
  }
  else{
    heatmap.4(tmp.dat,Rowv=as.dendrogram(gene.hc), Colv=as.dendrogram(hc), col=col, trace="none", dendrogram="none", cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors,breaks=breaks,main=main)
    cells.order=colnames(tmp.dat)[hc$order]
  }
  
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
 cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
 leg.clms = ceiling(length(unique(cl))/8)

 
  legend(x=0,y=0.67, 
         legend=c("Female", "Male"), 
         fill=c("hotpink1", "dodgerblue3"), 
         cex=.7, 
         title="Gender",bty="n")
  
  leg.region=as.matrix(legend.dat[[2]])
  region.col=ColSideColors["Region",]
  leg.region.u =leg.region[leg.region %in% region.col,]
  leg.region.u=as.data.frame(leg.region.u)
  leg.cols=factor(leg.region.u$leg.region.u, levels=unique(leg.region.u$leg.region.u))
  #leg.region.u= leg.region.u %>% rownames_to_column("Region")
  leg.cols= as.character(leg.cols)
  reg.clms=ceiling(length(leg.cols)/8)
  
  legend(x=0.4,y=1.0, 
         legend=rownames(leg.region.u), 
         fill=leg.cols, 
         cex=.7, 
         title="Region",
         ncol=reg.clms)
  
  #jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  cl.col = rainbow(length(unique(cl)))[as.factor(cl)]
  leg.clms = ceiling(length(unique(cl))/15)
  
  legend(x=0,y=0.6,       
         legend = unique(cl),
         fill=unique(cl.col),
         cex=.7,
         title="Cluster",
         ncol= leg.clms,
         bty="n"
  )
  
  if(!is.null(prefix)){
    dev.off()
  }
  return(cells.order)
}





display_cl4 <- function(cl, norm.dat,prefix=NULL, plot=!is.null(prefix), col=NULL, max.cl.size=NULL,markers=NULL,de.genes=NULL, main="",height=13, width=9, min.sep=10, legend.dat=NULL, ...)
{
  select.cells=names(cl)        
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if(!is.null(max.cl.size)){    
    tmp.cells = sample_cells(cl,  max.cl.size)
    cl = cl[tmp.cells]
  }
  select.cells=names(cl)
  cl.col = rainbow(length(unique(cl)))[as.factor(cl)]
  tmp.col =t(as.matrix(cl.col, ncol=1))
  colnames(tmp.col)= select.cells
  rownames(tmp.col)= c("cluster")
  if(!is.null(col)){
    tmp.col = rbind(tmp.col, col[,select.cells])
  }
  tmp.dat = as.matrix(norm.dat[,names(cl)])
  if(is.null(markers)){
    tmp = select_markers(tmp.dat,cl, de.genes=de.genes, ...)
    markers = tmp$markers
    de.genes=tmp$de.genes
  }
  cells_order=NULL
  if(plot & !is.null(markers)){
    cells_order=plot_cl_heatmap4(tmp.dat, cl, markers, ColSideColors=tmp.col,legend.dat=legend.dat, prefix=prefix, labels=NULL, by.cl=TRUE,min.sep=min.sep,main=main, height=height, width=width)
  }
  return(list(markers=markers,de.genes=de.genes, cells_order= cells_order))
}
