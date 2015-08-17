source('http://bioconductor.org/biocLite.R')
biocLite(c('affy','arrayQualityMetrics','limma','primeviewprobe'))
library(affy)
library(arrayQualityMetrics)
library(limma)
#library(primeviewprobe)

setwd("~/Syncplicity/Vit D Arrays")

##read in CEL files
fil<-readLines("vitDfiles.txt")
fil

Data<-ReadAffy(filenames=fil)
dim(Data)
dim(exprs(Data)) 

##QC pre-normalisation
#arrayQualityMetrics(expressionset=Data,outdir="QCreport1",force=T,do.logtransform=T) 

##normalise data: RMA
dat<-rma(Data,background=F)
dim(dat)
mat.rma<-exprs(dat)
mat.rma[1:4,1:4]

##QC post-normalisation
#arrayQualityMetrics(expressionset=dat,outdir="QCreport2",force=T,do.logtransform=T) 

##primeview symbols
load("~/Syncplicity/Vit D Arrays/primeviewSYMBOL.RData")
sym<-AnnotationDbi::as.list(primeviewSYMBOL) 
length(sym) 

##numbers of probes, genes 
sum(unlist(lapply(sym,length))>1)  
sum(unlist(lapply(sym,length))==1) 
length(sym) 
length(unique(sym)) 
unlistGenes<-unlist(unique(sym))
length(unlistGenes) 

sym_uni_probe<-sym[lapply(sym,length)==1]
sym_uni_probe[1:5]  
length(sym_uni_probe) 
length(unique(sym_uni_probe))  
uni2<-(unique(sym_uni_probe))
uni2[1:4] 
mat.rma<-na.omit(mat.rma[match(rownames(mat.rma), names(sym_uni_probe)),])
dim(mat.rma)  
mat.rma[1:5,1:5]

##obtain one signal per probe-set/gene
biocLite(c("WGCNA"))
library(WGCNA)

keepProbes<-intersect(rownames(mat.rma),names(sym_uni_probe)) 
length(keepProbes) 
datET<-mat.rma[match(keepProbes,rownames(mat.rma)),]
dim(datET) 
rowID<-rownames(datET)
rowGrp<-unlist(sym_uni_probe[match(rowID,names(sym_uni_probe))]) 
length(unique(rowGrp))  

collapse.object<-collapseRows(datET=datET,rowGroup=rowGrp,rowID=rowID) 
dim(collapse.object$datETcollapsed) 
collapse.object$datETcollapsed[1:10,1:5] 

datCol<-collapse.object$datETcollapsed
datCol<-datCol[-1,] 
colnames(datCol)<-gsub(".CEL","",unlist(lapply(strsplit(as.vector(colnames(datCol)),"_"),
                                               function(x) x[1]))) 

mat.rmaCollapsed<-data.frame(collapse.object$group2row,collapse.object$datETcollapsed)
mat.rmaCollapsed[1:5,1:5]
mat.rmaCollapsed<-mat.rmaCollapsed[-1,] 
mat.rmaCollapsed[1:5,1:5]
dim(mat.rmaCollapsed) 
mat.rma1<-mat.rmaCollapsed[,3:112]  
dim(mat.rma1) 
mat.rma1[1:5,1:5]

##sorting patient info for analysis
patno.celno<-read.csv("Patno_Cell file_grp matchA.csv") 

expt.dat<-data.frame(PatID=unlist(lapply(strsplit(as.vector(patno.celno[,1])," "),
                                         function(x) x[1])), 
                       Tissue=unlist(lapply(strsplit(as.vector(patno.celno[,1])," "),
                                            function(x) x[2])),
                       ArrayID=patno.celno[,2],
                       Group=patno.celno[,3])

arrayNames<-gsub("_2","",gsub(".CEL","",gsub("_(PrimeView).CEL","",
                                             colnames(datCol),fixed=T)))

edSort<-expt.dat[as.vector(na.omit(match(arrayNames,expt.dat$ArrayID))),]  

## Normal Tissue: Group 1 vs Group 2 (analysis performed pre-blinding; 
#Grp 1= placebo)
  
NBG<-edSort[na.omit(grep(rep("N",109),(edSort[,2]))),]

grp<-factor(NBG$Group,levels=c("1","2"))
design<-model.matrix(~grp)

arrayNames_A<-gsub(".CEL","",gsub("_2.CEL","",gsub("_.PrimeView..CEL","",
                                                   colnames(mat.rma1),fixed=T)))
nSort<-mat.rma1[match(NBG$ArrayID,arrayNames_A)]

cbind(NBG,colnames(nSort))
  
normFit<-eBayes(lmFit(nSort,design)) 
normTT<-topTable(normFit,coef=2,adjust="BH",n=nrow(nSort))
normTT[1:20,]

sum(normTT$adj.P.Val<0.05) 

biocLite(c("gtools"))
library("gtools")

biocLite("gplots")
library(gplots)

##heatmaps with more than one colour bar
heatmap.mik<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                       distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                          "row", 
                                                                          "column", 
                                                                          "none"), 
                       symm = FALSE, scale = c("none", 
                                                                                                                            "row", "column"), na.rm = TRUE, revC = identical(Colv, 
                                                                                                                                                                             "Rowv"), add.expr, breaks, col = "heat.colors", colsep, 
                       rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
                       notecex = 1, notecol = "cyan", na.color = par("bg"), trace = 
                         c("column", 
                                                                                      
                           "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                       vline = median(breaks), linecol = tracecol, margins = c(5, 
                                                                               5), 
                       ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                       cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
                       key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                                                   "density", 
                                                                   "none"), denscol = 
                         tracecol, symkey = min(x < 
                                                                                                                          0, na.rm = TRUE), densadj = 0.25, main = NULL, xlab = NULL, 
                       ylab = NULL, ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider 
            using only one or the other.")
  if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
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
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
    sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
    sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) 
    if (missing(col)) 
      breaks <- 16
  else breaks <- length(col) + 1
  if (length(breaks) == 1) {
    breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                  length = breaks)
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  else if (is.character(col) && length(col) == 1) 
    col <- do.call(col, list(ncol))
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x < min.breaks, min.breaks, x)
  x[] <- ifelse(x > max.breaks, max.breaks, x)
  lmat <- rbind(4:3, 2:1)
  lhei <- lwid <- c(keysize, 4)
  if (!missing(ColSideColors)) {
    #        if (!is.character(ColSideColors) || length(ColSideColors) != 
    #            nc) 
    #            stop("'ColSideColors' must be a character vector of length ncol(x)")
    if(is.null(nrow(ColSideColors))){
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    else{
      lmat <- rbind(lmat[1, ] + nrow(ColSideColors),
                    t(matrix(as.numeric(unlist(strsplit(paste("NA",1:nrow(ColSideColors)),
                                                        " "))),2,nrow(ColSideColors))),
                    lmat[2, ] + nrow(ColSideColors))
      #          lhei <- c(lhei[1], 0.2*nrow(ColSideColors), lhei[2])
      lhei <- c(lhei[1], rep(0.2,nrow(ColSideColors)), lhei[2])
    } 
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != 
          nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }
  lmat[is.na(lmat)] <- 0
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    if(is.null(nrow(ColSideColors))) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    else{
      for(kk in 1:nrow(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[kk,colInd], axes = FALSE)
        axis(4, 0, labels = rownames(ColSideColors)[kk], las = 2,
             line = -0.5, tick = 0, cex.axis = cexCol)
      }
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
    cellnote <- t(cellnote)
  }
  if (revC) {
    iy <- nr:1
    ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
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
                                                      1 - rsep) - 0.5, xright = 
                                ncol(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    for (i in colInd) {
      if (!is.null(vline)) {
        vline.vals <- scale01(vline, min.scale, max.scale)
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
    for (i in rowInd) {
      if (!is.null(hline)) {
        hline.vals <- scale01(hline, min.scale, max.scale)
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
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    if (symkey) {
      max.raw <- max(abs(x), na.rm = TRUE)
      min.raw <- -max.raw
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  invisible(list(rowInd = rowInd, colInd = colInd))
}

##heatmap of top 100 DEG- output sorted for thesis
KPnA<-rownames(normTT)[1:100]
yA<-t(apply(nSort[match(KPnA,rownames(nSort)),],
            1,function(x)(x-mean(x))/sd(x)))

yA[yA< -3]<- -3
yA[yA> 3]<-3
yA[1:5,1:5]

cA<-ifelse(design[,2]==1,"black",gray(0.6))
cutree<-cutree(hclust(dist(t(yA))),5)
rowLabKPnA<-as.vector(KPnA)
col_matrix_cA<-rbind(cA,greenred(5)[cutree])
rownames(col_matrix_cA)<-c("Treatmt/Control","Cluster")

heatmap.mik(yA,trace='none',scale='none',dendrogram="both",
            ColSideColors=col_matrix_cA,
            col=greenred(50),
            labRow=paste("Gene",rowLabKPnA),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(yA))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

##read in VDRE gene list
vdre<-readLines("VDRE_comb_080615.txt") # 3136; including orig list, 
#Meyer & Wang lists

##matching VDRE genes with primeview symbols/probes
vdre.probes<-list()
for(i in 1:length(vdre))vdre.probes[[i]]<-names(sym)[sym==vdre[i]]

names(vdre.probes)<-vdre 
allVDREprobes<-unlist(vdre.probes) # 5958 
kpComb<-intersect(allVDREprobes,rownames(mat.rma))   
length(kpComb) # 5851 

datET1_Comb<-na.omit(mat.rma[match(kpComb,rownames(mat.rma)),])

sum(unlist(lapply(allVDREprobes,length))>1) # 0
sum(unlist(lapply(allVDREprobes,length))==1) # 5958

##differentially expressed VDRE genes in normal tissue between grps
datVDRE<-na.omit(normTT[match(vdre,rownames(normTT)),])

VDREpVal<-datVDRE$P.Value
VDREp.adj<-round(p.adjust(sort(VDREpVal),"BH"),4)
VDREp.adj[1:10] # all 0.9999

##GSA
biocLite(c('org.Hs.eg.db','multtest','safe'))
library('org.Hs.eg.db')
library('multtest')
library('safe')

KEGG.list<-AnnotationDbi::as.list(org.Hs.egPATH)
GO.list<-AnnotationDbi::as.list(org.Hs.egGO)
GO.list2<-list()
for(i in 1:length(GO.list))GO.list2[[i]]<-names(GO.list[[i]])
GO.listx<-list()
for(i in 1:length(GO.list))GO.listx[[i]]<-GO.list[[i]][[1]][[1]]
names(GO.listx)<-names(GO.list)
SYMBOL.list<-AnnotationDbi::as.list(org.Hs.egSYMBOL)
names(KEGG.list)<-unlist(SYMBOL.list)
names(GO.list2)<-SYMBOL.list[1:length(GO.list2)]

##KEGG
biocLite(c('limma'))
library(limma)

C.matrixN<-getCmatrix(gene.list=KEGG.list,present.genes=rownames(normTT),
                      min.size=10,as.matrix=T)  
geneSetTest<-c()
for(i in 1:ncol(C.matrixN)){
  geneSetTest[i]<-geneSetTest(seq(1:nrow(C.matrixN))[C.matrixN[,i]==1],
                              statistics=as.numeric(normTT[,3]),ranks.only=T)
}
names(geneSetTest)<-colnames(C.matrixN)
geneSetTestadj<-p.adjust(geneSetTest,method="fdr")

biocLite(c('KEGG.db'))
library('KEGG.db')

path<-AnnotationDbi::as.list(KEGGPATHID2NAME)
names(geneSetTestadj)<-path[match(names(geneSetTestadj),names(path))]
head(sort(geneSetTestadj),n=20) 
HD<-head(sort(geneSetTestadj),n=20)

pway_size<-apply(C.matrixN,2,sum) 
pway_table<-cbind(names(pway_size),
                  names(geneSetTestadj),
                  geneSetTest,
                  geneSetTestadj,
                  pway_size)
colnames(pway_table)<-c("KEGG_id","pathway","p-value","adj p-value(fdr)",
                        "pathway size")
pway_table<-pway_table[order(as.numeric(pway_table[,4])),]
head(pway_table)

pway_table[1:12,] # 25 for p-values <0.01

##'met of xeno by cyt P450' p/way- barcode + mg plots  
match('Metabolism of xenobiotics by cytochrome P450',names(geneSetTestadj))
sum(C.matrixN[,65]) # 60
barcodeplot(normTT$t,C.matrixN[,65]==1) 

C.matrixN_HD<-C.matrixN[,match(names(HD),path[match(colnames(C.matrixN),
                                                    names(path))])] 
colnames(C.matrixN_HD)<-names(HD)
HD_Genes_KEGGmat<-apply(C.matrixN_HD,2,function(x)rownames(C.matrixN_HD)[x==1])
metXeno_P450<-HD_Genes_KEGGmat$`Metabolism of xenobiotics by cytochrome P450` 

vdre[na.omit(match(metXeno_P450,vdre))] # 19: "AKR1C2"  "CYP2B6"  "CYP2C9"  
#"CYP3A4"  "AKR1C1"  "CYP2C19" "CYP3A7"  "CYP2S1"  "GSTK1"   "EPHX1"  
#"CYP3A5"  "CYP3A43" "UGT1A1"  "MGST1"   "GSTM3"   "GSTP1"   "CYP1A1"  
#"AKR1C4"  "GSTA4"

##VDRE genes involved in p/ways
#drug<-HD_Genes_KEGGmat$`Drug metabolism - cytochrome P450`
#vdre[na.omit(match(drug,vdre))] # 15; "FMO2"    "CYP2B6"  "CYP2C9"  "CYP3A4"  
#"CYP2C19" "MAOA"    "CYP3A7"  "GSTK1"   "CYP3A5"  "CYP3A43"
#"UGT1A1"  "MGST1"   "GSTM3"   "GSTP1"   "GSTA4"
#vit<-HD_Genes_KEGGmat$`Vitamin digestion and absorption`
#vdre[na.omit(match(vit,vdre))] # "SLC19A3" "APOB"  "LMBRD1"

y1D<-t(apply(nSort[match(metXeno_P450,rownames(nSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1D[y1D< -3]<- -3
y1D[y1D> 3]<-3
y1D[1:5,1:5]

mgD<-svd(y1D)$v[,1]

NBG1<-factor(NBG$Group,levels=c(1,2),labels=c("Control","Treatment"))
boxplot(-mgD~NBG1,ylab="")
t.test(mgD~NBG1) # p-value= 0.6981
legend(0.45,0.25,c("p-value: 0.70"),box.col="white")

##FA met p/way- barcode + mg plots 
match('Fatty acid metabolism',names(geneSetTestadj))
sum(C.matrixN[,8]) # 40
barcodeplot(normTT$t,C.matrixN[,8]==1) 

HD_Genes_KEGGmat<-apply(C.matrixN_HD,2,function(x)rownames(C.matrixN_HD)
                        [x==1]) 
FAmet<-HD_Genes_KEGGmat$`Fatty acid metabolism`

vdre[na.omit(match(FAmet,vdre))] # "ACADM" "ACSL1"

y1B<-t(apply(nSort[match(FAmet,rownames(nSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1B[y1B< -3]<- -3
y1B[y1B> 3]<-3
y1B[1:5,1:5]

mgB<-svd(y1B)$v[,1]

#following (#) unnecessary as know in order
#NdatSort<-gsub(".CEL","",gsub("_.PrimeView..CEL","",colnames(nSort),
#fixed=T))
#NBGpair<-NBG[match(NdatSort,NBG$ArrayID),]
#cbind(NBG,NdatSort)

boxplot(-mgB~NBG1,ylab="")
t.test(mgB~NBG1) # p-value= 0.7286
legend(0.50,-0.38,c("p-value: 0.73"),box.col="white")

##GO
G.matrixN<-getCmatrix(gene.list=GO.list2,present.genes=rownames(normTT),
                      min.size=10,as.matrix=T) 
geneSetTestN<-c()
for(i in 1:ncol(G.matrixN)){
  geneSetTestN[i]<-geneSetTest(seq(1:nrow(G.matrixN))[G.matrixN[,i]==1],
                               statistic=as.numeric(normTT[,3]))
}
names(geneSetTestN)<-colnames(G.matrixN)
geneSetTestNadj<-p.adjust(geneSetTestN,method="fdr")
names(geneSetTestNadj)<-colnames(G.matrixN)

biocLite(c('GO.db'))
library(GO.db)

path2<-AnnotationDbi::as.list(GOTERM)
pathx<-list()
for(i in 1:length(path2)) pathx[[i]]<-path2[[i]]@Term
names(pathx)<-names(path2)
names(geneSetTestNadj)<-pathx[match(names(geneSetTestNadj),names(path2))]
head(sort(geneSetTestNadj),n=20) 
HD1<-head(sort(geneSetTestNadj),n=20)

pway_sizeA<-apply(G.matrixN,2,sum) 
pway_tableA<-cbind(names(pway_sizeA),
                     names(geneSetTestNadj),
                     geneSetTestN,
                     geneSetTestNadj,
                     pway_sizeA)
colnames(pway_tableA)<-c("GO_id","pathway","p-value","adj p-value (fdr)",
                         "pathway size")
pway_tableA<-pway_tableA[order(as.numeric(pway_tableA[,4])),]
head(pway_tableA)

pway_tableA[1:20,]

##'neg reg of growth' p/way- barcode + mg plots
match('negative regulation of growth',names(geneSetTestNadj)) 
sum(G.matrixN[,2225]) # 17
barcodeplot(normTT$t,G.matrixN[,2225]==1)

G.matrixN_HD<-G.matrixN[,match(names(HD1),pathx[match(colnames(G.matrixN),names(pathx))])] 
colnames(G.matrixN_HD)<-names(HD1)
HD1_Genes_KEGGmat<-apply(G.matrixN_HD,2,function(x)rownames(G.matrixN_HD)[x==1])
NRg<-HD1_Genes_KEGGmat$`negative regulation of growth` 

vdre[na.omit(match(NRg,vdre))] # "PNPT1" "HIF1A"

y1Fa<-t(apply(nSort[match(NRg,rownames(nSort)),],
              1,function(x)(x-mean(x))/sd(x)))

y1Fa[y1Fa< -3]<- -3
y1Fa[y1Fa> 3]<-3
y1Fa[1:5,1:5]

mgFa<-svd(y1Fa)$v[,1]

boxplot(-mgFa~NBG1,ylab="")
t.test(mgFa~NBG1) # p-value= 0.1837
legend(0.4,-0.18,c("p-value: 0.18"),box.col="white")

## Tumour Tissue: Group 1 vs Group 2

TBG<-edSort[na.omit(grep(rep("T",109),(edSort[,2]))),]

grp1<-factor(TBG$Group,levels=c("1","2"))
design1<-model.matrix(~grp1)

arrayNamesA<-gsub(".CEL","",gsub("_2.CEL","",gsub("_.PrimeView..CEL","",
                                                  colnames(mat.rma1),fixed=T)))
tSort<-mat.rma1[match(TBG$ArrayID,arrayNamesA)]

cbind(TBG,colnames(tSort))

tumFit<-eBayes(lmFit(tSort,design1)) 
tumTT<-topTable(tumFit,coef=2,adjust="BH",n=nrow(tSort))
tumTT[1:20,] 

sum(tumTT$adj.P.Val<0.05)

KPt<-rownames(tumTT)[1:100]
yC<-t(apply(tSort[match(KPt,rownames(tSort)),],
            1,function(x)(x-mean(x))/sd(x)))

yC[yC< -3]<- -3
yC[yC> 3]<-3
yC[1:5,1:5]

cC<-ifelse(design1[,2]==1,"black",gray(0.6))
cutreeC<-cutree(hclust(dist(t(yC))),2)
rowLabKPt<-as.vector(KPt)
col_matrix_cC<-rbind(cC,greenred(2)[cutreeC])
rownames(col_matrix_cC)<-c("Treatmt/Control","Cluster")

heatmap.mik(yC,trace='none',scale='none',dendrogram="both",
            ColSideColors=col_matrix_cC,
            col=greenred(50),
            labRow=paste("Gene",rowLabKPt),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(yC))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

datVDRE1<-na.omit(tumTT[match(vdre,rownames(tumTT)),])

VDREpVal1<-datVDRE1$P.Value
VDREp.adj1<-round(p.adjust(sort(VDREpVal1),"BH"),4)
VDREp.adj1[1:10] # most 0.9780 & others up to 0.9995

##KEGG
C.matrixT<-getCmatrix(gene.list=KEGG.list,present.genes=rownames(tumTT),
                      min.size=10,as.matrix=T) 
geneSetTestT<-c()
for(i in 1:ncol(C.matrixT)){
  geneSetTestT[i]<-geneSetTest(seq(1:nrow(C.matrixT))[C.matrixT[,i]==1],
                               statistics=as.numeric(tumTT[,3]),ranks.only=T)
}
names(geneSetTestT)<-colnames(C.matrixT)
geneSetTestTadj<-p.adjust(geneSetTestT,method="fdr")

names(geneSetTestTadj)<-path[match(names(geneSetTestTadj),names(path))]
head(sort(geneSetTestTadj),n=20) 
HD4<-head(sort(geneSetTestTadj),n=20)

pway_sizeB<-apply(C.matrixT,2,sum) 
pway_tableB<-cbind(names(pway_sizeB),
                     names(geneSetTestTadj),
                     geneSetTestT,
                     geneSetTestTadj,
                     pway_sizeB)
colnames(pway_tableB)<-c("KEGG_id","pathway","p-value","adj p-value (fdr)",
                         "pathway size")
pway_tableB<-pway_tableB[order(as.numeric(pway_tableB[,4])),]
head(pway_tableB)

pway_tableB[1:20,]

##FA met (T)
match('Fatty acid metabolism',names(geneSetTestTadj)) 
sum(C.matrixT[,8]) # 40
barcodeplot(tumTT$t,C.matrixT[,8]==1)
 
C.matrixT_HD<-C.matrixT[,match(names(HD4),path[match(colnames(C.matrixT),
                                                     names(path))])] 
colnames(C.matrixT_HD)<-names(HD4)
HD4_Genes_KEGGmat<-apply(C.matrixT_HD,2,function(x)rownames(C.matrixT_HD)[x==1]) 
FAmet1<-HD4_Genes_KEGGmat$`Fatty acid metabolism`

vdre[na.omit(match(FAmet1,vdre))] # "ACADM" "ACSL1"

y1G<-t(apply(tSort[match(FAmet1,rownames(tSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1G[y1G< -3]<- -3
y1G[y1G> 3]<-3
y1G[1:5,1:5]

mgG<-(-1)*svd(y1G)$v[,1]
mg.rnkG<-rank(-mgG)/length(mgG)
mg.colG<-greenred(length(mgG))[rank(-mgG)]

TBG1<-factor(TBG$Group,levels=c(1,2),labels=c("Control","Treatment"))
boxplot(-mgG~TBG1,ylab="")
t.test(mgG~TBG1) # p-value= 0.0320
legend(0.4,-0.26,c("p-value: 0.03"),box.col="white")

CG<-ifelse(design1[,2]==1,"black",gray(0.6))
rowLabG<-as.vector(FAmet1)
col_matrixG<-rbind(mg.colG,CG)
rownames(col_matrixG)<-c("Metagene","Treatmt/Control")

heatmap.mik(y1G,trace='none',scale='none',dendrogram="both",
            ColSideColors=col_matrixG,
            col=greenred(50),
            labRow=paste("Gene",rowLabG),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(y1G))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

##Oxidative phos
match("Oxidative phosphorylation",names(geneSetTestTadj)) 
sum(C.matrixT[,12]) # 120
barcodeplot(tumTT$t,C.matrixT[,12]==1)
  
HD4_Genes_KEGGmat<-apply(C.matrixT_HD,2,function(x)rownames(C.matrixT_HD)[x==1]) 
OP<-HD4_Genes_KEGGmat$`Oxidative phosphorylation`

vdre[na.omit(match(OP,vdre))] # 11; "ATP5O"   "COX4I1"  "NDUFS8"  "COX5A"   "UQCRFS1" 
#"UQCRC1"  "NDUFB10" "COX8A"   "ATP5G2"  "ATP5L"  "TCIRG1"

y1L<-t(apply(tSort[match(OP,rownames(tSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1L[y1L< -3]<- -3
y1L[y1L> 3]<-3
y1L[1:5,1:5]

mgL<-svd(y1L)$v[,1]

boxplot(-mgL~TBG1,ylab="")
t.test(mgL~TBG1) # p-value= 0.1267
legend(0.45,-0.35,c("p-value: 0.13"),box.col="white")

##GO
G.matrixT<-getCmatrix(gene.list=GO.list2,present.genes=rownames(tumTT),
                      min.size=10,as.matrix=T) 
geneSetTestT1<-c()
for(i in 1:ncol(G.matrixT)){
  geneSetTestT1[i]<-geneSetTest(seq(1:nrow(G.matrixT))[G.matrixT[,i]==1],
                                statistic=as.numeric(tumTT[,3]))
}
names(geneSetTestT1)<-colnames(G.matrixT)
geneSetTestT1adj<-p.adjust(geneSetTestT1,method="fdr")
names(geneSetTestT1adj)<-colnames(G.matrixT)

names(geneSetTestT1adj)<-pathx[match(names(geneSetTestT1adj),names(path2))]
head(sort(geneSetTestT1adj),n=20)
HD5<-head(sort(geneSetTestT1adj),n=20)

pway_sizeC<-apply(G.matrixT,2,sum) 
pway_tableC<-cbind(names(pway_sizeC),
                     names(geneSetTestT1adj),
                     geneSetTestT1,
                     geneSetTestT1adj,
                     pway_sizeC)
colnames(pway_tableC)<-c("GO_id","pathway","p-value","adj p-value (fdr)",
                         "pathway size")
pway_tableC<-pway_tableC[order(as.numeric(pway_tableC[,4])),]
head(pway_tableC)

pway_tableC[1:20,]

##FA beta-oxid
match("fatty acid beta-oxidation",names(geneSetTestT1adj)) 
sum(G.matrixT[,683]) # 33
barcodeplot(tumTT$t,G.matrixT[,683]==1)
 
G.matrixT_HD<-G.matrixT[,match(names(HD5),pathx[match(colnames(G.matrixT),
                                                      names(pathx))])] 
colnames(G.matrixT_HD)<-names(HD5)
HD5_Genes_KEGGmat<-apply(G.matrixT_HD,2,function(x)rownames(G.matrixT_HD)
                         [x==1]) 
FAB<-HD5_Genes_KEGGmat$`fatty acid beta-oxidation`

vdre[na.omit(match(FAB,vdre))] # "ACADM"    "ABCD1"    "SLC25A17" "MUT"      
#"ABCD2"    "PEX2"

y1R<-t(apply(tSort[match(FAB,rownames(tSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1R[y1R< -3]<- -3
y1R[y1R> 3]<-3
y1R[1:5,1:5]

mgR<-(-1)*svd(y1R)$v[,1]
mg.rnkR<-rank(-mgR)/length(mgR)
mg.colR<-greenred(length(mgR))[rank(-mgR)]

boxplot(-mgR~TBG1,ylab="")
t.test(mgR~TBG1) # p-value= 0.02792
legend(0.35,-0.27,c("p-value: 0.028"),box.col="white")

CR<-ifelse(design1[,2]==1,"black",gray(0.6))
rowLabR<-as.vector(FAB)
col_matrixR<-rbind(mg.colR,CR)
rownames(col_matrixR)<-c("Metagene","Treatmt/Control")

heatmap.mik(y1R,trace='none',scale='none',dendrogram="both",
            ColSideColors=col_matrixR,
            col=greenred(50),
            labRow=paste("Gene",rowLabR),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(y1R))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

##resp e- trans chain
match("respiratory electron transport chain",names(geneSetTestT1adj)) 
sum(G.matrixT[,1405]) # 86
barcodeplot(tumTT$t,G.matrixT[,1405]==1)
 
HD5_Genes_KEGGmat<-apply(G.matrixT_HD,2,function(x)rownames(G.matrixT_HD)[x==1]) 
RETC1<-HD5_Genes_KEGGmat$`respiratory electron transport chain`

vdre[na.omit(match(RETC1,vdre))] # 9: "ATP5O"   "COX4I1"  "NDUFS8"  "COX5A"   
#"UQCRFS1" "UQCRC1"  "NDUFB10" "COX8A"   "ATP5L" 

y1T<-t(apply(tSort[match(RETC1,rownames(tSort)),],
             1,function(x)(x-mean(x))/sd(x)))

y1T[y1T< -3]<- -3
y1T[y1T> 3]<-3
y1T[1:5,1:5]

mgT<-svd(y1T)$v[,1]

boxplot(-mgT~TBG1,ylab="")
t.test(mgT~TBG1) # p-value= 0.0894
legend(0.42,-0.32,c("p-value: 0.09"),box.col="white")

## Paired Tissue: Group 1 vs Group 2
#sorting data 
Pairing<-gsub("VDS","",expt.dat[,1])
Pairing<-data.frame(Pairing)
pair<-data.frame(PatID=unlist(lapply(strsplit(as.vector(patno.celno[,1])," "),
                                     function(x) x[1])),
                   Tissue=unlist(lapply(strsplit(as.vector(patno.celno[,1])," "),
                                        function(x) x[2])),
                   Pairing=Pairing,
                   ArrayID=patno.celno[,2],
                   Group=patno.celno[,3])

arrayNames<-gsub("_2","",gsub(".CEL","",gsub("_(PrimeView).CEL","",
                                             colnames(datCol),fixed=T)))
pairSort<-pair[as.vector(na.omit(match(arrayNames,pair$ArrayID))),]
pairSort<-pairSort[-match("VD56",pairSort[,4]),] 
pairSort<-pairSort[-match("VD85",pairSort[,4]),]
pairSort<-pairSort[-match("VD104",pairSort[,4]),]
pairSort$Tissue<-gsub("T1","T",pairSort[,2])

drp<-colnames(datCol)[is.na(match(colnames(datCol),pairSort$ArrayID))] 
drp 
datColx<-datCol[,-match(drp,colnames(datCol))] 
dim(datColx)

cat(match(colnames(datColx),pairSort$ArrayID),"\n") 

dim(pairSort) 

pairSort[match(names(table(pairSort$PatID))[table(pairSort$PatID)==1],pairSort[,1]),]  
drp1<-as.vector(pairSort[match(names(table(pairSort$PatID))[table(pairSort$PatID)==1],
                               pairSort[,1]),4]) 
drp1 

datColxx<-datColx[,-match(drp1,colnames(datColx))] 
pairSortxx<-pairSort[-match(drp1,as.vector(pairSort$ArrayID)),] 
dim(datColxx)
dim(pairSortxx) 

normDat3<-datColxx[,grep("N",pairSortxx[,2])]
tumDat1<-datColxx[,grep("T",pairSortxx[,2])] 
colnames(normDat3)<-unique(pairSortxx[,1]) 
colnames(tumDat1)<-unique(pairSortxx[,1]) 
diff<-normDat3-tumDat1 
diff[1:5,1:5]

arrayPair<-pairSortxx[as.vector(na.omit(match(colnames(diff),pairSortxx$PatID))),]
dim(arrayPair)

cbind(arrayPair,colnames(diff)) 

GRP<-arrayPair[,5] 
design2<-model.matrix(~GRP)
rownames(design2)<-colnames(diff)

grpFit<-eBayes(lmFit(diff, design2))
grpTT<-topTable(grpFit,coef=2,adjust="BH",n=nrow(diff)) 
grpTT[1:20,]

sum(grpTT$adj.P.Val<0.05)  

##Examining the N-T intercept for comparison with the paired tissue findings
grpTTA<-topTable(grpFit,coef=1,adjust="BH",n=nrow(diff))
grpTTA[1:20,]

##KEGG- for N-T intercept
C.matrixGA<-getCmatrix(gene.list=KEGG.list,present.genes=rownames(grpTTA),
                       min.size=10,as.matrix=T) 
geneSetTestGA<-c()
for(i in 1:ncol(C.matrixGA)){
  geneSetTestGA[i]<-geneSetTest(seq(1:nrow(C.matrixGA))[C.matrixGA[,i]==1],
                               statistics=as.numeric(grpTTA[,3]),ranks.only=T)
}
names(geneSetTestGA)<-colnames(C.matrixGA)
geneSetTestGAadj<-p.adjust(geneSetTestGA,method="fdr")

names(geneSetTestGAadj)<-path[match(names(geneSetTestGAadj),names(path))]
head(sort(geneSetTestGAadj),n=10)
HD7A<-head(sort(geneSetTestGAadj),n=10) 

pway_sizeDA<-apply(C.matrixGA,2,sum) 
pway_tableDA<-cbind(names(pway_sizeDA),
                     names(geneSetTestGAadj),
                     geneSetTestGA,
                     geneSetTestGAadj,
                     pway_sizeDA)
colnames(pway_tableDA)<-c("KEGG_id","pathway","p-value","adj p-value(fdr)",
                          "pathway size")
pway_tableDA<-pway_tableDA[order(as.numeric(pway_tableDA[,4])),]
head(pway_tableDA)

pway_tableDA[1:14,] 

##GO- for N-T intercept
G.matrixGA<-getCmatrix(gene.list=GO.list2,present.genes=rownames(grpTTA),
                       min.size=10,as.matrix=T) 
geneSetTestGA1<-c()
for(i in 1:ncol(G.matrixGA)){
  geneSetTestGA1[i]<-geneSetTest(seq(1:nrow(G.matrixGA))[G.matrixGA[,i]==1],
                                 statistic=as.numeric(grpTTA[,3]))
} 
names(geneSetTestGA1)<-colnames(G.matrixGA)
geneSetTestGA1adj<-p.adjust(geneSetTestGA1,method="fdr")
names(geneSetTestGA1adj)<-colnames(G.matrixGA)

names(geneSetTestGA1adj)<-pathx[match(names(geneSetTestGA1adj),names(path2))]
head(sort(geneSetTestGA1adj),n=20)
HD50A<-head(sort(geneSetTestGA1adj),n=20)

pway_sizeEA<-apply(G.matrixGA,2,sum) 
pway_tableEA<-cbind(names(pway_sizeEA),
                     names(geneSetTestGA1adj),
                     geneSetTestGA1,
                     geneSetTestGA1adj,
                     pway_sizeEA)
colnames(pway_tableEA)<-c("GO_id","pathway","p-value","adj p-value (fdr)",
                          "pathway size")
pway_tableEA<-pway_tableEA[order(as.numeric(pway_tableEA[,4])),]
head(pway_tableEA)

pway_tableEA[1:52,]
#pway_tableEA[1:15,] # for thesis

##heatmap of top 100 genes for paired analysis
KPdg<-rownames(grpTT)[1:100]
yE<-t(apply(diff[match(KPdg,rownames(diff)),],
            1,function(x)(x-mean(x))/sd(x)))

yE[yE< -3]<- -3
yE[yE> 3]<-3
yE[1:5,1:5]

cE<-ifelse(design2[,2]==1,"black",gray(0.6))
cutreeE<-cutree(hclust(dist(t(yE))),5)
rowLabKPdg<-as.vector(KPdg)
col_matrix_cE<-rbind(cE,greenred(5)[cutreeE])
rownames(col_matrix_cE)<-c("Treatmt/Control","Cluster")

heatmap.mik(yE,trace='none',scale='none',dendrogram="both",
            ColSideColors=col_matrix_cE,
            col=greenred(50),
            labRow=paste("Gene",rowLabKPdg),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(yE))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

##VDRE gene list analysis
datVDRE2<-na.omit(grpTT[match(vdre,rownames(grpTT)),])

VDREpVal2<-datVDRE2$P.Value
VDREp.adj2<-round(p.adjust(sort(VDREpVal2),"BH"),4)
VDREp.adj2[1:10] # lots 0.9551 & other up to 0.9994

##KEGG
C.matrixG<-getCmatrix(gene.list=KEGG.list,present.genes=rownames(grpTT),
                      min.size=10,as.matrix=T) 
geneSetTestG<-c()
for(i in 1:ncol(C.matrixG)){
  geneSetTestG[i]<-geneSetTest(seq(1:nrow(C.matrixG))[C.matrixG[,i]==1],
                               statistics=as.numeric(grpTT[,3]),ranks.only=T)
}
names(geneSetTestG)<-colnames(C.matrixG)
geneSetTestGadj<-p.adjust(geneSetTestG,method="fdr")

names(geneSetTestGadj)<-path[match(names(geneSetTestGadj),names(path))]
head(sort(geneSetTestGadj),n=10)
HD7<-head(sort(geneSetTestGadj),n=10)
 
pway_sizeD<-apply(C.matrixG,2,sum) 
pway_tableD<-cbind(names(pway_sizeD),
                     names(geneSetTestGadj),
                     geneSetTestG,
                     geneSetTestGadj,
                     pway_sizeD)
colnames(pway_tableD)<-c("KEGG_id","pathway","p-value","adj p-value (fdr)",
                         "pathway size")
pway_tableD<-pway_tableD[order(as.numeric(pway_tableD[,4])),]
head(pway_tableD)

pway_tableD[1:20,]

##ribosome
match('Ribosome',names(geneSetTestGadj)) 
sum(C.matrixG[,72]) # 86
barcodeplot(grpTT$t,C.matrixG[,72]==1)
  
C.matrixG_HD<-C.matrixG[,match(names(HD7),path[match(colnames(C.matrixG),
                                                     names(path))])] 
colnames(C.matrixG_HD)<-names(HD7)
HD7_Genes_KEGGmat<-apply(C.matrixG_HD,2,function(x)rownames(C.matrixG_HD)
                         [x==1]) 
rib<-HD7_Genes_KEGGmat$`Ribosome`

vdre[na.omit(match(rib,vdre))] # 7; "RPL22L1" "RPLP0"   "RPS19"   "RPL10A"  
#"RPL3L" "RPL37A"  "RPS29"

y1V<-t(apply(diff[match(rib,rownames(diff)),],
             1,function(x)(x-mean(x))/sd(x)))

y1V[y1V< -3]<- -3
y1V[y1V> 3]<-3
y1V[1:5,1:5]

mgV<-svd(y1V)$v[,1]

arrayPair<-pairSortxx[as.vector(na.omit(match(colnames(diff),pairSortxx$PatID))),]
arrayPair$Group<-factor(arrayPair$Group,levels=c("1","2"),labels=c("Control","Treatment"))
boxplot(mgV~arrayPair$Group,ylab="")
t.test(mgV~arrayPair$Group) # p-value= 0.2657
legend(0.37,0.4,c("p-value: 0.27"),box.col="white")

rowLabV<-as.vector(rib)
CV<-ifelse(design2[,2]==1,"black",gray(0.6))

heatmap.mik(y1V,trace='none',scale='none',dendrogram="both",
            ColSideColors=CV,
            col=greenred(50),
            labRow=paste("Gene",rowLabV),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(y1V))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))

##GO
G.matrixG<-getCmatrix(gene.list=GO.list2,present.genes=rownames(grpTT),
                      min.size=10,as.matrix=T) 
geneSetTestG1<-c()
for(i in 1:ncol(G.matrixG)){
  geneSetTestG1[i]<-geneSetTest(seq(1:nrow(G.matrixG))[G.matrixG[,i]==1],
                                statistic=as.numeric(grpTT[,3]))
} 
names(geneSetTestG1)<-colnames(G.matrixG)
geneSetTestG1adj<-p.adjust(geneSetTestG1,method="fdr")
names(geneSetTestG1adj)<-colnames(G.matrixG)

names(geneSetTestG1adj)<-pathx[match(names(geneSetTestG1adj),names(path2))]
head(sort(geneSetTestG1adj),n=20)
HD50<-head(sort(geneSetTestG1adj),n=20)

pway_sizeE<-apply(G.matrixG,2,sum) 
pway_tableE<-cbind(names(pway_sizeE),
                     names(geneSetTestG1adj),
                     geneSetTestG1,
                     geneSetTestG1adj,
                     pway_sizeE)
colnames(pway_tableE)<-c("GO_id","pathway","p-value","adj p-value (fdr)",
                         "pathway size")
pway_tableE<-pway_tableE[order(as.numeric(pway_tableE[,4])),]
head(pway_tableE)

pway_tableE[1:20,]

##trans termination
match('translational termination',names(geneSetTestG1adj)) 
sum(G.matrixG[,645]) # 83
barcodeplot(grpTT$t,G.matrixG[,645]==1)

G.matrixG_HD<-G.matrixG[,match(names(HD50),pathx[match(colnames(G.matrixG),
                                                       names(pathx))])] 
colnames(G.matrixG_HD)<-names(HD50)
HD50_Genes_KEGGmat<-apply(G.matrixG_HD,2,function(x)rownames(G.matrixG_HD)
                          [x==1]) 
TT<-HD50_Genes_KEGGmat$`translational termination`

vdre[na.omit(match(TT,vdre))] # 8; "ETF1"   "RPLP0"  "RPS19"  "RPL10A" "RPL3L"  
#"RPL37A" "GSPT2"  "RPS29" 

yTT<-t(apply(diff[match(TT,rownames(diff)),],
             1,function(x)(x-mean(x))/sd(x)))

yTT[yTT< -3]<- -3
yTT[yTT> 3]<-3
yTT[1:5,1:5]

mgTT<-svd(yTT)$v[,1]

boxplot(mgTT~arrayPair$Group,ylab="")
t.test(mgTT~arrayPair$Group) # p-value= 0.2826
legend(0.4,-0.3,c("p-value: 0.28"),box.col="white")

rowLabTT<-as.vector(TT)
Ctt<-ifelse(design2[,2]==1,"black",gray(0.6))

heatmap.mik(yTT,trace='none',scale='none',dendrogram="both",
            ColSideColors=Ctt,
            col=greenred(50),
            labRow=paste("Gene",rowLabTT),
            labCol=c(paste("Sample  ",1:9),paste("Sample",10:ncol(yTT))),
            adjCol=(c(1,0.5)),
            cexRow=1.5,cexCol=1.5,
            margins=c(10,10))  

##proliferation mg across grps
prolif<-readLines("Updated prolif gene names.txt") #38
prolMat<-na.omit(tSort[match(prolif,rownames(tSort)),]) #35 54

yyP<-t(apply(prolMat,1,function(x)(x-mean(x))/sd(x)))

yyP[yyP< -3]<- -3
yyP[yyP> 3]<-3
yyP[1:5,1:5]

mgProl<-svd(yyP)$v[,1]

prolMGsort<-factor(TBG$Group,levels=c("1","2"),labels=c("Control","Treatment"))
boxplot(mgProl~prolMGsort,ylab="Proliferation Metagene")
t.test(mgProl~prolMGsort) # 0.9236
legend(0.4,0.2,c("p-value: 0.92"),box.col="white")
